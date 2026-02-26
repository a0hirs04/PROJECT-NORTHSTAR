from __future__ import annotations

import struct
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Callable, Dict

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import SimulationMetrics


def _write_mat_v4(path: Path, matrices: Dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as f:
        for name, matrix in matrices.items():
            arr = np.asarray(matrix, dtype=np.float64)
            if arr.ndim != 2:
                raise ValueError(f"MAT-v4 writer expects 2D array for {name}; got shape {arr.shape}")

            name_bytes = name.encode("utf-8") + b"\x00"
            header = struct.pack(
                "<iiiii",
                0,  # type_code: double
                int(arr.shape[0]),
                int(arr.shape[1]),
                0,  # imagf
                len(name_bytes),
            )
            f.write(header)
            f.write(name_bytes)
            f.write(np.asfortranarray(arr).astype("<f8", copy=False).tobytes(order="F"))


def _build_cell_matrix(
    *,
    dead_second_tumor: bool,
    stromal_activation: tuple[float, float],
    tumor_hif1a: tuple[float, float],
    tumor_drug_sensitivity: tuple[float, float],
) -> np.ndarray:
    # Rows: [cell_type, dead, current_death_model, is_activated, is_mesenchymal,
    # drug_sensitivity, HIF1A, position(x,y,z)]
    arr = np.zeros((10, 4), dtype=np.float64)
    arr[0, :] = [0, 0, 1, 1]  # 2 tumor, 2 stromal
    arr[1, :] = [0, 1 if dead_second_tumor else 0, 0, 0]
    arr[2, :] = [0, 100 if dead_second_tumor else 0, 0, 0]
    arr[3, :] = [0, 0, stromal_activation[0], stromal_activation[1]]
    arr[4, :] = [0.9, 0.1, 0, 0]  # one mesenchymal tumor cell
    arr[5, :] = [tumor_drug_sensitivity[0], tumor_drug_sensitivity[1], 0.2, 0.2]
    arr[6, :] = [tumor_hif1a[0], tumor_hif1a[1], 0.1, 0.1]
    # Tumor cells at (0,0) and (100,0) -> extent=100
    arr[7, :] = [0, 100, 150, 200]
    arr[8, :] = [0, 0, 100, 200]
    arr[9, :] = [0, 0, 0, 0]
    return arr


def _build_micro_matrix(ecm_values: list[float], drug_values: list[float]) -> np.ndarray:
    # 6 voxels. Rows:
    # 0:x, 1:y, 2:z, 3:voxel_volume, 4:oxygen, 5:ecm_density, 6:drug
    x = np.array([0, 120, 150, 250, -100, 400], dtype=np.float64)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    voxel_volume = np.ones_like(x)
    oxygen = np.array([38, 25, 20, 15, 10, 30], dtype=np.float64)
    ecm = np.array(ecm_values, dtype=np.float64)
    drug = np.array(drug_values, dtype=np.float64)
    return np.vstack([x, y, z, voxel_volume, oxygen, ecm, drug])


def _write_snapshot(
    output_dir: Path,
    *,
    step: int,
    time_min: float,
    cell_matrix: np.ndarray,
    micro_matrix: np.ndarray,
) -> None:
    xml_name = f"output{step:08d}.xml"
    cell_mat_name = f"cells_{step:08d}.mat"
    micro_mat_name = f"micro_{step:08d}.mat"

    _write_mat_v4(output_dir / cell_mat_name, {"cells": cell_matrix})
    _write_mat_v4(output_dir / micro_mat_name, {"multiscale_microenvironment": micro_matrix})

    root = ET.Element("MultiCellDS")
    metadata = ET.SubElement(root, "metadata")
    ET.SubElement(metadata, "current_time").text = f"{time_min}"

    cellular_information = ET.SubElement(root, "cellular_information")
    cell_types = ET.SubElement(cellular_information, "cell_types")
    ET.SubElement(cell_types, "type", {"type": "0"}).text = "tumor_cell"
    ET.SubElement(cell_types, "type", {"type": "1"}).text = "stromal_cell"

    simplified = ET.SubElement(cellular_information, "simplified_data")
    labels = ET.SubElement(simplified, "labels")

    label_specs = [
        (0, 1, "cell_type"),
        (1, 1, "dead"),
        (2, 1, "current_death_model"),
        (3, 1, "is_activated"),
        (4, 1, "is_mesenchymal"),
        (5, 1, "drug_sensitivity"),
        (6, 1, "HIF1A"),
        (7, 3, "position"),
    ]
    for idx, size, name in label_specs:
        ET.SubElement(labels, "label", {"index": str(idx), "size": str(size)}).text = name
    ET.SubElement(simplified, "filename").text = cell_mat_name

    micro = ET.SubElement(root, "microenvironment")
    domain = ET.SubElement(micro, "domain")
    variables = ET.SubElement(domain, "variables")
    ET.SubElement(variables, "variable", {"ID": "0", "name": "oxygen"})
    ET.SubElement(variables, "variable", {"ID": "1", "name": "ecm_density"})
    ET.SubElement(variables, "variable", {"ID": "2", "name": "drug"})
    data = ET.SubElement(domain, "data")
    ET.SubElement(data, "filename").text = micro_mat_name

    tree = ET.ElementTree(root)
    tree.write(output_dir / xml_name, encoding="utf-8", xml_declaration=True)


@pytest.fixture
def metrics_factory() -> Callable[..., SimulationMetrics]:
    def _factory(**overrides) -> SimulationMetrics:
        base = {
            "total_tumor_cells": 50,
            "total_stromal_cells": 200,
            "live_tumor_cells": 25,
            "activated_cafs": 100,
            "mesenchymal_tumor_cells": 10,
            "mean_ecm_density": 0.5,
            "max_ecm_density": 0.8,
            "mean_tumor_drug_sensitivity": 0.7,
            "tumor_extent": 400.0,
            "stroma_barrier_score": 0.5,
            "drug_penetration": 0.5,
            "hypoxic_fraction": 0.2,
            "time": 0.0,
            "source_file": "",
        }
        base.update(overrides)
        return SimulationMetrics(**base)

    return _factory


@pytest.fixture
def mock_physicell_output(tmp_path: Path) -> Path:
    output_dir = tmp_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    step0_cells = _build_cell_matrix(
        dead_second_tumor=False,
        stromal_activation=(0.8, 0.3),
        tumor_hif1a=(0.7, 0.2),
        tumor_drug_sensitivity=(0.8, 0.6),
    )
    step1_cells = _build_cell_matrix(
        dead_second_tumor=True,
        stromal_activation=(0.9, 0.8),
        tumor_hif1a=(0.8, 0.6),
        tumor_drug_sensitivity=(0.7, 0.2),
    )

    step0_micro = _build_micro_matrix(
        ecm_values=[0.1, 0.4, 0.5, 0.6, 0.2, 0.3],
        drug_values=[0.0, 0.2, 0.2, 0.1, 0.1, 0.3],
    )
    step1_micro = _build_micro_matrix(
        ecm_values=[0.2, 0.6, 0.7, 0.8, 0.3, 0.5],
        drug_values=[0.1, 0.5, 0.4, 0.3, 0.2, 0.6],
    )

    _write_snapshot(output_dir, step=0, time_min=0.0, cell_matrix=step0_cells, micro_matrix=step0_micro)
    _write_snapshot(output_dir, step=1, time_min=60.0, cell_matrix=step1_cells, micro_matrix=step1_micro)

    return output_dir


@pytest.fixture
def ea_test_paths(tmp_path: Path) -> dict:
    fake_binary = tmp_path / "fake_physicell.sh"
    fake_binary.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    fake_binary.chmod(fake_binary.stat().st_mode | 0o111)

    config_xml = tmp_path / "PhysiCell_settings.xml"
    config_xml.write_text(
        "<PhysiCell_settings><save><folder>output</folder></save></PhysiCell_settings>\n",
        encoding="utf-8",
    )

    return {
        "binary": fake_binary,
        "config": config_xml,
        "sim_output_dir": tmp_path / "sim_runs",
        "stats_csv": tmp_path / "ea_stats.csv",
    }
 