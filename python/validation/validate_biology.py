from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import shutil
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
import xml.etree.ElementTree as ET

# Ensure project-root imports work when running as a script.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

try:
    from ..ea.fitness import compute_fitness
    from ..wrapper.output_parser import OutputParser, SimulationMetrics
    from ..wrapper.physicell_runner import PhysiCellRunner
except Exception:  # pragma: no cover
    from python.ea.fitness import compute_fitness  # type: ignore
    from python.wrapper.output_parser import OutputParser, SimulationMetrics  # type: ignore
    from python.wrapper.physicell_runner import PhysiCellRunner  # type: ignore


LOGGER = logging.getLogger(__name__)


@dataclass
class ScenarioResult:
    anchor_id: int
    scenario: str
    success: bool
    criterion: str
    details: str
    dependencies: List[str]
    run_dir: str
    exit_code: int
    wall_time: float
    fitness: float
    metrics: Dict[str, Any]
    extras: Dict[str, Any]


@dataclass
class ValidationSummary:
    total: int
    passed: int
    failed: int
    scenario_results: List[ScenarioResult]


@dataclass
class ScenarioSpec:
    anchor_id: int
    name: str
    description: str
    criterion: str
    intervention_payload: Dict[str, Any]
    user_parameter_overrides: Dict[str, Any]
    variable_overrides: Dict[str, Dict[str, Any]]
    dependencies: List[str]
    evaluator: Callable[[SimulationMetrics, Dict[str, Any], Dict[str, Any]], Tuple[bool, str]]


class BiologyValidator:
    def __init__(
        self,
        binary_path: Path | str,
        config_path: Path | str,
        output_dir: Path | str,
        timeout_seconds: int = 7200,
        sim_max_time: Optional[float] = None,
    ):
        self.binary_path = Path(binary_path).expanduser().resolve()
        self.base_config_path = Path(config_path).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.timeout_seconds = int(timeout_seconds)
        self.sim_max_time: Optional[float] = sim_max_time

        if not self.binary_path.exists():
            raise FileNotFoundError(f"PhysiCell binary not found: {self.binary_path}")
        if not self.base_config_path.exists():
            raise FileNotFoundError(f"PhysiCell config not found: {self.base_config_path}")

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.initial_tumor_cells = self._read_user_param_int(self.base_config_path, "number_of_tumor_cells", default=50)
        self._baseline_cache: Dict[str, Any] = {}

    def run_all(self) -> ValidationSummary:
        specs = self._scenario_specs()
        scenario_results: List[ScenarioResult] = []
        scenario_by_name: Dict[str, ScenarioResult] = {}
        for spec in specs:
            failed_dependencies = [
                dep
                for dep in spec.dependencies
                if dep not in scenario_by_name or not scenario_by_name[dep].success
            ]
            if failed_dependencies:
                dep_msg = ", ".join(failed_dependencies)
                LOGGER.warning(
                    "Skipping validation scenario %s (anchor %d): dependency failure: %s",
                    spec.name,
                    spec.anchor_id,
                    dep_msg,
                )
                skipped = ScenarioResult(
                    anchor_id=spec.anchor_id,
                    scenario=spec.name,
                    success=False,
                    criterion=spec.criterion,
                    details=f"skipped due to failed dependencies: {dep_msg}",
                    dependencies=list(spec.dependencies),
                    run_dir="",
                    exit_code=-2,
                    wall_time=0.0,
                    fitness=0.0,
                    metrics={},
                    extras={"skipped": True, "failed_dependencies": failed_dependencies},
                )
                scenario_results.append(skipped)
                scenario_by_name[spec.name] = skipped
                self._baseline_cache[spec.name] = {
                    "fitness": skipped.fitness,
                    "metrics": skipped.metrics,
                    "extras": skipped.extras,
                }
                continue

            LOGGER.info("Running validation scenario: %s", spec.name)
            result = self._run_scenario(spec)
            scenario_results.append(result)
            scenario_by_name[spec.name] = result
            self._baseline_cache[spec.name] = {
                "fitness": result.fitness,
                "metrics": result.metrics,
                "extras": result.extras,
            }

        passed = sum(1 for r in scenario_results if r.success)
        failed = len(scenario_results) - passed
        return ValidationSummary(
            total=len(scenario_results),
            passed=passed,
            failed=failed,
            scenario_results=scenario_results,
        )

    def save_summary(self, summary: ValidationSummary, output_prefix: Path | str) -> Tuple[Path, Path]:
        output_prefix = Path(output_prefix).expanduser().resolve()
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        json_path = output_prefix.with_suffix(".json")
        csv_path = output_prefix.with_suffix(".csv")

        payload = {
            "total": summary.total,
            "passed": summary.passed,
            "failed": summary.failed,
            "scenario_results": [asdict(r) for r in summary.scenario_results],
        }
        json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        fieldnames = [
            "anchor_id",
            "scenario",
            "success",
            "criterion",
            "details",
            "dependencies",
            "run_dir",
            "exit_code",
            "wall_time",
            "fitness",
        ]
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for result in summary.scenario_results:
                writer.writerow(
                    {
                        "anchor_id": result.anchor_id,
                        "scenario": result.scenario,
                        "success": result.success,
                        "criterion": result.criterion,
                        "details": result.details,
                        "dependencies": ";".join(result.dependencies),
                        "run_dir": result.run_dir,
                        "exit_code": result.exit_code,
                        "wall_time": f"{result.wall_time:.3f}",
                        "fitness": f"{result.fitness:.6f}",
                    }
                )
        return json_path, csv_path

    def _run_scenario(self, spec: ScenarioSpec) -> ScenarioResult:
        scenario_dir = self.output_dir / f"{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}_{spec.name.lower()}"
        scenario_dir.mkdir(parents=True, exist_ok=True)

        config_copy = scenario_dir / "scenario_config.xml"
        shutil.copy2(self.base_config_path, config_copy)
        self._patch_config(
            config_copy,
            user_parameter_overrides=spec.user_parameter_overrides,
            variable_overrides=spec.variable_overrides,
        )

        runner = PhysiCellRunner(
            binary_path=self.binary_path,
            config_path=config_copy,
            output_dir=scenario_dir / "sim_runs",
            timeout_seconds=self.timeout_seconds,
        )

        sim_result = runner.run(spec.intervention_payload)
        if not sim_result.success:
            return ScenarioResult(
                anchor_id=spec.anchor_id,
                scenario=spec.name,
                success=False,
                criterion=spec.criterion,
                details=f"simulation failed (exit_code={sim_result.exit_code})",
                dependencies=list(spec.dependencies),
                run_dir=str(sim_result.run_dir),
                exit_code=sim_result.exit_code,
                wall_time=sim_result.wall_time,
                fitness=0.0,
                metrics={},
                extras={},
            )

        parser = OutputParser(sim_result.run_dir / "output")
        metrics = parser.parse_final_state()
        fitness = float(compute_fitness(metrics))
        extras = self._compute_extras(parser, metrics)

        ok, details = spec.evaluator(metrics, extras, self._baseline_cache)
        return ScenarioResult(
            anchor_id=spec.anchor_id,
            scenario=spec.name,
            success=ok,
            criterion=spec.criterion,
            details=details,
            dependencies=list(spec.dependencies),
            run_dir=str(sim_result.run_dir),
            exit_code=sim_result.exit_code,
            wall_time=sim_result.wall_time,
            fitness=fitness,
            metrics=asdict(metrics),
            extras=extras,
        )

    def _compute_extras(self, parser: OutputParser, metrics: SimulationMetrics) -> Dict[str, Any]:
        extras: Dict[str, Any] = {}

        try:
            df = parser.parse_timeseries()
            extras["timeseries_rows"] = int(len(df))
            if len(df) >= 2:
                extras["ecm_delta"] = float(df["mean_ecm"].iloc[-1] - df["mean_ecm"].iloc[0])
                extras["activated_cafs_delta"] = float(df["activated_cafs"].iloc[-1] - df["activated_cafs"].iloc[0])
                extras["tumor_count_delta"] = float(df["tumor_count"].iloc[-1] - df["tumor_count"].iloc[0])
                extras["tumor_count_min"] = float(df["tumor_count"].min())
                extras["tumor_count_end"] = float(df["tumor_count"].iloc[-1])
                extras["tumor_regrowth_after_nadir"] = bool(
                    extras["tumor_count_end"] > extras["tumor_count_min"] * 1.05
                )
                extras.update(self._timeseries_snapshots(df))
            else:
                extras["ecm_delta"] = math.nan
                extras["activated_cafs_delta"] = math.nan
                extras["tumor_count_delta"] = math.nan
                extras["tumor_count_min"] = math.nan
                extras["tumor_count_end"] = math.nan
                extras["tumor_regrowth_after_nadir"] = False
        except Exception as exc:  # pragma: no cover
            LOGGER.warning("Timeseries parsing failed: %s", exc)
            extras["timeseries_rows"] = 0

        # Gene-level summary for resistance scenario (NRF2 / ABCB1 in surviving tumor cells).
        extras["mean_nrf2_surviving_tumor"] = self._mean_gene_in_surviving_tumor(parser, "NRF2")
        extras["mean_abcb1_surviving_tumor"] = self._mean_gene_in_surviving_tumor(parser, "ABCB1")
        extras["activated_caf_fraction"] = (
            float(metrics.activated_cafs) / float(metrics.total_stromal_cells)
            if metrics.total_stromal_cells > 0
            else math.nan
        )
        extras.update(self._compute_spatial_emt_metrics(parser))
        extras.update(self._compute_sanctuary_metrics(parser))
        return extras

    @staticmethod
    def _timeseries_snapshots(df) -> Dict[str, float]:
        out: Dict[str, float] = {}
        if len(df) == 0:
            return out

        idx_map = {
            "early": 0,
            "mid": len(df) // 2,
            "late": len(df) - 1,
        }
        columns = [
            "mean_ecm",
            "drug_penetration",
            "hypoxic_fraction",
            "tumor_count",
            "tumor_extent",
            "mesenchymal_tumor_cells",
            "activated_cafs",
        ]
        for col in columns:
            if col not in df.columns:
                continue
            for stage, idx in idx_map.items():
                out[f"{col}_{stage}"] = float(df[col].iloc[idx])
        return out

    def _compute_spatial_emt_metrics(self, parser: OutputParser) -> Dict[str, Any]:
        result: Dict[str, Any] = {
            "core_mesenchymal_fraction": math.nan,
            "peripheral_mesenchymal_fraction": math.nan,
            "emt_peripheral_minus_core": math.nan,
        }
        try:
            final_xml = parser._find_final_snapshot_xml()  # noqa: SLF001
            snapshot = parser._read_physicell_xml(final_xml)  # noqa: SLF001
            matrix = snapshot["cell_matrix"]
            label_name_map = snapshot["label_name_map"]
        except Exception:  # pragma: no cover
            return result

        cell_type = self._row_by_label(matrix, label_name_map, "cell_type")
        dead = self._row_by_label(matrix, label_name_map, "dead")
        death_model = self._row_by_label(matrix, label_name_map, "current_death_model")
        is_mesenchymal = self._row_by_label(matrix, label_name_map, "is_mesenchymal")
        positions = self._positions_by_label(matrix, label_name_map, "position")

        if cell_type is None or is_mesenchymal is None or positions.size == 0:
            return result

        import numpy as np

        tumor_mask = np_rint(cell_type) == 0
        if dead is not None:
            tumor_mask &= dead <= 0.5
        if death_model is not None:
            tumor_mask &= np_rint(death_model) != 100

        if int(np.sum(tumor_mask)) < 10:
            return result

        tumor_positions = positions[tumor_mask, :]
        tumor_mesenchymal = is_mesenchymal[tumor_mask] > 0.5
        centroid = np.nanmean(tumor_positions, axis=0)
        radial = np.linalg.norm(tumor_positions - centroid, axis=1)
        core_cutoff = float(np.quantile(radial, 0.40))
        peripheral_cutoff = float(np.quantile(radial, 0.60))
        core_mask = radial <= core_cutoff
        peripheral_mask = radial >= peripheral_cutoff

        if not core_mask.any() or not peripheral_mask.any():
            return result

        core_frac = float(np.mean(tumor_mesenchymal[core_mask]))
        peripheral_frac = float(np.mean(tumor_mesenchymal[peripheral_mask]))
        result["core_mesenchymal_fraction"] = core_frac
        result["peripheral_mesenchymal_fraction"] = peripheral_frac
        result["emt_peripheral_minus_core"] = float(peripheral_frac - core_frac)
        return result

    def _compute_sanctuary_metrics(self, parser: OutputParser) -> Dict[str, Any]:
        result: Dict[str, Any] = {
            "live_mean_local_ecm": math.nan,
            "all_mean_local_ecm": math.nan,
            "live_mean_local_drug": math.nan,
            "all_mean_local_drug": math.nan,
            "live_fraction_in_top_ecm_quartile": math.nan,
            "live_tumor_count_snapshot": 0,
        }
        try:
            final_xml = parser._find_final_snapshot_xml()  # noqa: SLF001
            snapshot = parser._read_physicell_xml(final_xml)  # noqa: SLF001
            matrix = snapshot["cell_matrix"]
            label_name_map = snapshot["label_name_map"]
            micro_coords = snapshot["micro_coords"]
            micro_values = snapshot["micro_values"]
        except Exception:  # pragma: no cover
            return result

        cell_type = self._row_by_label(matrix, label_name_map, "cell_type")
        dead = self._row_by_label(matrix, label_name_map, "dead")
        death_model = self._row_by_label(matrix, label_name_map, "current_death_model")
        positions = self._positions_by_label(matrix, label_name_map, "position")
        ecm_vals = micro_values.get("ecm_density")
        drug_vals = micro_values.get("drug")

        if cell_type is None or positions.size == 0:
            return result
        if ecm_vals is None or drug_vals is None or micro_coords.size == 0:
            return result

        import numpy as np

        tumor_mask = np_rint(cell_type) == 0
        live_tumor_mask = tumor_mask.copy()
        if dead is not None:
            live_tumor_mask &= dead <= 0.5
        if death_model is not None:
            live_tumor_mask &= np_rint(death_model) != 100

        if not tumor_mask.any():
            return result

        tumor_positions = positions[tumor_mask, :]
        live_positions = positions[live_tumor_mask, :]
        result["live_tumor_count_snapshot"] = int(np.sum(live_tumor_mask))

        tumor_ecm = self._sample_nearest_voxel_values(tumor_positions, micro_coords, ecm_vals)
        tumor_drug = self._sample_nearest_voxel_values(tumor_positions, micro_coords, drug_vals)
        result["all_mean_local_ecm"] = float(np.nanmean(tumor_ecm))
        result["all_mean_local_drug"] = float(np.nanmean(tumor_drug))

        if live_positions.size == 0:
            return result

        live_ecm = self._sample_nearest_voxel_values(live_positions, micro_coords, ecm_vals)
        live_drug = self._sample_nearest_voxel_values(live_positions, micro_coords, drug_vals)
        result["live_mean_local_ecm"] = float(np.nanmean(live_ecm))
        result["live_mean_local_drug"] = float(np.nanmean(live_drug))
        top_ecm_cutoff = float(np.quantile(tumor_ecm, 0.75))
        result["live_fraction_in_top_ecm_quartile"] = float(np.mean(live_ecm >= top_ecm_cutoff))
        return result

    @staticmethod
    def _positions_by_label(matrix, label_name_map: Dict[str, Dict[str, Any]], label: str):
        import numpy as np

        entry = label_name_map.get(label)
        if entry is None:
            return np.empty((0, 3), dtype=float)
        idx = int(entry.get("index", -1))
        size = int(entry.get("size", 0))
        if idx < 0 or size < 3:
            return np.empty((0, 3), dtype=float)
        if idx + 3 > matrix.shape[0]:
            return np.empty((0, 3), dtype=float)
        return matrix[idx : idx + 3, :].T

    @staticmethod
    def _sample_nearest_voxel_values(points, voxels, values):
        import numpy as np

        if points.size == 0:
            return np.asarray([], dtype=float)
        d2 = np.sum((points[:, None, :] - voxels[None, :, :]) ** 2, axis=2)
        nearest_idx = np.argmin(d2, axis=1)
        return np.asarray(values[nearest_idx], dtype=float)

    def _mean_gene_in_surviving_tumor(self, parser: OutputParser, gene_name: str) -> float:
        try:
            final_xml = parser._find_final_snapshot_xml()  # noqa: SLF001
            snapshot = parser._read_physicell_xml(final_xml)  # noqa: SLF001
            matrix = snapshot["cell_matrix"]
            label_name_map = snapshot["label_name_map"]
        except Exception:  # pragma: no cover
            return math.nan

        gene_row = self._row_by_label(matrix, label_name_map, gene_name)
        cell_type = self._row_by_label(matrix, label_name_map, "cell_type")
        dead = self._row_by_label(matrix, label_name_map, "dead")
        death_model = self._row_by_label(matrix, label_name_map, "current_death_model")

        if gene_row is None or cell_type is None:
            return math.nan

        tumor_mask = np_rint(cell_type) == 0
        if dead is not None:
            tumor_mask &= dead <= 0.5
        if death_model is not None:
            tumor_mask &= np_rint(death_model) != 100

        if not tumor_mask.any():
            return math.nan
        return float(gene_row[tumor_mask].mean())

    @staticmethod
    def _row_by_label(matrix, label_name_map: Dict[str, Dict[str, Any]], label: str):
        entry = label_name_map.get(label)
        if entry is None:
            return None
        idx = int(entry.get("index", -1))
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    def _patch_config(
        self,
        config_path: Path,
        user_parameter_overrides: Dict[str, Any],
        variable_overrides: Dict[str, Dict[str, Any]],
    ) -> None:
        tree = ET.parse(config_path)
        root = tree.getroot()

        if self.sim_max_time is not None:
            for tag in ("max_time", "overall_time"):
                node = root.find(f".//{tag}")
                if node is not None:
                    node.text = str(self.sim_max_time)

        for key, value in user_parameter_overrides.items():
            node = root.find(f".//user_parameters/{key}")
            if node is None:
                continue
            node.text = str(value)

        for variable_name, overrides in variable_overrides.items():
            var_node = root.find(f".//microenvironment_setup/variable[@name='{variable_name}']")
            if var_node is None:
                continue

            if "Dirichlet_boundary_condition" in overrides:
                dbc = var_node.find("./Dirichlet_boundary_condition")
                if dbc is not None:
                    cfg = overrides["Dirichlet_boundary_condition"]
                    if "value" in cfg:
                        dbc.text = str(cfg["value"])
                    if "enabled" in cfg:
                        dbc.set("enabled", "true" if bool(cfg["enabled"]) else "false")

            if "Dirichlet_options" in overrides:
                bcfg = overrides["Dirichlet_options"]
                for bnode in var_node.findall("./Dirichlet_options/boundary_value"):
                    if "value" in bcfg:
                        bnode.text = str(bcfg["value"])
                    if "enabled" in bcfg:
                        bnode.set("enabled", "true" if bool(bcfg["enabled"]) else "false")

            if "initial_condition" in overrides:
                inode = var_node.find("./initial_condition")
                if inode is not None:
                    inode.text = str(overrides["initial_condition"])

        tree.write(config_path, encoding="utf-8", xml_declaration=True)

    @staticmethod
    def _read_user_param_int(config_path: Path, name: str, default: int) -> int:
        try:
            tree = ET.parse(config_path)
            root = tree.getroot()
            node = root.find(f".//user_parameters/{name}")
            if node is not None and node.text is not None:
                return int(float(node.text.strip()))
        except Exception:  # pragma: no cover
            pass
        return int(default)

    def _scenario_specs(self) -> List[ScenarioSpec]:
        no_cytotoxic_user = {"drug_start_time": 1.0e9, "drug_concentration": 0.0}
        no_cytotoxic_var = {
            "drug": {
                "Dirichlet_boundary_condition": {"enabled": False, "value": 0.0},
                "Dirichlet_options": {"enabled": False, "value": 0.0},
                "initial_condition": 0.0,
            }
        }
        cytotoxic_user = {"drug_start_time": 0.0, "drug_concentration": 1.0}
        cytotoxic_var = {
            "drug": {
                "Dirichlet_boundary_condition": {"enabled": True, "value": 1.0},
                "Dirichlet_options": {"enabled": True, "value": 1.0},
                "initial_condition": 1.0,
            }
        }
        high_tgfb_var = {
            "tgfb": {
                "Dirichlet_boundary_condition": {"enabled": True, "value": 1.0},
                "Dirichlet_options": {"enabled": True, "value": 1.0},
                "initial_condition": 1.0,
            }
        }
        high_stroma_user = {"number_of_stromal_cells": 400}
        low_compaction_user = {"mechanical_compaction_strength": 0.1}
        shh_off = [
            {"knob": "shh_secretion_rate", "effect": "INHIBIT", "strength": 1.0, "name": "SHH_OFF"}
        ]

        specs: List[ScenarioSpec] = [
            ScenarioSpec(
                anchor_id=1,
                name="ANCHOR_1_SELF_ASSEMBLY",
                description="Barrier self-assembles from tumor->stroma paracrine signaling.",
                criterion="CAF activation rises and ECM accumulates without external instruction",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=[],
                evaluator=self._eval_anchor1_self_assembly,
            ),
            ScenarioSpec(
                anchor_id=3,
                name="ANCHOR_3_SHH_PARADOX",
                description="SHH inhibition alone should worsen tumor outcome vs baseline.",
                criterion="tumor burden and extent increase vs Anchor 1 without cytotoxic",
                intervention_payload={
                    "calibration_profile": "AsPC-1",
                    "knob_interventions": shh_off,
                },
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=["ANCHOR_1_SELF_ASSEMBLY"],
                evaluator=self._eval_anchor3_shh_paradox,
            ),
            ScenarioSpec(
                anchor_id=2,
                name="ANCHOR_2_DRUG_PENETRATION_MATURITY",
                description="Drug penetration declines as ECM barrier matures over time.",
                criterion="ECM rises early->mid->late while tumor drug penetration drops monotonically",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=cytotoxic_user,
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_1_SELF_ASSEMBLY"],
                evaluator=self._eval_anchor2_penetration_maturity,
            ),
            ScenarioSpec(
                anchor_id=4,
                name="ANCHOR_4_SMAD4_ASYMMETRY",
                description="SMAD4-intact profile should show stronger TGF-beta growth braking.",
                criterion="PANC-1 (high knob2) grows slower than AsPC-1 baseline under high TGF-beta",
                intervention_payload={
                    "calibration_profile": "PANC-1",
                    "knob_interventions": [],
                },
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides={**no_cytotoxic_var, **high_tgfb_var},
                dependencies=["ANCHOR_1_SELF_ASSEMBLY"],
                evaluator=self._eval_anchor4_smad4_asymmetry,
            ),
            ScenarioSpec(
                anchor_id=5,
                name="ANCHOR_5_CENTRAL_HYPOXIA",
                description="Hypoxia should emerge in the interior as barrier density increases.",
                criterion="hypoxic fraction rises with ECM maturation over time",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=["ANCHOR_1_SELF_ASSEMBLY"],
                evaluator=self._eval_anchor5_central_hypoxia,
            ),
            ScenarioSpec(
                anchor_id=6,
                name="ANCHOR_6_PERIPHERAL_EMT",
                description="EMT should be spatially enriched at tumor periphery.",
                criterion="peripheral mesenchymal fraction is higher than core fraction",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=["ANCHOR_5_CENTRAL_HYPOXIA"],
                evaluator=self._eval_anchor6_peripheral_emt,
            ),
            ScenarioSpec(
                anchor_id=7,
                name="ANCHOR_7_ECM_DEGRADATION_LIMITS",
                description="Barrier opening without cytotoxic is insufficient; with cytotoxic it improves kill.",
                criterion="SHH-off + cytotoxic improves kill vs SHH-off alone while residual tumor persists",
                intervention_payload={
                    "calibration_profile": "AsPC-1",
                    "knob_interventions": shh_off,
                },
                user_parameter_overrides=cytotoxic_user,
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_3_SHH_PARADOX", "ANCHOR_2_DRUG_PENETRATION_MATURITY"],
                evaluator=self._eval_anchor7_ecm_degradation_limits,
            ),
            ScenarioSpec(
                anchor_id=8,
                name="ANCHOR_8_TWO_COMPONENT_BARRIER",
                description="Mechanical and diffusion barrier effects should be distinguishable.",
                criterion="collagen-like mechanical depletion proxy increases spread with limited penetration change",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides={**cytotoxic_user, **low_compaction_user},
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_7_ECM_DEGRADATION_LIMITS"],
                evaluator=self._eval_anchor8_two_component_barrier,
            ),
            ScenarioSpec(
                anchor_id=9,
                name="ANCHOR_9_BARRIER_DENSITY_PROGNOSTIC",
                description="Higher PSC density should worsen response for same drug regimen.",
                criterion="higher stromal seeding yields lower penetration and higher residual tumor",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides={**cytotoxic_user, **high_stroma_user},
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_2_DRUG_PENETRATION_MATURITY"],
                evaluator=self._eval_anchor9_density_prognostic,
            ),
            ScenarioSpec(
                anchor_id=10,
                name="ANCHOR_10_SPATIAL_SANCTUARY",
                description="Residual survivors should cluster in dense ECM sanctuary under strong drug pressure.",
                criterion="survivors localize to high-ECM/low-drug niches and show regrowth tendency",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides={**cytotoxic_user, **high_stroma_user},
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_9_BARRIER_DENSITY_PROGNOSTIC"],
                evaluator=self._eval_anchor10_spatial_sanctuary,
            ),
        ]
        return specs

    @staticmethod
    def _eval_anchor1_self_assembly(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        activated_caf_fraction = float(extras.get("activated_caf_fraction", math.nan))
        ecm_delta = float(extras.get("ecm_delta", math.nan))
        ok = (
            metrics.activated_cafs > 0
            and activated_caf_fraction >= 0.05
            and metrics.mean_ecm_density > 0.02
            and math.isfinite(ecm_delta)
            and ecm_delta > 0.0
        )
        return ok, (
            f"activated_caf_fraction={activated_caf_fraction:.3f}, "
            f"mean_ecm={metrics.mean_ecm_density:.3f}, "
            f"ecm_delta={ecm_delta:.3f}"
        )

    @staticmethod
    def _eval_anchor3_shh_paradox(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ref = baseline.get("ANCHOR_1_SELF_ASSEMBLY", {}).get("metrics", {})
        if not ref:
            return False, "missing ANCHOR_1_SELF_ASSEMBLY reference"

        ref_tumor = float(ref.get("total_tumor_cells", math.nan))
        ref_extent = float(ref.get("tumor_extent", math.nan))
        ref_ecm = float(ref.get("mean_ecm_density", math.nan))

        worse_burden = metrics.total_tumor_cells > ref_tumor * 1.05
        larger_extent = metrics.tumor_extent > ref_extent * 1.02
        thinner_barrier = metrics.mean_ecm_density < ref_ecm
        ok = worse_burden and larger_extent and thinner_barrier
        details = (
            f"tumor={metrics.total_tumor_cells} vs a1={ref_tumor:.1f}, "
            f"extent={metrics.tumor_extent:.1f} vs a1={ref_extent:.1f}, "
            f"mean_ecm={metrics.mean_ecm_density:.3f} vs a1={ref_ecm:.3f}"
        )
        return ok, details

    @staticmethod
    def _eval_anchor2_penetration_maturity(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ecm_early = float(extras.get("mean_ecm_early", math.nan))
        ecm_mid = float(extras.get("mean_ecm_mid", math.nan))
        ecm_late = float(extras.get("mean_ecm_late", math.nan))
        pen_early = float(extras.get("drug_penetration_early", math.nan))
        pen_mid = float(extras.get("drug_penetration_mid", math.nan))
        pen_late = float(extras.get("drug_penetration_late", math.nan))

        ecm_matures = ecm_early < ecm_mid < ecm_late
        penetration_drops = pen_early > pen_mid > pen_late
        drop_magnitude = pen_late <= pen_early * 0.9 if math.isfinite(pen_early) else False
        ok = ecm_matures and penetration_drops and drop_magnitude
        details = (
            f"ecm: {ecm_early:.3f}->{ecm_mid:.3f}->{ecm_late:.3f}; "
            f"penetration: {pen_early:.3f}->{pen_mid:.3f}->{pen_late:.3f}"
        )
        return ok, details

    @staticmethod
    def _eval_anchor4_smad4_asymmetry(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        aspc1 = baseline.get("ANCHOR_1_SELF_ASSEMBLY", {}).get("metrics", {})
        if not aspc1:
            return False, "missing ANCHOR_1_SELF_ASSEMBLY reference"

        aspc1_tumor = float(aspc1.get("total_tumor_cells", math.nan))
        aspc1_ecm = float(aspc1.get("mean_ecm_density", math.nan))
        slower_growth = metrics.total_tumor_cells < aspc1_tumor * 0.95
        thinner_or_similar_barrier = metrics.mean_ecm_density <= aspc1_ecm * 1.1
        ok = slower_growth and thinner_or_similar_barrier
        details = (
            f"panc1_tumor={metrics.total_tumor_cells} vs aspc1={aspc1_tumor:.1f}, "
            f"panc1_ecm={metrics.mean_ecm_density:.3f} vs aspc1={aspc1_ecm:.3f}"
        )
        return ok, details

    @staticmethod
    def _eval_anchor5_central_hypoxia(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        h_early = float(extras.get("hypoxic_fraction_early", math.nan))
        h_mid = float(extras.get("hypoxic_fraction_mid", math.nan))
        h_late = float(extras.get("hypoxic_fraction_late", math.nan))
        e_early = float(extras.get("mean_ecm_early", math.nan))
        e_late = float(extras.get("mean_ecm_late", math.nan))
        ok = (h_early <= h_mid <= h_late) and (h_late >= 0.1) and (e_late > e_early)
        details = (
            f"hypoxia={h_early:.3f}->{h_mid:.3f}->{h_late:.3f}, "
            f"ecm={e_early:.3f}->{e_late:.3f}"
        )
        return ok, details

    @staticmethod
    def _eval_anchor6_peripheral_emt(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        core = float(extras.get("core_mesenchymal_fraction", math.nan))
        peripheral = float(extras.get("peripheral_mesenchymal_fraction", math.nan))
        delta = float(extras.get("emt_peripheral_minus_core", math.nan))
        ok = math.isfinite(delta) and delta > 0.05 and peripheral > core
        details = f"core_mes={core:.3f}, peripheral_mes={peripheral:.3f}, delta={delta:.3f}"
        return ok, details

    @staticmethod
    def _eval_anchor7_ecm_degradation_limits(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        shh_no_cyto = baseline.get("ANCHOR_3_SHH_PARADOX", {}).get("metrics", {})
        drug_only = baseline.get("ANCHOR_2_DRUG_PENETRATION_MATURITY", {}).get("metrics", {})
        if not shh_no_cyto or not drug_only:
            return False, "missing ANCHOR_3_SHH_PARADOX or ANCHOR_2_DRUG_PENETRATION_MATURITY reference"

        shh_no_cyto_tumor = float(shh_no_cyto.get("total_tumor_cells", math.nan))
        drug_only_live = float(drug_only.get("live_tumor_cells", math.nan))
        drug_only_pen = float(drug_only.get("drug_penetration", math.nan))

        improved_vs_shh_alone = metrics.total_tumor_cells < shh_no_cyto_tumor
        improved_penetration = metrics.drug_penetration > drug_only_pen
        residual_persists = metrics.live_tumor_cells > 0
        better_than_drug_only = metrics.live_tumor_cells <= drug_only_live

        ok = improved_vs_shh_alone and improved_penetration and residual_persists and better_than_drug_only
        details = (
            f"tumor={metrics.total_tumor_cells} vs shh_alone={shh_no_cyto_tumor:.1f}, "
            f"live={metrics.live_tumor_cells} vs drug_only={drug_only_live:.1f}, "
            f"penetration={metrics.drug_penetration:.3f} vs drug_only={drug_only_pen:.3f}"
        )
        return ok, details

    @staticmethod
    def _eval_anchor8_two_component_barrier(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ref = baseline.get("ANCHOR_7_ECM_DEGRADATION_LIMITS", {}).get("metrics", {})
        if not ref:
            return False, "missing ANCHOR_7_ECM_DEGRADATION_LIMITS reference"

        ref_extent = float(ref.get("tumor_extent", math.nan))
        ref_pen = float(ref.get("drug_penetration", math.nan))
        extent_ratio = metrics.tumor_extent / max(ref_extent, 1e-9)
        penetration_ratio = metrics.drug_penetration / max(ref_pen, 1e-9)

        # Proxy interpretation:
        # low compaction (collagen-like depletion) should release spatial restraint
        # more strongly than it improves penetration.
        ok = extent_ratio > 1.05 and penetration_ratio < 1.15
        details = (
            f"extent_ratio={extent_ratio:.3f}, penetration_ratio={penetration_ratio:.3f} "
            f"(ref extent={ref_extent:.1f}, ref pen={ref_pen:.3f})"
        )
        return ok, details

    @staticmethod
    def _eval_anchor9_density_prognostic(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ref = baseline.get("ANCHOR_2_DRUG_PENETRATION_MATURITY", {}).get("metrics", {})
        if not ref:
            return False, "missing ANCHOR_2_DRUG_PENETRATION_MATURITY reference"

        ref_pen = float(ref.get("drug_penetration", math.nan))
        ref_live = float(ref.get("live_tumor_cells", math.nan))
        ref_stroma = float(ref.get("total_stromal_cells", math.nan))
        ok = (
            metrics.total_stromal_cells > ref_stroma
            and metrics.drug_penetration < ref_pen
            and metrics.live_tumor_cells > ref_live
        )
        details = (
            f"stroma={metrics.total_stromal_cells} vs ref={ref_stroma:.0f}, "
            f"penetration={metrics.drug_penetration:.3f} vs ref={ref_pen:.3f}, "
            f"live_tumor={metrics.live_tumor_cells} vs ref={ref_live:.1f}"
        )
        return ok, details

    @staticmethod
    def _eval_anchor10_spatial_sanctuary(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        live_local_ecm = float(extras.get("live_mean_local_ecm", math.nan))
        all_local_ecm = float(extras.get("all_mean_local_ecm", math.nan))
        live_local_drug = float(extras.get("live_mean_local_drug", math.nan))
        all_local_drug = float(extras.get("all_mean_local_drug", math.nan))
        top_quartile_frac = float(extras.get("live_fraction_in_top_ecm_quartile", math.nan))
        regrowth_proxy = bool(extras.get("tumor_regrowth_after_nadir", False))

        has_residual = metrics.live_tumor_cells > 0
        ecm_sanctuary = math.isfinite(live_local_ecm) and math.isfinite(all_local_ecm) and live_local_ecm > all_local_ecm
        drug_shadow = math.isfinite(live_local_drug) and math.isfinite(all_local_drug) and live_local_drug < all_local_drug
        concentrated_in_dense_stroma = math.isfinite(top_quartile_frac) and top_quartile_frac >= 0.5

        ok = has_residual and ecm_sanctuary and drug_shadow and concentrated_in_dense_stroma and regrowth_proxy
        details = (
            f"live_tumor={metrics.live_tumor_cells}, "
            f"live_ecm={live_local_ecm:.3f} vs all_ecm={all_local_ecm:.3f}, "
            f"live_drug={live_local_drug:.3f} vs all_drug={all_local_drug:.3f}, "
            f"top_ecm_frac={top_quartile_frac:.3f}, regrowth_proxy={regrowth_proxy}"
        )
        return ok, details


def np_rint(values):
    # Local helper to avoid importing numpy at module import-time in minimal environments.
    import numpy as np

    return np.rint(values).astype(int)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run biological validation scenarios for Stroma World.")
    parser.add_argument("--physicell-binary", required=True, help="Path to stroma_world executable.")
    parser.add_argument("--physicell-config", "--base-config", dest="physicell_config", required=True, help="Path to PhysiCell_settings.xml.")
    parser.add_argument("--output-dir", required=True, help="Directory for validation outputs.")
    parser.add_argument("--timeout-seconds", type=int, default=7200, help="Per-simulation timeout in seconds.")
    parser.add_argument(
        "--sim-max-time",
        type=float,
        default=None,
        help="Override max_time (minutes) in the PhysiCell config for each validation run. "
             "Use a short value (e.g. 1440) to limit output size.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level.",
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)

    validator = BiologyValidator(
        binary_path=Path(args.physicell_binary),
        config_path=Path(args.physicell_config),
        output_dir=Path(args.output_dir),
        timeout_seconds=args.timeout_seconds,
        sim_max_time=args.sim_max_time,
    )
    summary = validator.run_all()

    prefix = Path(args.output_dir).expanduser().resolve() / "biology_validation_summary"
    json_path, csv_path = validator.save_summary(summary, prefix)

    print(f"Biology validation complete: {summary.passed}/{summary.total} passed")
    print(f"JSON summary: {json_path}")
    print(f"CSV summary: {csv_path}")
    for scenario in summary.scenario_results:
        status = "PASS" if scenario.success else "FAIL"
        print(f"- [{status}] Anchor {scenario.anchor_id} {scenario.scenario}: {scenario.details}")
    return 0 if summary.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
