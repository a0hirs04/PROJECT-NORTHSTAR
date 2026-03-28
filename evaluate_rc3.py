#!/usr/bin/env python3
"""Reality Check 3 evaluator: vismodegib paradox."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))
from python.wrapper.output_parser import OutputParser

SAVE_INTERVAL = 360.0
RC3_ENDPOINT_TIME = 80640.0
RC2_TREAT_END_TIME = 40320.0
PERITUMORAL_SHELL_WIDTH = 100.0
ARM_ORDER = ["Arm_A_Control", "Arm_B_SHH_Only", "Arm_C_Combo"]
REQUIRED_SEEDS = [42, 43, 44, 45, 46]
ARM_SHORT = {
    "Arm_A_Control": "A",
    "Arm_B_SHH_Only": "B",
    "Arm_C_Combo": "C",
}
SEED_RE = re.compile(r"seed[_-]?(\d+)")


def _row(matrix: np.ndarray, labels: Dict[str, Dict[str, Any]], name: str):
    entry = labels.get(name)
    if entry is None:
        return None
    idx = int(entry["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _sample_field_at_positions(positions: np.ndarray, micro_coords: np.ndarray, field_vals: np.ndarray) -> np.ndarray:
    if positions.size == 0 or micro_coords.size == 0 or field_vals.size == 0:
        return np.array([], dtype=float)
    out = []
    vox2 = micro_coords[:, :2] if micro_coords.shape[1] > 2 else micro_coords
    for i in range(0, positions.shape[0], 500):
        p = positions[i : i + 500, :2]
        d2 = np.sum((p[:, None, :] - vox2[None, :, :]) ** 2, axis=2)
        out.append(field_vals[np.argmin(d2, axis=1)])
    return np.concatenate(out) if out else np.array([], dtype=float)


def _fmt(v: float, nd: int = 2) -> str:
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "nan"
    return f"{float(v):.{nd}f}"


def _mean_sd(values: Sequence[float]) -> tuple[float, float, int]:
    finite = np.array([float(v) for v in values if math.isfinite(float(v))], dtype=float)
    if finite.size == 0:
        return math.nan, math.nan, 0
    mean = float(np.mean(finite))
    sd = float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0
    return mean, sd, int(finite.size)


def _median(values: Sequence[float]) -> float:
    finite = np.array([float(v) for v in values if math.isfinite(float(v))], dtype=float)
    if finite.size == 0:
        return math.nan
    return float(np.median(finite))


def _infer_seed(path: Path) -> int | None:
    m = SEED_RE.search(path.as_posix())
    return int(m.group(1)) if m else None


def validate_cell_type_mapping(snapshot: Dict[str, Any]) -> tuple[int, int]:
    cell_type_names = snapshot.get("cell_type_names", {})
    if isinstance(cell_type_names, dict) and cell_type_names:
        tumor_ids = [cid for cid, name in cell_type_names.items() if name == "tumor_cell"]
        stroma_ids = [cid for cid, name in cell_type_names.items() if name == "stromal_cell"]
        if len(tumor_ids) == 1 and len(stroma_ids) == 1:
            return int(tumor_ids[0]), int(stroma_ids[0])
        raise ValueError(
            "Could not validate tumor/stromal cell type IDs from snapshot metadata: "
            f"{cell_type_names}"
        )
    raise ValueError(
        "Snapshot is missing cell_type_names metadata; refusing to silently assume tumor=0 and stroma=1."
    )


def _load_snapshot(parser: OutputParser, out_dir: Path, target_time_min: float):
    idx = int(round(target_time_min / SAVE_INTERVAL))
    for i in range(idx, -1, -1):
        f = out_dir / f"output{i:08d}.xml"
        if f.exists() and f.stat().st_size > 0:
            return parser._read_physicell_xml(f), f
    final_xml = out_dir / "final.xml"
    if final_xml.exists() and final_xml.stat().st_size > 0:
        return parser._read_physicell_xml(final_xml), final_xml
    raise FileNotFoundError(f"No readable snapshot at or before t={target_time_min} in {out_dir}")


def _parse_tumor_metrics(snapshot: Dict[str, Any], parser: OutputParser) -> Dict[str, Any]:
    matrix = snapshot["cell_matrix"]
    labels = snapshot["label_name_map"]
    micro_coords = snapshot["micro_coords"]
    micro_values = snapshot["micro_values"]

    cell_type = _row(matrix, labels, "cell_type")
    dead = _row(matrix, labels, "dead")
    death_model = _row(matrix, labels, "current_death_model")
    n_cells = matrix.shape[1]

    live_mask = np.ones(n_cells, dtype=bool)
    if dead is not None:
        live_mask &= dead <= 0.5
    if death_model is not None:
        live_mask &= np.rint(death_model).astype(int) != 100

    tumor_type_id, _ = validate_cell_type_mapping(snapshot)
    ctype = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1)
    tumor_mask = ctype == tumor_type_id
    live_tumor = live_mask & tumor_mask

    positions = parser._get_positions(matrix, labels)
    tumor_pos = positions[live_tumor] if positions.size and np.any(live_tumor) else np.empty((0, 3))

    if tumor_pos.shape[0] > 0:
        centroid = np.mean(tumor_pos, axis=0)
        radii = np.linalg.norm(tumor_pos - centroid, axis=1)
        radius_max = float(np.max(radii)) if radii.size else 0.0
        radius_p95 = float(np.percentile(radii, 95)) if radii.size else 0.0
    else:
        centroid = np.array([math.nan, math.nan, math.nan])
        radius_max = 0.0
        radius_p95 = 0.0

    peri_ecm = math.nan
    ecm_field = micro_values.get("ecm_density")
    if ecm_field is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        shell_mask = np.zeros(micro_coords.shape[0], dtype=bool)
        for i in range(0, micro_coords.shape[0], 5000):
            v = micro_coords[i : i + 5000]
            d2 = np.sum((v[:, None, :] - tumor_pos[None, :, :]) ** 2, axis=2)
            dmin = np.min(d2, axis=1)
            shell_mask[i : i + 5000] = (dmin > 0.0) & (dmin <= PERITUMORAL_SHELL_WIDTH**2)
        if np.any(shell_mask):
            peri_ecm = float(np.nanmean(ecm_field[shell_mask]))
    if ecm_field is not None and ecm_field.size > 0:
        global_ecm_max = float(np.nanmax(ecm_field))
        if global_ecm_max == 0.0 and math.isfinite(peri_ecm) and peri_ecm != 0.0:
            raise ValueError(
                "ECM consistency check failed: global ecm_density max is 0 but periECM is "
                f"{peri_ecm} in {snapshot['filepath']}"
            )

    return {
        "time_min": float(snapshot["time"]),
        "live_tumor": int(np.sum(live_tumor)),
        "centroid": centroid,
        "radius_max": radius_max,
        "radius_p95": radius_p95,
        "peri_ecm": peri_ecm,
    }


def _discover_output_dirs(base_path: Path) -> List[Path]:
    def is_output_dir(path: Path) -> bool:
        return path.is_dir() and ((path / "final.xml").exists() or any(path.glob("output*.xml")))

    if is_output_dir(base_path):
        return [base_path]
    nested_output = base_path / "output"
    if is_output_dir(nested_output):
        return [nested_output]

    discovered = []
    for output_dir in sorted(base_path.rglob("output")):
        if is_output_dir(output_dir):
            discovered.append(output_dir)
    return discovered


def _read_xml_float(config_path: Path, xpath: str) -> float:
    import xml.etree.ElementTree as ET

    root = ET.parse(config_path).getroot()
    node = root.find(xpath)
    if node is None or node.text is None:
        raise ValueError(f"Missing XML node {xpath} in {config_path}")
    return float(node.text)


def _validate_rc2_benchmark_run(rep_dir: Path) -> None:
    config_path = rep_dir / "config.xml"
    if not config_path.exists():
        raise FileNotFoundError(f"RC2 benchmark replicate is missing config.xml: {rep_dir}")

    drug_start = _read_xml_float(config_path, ".//user_parameters/drug_start_time")
    drug_end = _read_xml_float(config_path, ".//user_parameters/drug_end_time")
    drug_conc = _read_xml_float(config_path, ".//user_parameters/drug_concentration")
    shh_strength = _read_xml_float(config_path, ".//user_parameters/shh_inhibition_strength")
    drug_kill_multiplier = _read_xml_float(config_path, ".//user_parameters/drug_kill_multiplier")

    if not math.isclose(drug_start, 20160.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(f"RC2 benchmark run has unexpected drug_start_time={drug_start}: {rep_dir}")
    if not math.isclose(drug_end, 40320.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(f"RC2 benchmark run has unexpected drug_end_time={drug_end}: {rep_dir}")
    if drug_conc <= 0.0:
        raise ValueError(f"RC2 benchmark run is not drug-on: drug_concentration={drug_conc} in {rep_dir}")
    if not math.isclose(shh_strength, 0.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(
            f"RC2 benchmark run is not drug-alone: shh_inhibition_strength={shh_strength} in {rep_dir}"
        )
    if not math.isclose(drug_kill_multiplier, 0.01, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(
            "RC2 benchmark run is not the frozen km=0.01 gold comparator: "
            f"drug_kill_multiplier={drug_kill_multiplier} in {rep_dir}"
        )


def _clean_benchmark_replicate_dirs(path: Path) -> List[Path]:
    direct_dirs = sorted([p for p in path.iterdir() if p.is_dir()])
    replicate_dirs = [p for p in direct_dirs if re.fullmatch(r"replicate_\d+_seed\d+", p.name)]
    if not replicate_dirs:
        raise ValueError(
            "RC2 benchmark directory must contain replicate_* drug-alone runs only. "
            f"No replicate_* directories found under {path}"
        )
    unexpected_dirs = [p for p in direct_dirs if p not in replicate_dirs]
    if unexpected_dirs:
        names = ", ".join(p.name for p in unexpected_dirs)
        raise ValueError(
            "RC2 benchmark directory is not clean. Expected only replicate_* drug-alone runs, "
            f"found extra directories: {names}"
        )
    for rep_dir in replicate_dirs:
        _validate_rc2_benchmark_run(rep_dir)
    return replicate_dirs


def _load_rc2_benchmark_from_dir(path: Path) -> Dict[str, Any]:
    replicate_dirs = _clean_benchmark_replicate_dirs(path)

    runs = []
    for rep_dir in replicate_dirs:
        out_dir = rep_dir / "output"
        if not out_dir.is_dir():
            raise FileNotFoundError(f"RC2 benchmark replicate is missing output/: {rep_dir}")
        parser = OutputParser(out_dir)
        snapshot, _ = _load_snapshot(parser, Path(parser.output_dir), RC2_TREAT_END_TIME)
        metrics = _parse_tumor_metrics(snapshot, parser)
        seed = _infer_seed(rep_dir)
        if seed is None:
            raise ValueError(f"Could not infer seed from benchmark replicate path: {rep_dir}")
        runs.append(
            {
                "seed": seed,
                "day28_live_tumor": metrics["live_tumor"],
                "source": str(out_dir),
                "replicate": rep_dir.name,
            }
        )

    values = [float(run["day28_live_tumor"]) for run in runs]
    mean, sd, n = _mean_sd(values)
    by_seed = {int(run["seed"]): float(run["day28_live_tumor"]) for run in runs if run["seed"] is not None}
    return {
        "mode": "directory",
        "runs": runs,
        "by_seed": by_seed,
        "mean": mean,
        "sd": sd,
        "n": n,
        "source": str(path),
    }


def _load_rc2_benchmark_from_json(path: Path) -> Dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    runs = []
    by_seed = {}

    if isinstance(data, dict) and isinstance(data.get("runs"), list):
        for entry in data["runs"]:
            value = entry.get("day28_live_tumor", entry.get("endpoint_live_tumor"))
            if value is None:
                continue
            seed = entry.get("seed")
            value = float(value)
            runs.append({"seed": seed, "day28_live_tumor": value, "source": str(path)})
            if seed is not None:
                by_seed[int(seed)] = value
    elif isinstance(data, dict):
        mean = data.get("mean_day28_live_tumor", data.get("mean_endpoint_live_tumor"))
        if mean is not None:
            mean = float(mean)
            return {
                "mode": "json-summary",
                "runs": [],
                "by_seed": {},
                "mean": mean,
                "sd": float(data.get("sd_day28_live_tumor", data.get("sd_endpoint_live_tumor", math.nan))),
                "n": int(data.get("n", 1)),
                "source": str(path),
            }

    if not runs:
        raise ValueError(f"Unsupported RC2 benchmark JSON format: {path}")

    values = [float(run["day28_live_tumor"]) for run in runs]
    mean, sd, n = _mean_sd(values)
    return {
        "mode": "json-runs",
        "runs": runs,
        "by_seed": by_seed,
        "mean": mean,
        "sd": sd,
        "n": n,
        "source": str(path),
    }


def _load_rc2_benchmark_from_text(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    matches = [float(x) for x in re.findall(r"Day 28 live tumor:\s*(\d+)", text)]
    if not matches:
        matches = [float(x) for x in re.findall(r"d28\s*=\s*(\d+)", text)]
    if not matches:
        raise ValueError(f"Could not parse RC2 benchmark from text file: {path}")
    mean, sd, n = _mean_sd(matches)
    return {
        "mode": "text-summary",
        "runs": [],
        "by_seed": {},
        "mean": mean,
        "sd": sd,
        "n": n,
        "source": str(path),
    }


def load_rc2_benchmark(path_str: str) -> Dict[str, Any]:
    path = Path(path_str).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"RC2 benchmark path does not exist: {path}")
    if path.is_dir():
        return _load_rc2_benchmark_from_dir(path)
    if path.suffix.lower() == ".json":
        return _load_rc2_benchmark_from_json(path)
    return _load_rc2_benchmark_from_text(path)


def discover_rc3_runs(base_dir: Path | None, arm_dirs: Dict[str, Path | None]) -> List[Dict[str, Any]]:
    resolved_arm_dirs: Dict[str, Path] = {}
    for arm in ARM_ORDER:
        override = arm_dirs.get(arm)
        if override is not None:
            resolved_arm_dirs[arm] = override
            continue
        if base_dir is None:
            raise FileNotFoundError(
                f"Missing RC3 arm directory for {arm}. Provide --base-dir or explicit --arm-*-dir overrides."
            )
        resolved_arm_dirs[arm] = base_dir / arm

    missing_arms = [arm for arm, arm_dir in resolved_arm_dirs.items() if not arm_dir.is_dir()]
    if missing_arms:
        raise FileNotFoundError(
            "RC3 arm set is incomplete. Missing arm directories: " + ", ".join(missing_arms)
        )

    runs = []
    for arm in ARM_ORDER:
        arm_dir = resolved_arm_dirs[arm]
        for seed in REQUIRED_SEEDS:
            seed_dir = arm_dir / f"seed_{seed}"
            output_dir = seed_dir / "output"
            if not seed_dir.is_dir():
                raise FileNotFoundError(f"Missing RC3 run directory: {seed_dir}")
            if not output_dir.is_dir():
                raise FileNotFoundError(f"Missing RC3 output directory: {output_dir}")
            runs.append(
                {
                    "arm": arm,
                    "seed": seed,
                    "run_dir": seed_dir,
                    "output_dir": output_dir,
                }
            )
    return runs


def evaluate_runs(base_dir: Path) -> List[Dict[str, Any]]:
    evaluated = []
    for run in discover_rc3_runs(base_dir, {}):
        parser = OutputParser(run["output_dir"])
        snapshot, snapshot_file = _load_snapshot(parser, Path(parser.output_dir), RC3_ENDPOINT_TIME)
        metrics = _parse_tumor_metrics(snapshot, parser)
        metrics.update(
            {
                "arm": run["arm"],
                "seed": run["seed"],
                "snapshot_file": str(snapshot_file),
                "output_dir": str(run["output_dir"]),
            }
        )
        evaluated.append(metrics)
    return sorted(evaluated, key=lambda x: (ARM_ORDER.index(x["arm"]), x["seed"]))


def summarize_by_arm(runs: Sequence[Dict[str, Any]]) -> Dict[str, Dict[str, float]]:
    summary: Dict[str, Dict[str, float]] = {}
    for arm in ARM_ORDER:
        arm_runs = [run for run in runs if run["arm"] == arm]
        summary[arm] = {
            "n_runs": len(arm_runs),
        }
        for field in ("live_tumor", "radius_max", "radius_p95", "peri_ecm"):
            values = [float(run[field]) for run in arm_runs]
            mean, sd, n = _mean_sd(values)
            summary[arm][f"{field}_mean"] = mean
            summary[arm][f"{field}_sd"] = sd
            summary[arm][f"{field}_n"] = n
            summary[arm][f"{field}_median"] = _median(values)
    return summary


def choose_benchmark_value(benchmark: Dict[str, Any], seeds: Sequence[int]) -> tuple[float, str]:
    by_seed = benchmark.get("by_seed", {})
    if by_seed and all(seed in by_seed for seed in seeds):
        values = [float(by_seed[seed]) for seed in seeds]
        return _median(values), f"matched-seed median over seeds {', '.join(str(s) for s in seeds)}"
    if math.isfinite(float(benchmark.get("mean", math.nan))):
        runs = benchmark.get("runs", [])
        if runs:
            available_seeds = sorted(int(seed) for seed in by_seed.keys())
            return _median([float(run["day28_live_tumor"]) for run in runs]), (
                "benchmark-set median from "
                f"{benchmark.get('n', 0)} run(s); PROVISIONAL because matched seeds "
                f"{', '.join(str(s) for s in seeds)} are not all available"
                + (f" (available: {', '.join(str(s) for s in available_seeds)})" if available_seeds else "")
            )
        median_value = benchmark.get("median", math.nan)
        if math.isfinite(float(median_value)):
            return float(median_value), (
                f"benchmark summary median from {benchmark.get('n', 0)} run(s); PROVISIONAL because matched seeds are unavailable"
            )
        return float(benchmark["mean"]), (
            f"benchmark summary mean fallback from {benchmark.get('n', 0)} run(s); PROVISIONAL because matched seeds are unavailable"
        )
    raise ValueError("RC2 benchmark did not yield a usable comparison value for C3.5")


def evaluate_checks(summary: Dict[str, Dict[str, float]], benchmark_value: float) -> Dict[str, Dict[str, Any]]:
    a = summary["Arm_A_Control"]
    b = summary["Arm_B_SHH_Only"]
    c = summary["Arm_C_Combo"]

    checks = {
        "C3.1": {
            "passed": b["live_tumor_median"] > a["live_tumor_median"],
            "detail": f"median_tumor_B={_fmt(b['live_tumor_median'])} > median_tumor_A={_fmt(a['live_tumor_median'])}",
        },
        "C3.2": {
            "passed": b["radius_max_median"] > a["radius_max_median"],
            "detail": f"median_radius_B={_fmt(b['radius_max_median'])} > median_radius_A={_fmt(a['radius_max_median'])}",
        },
        "C3.3": {
            "passed": b["peri_ecm_median"] < a["peri_ecm_median"],
            "detail": f"median_periECM_B={_fmt(b['peri_ecm_median'], 4)} < median_periECM_A={_fmt(a['peri_ecm_median'], 4)}",
        },
        "C3.4": {
            "passed": c["live_tumor_median"] < a["live_tumor_median"],
            "detail": f"median_tumor_C={_fmt(c['live_tumor_median'])} < median_tumor_A={_fmt(a['live_tumor_median'])}",
        },
        "C3.5": {
            "passed": c["live_tumor_median"] < benchmark_value,
            "detail": f"median_tumor_C={_fmt(c['live_tumor_median'])} < rc2_benchmark={_fmt(benchmark_value)}",
        },
        "C3.6": {
            "passed": c["live_tumor_median"] < a["live_tumor_median"] < b["live_tumor_median"],
            "detail": (
                f"rank order tumor medians: C={_fmt(c['live_tumor_median'])} < "
                f"A={_fmt(a['live_tumor_median'])} < B={_fmt(b['live_tumor_median'])}"
            ),
        },
    }
    return checks


def interpret_failures(checks: Dict[str, Dict[str, Any]]) -> List[str]:
    messages = []
    if not checks["C3.1"]["passed"]:
        messages.append("B <= A: likely mechanical confinement is too weak.")
    if not checks["C3.3"]["passed"]:
        messages.append("periECM(B) not lower than periECM(A): SHH->GLI1->ECM wiring is broken.")
    if not checks["C3.4"]["passed"]:
        messages.append("C not better than A: thinner stroma is not improving drug access enough.")
    if not checks["C3.5"]["passed"]:
        messages.append("C not better than the RC2 drug-alone benchmark: combo benefit is too weak.")
    if not checks["C3.6"]["passed"]:
        messages.append("C < A < B rank order is broken at the arm-mean level.")
    return messages


def _matched_seed_rows(runs: Sequence[Dict[str, Any]]) -> List[str]:
    by_seed_arm: Dict[tuple[int, str], Dict[str, Any]] = {}
    seeds = sorted({int(run["seed"]) for run in runs})
    for run in runs:
        by_seed_arm[(int(run["seed"]), run["arm"])] = run

    rows = []
    rows.append(
        "Seed  "
        "A_live  A_rad  A_peri   "
        "B_live  B_rad  B_peri   "
        "C_live  C_rad  C_peri"
    )
    for seed in seeds:
        a = by_seed_arm.get((seed, "Arm_A_Control"))
        b = by_seed_arm.get((seed, "Arm_B_SHH_Only"))
        c = by_seed_arm.get((seed, "Arm_C_Combo"))
        rows.append(
            f"{seed:<5} "
            f"{(a['live_tumor'] if a else 'na')!s:<6} {_fmt(a['radius_max']) if a else 'na':<6} {_fmt(a['peri_ecm'], 4) if a else 'na':<8} "
            f"{(b['live_tumor'] if b else 'na')!s:<6} {_fmt(b['radius_max']) if b else 'na':<6} {_fmt(b['peri_ecm'], 4) if b else 'na':<8} "
            f"{(c['live_tumor'] if c else 'na')!s:<6} {_fmt(c['radius_max']) if c else 'na':<6} {_fmt(c['peri_ecm'], 4) if c else 'na':<8}"
        )
    return rows


def build_report(base_dir: Path, runs: Sequence[Dict[str, Any]], summary: Dict[str, Dict[str, float]], benchmark: Dict[str, Any], benchmark_value: float, benchmark_note: str, checks: Dict[str, Dict[str, Any]]) -> str:
    lines: List[str] = []
    lines.append("=" * 100)
    lines.append("RC3 VISMODEGIB PARADOX REPORT")
    lines.append(f"Base: {base_dir}")
    lines.append("=" * 100)
    lines.append("")
    lines.append("Per-Seed Endpoint Table")
    lines.extend(_matched_seed_rows(runs))
    lines.append("")
    lines.append("Per-Run Endpoint Diagnostics")
    lines.append("Arm  Seed  LiveTumor  RadiusMax  RadiusP95  PeriECM    Snapshot")
    for run in runs:
        lines.append(
            f"{ARM_SHORT[run['arm']]:<4}"
            f" {run['seed']:<5}"
            f" {run['live_tumor']:<10}"
            f" {_fmt(run['radius_max']):<10}"
            f" {_fmt(run['radius_p95']):<10}"
            f" {_fmt(run['peri_ecm'], 4):<10}"
            f" {Path(run['snapshot_file']).name}"
        )
    lines.append("")
    lines.append("Arm Summary Table")
    lines.append("Arm  n  LiveTumor(mean±sd)  RadiusMax(mean±sd)  RadiusP95(mean±sd)  PeriECM(mean±sd)")
    for arm in ARM_ORDER:
        s = summary[arm]
        lines.append(
            f"{ARM_SHORT[arm]:<4} {int(s['n_runs']):<2}"
            f" {_fmt(s['live_tumor_mean'])}±{_fmt(s['live_tumor_sd'])}"
            f"    {_fmt(s['radius_max_mean'])}±{_fmt(s['radius_max_sd'])}"
            f"    {_fmt(s['radius_p95_mean'])}±{_fmt(s['radius_p95_sd'])}"
            f"    {_fmt(s['peri_ecm_mean'], 4)}±{_fmt(s['peri_ecm_sd'], 4)}"
        )
    lines.append("")
    lines.append(f"RC2 Benchmark: {benchmark.get('source')} ({benchmark_note})")
    lines.append(f"Benchmark value for C3.5: {_fmt(benchmark_value)}")
    benchmark_runs = benchmark.get("runs", [])
    if benchmark_runs:
        lines.append("Benchmark runs loaded for C3.5:")
        for run in benchmark_runs:
            rep = run.get("replicate", "summary")
            seed = run.get("seed", "na")
            val = run.get("day28_live_tumor", math.nan)
            src = run.get("source", "unknown")
            lines.append(f"- {rep} seed={seed} day28_live_tumor={_fmt(val)} source={src}")
    lines.append("")
    lines.append("RC3 Checks")
    for name in ("C3.1", "C3.2", "C3.3", "C3.4", "C3.5", "C3.6"):
        verdict = "PASS" if checks[name]["passed"] else "FAIL"
        lines.append(f"{name}: {verdict} | {checks[name]['detail']}")
    lines.append("")
    overall = all(bool(check["passed"]) for check in checks.values())
    lines.append(f"Overall Verdict: {'PASS' if overall else 'FAIL'}")
    if overall:
        lines.append("Interpretation: SHH-only worsens outcome and expands the tumor while lowering periECM, while combo outperforms both control and RC2 drug-alone.")
    else:
        lines.append("Failure Interpretation:")
        for message in interpret_failures(checks):
            lines.append(f"- {message}")
    lines.append("")
    return "\n".join(lines)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate RC3 runs.")
    parser.add_argument("--base-dir", help="Base path to rc3_vismodegib/")
    parser.add_argument("--arm-a-dir", help="Optional explicit Arm_A_Control directory to use for evaluation.")
    parser.add_argument("--arm-b-dir", help="Optional explicit Arm_B_SHH_Only directory to use for evaluation.")
    parser.add_argument("--arm-c-dir", help="Optional explicit Arm_C_Combo directory to use for evaluation.")
    parser.add_argument(
        "--rc2-benchmark",
        required=True,
        help="Path to a clean RC2 drug-alone benchmark set or an explicit benchmark summary file.",
    )
    parser.add_argument("--summary-out", help="Optional path to write the full text summary.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    base_dir = Path(args.base_dir).expanduser().resolve() if args.base_dir else None
    if base_dir is not None and not base_dir.exists():
        print(f"ERROR: base dir does not exist: {base_dir}", file=sys.stderr)
        return 1

    arm_dirs = {
        "Arm_A_Control": Path(args.arm_a_dir).expanduser().resolve() if args.arm_a_dir else None,
        "Arm_B_SHH_Only": Path(args.arm_b_dir).expanduser().resolve() if args.arm_b_dir else None,
        "Arm_C_Combo": Path(args.arm_c_dir).expanduser().resolve() if args.arm_c_dir else None,
    }
    if base_dir is None and not all(arm_dirs.values()):
        print(
            "ERROR: provide --base-dir or explicit --arm-a-dir/--arm-b-dir/--arm-c-dir overrides",
            file=sys.stderr,
        )
        return 1

    try:
        runs = []
        for run in discover_rc3_runs(base_dir, arm_dirs):
            parser = OutputParser(run["output_dir"])
            snapshot, snapshot_file = _load_snapshot(parser, Path(parser.output_dir), RC3_ENDPOINT_TIME)
            metrics = _parse_tumor_metrics(snapshot, parser)
            metrics.update(
                {
                    "arm": run["arm"],
                    "seed": run["seed"],
                    "snapshot_file": str(snapshot_file),
                    "output_dir": str(run["output_dir"]),
                }
            )
            runs.append(metrics)
        runs = sorted(runs, key=lambda x: (ARM_ORDER.index(x["arm"]), x["seed"]))
        summary = summarize_by_arm(runs)
        benchmark = load_rc2_benchmark(args.rc2_benchmark)
        combo_seeds = [int(run["seed"]) for run in runs if run["arm"] == "Arm_C_Combo"]
        benchmark_value, benchmark_note = choose_benchmark_value(benchmark, combo_seeds)
        checks = evaluate_checks(summary, benchmark_value)
        report_base = base_dir if base_dir is not None else Path("arm-overrides")
        report = build_report(report_base, runs, summary, benchmark, benchmark_value, benchmark_note, checks)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(report)
    if args.summary_out:
        out_path = Path(args.summary_out).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(report + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
