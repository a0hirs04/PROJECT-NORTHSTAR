#!/usr/bin/env python3
"""
Reality Check 3 launcher: vismodegib paradox.

Creates RC3 run directories from the frozen RC2 gold-master config, with
 either a full 15-run layout or a 2-run smoke test layout. By default the
script prepares configs and SLURM scripts but does not submit anything.
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
import textwrap
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from python.wrapper.output_parser import OutputParser


PROJECT_ROOT = Path(__file__).resolve().parent
BINARY = PROJECT_ROOT / "stroma_world"
DEFAULT_PARENT_CANDIDATES = [
    PROJECT_ROOT / "build" / "config" / "PhysiCell_settings.xml",
    PROJECT_ROOT / "artifacts" / "rc2_true_biological_pass_km0p01_seed42" / "winning_config.xml",
    PROJECT_ROOT / "config" / "PhysiCell_settings.xml",
]

DEFAULT_FULL_ROOT = PROJECT_ROOT / "build" / "rc3_vismodegib"
DEFAULT_SMOKE_ROOT = PROJECT_ROOT / "build" / "rc3_vismodegib_smoke"
DEFAULT_SEEDS = [42, 43, 44, 45, 46]

SLURM_PARTITION = "cpu384g"
SLURM_CPUS = 32
SLURM_MEM = "128G"
SLURM_TIME = "06:00:00"
RC3_MAX_TIME = 80640.0
SMOKE_SIGNAL_CHECK_TIME = 21600.0
SAVE_INTERVAL = 360.0
TERMINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
    "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
    "BOOT_FAIL", "DEADLINE",
}


@dataclass(frozen=True)
class ArmSpec:
    name: str
    drug_concentration: float
    shh_inhibition_strength: float
    shh_inhibition_start_time: float
    drug_end_time: float | None = None


@dataclass(frozen=True)
class RunSpec:
    arm: ArmSpec
    seed: int
    run_dir: Path
    output_dir: Path
    config_path: Path
    slurm_script: Path


@dataclass(frozen=True)
class ParentConfigSummary:
    drug_start_time: float
    drug_end_time: float
    drug_concentration: float
    drug_kill_multiplier: float
    shh_inhibition_start_time: float
    shh_inhibition_strength: float


@dataclass(frozen=True)
class SmokeSnapshot:
    time_min: float
    stromal_gli1_mean: float
    stromal_gli1_positive_frac: float
    shh_at_stroma_mean: float
    peri_ecm: float


ARM_A = ArmSpec(
    name="Arm_A_Control",
    drug_concentration=0.0,
    shh_inhibition_strength=0.0,
    shh_inhibition_start_time=1e18,
)
ARM_B = ArmSpec(
    name="Arm_B_SHH_Only",
    drug_concentration=0.0,
    shh_inhibition_strength=1.0,
    shh_inhibition_start_time=20160.0,
)
ARM_C = ArmSpec(
    name="Arm_C_Combo",
    drug_concentration=1.0,
    shh_inhibition_strength=1.0,
    shh_inhibition_start_time=20160.0,
    drug_end_time=40320.0,
)
ALL_ARMS = [ARM_A, ARM_B, ARM_C]
ARM_BY_NAME = {arm.name: arm for arm in ALL_ARMS}


def resolve_parent_config(explicit: str | None) -> tuple[Path, str]:
    if explicit:
        candidate = Path(explicit).expanduser().resolve()
        if not candidate.exists():
            raise FileNotFoundError(f"Requested parent config does not exist: {candidate}")
        return candidate, "explicit"

    for candidate in DEFAULT_PARENT_CANDIDATES:
        resolved = candidate.expanduser().resolve()
        if resolved.exists():
            if candidate == DEFAULT_PARENT_CANDIDATES[0]:
                return resolved, "default build/config"
            if candidate == DEFAULT_PARENT_CANDIDATES[1]:
                return resolved, "frozen RC2 gold-master artifact"
            return resolved, "repo config fallback"

    tried = "\n".join(str(p) for p in DEFAULT_PARENT_CANDIDATES)
    raise FileNotFoundError(f"No usable parent config found. Tried:\n{tried}")


def read_parent_summary(parent_config: Path) -> ParentConfigSummary:
    tree = ET.parse(parent_config)
    root = tree.getroot()

    def read_float(xpath: str) -> float:
        node = root.find(xpath)
        if node is None or node.text is None:
            raise ValueError(f"Missing XML node for xpath: {xpath}")
        return float(node.text)

    return ParentConfigSummary(
        drug_start_time=read_float(".//user_parameters/drug_start_time"),
        drug_end_time=read_float(".//user_parameters/drug_end_time"),
        drug_concentration=read_float(".//user_parameters/drug_concentration"),
        drug_kill_multiplier=read_float(".//user_parameters/drug_kill_multiplier"),
        shh_inhibition_start_time=read_float(".//user_parameters/shh_inhibition_start_time"),
        shh_inhibition_strength=read_float(".//user_parameters/shh_inhibition_strength"),
    )


def validate_parent_summary(summary: ParentConfigSummary) -> None:
    expected = {
        "drug_start_time": 20160.0,
        "drug_end_time": 40320.0,
        "drug_concentration": 1.0,
        "drug_kill_multiplier": 0.01,
    }
    for name, expected_value in expected.items():
        actual = getattr(summary, name)
        if not math.isclose(actual, expected_value, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(
                f"Parent config is not the frozen RC2 gold-master: "
                f"{name}={actual} (expected {expected_value})"
            )


def patch_config(parent_config: Path, run: RunSpec) -> None:
    tree = ET.parse(parent_config)
    root = tree.getroot()

    def set_text(xpath: str, value: float | int | str) -> None:
        node = root.find(xpath)
        if node is None:
            raise ValueError(f"Missing XML node for xpath: {xpath}")
        node.text = str(value)

    set_text("./save/folder", str(run.output_dir))
    set_text(".//overall/max_time", RC3_MAX_TIME)
    set_text(".//options/random_seed", run.seed)
    set_text(".//parallel/omp_num_threads", SLURM_CPUS)
    set_text(".//user_parameters/drug_concentration", run.arm.drug_concentration)
    set_text(".//user_parameters/shh_inhibition_strength", run.arm.shh_inhibition_strength)
    set_text(".//user_parameters/shh_inhibition_start_time", run.arm.shh_inhibition_start_time)
    if run.arm.drug_end_time is not None:
        set_text(".//user_parameters/drug_end_time", run.arm.drug_end_time)

    tree.write(run.config_path, encoding="utf-8", xml_declaration=True)


def _query_job(job_id: str) -> tuple[str, int, bool]:
    try:
        result = subprocess.run(
            ["sacct", "-j", job_id, "--format=State,ExitCode", "--noheader", "--parsable2"],
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
        for line in result.stdout.strip().splitlines():
            parts = line.strip().split("|")
            if len(parts) >= 2:
                state = parts[0].strip().split()[0] if parts[0].strip() else "UNKNOWN"
                try:
                    exit_code = int(parts[1].split(":")[0])
                except (ValueError, IndexError):
                    exit_code = -1
                return state, exit_code, state in TERMINAL_STATES
    except Exception:
        pass

    try:
        result = subprocess.run(
            ["squeue", "-j", job_id, "--format=%T", "--noheader"],
            capture_output=True,
            text=True,
            check=False,
            timeout=15,
        )
        state = result.stdout.strip() or "COMPLETED"
        return state, 0, state in TERMINAL_STATES
    except Exception:
        return "UNKNOWN", -1, False


def wait_for_jobs(submitted: Sequence[tuple[RunSpec, str]], poll_seconds: int = 20) -> Dict[str, tuple[str, int]]:
    pending = {job_id: run for run, job_id in submitted}
    final: Dict[str, tuple[str, int]] = {}
    while pending:
        for job_id in list(pending.keys()):
            state, exit_code, terminal = _query_job(job_id)
            if terminal:
                final[job_id] = (state, exit_code)
                pending.pop(job_id, None)
        if pending:
            time.sleep(poll_seconds)
    return final


def _load_snapshot(output_dir: Path, target_time_min: float) -> tuple[OutputParser, Dict[str, Any]]:
    parser = OutputParser(output_dir)
    idx = int(round(target_time_min / SAVE_INTERVAL))
    for i in range(idx, -1, -1):
        candidate = Path(parser.output_dir) / f"output{i:08d}.xml"
        if candidate.exists() and candidate.stat().st_size > 0:
            return parser, parser._read_physicell_xml(candidate)
    final_xml = Path(parser.output_dir) / "final.xml"
    if final_xml.exists() and final_xml.stat().st_size > 0:
        return parser, parser._read_physicell_xml(final_xml)
    raise FileNotFoundError(f"No readable snapshot at or before t={target_time_min} in {output_dir}")


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


def parse_smoke_snapshot(output_dir: Path, target_time_min: float) -> SmokeSnapshot:
    parser, snapshot = _load_snapshot(output_dir, target_time_min)
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

    tumor_type_id, stromal_type_id = validate_cell_type_mapping(snapshot)
    ctype = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1, dtype=int)
    live_tumor = live_mask & (ctype == tumor_type_id)
    live_stroma = live_mask & (ctype == stromal_type_id)

    positions = parser._get_positions(matrix, labels)
    tumor_pos = positions[live_tumor] if positions.size and np.any(live_tumor) else np.empty((0, 3))
    stroma_pos = positions[live_stroma] if positions.size and np.any(live_stroma) else np.empty((0, 3))

    gli1 = _row(matrix, labels, "gli1_active")
    if gli1 is None:
        raise ValueError(
            "Snapshot is missing the real stromal 'gli1_active' label. "
            "RC3 smoke validation refuses to substitute tumor markers for stromal GLI1."
        )
    gli1_vals = gli1[live_stroma] if np.any(live_stroma) else np.array([], dtype=float)
    gli1_mean = float(np.mean(gli1_vals)) if gli1_vals.size else math.nan
    gli1_frac = float(np.mean(gli1_vals > 0.5)) if gli1_vals.size else math.nan

    shh_field = micro_values.get("shh")
    shh_at_stroma = math.nan
    if shh_field is not None and stroma_pos.shape[0] > 0 and micro_coords.size > 0:
        sampled = _sample_field_at_positions(stroma_pos, micro_coords, shh_field)
        if sampled.size > 0:
            shh_at_stroma = float(np.mean(sampled))

    peri_ecm = math.nan
    ecm_field = micro_values.get("ecm_density")
    if ecm_field is not None and tumor_pos.shape[0] > 0 and micro_coords.size > 0:
        shell_mask = np.zeros(micro_coords.shape[0], dtype=bool)
        for i in range(0, micro_coords.shape[0], 5000):
            v = micro_coords[i : i + 5000]
            d2 = np.sum((v[:, None, :] - tumor_pos[None, :, :]) ** 2, axis=2)
            dmin = np.min(d2, axis=1)
            shell_mask[i : i + 5000] = (dmin > 0.0) & (dmin <= 100.0**2)
        if np.any(shell_mask):
            peri_ecm = float(np.nanmean(ecm_field[shell_mask]))
    if ecm_field is not None and ecm_field.size > 0:
        global_ecm_max = float(np.nanmax(ecm_field))
        if global_ecm_max == 0.0 and math.isfinite(peri_ecm) and peri_ecm != 0.0:
            raise ValueError(
                "ECM consistency check failed: global ecm_density max is 0 but periECM is "
                f"{peri_ecm} in {snapshot['filepath']}"
            )

    return SmokeSnapshot(
        time_min=float(snapshot["time"]),
        stromal_gli1_mean=gli1_mean,
        stromal_gli1_positive_frac=gli1_frac,
        shh_at_stroma_mean=shh_at_stroma,
        peri_ecm=peri_ecm,
    )


def run_smoke_report(root_dir: Path) -> tuple[bool, str]:
    arm_a_output = root_dir / ARM_A.name / "seed_42" / "output"
    arm_b_output = root_dir / ARM_B.name / "seed_42" / "output"
    if not arm_a_output.exists() or not arm_b_output.exists():
        raise FileNotFoundError("Smoke-test outputs for Arm A and Arm B were not found.")

    a_pre = parse_smoke_snapshot(arm_a_output, 20160.0)
    b_pre = parse_smoke_snapshot(arm_b_output, 20160.0)
    a_post = parse_smoke_snapshot(arm_a_output, SMOKE_SIGNAL_CHECK_TIME)
    b_post = parse_smoke_snapshot(arm_b_output, SMOKE_SIGNAL_CHECK_TIME)
    a_final = parse_smoke_snapshot(arm_a_output, RC3_MAX_TIME)
    b_final = parse_smoke_snapshot(arm_b_output, RC3_MAX_TIME)

    checks = {
        "gli1_drop_vs_control": (
            math.isfinite(b_post.stromal_gli1_mean)
            and math.isfinite(a_post.stromal_gli1_mean)
            and b_post.stromal_gli1_mean < a_post.stromal_gli1_mean
        ),
        "gli1_drop_vs_pre": (
            math.isfinite(b_post.stromal_gli1_mean)
            and math.isfinite(b_pre.stromal_gli1_mean)
            and b_post.stromal_gli1_mean < b_pre.stromal_gli1_mean
        ),
        "gli1_positive_drop": (
            math.isfinite(b_post.stromal_gli1_positive_frac)
            and math.isfinite(a_post.stromal_gli1_positive_frac)
            and b_post.stromal_gli1_positive_frac < a_post.stromal_gli1_positive_frac
        ),
        "peri_ecm_drop": (
            math.isfinite(b_final.peri_ecm)
            and math.isfinite(a_final.peri_ecm)
            and b_final.peri_ecm < a_final.peri_ecm
        ),
    }
    overall = all(checks.values())

    lines = [
        "RC3 SHH smoke test",
        f"Signal check time: t={SMOKE_SIGNAL_CHECK_TIME:.0f} min",
        "",
        "Arm A / Arm B comparison",
        f"  pre GLI1 mean:        A={a_pre.stromal_gli1_mean:.4f}  B={b_pre.stromal_gli1_mean:.4f}",
        f"  post GLI1 mean:       A={a_post.stromal_gli1_mean:.4f}  B={b_post.stromal_gli1_mean:.4f}",
        f"  post GLI1+ frac:      A={a_post.stromal_gli1_positive_frac:.4f}  B={b_post.stromal_gli1_positive_frac:.4f}",
        f"  post SHH@stroma mean: A={a_post.shh_at_stroma_mean:.4f}  B={b_post.shh_at_stroma_mean:.4f}",
        f"  final periECM:        A={a_final.peri_ecm:.4f}  B={b_final.peri_ecm:.4f}",
        "",
        f"  gli1_drop_vs_control: {'PASS' if checks['gli1_drop_vs_control'] else 'FAIL'}",
        f"  gli1_drop_vs_pre:     {'PASS' if checks['gli1_drop_vs_pre'] else 'FAIL'}",
        f"  gli1_positive_drop:   {'PASS' if checks['gli1_positive_drop'] else 'FAIL'}",
        f"  peri_ecm_drop:        {'PASS' if checks['peri_ecm_drop'] else 'FAIL'}",
        "",
        "Note: the current intervention attenuates stromal SHH sensing, so the decisive smoke metric is the stromal GLI1 readout rather than the extracellular SHH field itself.",
        "",
        f"Overall smoke verdict: {'PASS' if overall else 'FAIL'}",
    ]
    return overall, "\n".join(lines)


def copy_local_config_bundle(run: RunSpec) -> None:
    local_cfg = run.run_dir / "config"
    local_cfg.mkdir(parents=True, exist_ok=True)
    for name in ("tumor_calibration_knobs.json", "gene_params_default.json"):
        src = PROJECT_ROOT / "config" / name
        if src.exists():
            shutil.copy2(src, local_cfg / name)


def write_slurm_script(run: RunSpec) -> None:
    job_name = f"rc3_{run.arm.name.lower()}_s{run.seed}".replace("arm_", "").replace("__", "_")
    run.slurm_script.write_text(
        textwrap.dedent(
            f"""\
            #!/bin/bash
            #SBATCH --job-name={job_name}
            #SBATCH --partition={SLURM_PARTITION}
            #SBATCH --nodes=1
            #SBATCH --ntasks-per-node=1
            #SBATCH --cpus-per-task={SLURM_CPUS}
            #SBATCH --mem={SLURM_MEM}
            #SBATCH --time={SLURM_TIME}
            #SBATCH --output={run.run_dir}/slurm_%j.out
            #SBATCH --error={run.run_dir}/slurm_%j.err

            set -euo pipefail
            export OMP_NUM_THREADS={SLURM_CPUS}
            export OMP_PROC_BIND=spread
            export OMP_PLACES=cores

            cd {PROJECT_ROOT}
            echo "=== RC3 {run.arm.name} seed {run.seed} ==="
            echo "Start: $(date)"
            {BINARY} {run.config_path}
            echo "End:   $(date)"
            """
        ),
        encoding="utf-8",
    )
    run.slurm_script.chmod(0o755)


def build_run_specs(root_dir: Path, arms: Sequence[ArmSpec], seeds: Sequence[int]) -> List[RunSpec]:
    specs: List[RunSpec] = []
    for arm in arms:
        for seed in seeds:
            run_dir = root_dir / arm.name / f"seed_{seed}"
            specs.append(
                RunSpec(
                    arm=arm,
                    seed=seed,
                    run_dir=run_dir,
                    output_dir=run_dir / "output",
                    config_path=run_dir / "config.xml",
                    slurm_script=run_dir / "run.slurm.sh",
                )
            )
    return specs


def prepare_run_layout(root_dir: Path, runs: Sequence[RunSpec], parent_config: Path, force_clean: bool) -> None:
    if root_dir.exists():
        if not force_clean:
            raise FileExistsError(
                f"Target root already exists: {root_dir}\n"
                "Use --force-clean to remove and recreate it."
            )
        shutil.rmtree(root_dir)

    for run in runs:
        run.output_dir.mkdir(parents=True, exist_ok=True)
        patch_config(parent_config, run)
        copy_local_config_bundle(run)
        write_slurm_script(run)


def submit_runs(runs: Sequence[RunSpec]) -> List[tuple[RunSpec, str]]:
    submitted: List[tuple[RunSpec, str]] = []
    for run in runs:
        cmd = ["sbatch", "--parsable", str(run.slurm_script)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(
                f"sbatch failed for {run.arm.name} seed {run.seed}: {result.stderr.strip()}"
            )
        job_id = result.stdout.strip().split(";")[0]
        submitted.append((run, job_id))
    return submitted


def print_plan(
    *,
    parent_config: Path,
    parent_source: str,
    parent_summary: ParentConfigSummary,
    root_dir: Path,
    runs: Sequence[RunSpec],
    submitted: Sequence[tuple[RunSpec, str]] | None = None,
) -> None:
    print(f"RC3 parent config: {parent_config} ({parent_source})")
    print(
        "Parent regimen: "
        f"drug_start={parent_summary.drug_start_time:.0f}, "
        f"drug_end={parent_summary.drug_end_time:.0f}, "
        f"drug_concentration={parent_summary.drug_concentration:.2f}, "
        f"drug_kill_multiplier={parent_summary.drug_kill_multiplier:.4f}, "
        f"shh_inhibition_start={parent_summary.shh_inhibition_start_time:.0f}, "
        f"shh_inhibition_strength={parent_summary.shh_inhibition_strength:.2f}"
    )
    print(f"RC3 root: {root_dir}")
    print(f"Total runs: {len(runs)}")
    print("")

    print("Arm parameter diffs from parent:")
    selected_arms: Dict[str, ArmSpec] = {}
    for run in runs:
        selected_arms[run.arm.name] = run.arm
    for arm_name in [arm.name for arm in ALL_ARMS if arm.name in selected_arms]:
        arm = selected_arms[arm_name]
        print(
            f"  {arm.name}: "
            f"max_time={RC3_MAX_TIME:.0f}, "
            f"drug_concentration={arm.drug_concentration:.2f}, "
            f"drug_end_time={(arm.drug_end_time if arm.drug_end_time is not None else parent_summary.drug_end_time):.0f}, "
            f"shh_inhibition_strength={arm.shh_inhibition_strength:.2f}, "
            f"shh_inhibition_start_time={arm.shh_inhibition_start_time:.0f}"
        )
    print("")

    submitted_map = {id(run): job_id for run, job_id in (submitted or [])}
    for run in runs:
        launch_cmd = f"sbatch --parsable {run.slurm_script}"
        print(f"[{run.arm.name}] seed {run.seed}")
        print(f"  config:  {run.config_path}")
        print(f"  output:  {run.output_dir}")
        if id(run) in submitted_map:
            job_id = submitted_map[id(run)]
            print(f"  logs:    {run.run_dir / f'slurm_{job_id}.out'}")
            print(f"           {run.run_dir / f'slurm_{job_id}.err'}")
            print(f"  launch:  {launch_cmd}")
            print(f"  job_id:  {job_id}")
        else:
            print(f"  logs:    {run.run_dir / 'slurm_%j.out'}")
            print(f"           {run.run_dir / 'slurm_%j.err'}")
            print(f"  launch:  {launch_cmd}")
        print("")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare or launch RC3 runs.")
    parser.add_argument(
        "--parent-config",
        help="Optional explicit parent config. Defaults to the frozen RC2 gold-master artifact if available.",
    )
    parser.add_argument(
        "--root-dir",
        default=str(DEFAULT_FULL_ROOT),
        help="Root directory for the full 15-run RC3 layout.",
    )
    parser.add_argument(
        "--smoke-root-dir",
        default=str(DEFAULT_SMOKE_ROOT),
        help="Root directory for the 2-run RC3 smoke test layout.",
    )
    parser.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=DEFAULT_SEEDS,
        help="Seed list for the full RC3 sweep. Default: 42 43 44 45 46",
    )
    parser.add_argument(
        "--arms",
        nargs="+",
        choices=[arm.name for arm in ALL_ARMS],
        help="Optional subset of RC3 arms to prepare or launch. Defaults to all three arms.",
    )
    parser.add_argument(
        "--arm-c-drug-end-time",
        type=float,
        help="Optional override for Arm_C_Combo drug_end_time, used for schedule-only reruns.",
    )
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Prepare or launch only the Arm A / Arm B seed-42 smoke test.",
    )
    parser.add_argument(
        "--launch",
        "--submit",
        dest="launch",
        action="store_true",
        help="Submit the prepared SLURM scripts after writing configs.",
    )
    parser.add_argument(
        "--wait",
        action="store_true",
        help="Wait for submitted jobs to finish. In smoke-test mode this also prints the smoke report.",
    )
    parser.add_argument(
        "--force-clean",
        action="store_true",
        help="Remove the target root directory before recreating the RC3 layout.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])

    if not BINARY.exists():
        print(f"ERROR: binary not found: {BINARY}", file=sys.stderr)
        return 1
    if args.wait and not args.launch:
        print("ERROR: --wait requires --launch/--submit", file=sys.stderr)
        return 1

    try:
        parent_config, parent_source = resolve_parent_config(args.parent_config)
        parent_summary = read_parent_summary(parent_config)
        validate_parent_summary(parent_summary)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.smoke_test:
        root_dir = Path(args.smoke_root_dir).expanduser().resolve()
        arms = [ARM_A, ARM_B]
        seeds = [42]
    else:
        root_dir = Path(args.root_dir).expanduser().resolve()
        selected_names = args.arms if args.arms else [arm.name for arm in ALL_ARMS]
        arms = [ARM_BY_NAME[name] for name in selected_names]
        if args.arm_c_drug_end_time is not None:
            arms = [
                replace(arm, drug_end_time=args.arm_c_drug_end_time)
                if arm.name == "Arm_C_Combo" else arm
                for arm in arms
            ]
        seeds = list(args.seeds)

    runs = build_run_specs(root_dir, arms, seeds)

    try:
        prepare_run_layout(root_dir, runs, parent_config, args.force_clean)
        submitted = submit_runs(runs) if args.launch else None
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print_plan(
        parent_config=parent_config,
        parent_source=parent_source,
        parent_summary=parent_summary,
        root_dir=root_dir,
        runs=runs,
        submitted=submitted,
    )

    if args.launch and args.wait and submitted is not None:
        final_states = wait_for_jobs(submitted)
        print("Final job states:")
        for run, job_id in submitted:
            state, exit_code = final_states[job_id]
            print(f"  {run.arm.name} seed {run.seed}: {state} (exit={exit_code})")
        print("")

        if args.smoke_test:
            try:
                passed, smoke_report = run_smoke_report(root_dir)
            except Exception as exc:
                print(f"ERROR: smoke test evaluation failed: {exc}", file=sys.stderr)
                return 1
            print(smoke_report)
            return 0 if passed else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
