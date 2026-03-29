#!/usr/bin/env python3
"""
Prepare or launch a clean 5-seed RC2 gold benchmark set.

This pins the run to the frozen RC2 gold-master config and only varies the
random seed across replicates.
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import textwrap
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parent
BINARY = PROJECT_ROOT / "stroma_world"
DEFAULT_PARENT = PROJECT_ROOT / "artifacts" / "rc2_true_biological_pass_km0p01_seed42" / "winning_config.xml"
DEFAULT_ROOT = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_gold_benchmark_5seed")
DEFAULT_SEEDS = [42, 43, 44, 45, 46]

SLURM_PARTITION = "cpu384g"
SLURM_CPUS = 32
SLURM_MEM = "128G"
SLURM_TIME = "06:00:00"

EXPECTED = {
    "max_time": 60480.0,
    "drug_start_time": 20160.0,
    "drug_end_time": 40320.0,
    "drug_concentration": 1.0,
    "drug_kill_multiplier": 0.01,
    "shh_inhibition_strength": 0.0,
}


@dataclass(frozen=True)
class RunSpec:
    replicate_index: int
    seed: int
    run_dir: Path
    output_dir: Path
    config_path: Path
    slurm_script: Path


def _read_float(root: ET.Element, xpath: str) -> float:
    node = root.find(xpath)
    if node is None or node.text is None:
        raise ValueError(f"Missing XML node for xpath: {xpath}")
    return float(node.text)


def validate_parent_config(path: Path) -> None:
    root = ET.parse(path).getroot()
    actual = {
        "max_time": _read_float(root, ".//overall/max_time"),
        "drug_start_time": _read_float(root, ".//user_parameters/drug_start_time"),
        "drug_end_time": _read_float(root, ".//user_parameters/drug_end_time"),
        "drug_concentration": _read_float(root, ".//user_parameters/drug_concentration"),
        "drug_kill_multiplier": _read_float(root, ".//user_parameters/drug_kill_multiplier"),
        "shh_inhibition_strength": _read_float(root, ".//user_parameters/shh_inhibition_strength"),
    }
    for key, expected in EXPECTED.items():
        if not math.isclose(actual[key], expected, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(
                f"Parent config is not the frozen RC2 gold baseline: {key}={actual[key]} (expected {expected})"
            )


def build_run_specs(root_dir: Path, seeds: Sequence[int]) -> list[RunSpec]:
    specs = []
    for i, seed in enumerate(seeds, start=1):
        run_dir = root_dir / f"replicate_{i:02d}_seed{seed}"
        specs.append(
            RunSpec(
                replicate_index=i,
                seed=seed,
                run_dir=run_dir,
                output_dir=run_dir / "output",
                config_path=run_dir / "config.xml",
                slurm_script=run_dir / "run.slurm.sh",
            )
        )
    return specs


def patch_config(parent_config: Path, run: RunSpec) -> None:
    tree = ET.parse(parent_config)
    root = tree.getroot()

    def set_text(xpath: str, value: str | int | float) -> None:
        node = root.find(xpath)
        if node is None:
            raise ValueError(f"Missing XML node for xpath: {xpath}")
        node.text = str(value)

    set_text("./save/folder", str(run.output_dir))
    set_text(".//options/random_seed", run.seed)
    set_text(".//parallel/omp_num_threads", SLURM_CPUS)
    tree.write(run.config_path, encoding="utf-8", xml_declaration=True)


def copy_local_config_bundle(run: RunSpec) -> None:
    local_cfg = run.run_dir / "config"
    local_cfg.mkdir(parents=True, exist_ok=True)
    for name in ("tumor_calibration_knobs.json", "gene_params_default.json"):
        src = PROJECT_ROOT / "config" / name
        if src.exists():
            shutil.copy2(src, local_cfg / name)


def write_slurm_script(run: RunSpec) -> None:
    job_name = f"rc2_gold_r{run.replicate_index}_s{run.seed}"
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
            echo "=== RC2 Gold Benchmark replicate {run.replicate_index} seed {run.seed} ==="
            echo "Start: $(date)"
            {BINARY} {run.config_path}
            echo "End:   $(date)"
            """
        ),
        encoding="utf-8",
    )
    run.slurm_script.chmod(0o755)


def prepare_layout(root_dir: Path, parent_config: Path, runs: Sequence[RunSpec], force_clean: bool) -> None:
    if root_dir.exists():
        if not force_clean:
            raise FileExistsError(
                f"Target root already exists: {root_dir}\nUse --force-clean to remove and recreate it."
            )
        shutil.rmtree(root_dir)

    for run in runs:
        run.output_dir.mkdir(parents=True, exist_ok=True)
        patch_config(parent_config, run)
        copy_local_config_bundle(run)
        write_slurm_script(run)


def submit_runs(runs: Sequence[RunSpec]) -> list[tuple[RunSpec, str]]:
    submitted = []
    for run in runs:
        result = subprocess.run(
            ["sbatch", "--parsable", str(run.slurm_script)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"sbatch failed for replicate {run.replicate_index} seed {run.seed}: {result.stderr.strip()}"
            )
        submitted.append((run, result.stdout.strip().split(";")[0]))
    return submitted


def print_plan(parent_config: Path, root_dir: Path, seeds: Sequence[int], runs: Sequence[RunSpec], submitted: Sequence[tuple[RunSpec, str]] | None = None) -> None:
    print(f"RC2 gold parent config: {parent_config}")
    print(f"RC2 gold output root: {root_dir}")
    print(f"Seeds: {list(seeds)}")
    print("")

    submitted_map = {id(run): job_id for run, job_id in (submitted or [])}
    for run in runs:
        print(f"[replicate {run.replicate_index:02d}] seed {run.seed}")
        print(f"  config:  {run.config_path}")
        print(f"  output:  {run.output_dir}")
        if id(run) in submitted_map:
            job_id = submitted_map[id(run)]
            print(f"  logs:    {run.run_dir / f'slurm_{job_id}.out'}")
            print(f"           {run.run_dir / f'slurm_{job_id}.err'}")
            print(f"  launch:  sbatch --parsable {run.slurm_script}")
            print(f"  job_id:  {job_id}")
        else:
            print(f"  logs:    {run.run_dir / 'slurm_%j.out'}")
            print(f"           {run.run_dir / 'slurm_%j.err'}")
            print(f"  launch:  sbatch --parsable {run.slurm_script}")
        print("")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare or launch the 5-seed RC2 gold benchmark set.")
    parser.add_argument("--parent-config", default=str(DEFAULT_PARENT), help="Frozen RC2 gold parent config.")
    parser.add_argument("--root-dir", default=str(DEFAULT_ROOT), help="Output root for the benchmark set.")
    parser.add_argument("--seeds", nargs="+", type=int, default=DEFAULT_SEEDS, help="Seed list to benchmark.")
    parser.add_argument("--submit", action="store_true", help="Submit the prepared runs to SLURM.")
    parser.add_argument("--force-clean", action="store_true", help="Delete and recreate the output root.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    parent_config = Path(args.parent_config).expanduser().resolve()
    root_dir = Path(args.root_dir).expanduser().resolve()
    seeds = list(dict.fromkeys(args.seeds))

    if not BINARY.exists():
        print(f"ERROR: binary not found: {BINARY}")
        return 1
    if not parent_config.exists():
        print(f"ERROR: parent config not found: {parent_config}")
        return 1

    try:
        validate_parent_config(parent_config)
        runs = build_run_specs(root_dir, seeds)
        prepare_layout(root_dir, parent_config, runs, args.force_clean)
        submitted = submit_runs(runs) if args.submit else None
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1

    print_plan(parent_config, root_dir, seeds, runs, submitted)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
