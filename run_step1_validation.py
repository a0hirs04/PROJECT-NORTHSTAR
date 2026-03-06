#!/usr/bin/env python3
"""
Step 1 Validation — Submit RC1 baseline + RC2 probe for seed 42.

Run A: RC1-style baseline (14 days, no drug) — checks pre-drug tumor viability
Run B: RC2 probe (16 days, drug ON at day 14) — checks drug penetration

Both submitted as parallel SLURM jobs on compute partition.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
BINARY = PROJECT_ROOT / "stroma_world"
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"

SEED = 42
SLURM_PARTITION = "compute"
SLURM_CPUS = 128
SLURM_MEM = "128G"
SLURM_TIME = "04:00:00"
SAVE_INTERVAL = 360  # minutes


def _patch_config(src, dst, output_dir, seed, max_time, drug_start, drug_end, drug_conc):
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath, value):
        node = root.find(xpath)
        if node is not None:
            node.text = str(value)

    _set("./save/folder", str(output_dir))
    _set(".//overall/max_time", str(max_time))
    _set(".//options/random_seed", str(seed))
    _set(".//parallel/omp_num_threads", str(SLURM_CPUS))

    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INTERVAL)

    _set(".//user_parameters/drug_start_time", str(drug_start))
    _set(".//user_parameters/drug_end_time", str(drug_end))
    _set(".//user_parameters/drug_concentration", str(drug_conc))

    # Enforce Dirichlet BCs
    for var_name, val in [("oxygen", "38"), ("tgfb", "0"), ("shh", "0"),
                          ("drug", "0"), ("ecm_density", "0")]:
        var = root.find(f".//microenvironment_setup/variable[@name='{var_name}']")
        if var is None:
            continue
        dbc = var.find("./Dirichlet_boundary_condition")
        if dbc is not None:
            dbc.text = val
            dbc.set("enabled", "true")
        for bv in var.findall("./Dirichlet_options/boundary_value"):
            bv.text = val
            bv.set("enabled", "true")

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def _write_slurm(rep_dir, config_path, job_name):
    script = rep_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name={job_name}
        #SBATCH --partition={SLURM_PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={SLURM_CPUS}
        #SBATCH --mem={SLURM_MEM}
        #SBATCH --time={SLURM_TIME}
        #SBATCH --output={rep_dir}/slurm_%j.out
        #SBATCH --error={rep_dir}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={SLURM_CPUS}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores

        cd {PROJECT_ROOT}
        echo "=== {job_name} ==="
        echo "Start: $(date)"
        {BINARY} {config_path}
        echo "End:   $(date)"
    """))
    script.chmod(0o755)
    return script


def main():
    if not BINARY.exists():
        print(f"ERROR: binary not found: {BINARY}")
        return 1

    work_dir = PROJECT_ROOT / "build" / "step1_validation"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    # --- Run A: RC1 baseline (14 days, no drug) ---
    run_a = work_dir / "rc1_baseline_seed42"
    run_a.mkdir(parents=True)
    out_a = run_a / "output"
    out_a.mkdir()

    config_a = run_a / "config.xml"
    _patch_config(BASE_CONFIG, config_a, out_a, SEED,
                  max_time=20160.0,           # 14 days
                  drug_start=999999.0,        # no drug
                  drug_end=999999.0,
                  drug_conc=0.0)

    # Copy calibration knobs
    for run_dir in [run_a]:
        local_cfg = run_dir / "config"
        local_cfg.mkdir(exist_ok=True)
        for f in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
            src = PROJECT_ROOT / "config" / f
            if src.exists():
                shutil.copy(src, local_cfg / f)

    script_a = _write_slurm(run_a, config_a, "step1_rc1_s42")

    # --- Run B: RC2 probe (16 days, drug ON at day 14) ---
    run_b = work_dir / "rc2_probe_seed42"
    run_b.mkdir(parents=True)
    out_b = run_b / "output"
    out_b.mkdir()

    config_b = run_b / "config.xml"
    _patch_config(BASE_CONFIG, config_b, out_b, SEED,
                  max_time=23040.0,           # 16 days
                  drug_start=20160.0,         # drug ON at day 14
                  drug_end=40320.0,           # drug OFF at day 28 (beyond sim end)
                  drug_conc=1.0)

    local_cfg_b = run_b / "config"
    local_cfg_b.mkdir(exist_ok=True)
    for f in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
        src = PROJECT_ROOT / "config" / f
        if src.exists():
            shutil.copy(src, local_cfg_b / f)

    script_b = _write_slurm(run_b, config_b, "step1_rc2probe_s42")

    # --- Submit both ---
    jobs = []
    for label, script in [("RC1 baseline", script_a), ("RC2 probe", script_b)]:
        result = subprocess.run(
            ["sbatch", "--parsable", str(script)],
            capture_output=True, text=True, check=False,
        )
        if result.returncode != 0:
            print(f"ERROR: sbatch failed for {label}: {result.stderr.strip()}")
            return 1
        job_id = result.stdout.strip().split(";")[0]
        jobs.append((label, job_id))
        print(f"  Submitted {label} -> SLURM job {job_id}")

    print(f"\n  Both jobs submitted. Monitor with:")
    for label, jid in jobs:
        print(f"    sacct -j {jid} --format=JobID,State,Elapsed,ExitCode")
    print(f"\n  Output dirs:")
    print(f"    RC1: {out_a}")
    print(f"    RC2: {out_b}")
    print(f"\n  When done, run: python diagnose_step1.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
