#!/usr/bin/env python3
"""
Stage 2 RC2 Counter-Measurement Sweep
======================================
7 jobs covering all likely failure modes:

  Job  | Strategy                     | km   | ap     | thr  | uptake | conc | kill_coeff
  -----|------------------------------|------|--------|------|--------|------|----------
  A    | Strong drug + fast resist    | 2.0  | 0.002  | 0.03 | 0.10   | 0.10 | 0.05
  B    | Strong drug + med resist     | 2.0  | 0.001  | 0.03 | 0.10   | 0.10 | 0.05
  C    | Med drug + fast resist       | 1.5  | 0.002  | 0.03 | 0.10   | 0.10 | 0.05
  D    | Balanced middle              | 1.5  | 0.001  | 0.03 | 0.10   | 0.10 | 0.05
  E    | Low threshold (easy NRF2)    | 1.5  | 0.0005 | 0.01 | 0.10   | 0.10 | 0.05
  F    | Force uptake (2x rate)       | 2.0  | 0.001  | 0.02 | 0.25   | 0.10 | 0.05
  G    | Flood drug + low threshold   | 1.5  | 0.001  | 0.02 | 0.15   | 0.25 | 0.08

Failure modes covered:
  - Drug too weak:          A, B, F (higher km)
  - Drug never enters:      F (2.5x uptake), G (1.5x uptake + 2.5x concentration)
  - NRF2 never triggers:    E (threshold 0.01), F,G (threshold 0.02)
  - Resistance too slow:    A, C (ap=0.002)
  - Total wipeout:          C, D, E (moderate km with resistance support)
"""

import os, sys, json, subprocess, shutil
import xml.etree.ElementTree as ET

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
CONFIG_SRC   = os.path.join(PROJECT_ROOT, "config", "PhysiCell_settings.xml")
BINARY       = os.path.join(PROJECT_ROOT, "build", "bin", "release", "stroma_world")
SWEEP_BASE   = "/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_countermeasures"

# ── Variant definitions ─────────────────────────────────────────────
VARIANTS = {
    "A_strong_fast":  {"drug_kill_multiplier": 2.0, "abcb1_production_rate": 0.002,
                       "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
                       "drug_concentration": 0.10, "drug_kill_coefficient": 0.05},

    "B_strong_med":   {"drug_kill_multiplier": 2.0, "abcb1_production_rate": 0.001,
                       "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
                       "drug_concentration": 0.10, "drug_kill_coefficient": 0.05},

    "C_med_fast":     {"drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.002,
                       "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
                       "drug_concentration": 0.10, "drug_kill_coefficient": 0.05},

    "D_balanced":     {"drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
                       "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
                       "drug_concentration": 0.10, "drug_kill_coefficient": 0.05},

    "E_low_thresh":   {"drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.0005,
                       "drug_stress_threshold": 0.01, "drug_uptake_rate": 0.10,
                       "drug_concentration": 0.10, "drug_kill_coefficient": 0.05},

    "F_force_uptake": {"drug_kill_multiplier": 2.0, "abcb1_production_rate": 0.001,
                       "drug_stress_threshold": 0.02, "drug_uptake_rate": 0.25,
                       "drug_concentration": 0.10, "drug_kill_coefficient": 0.05},

    "G_flood_drug":   {"drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
                       "drug_stress_threshold": 0.02, "drug_uptake_rate": 0.15,
                       "drug_concentration": 0.25, "drug_kill_coefficient": 0.08},
}

# ── Stage 2 timing constants ───────────────────────────────────────
T_PRE       = 20160.0   # day 14 — drug ON
T_TREAT_END = 40320.0   # day 28 — drug OFF
T_POST      = 60480.0   # day 42 — end
SAVE_INT    = 360        # 6 h snapshots
SEED        = 42
CPUS        = 32

def patch_config(src, dst, out_dir, params):
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath, val):
        n = root.find(xpath)
        if n is not None:
            n.text = str(val)

    # Standard Stage 2 patches
    _set("./save/folder",                          out_dir)
    _set(".//overall/max_time",                    str(T_POST))
    _set(".//options/random_seed",                 str(SEED))
    _set(".//parallel/omp_num_threads",            str(CPUS))
    _set(".//user_parameters/drug_start_time",     str(T_PRE))
    _set(".//user_parameters/drug_end_time",       str(T_TREAT_END))

    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INT)

    # Dirichlet BCs
    for var_name, val in [("oxygen","38"),("tgfb","0"),("shh","0"),
                          ("drug","0"),("ecm_density","0")]:
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

    # Sweep-specific parameter patches
    for pname, pval in params.items():
        _set(f".//user_parameters/{pname}", str(pval))

    # Also patch drug_concentration in the microenvironment if provided
    if "drug_concentration" in params:
        _set(".//user_parameters/drug_concentration", str(params["drug_concentration"]))

    tree.write(dst, encoding="utf-8", xml_declaration=True)

def submit_variant(name, params):
    work_dir = os.path.join(SWEEP_BASE, name)
    out_dir  = os.path.join(work_dir, "output")
    run_dir  = work_dir
    os.makedirs(out_dir, exist_ok=True)

    # Copy auxiliary config files
    local_cfg = os.path.join(run_dir, "config")
    os.makedirs(local_cfg, exist_ok=True)
    for f in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
        src = os.path.join(PROJECT_ROOT, "config", f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(local_cfg, f))

    # Patch config
    config_dst = os.path.join(run_dir, "config.xml")
    patch_config(CONFIG_SRC, config_dst, out_dir, params)

    # SLURM script
    sbatch_script = os.path.join(work_dir, "run.sh")
    with open(sbatch_script, "w") as fh:
        fh.write(f"""#!/bin/bash
#SBATCH --job-name=rc2_{name[:8]}
#SBATCH --partition=cpu384g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output={work_dir}/slurm_%j.out
#SBATCH --error={work_dir}/slurm_%j.err

set -euo pipefail
module purge 2>/dev/null || true
module load gcc/12 2>/dev/null || true

export OMP_NUM_THREADS={CPUS}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

cd "{run_dir}"
echo "[rc2_{name}] start: $(date)"
"{BINARY}" "{config_dst}"
echo "[rc2_{name}] done:  $(date)"
""")

    result = subprocess.run(["sbatch", sbatch_script],
                            capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.returncode == 0 else "FAILED"
    return job_id

# ── Main ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    os.makedirs(SWEEP_BASE, exist_ok=True)

    manifest = {"sweep_type": "rc2_countermeasures", "variants": {}}
    print(f"{'Variant':<20} {'JobID':<8} {'km':>5} {'ap':>8} {'thr':>6} {'uptk':>6} {'conc':>6} {'kc':>6}")
    print("-" * 75)

    for name, params in VARIANTS.items():
        job_id = submit_variant(name, params)
        manifest["variants"][name] = {"params": params, "job_id": job_id,
                                       "work_dir": os.path.join(SWEEP_BASE, name)}
        print(f"  {name:<20} {job_id:<8} {params['drug_kill_multiplier']:>5.1f} "
              f"{params['abcb1_production_rate']:>8.4f} {params['drug_stress_threshold']:>6.3f} "
              f"{params['drug_uptake_rate']:>6.2f} {params['drug_concentration']:>6.2f} "
              f"{params['drug_kill_coefficient']:>6.3f}")

    manifest_path = os.path.join(SWEEP_BASE, "manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\n{'='*75}")
    print(f"  {len(VARIANTS)} RC2 counter-measurement jobs submitted")
    print(f"  + 1 baseline (job 10620: km=1.0, ap=0.0005, thr=0.03)")
    print(f"  = 8 total Stage 2 RC2 runs in parallel")
    print(f"  Manifest: {manifest_path}")
    print(f"  Monitor:  squeue -u $USER")
    print(f"{'='*75}")
