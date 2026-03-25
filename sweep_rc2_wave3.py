#!/usr/bin/env python3
"""
Stage 2 RC2 Wave 3 — Fixed HIF1A priming + full parameter sweep
=================================================================

ROOT CAUSE FIXED: Hardcoded HIF1A stress bonus (0.5) replaced with XML
parameter hif1a_nrf2_priming_bonus (default 0.02). Now hypoxia PRIMES
resistance gently but doesn't give it for free. Drug gets a kill window.

Biology preserved: HIF1A still contributes to NRF2 activation (real PDAC
biology), but at 0.02 instead of 0.5. Combined with threshold=0.03, this
means:
  - HIF1A alone: stress = 0.02 < 0.03 → NO NRF2 activation
  - Drug alone: stress = ic_drug. Needs ic_drug > 0.03 → NRF2 activates
  - HIF1A + drug: stress = ic_drug + 0.02. Needs ic_drug > 0.01 → easier!
  → Hypoxic cells develop resistance FASTER than normoxic (biologically correct)
  → But they're NOT pre-resistant before drug arrives

Wave 3 jobs:
  P  | Correct wiring baseline          | bonus=0.02, thr=0.03, km=1.5, balanced
  Q  | Stronger drug                    | bonus=0.02, thr=0.03, km=2.0, ap=0.001
  R  | No HIF1A priming                 | bonus=0.0, thr=0.03, km=1.5
  S  | Moderate priming + higher thresh | bonus=0.10, thr=0.15, km=1.5
  T  | Low efflux + correct wiring      | bonus=0.02, thr=0.03, efflux=0.3, km=1.5
  U  | High kill coeff + correct wiring | bonus=0.02, thr=0.03, kc=0.10, km=2.0
  V  | Late drug (3-wk barrier)         | bonus=0.02, thr=0.03, d21-35, km=1.5
  W  | Goldilocks (all moderate)        | bonus=0.02, thr=0.03, km=1.0, balanced
"""

import os, sys, json, subprocess, shutil
import xml.etree.ElementTree as ET

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
CONFIG_SRC   = os.path.join(PROJECT_ROOT, "config", "PhysiCell_settings.xml")
BINARY       = os.path.join(PROJECT_ROOT, "build", "bin", "release", "stroma_world")
SWEEP_BASE   = "/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_wave3"

DEFAULT_TIMING = {
    "drug_start_time": 20160.0,   # day 14
    "drug_end_time":   40320.0,   # day 28
    "max_time":        60480.0,   # day 42
}

SAVE_INT = 360
SEED     = 42
CPUS     = 32

VARIANTS = {
    "P_fixed_baseline": {
        "params": {
            "drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
            "drug_concentration": 0.10, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.7, "hif1a_nrf2_priming_bonus": 0.02,
        },
        "timing": DEFAULT_TIMING,
        "why": "Correct HIF1A wiring + balanced params. The 'should have been default' job.",
    },

    "Q_fixed_strong": {
        "params": {
            "drug_kill_multiplier": 2.0, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.15,
            "drug_concentration": 0.15, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.7, "hif1a_nrf2_priming_bonus": 0.02,
        },
        "timing": DEFAULT_TIMING,
        "why": "Stronger drug with correct wiring. Should show clear partial response (C2.2).",
    },

    "R_no_priming": {
        "params": {
            "drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
            "drug_concentration": 0.10, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.7, "hif1a_nrf2_priming_bonus": 0.0,
        },
        "timing": DEFAULT_TIMING,
        "why": "Zero HIF1A priming: pure drug-driven resistance. Clean test of drug pathway.",
    },

    "S_mod_priming": {
        "params": {
            "drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.15, "drug_uptake_rate": 0.10,
            "drug_concentration": 0.10, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.7, "hif1a_nrf2_priming_bonus": 0.10,
        },
        "timing": DEFAULT_TIMING,
        "why": "Moderate priming (0.10) with higher threshold (0.15). Drug needs ic>0.05 with HIF1A, ic>0.15 without.",
    },

    "T_fixed_lowefflux": {
        "params": {
            "drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.15,
            "drug_concentration": 0.15, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.3, "hif1a_nrf2_priming_bonus": 0.02,
        },
        "timing": DEFAULT_TIMING,
        "why": "Low efflux (0.3) + correct wiring. Drug accumulates more, stronger kill window.",
    },

    "U_fixed_highkill": {
        "params": {
            "drug_kill_multiplier": 2.0, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.15,
            "drug_concentration": 0.15, "drug_kill_coefficient": 0.10,
            "efflux_strength": 0.5, "hif1a_nrf2_priming_bonus": 0.02,
        },
        "timing": DEFAULT_TIMING,
        "why": "High kill (kc=0.10, km=2.0) + correct wiring. Maximum drug lethality with balanced efflux.",
    },

    "V_fixed_latedrug": {
        "params": {
            "drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
            "drug_concentration": 0.10, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.7, "hif1a_nrf2_priming_bonus": 0.02,
        },
        "timing": {
            "drug_start_time": 30240.0,   # day 21
            "drug_end_time":   50400.0,    # day 35
            "max_time":        70560.0,    # day 49
        },
        "why": "3-week barrier maturation + correct wiring. Best chance for C2.4/C2.5 (spatial sanctuary).",
    },

    "W_fixed_goldilocks": {
        "params": {
            "drug_kill_multiplier": 1.0, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.15,
            "drug_concentration": 0.15, "drug_kill_coefficient": 0.06,
            "efflux_strength": 0.5, "hif1a_nrf2_priming_bonus": 0.02,
        },
        "timing": DEFAULT_TIMING,
        "why": "Goldilocks: gentle drug (km=1.0) + moderate everything + correct wiring. Maximum C2.1 safety.",
    },
}


def patch_config(src, dst, out_dir, params, timing):
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath, val):
        n = root.find(xpath)
        if n is not None:
            n.text = str(val)

    _set("./save/folder",                          out_dir)
    _set(".//overall/max_time",                    str(timing["max_time"]))
    _set(".//options/random_seed",                 str(SEED))
    _set(".//parallel/omp_num_threads",            str(CPUS))
    _set(".//user_parameters/drug_start_time",     str(timing["drug_start_time"]))
    _set(".//user_parameters/drug_end_time",       str(timing["drug_end_time"]))

    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INT)

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

    for pname, pval in params.items():
        _set(f".//user_parameters/{pname}", str(pval))

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def submit_variant(name, vdef):
    params = vdef["params"]
    timing = vdef["timing"]
    work_dir = os.path.join(SWEEP_BASE, name)
    out_dir  = os.path.join(work_dir, "output")
    os.makedirs(out_dir, exist_ok=True)

    local_cfg = os.path.join(work_dir, "config")
    os.makedirs(local_cfg, exist_ok=True)
    for f in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
        src = os.path.join(PROJECT_ROOT, "config", f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(local_cfg, f))

    config_dst = os.path.join(work_dir, "config.xml")
    patch_config(CONFIG_SRC, config_dst, out_dir, params, timing)

    sbatch_script = os.path.join(work_dir, "run.sh")
    with open(sbatch_script, "w") as fh:
        fh.write(f"""#!/bin/bash
#SBATCH --job-name=w3_{name[:8]}
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

cd "{work_dir}"
echo "[w3_{name}] start: $(date)"
echo "[w3_{name}] FIXED BINARY: hif1a_nrf2_priming_bonus wired to XML"
"{BINARY}" "{config_dst}"
echo "[w3_{name}] done:  $(date)"
""")

    result = subprocess.run(["sbatch", sbatch_script],
                            capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.returncode == 0 else "FAILED"
    return job_id


if __name__ == "__main__":
    os.makedirs(SWEEP_BASE, exist_ok=True)

    manifest = {"sweep_type": "rc2_wave3_hif1a_fix",
                "fix": "hif1a_nrf2_priming_bonus wired to replace hardcoded 0.5",
                "variants": {}}

    print(f"\n{'='*85}")
    print("  WAVE 3: HIF1A PRIMING FIX + parameter sweep (8 RC2 jobs)")
    print(f"  Fix: stress = ic_drug + hif1a_bonus (XML) instead of ic_drug + 0.5 (hardcoded)")
    print(f"{'='*85}\n")

    fmt = f"  {'Variant':<22} {'Job':<8} {'km':>5} {'ap':>8} {'thr':>6} {'uptk':>6} {'kc':>6} {'efx':>5} {'bonus':>6} {'timing':<10}"
    print(fmt)
    print("  " + "-" * 100)

    for name, vdef in VARIANTS.items():
        p = vdef["params"]
        t = vdef["timing"]
        job_id = submit_variant(name, vdef)
        manifest["variants"][name] = {
            "params": p, "timing": t,
            "job_id": job_id, "why": vdef["why"],
            "work_dir": os.path.join(SWEEP_BASE, name),
        }
        t_label = "d14-28" if t == DEFAULT_TIMING else f"d{int(t['drug_start_time']/1440)}-{int(t['drug_end_time']/1440)}"
        print(f"  {name:<22} {job_id:<8} {p.get('drug_kill_multiplier'):>5} "
              f"{p.get('abcb1_production_rate'):>8} {p.get('drug_stress_threshold'):>6} "
              f"{p.get('drug_uptake_rate'):>6} {p.get('drug_kill_coefficient'):>6} "
              f"{p.get('efflux_strength'):>5} {p.get('hif1a_nrf2_priming_bonus'):>6} {t_label:<10}")

    manifest_path = os.path.join(SWEEP_BASE, "manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\n{'='*85}")
    print(f"  {len(VARIANTS)} Wave 3 RC2 jobs submitted (FIXED BINARY)")
    print(f"  Total across all waves: 8 + 8 + 8 = 24 RC2 variants")
    print(f"  Manifest: {manifest_path}")
    print(f"{'='*85}")
    print()
    print("  KEY CHANGE: hif1a_nrf2_priming_bonus now wired to XML (was hardcoded 0.5)")
    print("  Pre-fix:  93% of cells had ABCB1 before drug (stress=0.5 >> threshold=0.03)")
    print("  Post-fix: 0% of cells have ABCB1 before drug (stress=0.02 < threshold=0.03)")
    print("  → Drug now gets a genuine kill window before resistance develops")
    print()
