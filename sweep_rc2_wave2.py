#!/usr/bin/env python3
"""
Stage 2 RC2 Counter-Measurement Wave 2 — Cover EVERY remaining failure mode
=============================================================================

Wave 1 (jobs A-G) swept: km, ap, threshold, uptake, concentration, kill_coeff
Wave 2 (jobs H-O) covers everything else that could possibly go wrong:

  Job | Strategy                           | Critical new axis
  ----|------------------------------------|---------------------------------
  H   | Nuclear uptake + low efflux        | efflux=0.3 (drug can accumulate)
  I   | Ultra-gentle (guaranteed survivors) | km=0.5, kc=0.03
  J   | Max resistance + fast induction    | efflux=0.9, ap=0.003, delay=0.2
  K   | High kill coefficient              | kc=0.15 (10x effective_drug→death)
  L   | Late drug (3-wk barrier maturation)| drug d21-d35, sim d49
  M   | Short treatment (1 week only)      | drug d14-d21 only, more regrowth
  N   | Zero inheritance floor             | abcb1_inheritance_floor=0.0
  O   | Goldilocks (all moderate)          | balanced mid-range everything

RC2 failure modes covered across ALL 16 jobs (Wave 1 + Wave 2):
  C2.1 Drug eradicates:     I(ultra-gentle), M(short treat), J(max resist)
  C2.2 Drug has no effect:  H(force uptake), K(high kc), F(2.5x uptake)
  C2.3 No regrowth:         M(more regrowth time), L(mature barrier)
  C2.4 Barrier collapses:   L(thicker barrier at drug start)
  C2.5 No spatial gradient: L(mature barrier), O(balanced)
  C2.6 No resistance:       J(fast induction), N(clean inheritance)
  General balance:          O(Goldilocks), D(balanced from Wave 1)

Additional parameter axes covered:
  efflux_strength:          0.3(H), 0.5(K,N,O), 0.7(default), 0.9(J)
  efflux_induction_delay:   0.2(J), 0.3(O), 0.5(default)
  drug_kill_coefficient:    0.03(I), 0.05(def), 0.06(O), 0.08(N), 0.15(K)
  abcb1_inheritance_floor:  0.0(N), 0.6(default)
  drug timing:              d14-28(default), d21-35(L), d14-21(M)
  resistance_withdrawal_tau: 10080(default) — covered by L/M timing changes
"""

import os, sys, json, subprocess, shutil
import xml.etree.ElementTree as ET

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
CONFIG_SRC   = os.path.join(PROJECT_ROOT, "config", "PhysiCell_settings.xml")
BINARY       = os.path.join(PROJECT_ROOT, "build", "bin", "release", "stroma_world")
SWEEP_BASE   = "/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_wave2"

# Standard Stage 2 timing (overridden per-job where noted)
DEFAULT_TIMING = {
    "drug_start_time": 20160.0,   # day 14
    "drug_end_time":   40320.0,   # day 28
    "max_time":        60480.0,   # day 42
}

SAVE_INT = 360
SEED     = 42
CPUS     = 32

# ── Variant definitions ─────────────────────────────────────────────
VARIANTS = {
    "H_nuke_lowefflux": {
        "params": {
            "drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.02, "drug_uptake_rate": 0.50,
            "drug_concentration": 0.30, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.3,
        },
        "timing": DEFAULT_TIMING,
        "why": "Force drug accumulation: 5x uptake + low efflux. If drug still reads 0, pathway is structurally broken.",
    },

    "I_ultragentle": {
        "params": {
            "drug_kill_multiplier": 0.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
            "drug_concentration": 0.10, "drug_kill_coefficient": 0.03,
            "efflux_strength": 0.7,
        },
        "timing": DEFAULT_TIMING,
        "why": "Ultra-low lethality: km=0.5, kc=0.03. Guarantees C2.1 (survivors). Tests whether drug does ANYTHING at all (C2.2).",
    },

    "J_maxresist": {
        "params": {
            "drug_kill_multiplier": 2.0, "abcb1_production_rate": 0.003,
            "drug_stress_threshold": 0.01, "drug_uptake_rate": 0.10,
            "drug_concentration": 0.10, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.9, "efflux_induction_delay": 0.2,
        },
        "timing": DEFAULT_TIMING,
        "why": "Maximum resistance speed: 6x ABCB1 production, lowest threshold, fastest induction. Even strong drug can't wipe out.",
    },

    "K_highkc": {
        "params": {
            "drug_kill_multiplier": 1.0, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.02, "drug_uptake_rate": 0.15,
            "drug_concentration": 0.15, "drug_kill_coefficient": 0.15,
            "efflux_strength": 0.5,
        },
        "timing": DEFAULT_TIMING,
        "why": "3x kill coefficient: if drug enters cells but kill_coeff is bottleneck, this catches it. C2.2 safety net.",
    },

    "L_latedrug": {
        "params": {
            "drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.03, "drug_uptake_rate": 0.10,
            "drug_concentration": 0.10, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.7,
        },
        "timing": {
            "drug_start_time": 30240.0,   # day 21
            "drug_end_time":   50400.0,    # day 35
            "max_time":        70560.0,    # day 49
        },
        "why": "3-week barrier maturation before drug. Thicker ECM = better spatial gradient (C2.5) + barrier persistence (C2.4).",
    },

    "M_shorttreat": {
        "params": {
            "drug_kill_multiplier": 2.0, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.02, "drug_uptake_rate": 0.15,
            "drug_concentration": 0.15, "drug_kill_coefficient": 0.05,
            "efflux_strength": 0.7,
        },
        "timing": {
            "drug_start_time": 20160.0,   # day 14
            "drug_end_time":   30240.0,    # day 21 (1 week only!)
            "max_time":        60480.0,    # day 42 (3 weeks regrowth)
        },
        "why": "Short 1-week treatment: strong drug but limited exposure. Guarantees survivors (C2.1) + long regrowth window (C2.3).",
    },

    "N_nofloor": {
        "params": {
            "drug_kill_multiplier": 1.5, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.02, "drug_uptake_rate": 0.15,
            "drug_concentration": 0.15, "drug_kill_coefficient": 0.08,
            "efflux_strength": 0.5, "abcb1_inheritance_floor": 0.0,
        },
        "timing": DEFAULT_TIMING,
        "why": "Zero ABCB1 inheritance floor (XML default is 0.6!). Tests whether inheritance floor is masking true drug dynamics.",
    },

    "O_goldilocks": {
        "params": {
            "drug_kill_multiplier": 1.0, "abcb1_production_rate": 0.001,
            "drug_stress_threshold": 0.015, "drug_uptake_rate": 0.15,
            "drug_concentration": 0.15, "drug_kill_coefficient": 0.06,
            "efflux_strength": 0.5, "efflux_induction_delay": 0.3,
        },
        "timing": DEFAULT_TIMING,
        "why": "All parameters at moderate mid-range. Maximum balance. If nothing else passes, this is the best bet.",
    },
}


def patch_config(src, dst, out_dir, params, timing):
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath, val):
        n = root.find(xpath)
        if n is not None:
            n.text = str(val)

    # Standard Stage 2 patches
    _set("./save/folder",                          out_dir)
    _set(".//overall/max_time",                    str(timing["max_time"]))
    _set(".//options/random_seed",                 str(SEED))
    _set(".//parallel/omp_num_threads",            str(CPUS))
    _set(".//user_parameters/drug_start_time",     str(timing["drug_start_time"]))
    _set(".//user_parameters/drug_end_time",       str(timing["drug_end_time"]))

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

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def submit_variant(name, vdef):
    params = vdef["params"]
    timing = vdef["timing"]

    work_dir = os.path.join(SWEEP_BASE, name)
    out_dir  = os.path.join(work_dir, "output")
    os.makedirs(out_dir, exist_ok=True)

    # Copy auxiliary config files
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
#SBATCH --job-name=w2_{name[:8]}
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
echo "[w2_{name}] start: $(date)"
"{BINARY}" "{config_dst}"
echo "[w2_{name}] done:  $(date)"
""")

    result = subprocess.run(["sbatch", sbatch_script],
                            capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.returncode == 0 else "FAILED"
    return job_id


if __name__ == "__main__":
    os.makedirs(SWEEP_BASE, exist_ok=True)

    manifest = {"sweep_type": "rc2_wave2", "variants": {}}

    print(f"\n{'='*80}")
    print("  WAVE 2: Full failure-mode coverage (8 additional RC2 jobs)")
    print(f"{'='*80}\n")

    fmt = f"  {'Variant':<20} {'JobID':<8} {'km':>5} {'ap':>8} {'thr':>6} {'uptk':>6} {'conc':>6} {'kc':>6} {'efx':>5} {'timing':<12}"
    print(fmt)
    print("  " + "-" * 95)

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
        print(f"  {name:<20} {job_id:<8} {p.get('drug_kill_multiplier','-'):>5} "
              f"{p.get('abcb1_production_rate','-'):>8} {p.get('drug_stress_threshold','-'):>6} "
              f"{p.get('drug_uptake_rate','-'):>6} {p.get('drug_concentration','-'):>6} "
              f"{p.get('drug_kill_coefficient','-'):>6} {p.get('efflux_strength',0.7):>5} {t_label:<12}")

    manifest_path = os.path.join(SWEEP_BASE, "manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\n{'='*80}")
    print(f"  {len(VARIANTS)} Wave 2 RC2 jobs submitted")
    print(f"  Combined with Wave 1 (8 jobs) = 16 total RC2 variants")
    print(f"  Manifest: {manifest_path}")
    print(f"{'='*80}\n")

    print("  Failure mode coverage matrix:")
    print("  C2.1 (Drug eradicates):      I(ultra-gentle) M(short-treat) J(max-resist)")
    print("  C2.2 (Drug has no effect):    H(nuke-uptake) K(high-kc) + Wave1:F(force-uptake)")
    print("  C2.3 (No regrowth):           M(short-treat=3wk regrowth) L(mature-barrier)")
    print("  C2.4 (Barrier collapses):     L(3wk mature barrier at drug start)")
    print("  C2.5 (No spatial gradient):   L(thick barrier) O(goldilocks)")
    print("  C2.6 (No resistance):         J(fast-induction) N(clean-inheritance)")
    print("  General balance:              O(goldilocks) + Wave1:D(balanced)")
    print()
