#!/usr/bin/env python3
"""Mechanistic RC2 evaluator.

Produces a detailed, interpretation-first report for each RC2 run:
  1) Tumor timeline with growth rates
  2) Explained RC2 criteria breakdown
  3) Resistance dynamics (ABCB1 / NRF2)
  4) Mechanism check (selection vs loophole)
  5) Drug dynamics (intra + extra)
  6) EMT / ZEB1 state
  7) Spatial / ECM analysis with p95 radius
  8) Final verdict (true biological pass vs false pass)
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))
from python.wrapper.output_parser import OutputParser
from python.wrapper.workdir_utils import default_reality_check_dir


DEFAULT_OUT_DIR = (
    default_reality_check_dir(PROJECT_ROOT, "reality_check_2")
    / "replicate_01_seed42"
    / "output"
)
SAVE_INTERVAL = 360.0
BARRIER_MIN_PRE = 0.20
SANCTUARY_MIN_DELTA = 0.05
ABCB1_COLLAPSE_DROP = 0.50
MIN_TREATMENT_DRUG = 1e-6


def _row(matrix, labels, name):
    entry = labels.get(name)
    if entry is None:
        return None
    idx = int(entry["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _sample_field_at_positions(positions, micro_coords, field_vals):
    if positions.size == 0 or micro_coords.size == 0:
        return np.array([], dtype=float)
    out = []
    vox2 = micro_coords[:, :2] if micro_coords.shape[1] > 2 else micro_coords
    for i in range(0, positions.shape[0], 500):
        p = positions[i : i + 500, :2]
        d2 = np.sum((p[:, None, :] - vox2[None, :, :]) ** 2, axis=2)
        out.append(field_vals[np.argmin(d2, axis=1)])
    return np.concatenate(out) if out else np.array([], dtype=float)


def _fmt(v, nd=4):
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "nan"
    return f"{float(v):.{nd}f}"


def _infer_seed_label(out_dir: Path) -> str:
    m = re.search(r"seed(\d+)", out_dir.as_posix())
    return m.group(1) if m else "?"


def _load_snapshot(parser: OutputParser, out_dir: Path, target_time_min: float):
    idx = int(round(target_time_min / SAVE_INTERVAL))
    for i in range(idx, -1, -1):
        f = out_dir / f"output{i:08d}.xml"
        if f.exists() and f.stat().st_size > 0:
            try:
                return parser._read_physicell_xml(f), i
            except Exception:
                pass
    raise FileNotFoundError(f"No readable snapshot at or before t={target_time_min} (idx={idx})")


def _parse_snapshot(parser: OutputParser, snap: dict):
    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]
    micro_coords = snap["micro_coords"]
    micro_values = snap["micro_values"]

    cell_type = _row(matrix, labels, "cell_type")
    dead = _row(matrix, labels, "dead")
    death_model = _row(matrix, labels, "current_death_model")
    n_cells = matrix.shape[1]

    live_mask = np.ones(n_cells, dtype=bool)
    if dead is not None:
        live_mask &= dead <= 0.5
    if death_model is not None:
        live_mask &= np.rint(death_model).astype(int) != 100

    ctype = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1)
    tumor_mask = ctype == 0
    live_tumor = live_mask & tumor_mask
    dead_tumor = tumor_mask & ~live_mask

    pos = parser._get_positions(matrix, labels)
    tumor_pos = pos[live_tumor] if pos.size and np.any(live_tumor) else np.empty((0, 3))
    centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.shape[0] > 0 else np.zeros(3)
    if tumor_pos.shape[0] > 1:
        d = np.linalg.norm(tumor_pos - centroid, axis=1)
        radius_p95 = float(np.percentile(d, 95))
        radius_max = float(np.max(d))
    else:
        radius_p95 = 0.0
        radius_max = 0.0

    # Signals / states
    abcb1 = _row(matrix, labels, "abcb1_active")
    if abcb1 is None:
        abcb1 = _row(matrix, labels, "ABCB1")
    nrf2 = _row(matrix, labels, "nrf2_active")
    if nrf2 is None:
        nrf2 = _row(matrix, labels, "NRF2")
    zeb1 = _row(matrix, labels, "zeb1_active")
    if zeb1 is None:
        zeb1 = _row(matrix, labels, "ZEB1")
    ic_drug = _row(matrix, labels, "intracellular_drug")
    cycle_exit_rate = _row(matrix, labels, "current_cycle_phase_exit_rate")

    frac_abcb1 = float(np.mean(abcb1[live_tumor] > 0.5)) if (abcb1 is not None and np.any(live_tumor)) else math.nan
    mean_abcb1 = float(np.mean(abcb1[live_tumor])) if (abcb1 is not None and np.any(live_tumor)) else math.nan
    mean_nrf2 = float(np.mean(nrf2[live_tumor])) if (nrf2 is not None and np.any(live_tumor)) else math.nan
    frac_zeb1 = float(np.mean(zeb1[live_tumor] > 0.5)) if (zeb1 is not None and np.any(live_tumor)) else math.nan
    mean_ic_drug = float(np.mean(ic_drug[live_tumor])) if (ic_drug is not None and np.any(live_tumor)) else math.nan
    mean_cycle_rate = (
        float(np.mean(cycle_exit_rate[live_tumor]))
        if (cycle_exit_rate is not None and np.any(live_tumor))
        else math.nan
    )

    # Extracellular drug near live tumor
    ec_drug_field = micro_values.get("drug")
    mean_ec_drug = math.nan
    if ec_drug_field is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        sampled = _sample_field_at_positions(tumor_pos, micro_coords, ec_drug_field)
        if sampled.size > 0:
            mean_ec_drug = float(np.mean(sampled))

    # ECM metrics
    ecm_field = micro_values.get("ecm_density")
    peri_ecm = math.nan
    ecm_at_survivors = math.nan
    if ecm_field is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        shell = (vd >= radius_p95) & (vd <= radius_p95 + 100.0)
        if np.any(shell):
            peri_ecm = float(np.nanmean(ecm_field[shell]))
        sampled_ecm = _sample_field_at_positions(tumor_pos, micro_coords, ecm_field)
        if sampled_ecm.size > 0:
            ecm_at_survivors = float(np.mean(sampled_ecm))

    return {
        "time_min": float(snap["time"]),
        "n_live_tumor": int(np.sum(live_tumor)),
        "n_dead_tumor": int(np.sum(dead_tumor)),
        "frac_abcb1": frac_abcb1,
        "mean_abcb1": mean_abcb1,
        "mean_nrf2": mean_nrf2,
        "frac_zeb1": frac_zeb1,
        "mean_ic_drug": mean_ic_drug,
        "mean_cycle_rate": mean_cycle_rate,
        "mean_ec_drug": mean_ec_drug,
        "peri_ecm": peri_ecm,
        "ecm_at_survivors": ecm_at_survivors,
        "tumor_radius_p95": radius_p95,
        "tumor_radius_max": radius_max,
    }


def _criteria(pre, treat, post):
    reduction = 1.0 - treat["n_live_tumor"] / max(pre["n_live_tumor"], 1)
    c1 = reduction >= 0.10
    c2 = treat["n_live_tumor"] > 0

    peri_pre = pre["peri_ecm"]
    peri_tr = treat["peri_ecm"]
    barrier_valid = math.isfinite(peri_pre) and peri_pre >= BARRIER_MIN_PRE
    peri_thr = 0.80 * peri_pre if math.isfinite(peri_pre) else math.nan
    c3_pass = barrier_valid and math.isfinite(peri_tr) and peri_tr >= peri_thr
    c3_status = "INVALID" if not barrier_valid else ("PASS" if c3_pass else "FAIL")

    ecm_pre = pre["ecm_at_survivors"]
    ecm_tr = treat["ecm_at_survivors"]
    ecm_delta = ecm_tr - ecm_pre if (math.isfinite(ecm_pre) and math.isfinite(ecm_tr)) else math.nan
    c4 = math.isfinite(ecm_delta) and ecm_delta >= SANCTUARY_MIN_DELTA

    c5 = (
        math.isfinite(pre["frac_abcb1"]) and math.isfinite(treat["frac_abcb1"]) and
        treat["frac_abcb1"] > pre["frac_abcb1"]
    )
    c6 = post["n_live_tumor"] > treat["n_live_tumor"]

    return {
        "RC2-1": {"status": "PASS" if c1 else "FAIL", "hard_pass": c1, "reduction": reduction},
        "RC2-2": {"status": "PASS" if c2 else "FAIL", "hard_pass": c2, "survivors": treat["n_live_tumor"]},
        "RC2-3": {
            "status": c3_status,
            "hard_pass": c3_pass,
            "barrier_valid": barrier_valid,
            "peri_pre": peri_pre,
            "peri_treat": peri_tr,
            "peri_threshold": peri_thr,
        },
        "RC2-4": {
            "status": "PASS" if c4 else "FAIL",
            "hard_pass": c4,
            "ecm_pre": ecm_pre,
            "ecm_treat": ecm_tr,
            "ecm_delta": ecm_delta,
            "delta_threshold": SANCTUARY_MIN_DELTA,
        },
        "RC2-5": {
            "status": "PASS" if c5 else "FAIL",
            "soft_pass": c5,
            "ab_pre": pre["frac_abcb1"],
            "ab_treat": treat["frac_abcb1"],
        },
        "RC2-6": {
            "status": "PASS" if c6 else "FAIL",
            "hard_pass": c6,
            "n_treat": treat["n_live_tumor"],
            "n_post": post["n_live_tumor"],
        },
    }


def _report_one(run_name: str, out_dir: Path, timing: dict, fast_mode: bool = False):
    parser = OutputParser(out_dir)

    day_to_min = {
        14: float(timing["drug_start_time"]),
        21: 21.0 * 1440.0,
        28: float(timing["drug_end_time"]),
        31: 31.0 * 1440.0,
        35: 35.0 * 1440.0,
        42: float(timing["max_time"]),
    }
    days = [14, 21, 28, 31] if fast_mode else [14, 21, 28, 31, 35, 42]

    snaps = {}
    for d in days:
        snap, _ = _load_snapshot(parser, out_dir, day_to_min[d])
        snaps[d] = _parse_snapshot(parser, snap)
    for d in [29, 30]:
        snap, _ = _load_snapshot(parser, out_dir, d * 1440.0)
        snaps[d] = _parse_snapshot(parser, snap)

    pre = snaps[14]
    treat = snaps[28]
    post = snaps[42] if not fast_mode else None

    print("\n" + "=" * 100)
    print(f"RC2 MECHANISTIC REPORT :: {run_name}")
    print(f"Output: {out_dir}")
    print("=" * 100)

    # A) Timeline
    print("\nA) Timeline Summary")
    print("  Day   LiveTumor  DeadTumor   Δcells/day")
    prev_day = None
    prev_live = None
    for d in days:
        live = snaps[d]["n_live_tumor"]
        dead = snaps[d]["n_dead_tumor"]
        rate = math.nan if prev_day is None else (live - prev_live) / (d - prev_day)
        print(f"  {d:>3}   {live:>9}  {dead:>9}   {_fmt(rate, 2):>9}")
        prev_day = d
        prev_live = live

    shrink = treat["n_live_tumor"] < pre["n_live_tumor"]
    regrow = (post is not None) and (post["n_live_tumor"] > treat["n_live_tumor"])
    post_rate_1 = (snaps[31]["n_live_tumor"] - snaps[28]["n_live_tumor"]) / 3.0
    post_rate_2 = ((snaps[42]["n_live_tumor"] - snaps[31]["n_live_tumor"]) / 11.0) if not fast_mode else math.nan
    explosive = (not fast_mode) and regrow and post_rate_2 > max(1.0, 1.5 * post_rate_1)
    print(
        f"\n  Interpretation: treatment shrink={'YES' if shrink else 'NO'}; "
        f"withdrawal regrowth={'YES' if regrow else 'NO'}; "
        f"regrowth profile={'explosive' if explosive else ('N/A-fast-mode' if fast_mode else 'gradual/mixed')}."
    )

    # Shared strict checks
    reduction = 1.0 - treat["n_live_tumor"] / max(pre["n_live_tumor"], 1)
    c1 = reduction >= 0.10
    c2 = treat["n_live_tumor"] > 0
    peri_pre = pre["peri_ecm"]
    peri_treat = treat["peri_ecm"]
    barrier_valid = math.isfinite(peri_pre) and peri_pre >= BARRIER_MIN_PRE
    peri_thr = 0.80 * peri_pre if math.isfinite(peri_pre) else math.nan
    c3_status = "INVALID" if not barrier_valid else ("PASS" if (math.isfinite(peri_treat) and peri_treat >= peri_thr) else "FAIL")
    ecm_delta = (
        treat["ecm_at_survivors"] - pre["ecm_at_survivors"]
        if (math.isfinite(treat["ecm_at_survivors"]) and math.isfinite(pre["ecm_at_survivors"]))
        else math.nan
    )
    c4 = math.isfinite(ecm_delta) and ecm_delta >= SANCTUARY_MIN_DELTA
    c5 = math.isfinite(pre["frac_abcb1"]) and math.isfinite(treat["frac_abcb1"]) and (treat["frac_abcb1"] > pre["frac_abcb1"])

    # B) RC2 breakdown
    print("\nB) RC2 Breakdown")
    print(f"  RC2-1 Partial response: {'PASS' if c1 else 'FAIL'} | reduction={reduction:.2%} (d14 {pre['n_live_tumor']} -> d28 {treat['n_live_tumor']})")
    print(f"  RC2-2 Not eradicated: {'PASS' if c2 else 'FAIL'} | survivors@d28={treat['n_live_tumor']}")
    print(f"  RC2-3 Barrier persists: {c3_status} | peri_ecm d14={_fmt(peri_pre)} d28={_fmt(peri_treat)} need>={_fmt(peri_thr)}")
    print(f"  RC2-4 Spatial sanctuary: {'PASS' if c4 else 'FAIL'} | ecm@all_pre={_fmt(pre['ecm_at_survivors'])} ecm@survivors_d28={_fmt(treat['ecm_at_survivors'])} Δ={_fmt(ecm_delta)} need>={SANCTUARY_MIN_DELTA:.2f}")
    print(f"  RC2-5 ABCB1 emergence (soft): {'PASS' if c5 else 'FAIL'} | frac_abcb1 d14={_fmt(pre['frac_abcb1'])} d28={_fmt(treat['frac_abcb1'])}")
    if fast_mode:
        print("  RC2-6 Regrowth: N/A in --fast_mode (requires d42)")
    else:
        c6 = post["n_live_tumor"] > treat["n_live_tumor"]
        print(f"  RC2-6 Regrowth: {'PASS' if c6 else 'FAIL'} | d28={treat['n_live_tumor']} d42={post['n_live_tumor']}")

    # C) Resistance dynamics
    print("\nC) Resistance Dynamics")
    print("  Day   ABCB1+%   mean_ABCB1   mean_NRF2")
    ab_series = []
    for d in [28, 29, 30, 31]:
        s = snaps[d]
        ab_series.append(s["frac_abcb1"])
        print(f"  {d:>3}   {_fmt(100.0 * s['frac_abcb1'], 2):>7}%   {_fmt(s['mean_abcb1'], 4):>10}   {_fmt(s['mean_nrf2'], 4):>9}")
    peak = ab_series[0]
    min_post = min(ab_series[1:])
    collapse_ratio = ((peak - min_post) / max(peak, 1e-9)) if (math.isfinite(peak) and math.isfinite(min_post)) else math.nan
    resistance_flag = "COLLAPSE" if (math.isfinite(collapse_ratio) and collapse_ratio > ABCB1_COLLAPSE_DROP) else "STABLE"
    print(f"\n  collapse_ratio={_fmt(collapse_ratio, 3)} threshold={ABCB1_COLLAPSE_DROP:.2f} -> {resistance_flag}")

    # D) Drug dynamics
    print("\nD) Drug Dynamics")
    print("  Day   mean_intra_drug   mean_extra_drug")
    intra = []
    extra = []
    for d in [28, 29, 30, 31]:
        s = snaps[d]
        intra.append(s["mean_ic_drug"])
        extra.append(s["mean_ec_drug"])
        print(f"  {d:>3}      {_fmt(s['mean_ic_drug'], 6):>12}      {_fmt(s['mean_ec_drug'], 6):>12}")
    intra_clear = all(math.isfinite(intra[i]) and math.isfinite(intra[i + 1]) and intra[i + 1] <= intra[i] + 1e-12 for i in range(3))
    extra_clear = all(math.isfinite(extra[i]) and math.isfinite(extra[i + 1]) and extra[i + 1] <= extra[i] + 1e-12 for i in range(3))
    tx_intra = []
    for d in [21, 24, 27, 28]:
        snap, _ = _load_snapshot(parser, out_dir, d * 1440.0)
        tx_intra.append(_parse_snapshot(parser, snap)["mean_ic_drug"])
    treatment_intra = float(np.nanmax(tx_intra)) if tx_intra else math.nan
    drug_flag = "NO_EFFECTIVE_TREATMENT" if (not math.isfinite(treatment_intra) or treatment_intra < MIN_TREATMENT_DRUG) else "OK"
    print(
        f"\n  Interpretation: intracellular clearing={'YES' if intra_clear else 'NO/MIXED'}; "
        f"extracellular clearing={'YES' if extra_clear else 'NO/MIXED'}; "
        f"treatment intracellular peak={_fmt(treatment_intra, 6)} ({drug_flag})."
    )

    # E) EMT/ZEB1
    print("\nE) EMT / ZEB1")
    for d in [28, 31]:
        print(f"  Day {d}: ZEB1+ fraction = {_fmt(100.0 * snaps[d]['frac_zeb1'], 2)}%")
    if not fast_mode:
        print(f"  Day 42: ZEB1+ fraction = {_fmt(100.0 * snaps[42]['frac_zeb1'], 2)}%")

    # F) Spatial
    print("\nF) Spatial Analysis")
    print(f"  ECM at survivors d28: {_fmt(treat['ecm_at_survivors'])}")
    print(f"  ECM at all tumor d14: {_fmt(pre['ecm_at_survivors'])}")
    if fast_mode:
        print(f"  Tumor radius p95 d14/d28/d31: {_fmt(pre['tumor_radius_p95'],2)} / {_fmt(treat['tumor_radius_p95'],2)} / {_fmt(snaps[31]['tumor_radius_p95'],2)}")
    else:
        print(f"  Tumor radius p95 d14/d28/d42: {_fmt(pre['tumor_radius_p95'],2)} / {_fmt(treat['tumor_radius_p95'],2)} / {_fmt(post['tumor_radius_p95'],2)}")
    print(f"  Sanctuary zones: {'YES' if c4 else 'NO'}")
    print(f"  Barrier intact: {'INVALID' if c3_status == 'INVALID' else ('YES' if c3_status == 'PASS' else 'NO')}")

    # G) Mechanism interpretation
    print("\nG) Mechanism Interpretation")
    pref_surv = (math.isfinite(pre["frac_abcb1"]) and math.isfinite(treat["frac_abcb1"]) and (treat["frac_abcb1"] - pre["frac_abcb1"]) > 0.10)
    if (not fast_mode) and regrow:
        z42 = snaps[42]["frac_zeb1"]
        c28_rate = snaps[28]["mean_cycle_rate"]
        c42_rate = snaps[42]["mean_cycle_rate"]
        if math.isfinite(z42) and z42 < 0.5 and math.isfinite(c28_rate) and math.isfinite(c42_rate) and c42_rate >= 1.10 * max(c28_rate, 1e-12):
            mechanism = "epithelial_proliferation"
        elif math.isfinite(z42) and z42 >= 0.5:
            mechanism = "mesenchymal_persistence"
        else:
            mechanism = "ambiguous"
    else:
        mechanism = "ambiguous"
    print(f"  Preferential resistant survival: {'YES' if pref_surv else 'NO/WEAK'}")
    print(f"  Regrowth mechanism: {mechanism}")
    print(f"  Resistance status: {resistance_flag}")

    risk_flags = []
    if (not math.isfinite(peri_pre)) or (peri_pre < BARRIER_MIN_PRE):
        risk_flags.append("LOW_PRE_BARRIER")
        risk_flags.append("RC2_3_FAIL")
    elif c3_status == "FAIL":
        risk_flags.append("RC2_3_FAIL")
    if resistance_flag == "COLLAPSE":
        risk_flags.append("ABCB1_COLLAPSE")
    if drug_flag == "NO_EFFECTIVE_TREATMENT":
        risk_flags.append("WEAK_DRUG_EXPOSURE")
    if not c4:
        risk_flags.append("SANCTUARY_WEAK")

    # H) Final verdict
    print("\nH) FINAL VERDICT")
    if fast_mode:
        early_fail = (not c1) or (resistance_flag == "COLLAPSE")
        print("  Hard Score: N/A (fast mode)")
        print(f"  Classification: {'FALSE PASS' if early_fail else 'CONDITIONAL PASS'}")
        print(f"  Reason: {'Early rejection triggered.' if early_fail else 'Passed early rejection checks through day 31.'}")
        print(f"  Risk Flags: {', '.join(risk_flags) if risk_flags else 'NONE'}")
        return 0, ("FALSE PASS" if early_fail else "CONDITIONAL PASS")

    c6 = post["n_live_tumor"] > treat["n_live_tumor"]
    hard_pass = all([c1, c2, c3_status == "PASS", c4, c6])
    hard_score = sum([c1, c2, c3_status == "PASS", c4, c6])

    if hard_pass and resistance_flag == "STABLE" and peri_pre >= BARRIER_MIN_PRE:
        classification = "TRUE BIOLOGICAL PASS"
        reason = "Hard pass + stable resistance dynamics + valid pre-treatment barrier."
    elif hard_pass:
        classification = "CONDITIONAL PASS"
        reason = "Hard pass achieved but strict biological checks still flag risk."
    else:
        classification = "FALSE PASS"
        reason = "Hard pass not achieved under strict RC2 biological criteria."

    print(f"  Hard Score: {hard_score}/5")
    print(f"  Classification: {classification}")
    print(f"  Reason: {reason}")
    print(f"  Risk Flags: {', '.join(risk_flags) if risk_flags else 'NONE'}")
    return hard_score, classification


def _jobs_from_manifest(manifest_path: Path):
    m = json.loads(manifest_path.read_text())
    jobs = []
    for name, v in m.get("variants", {}).items():
        out = Path(v["work_dir"]) / "output"
        timing = v.get("timing", {
            "drug_start_time": 20160.0,
            "drug_end_time": 40320.0,
            "max_time": 60480.0,
        })
        jobs.append((name, out, timing))
    return jobs


def main():
    ap = argparse.ArgumentParser(description="Mechanistic RC2 evaluator")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Single RC2 output dir (contains output*.xml)")
    ap.add_argument("--manifest", type=Path, default=None,
                    help="Manifest JSON from a sweep launcher")
    ap.add_argument("--all-known", action="store_true",
                    help="Evaluate baseline + wave1 + wave2 + wave3 manifests")
    ap.add_argument("--fast_mode", action="store_true",
                    help="Early-reject mode: evaluate only through day 31")
    args = ap.parse_args()

    jobs = []

    if args.all_known:
        jobs.append((
            "baseline_10620",
            (PROJECT_ROOT / "build" / "rc2_full_seed42" / "replicate_01_seed42" / "output").resolve(),
            {"drug_start_time": 20160.0, "drug_end_time": 40320.0, "max_time": 60480.0},
        ))
        for mp in [
            Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_countermeasures/manifest.json"),
            Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_wave2/manifest.json"),
            Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_wave3/manifest.json"),
        ]:
            if mp.exists():
                jobs.extend(_jobs_from_manifest(mp))
    elif args.manifest is not None:
        jobs = _jobs_from_manifest(args.manifest.resolve())
    else:
        out_dir = args.out_dir.resolve() if args.out_dir is not None else DEFAULT_OUT_DIR.resolve()
        jobs = [(
            f"seed_{_infer_seed_label(out_dir)}",
            out_dir,
            {"drug_start_time": 20160.0, "drug_end_time": 40320.0, "max_time": 60480.0},
        )]

    if not jobs:
        print("No jobs found to evaluate.")
        return 1

    summary = []
    for name, out_dir, timing in jobs:
        if not out_dir.exists():
            print(f"\n[SKIP] {name}: output dir missing: {out_dir}")
            continue
        try:
            hard_score, classification = _report_one(name, out_dir, timing, fast_mode=args.fast_mode)
            summary.append((name, hard_score, classification))
        except Exception as e:
            print(f"\n[ERROR] {name}: {e}")

    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)
    for name, hard, cls in summary:
        hard_txt = "N/A" if args.fast_mode else f"{hard}/5"
        print(f"  {name:<28} hard={hard_txt}   biological={cls}")

    ranked = {"TRUE BIOLOGICAL PASS": [], "CONDITIONAL PASS": [], "FALSE PASS": []}
    for name, hard, cls in summary:
        ranked.setdefault(cls, []).append((name, hard))

    print("\nRANKED SUMMARY")
    if args.fast_mode:
        print(f"TRUE: {', '.join(n for n, _ in ranked.get('TRUE BIOLOGICAL PASS', [])) or 'none'}")
        print(f"CONDITIONAL: {', '.join(n for n, _ in ranked.get('CONDITIONAL PASS', [])) or 'none'}")
        print(f"FALSE: {', '.join(n for n, _ in ranked.get('FALSE PASS', [])) or 'none'}")
    else:
        print(f"TRUE: {', '.join(f'{n}({h}/5)' for n, h in ranked.get('TRUE BIOLOGICAL PASS', [])) or 'none'}")
        print(f"CONDITIONAL: {', '.join(f'{n}({h}/5)' for n, h in ranked.get('CONDITIONAL PASS', [])) or 'none'}")
        print(f"FALSE: {', '.join(f'{n}({h}/5)' for n, h in ranked.get('FALSE PASS', [])) or 'none'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
