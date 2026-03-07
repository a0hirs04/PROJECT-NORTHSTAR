#!/usr/bin/env python3
"""
RC2 Full Diagnostic — Evaluate 42-day run for seed 42.
Reports metrics at day 14, 28, 42 and evaluates RC2-1 through RC2-6.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import OutputParser

OUT_DIR = PROJECT_ROOT / "build" / "rc2_full_seed42" / "replicate_01_seed42" / "output"

CHECKPOINTS = {
    "day 14 (pre-drug)":   20160.0,
    "day 28 (treat-end)":  40320.0,
    "day 42 (post-withdrawal)": 60480.0,
}

STROMAL_TO_TUMOR_LABEL = {
    "acta2_active":   "ZEB1",
    "gli1_active":    "CDH1",
    "hif1a_active":   "TGFB1_expr",
    "mechanical_pressure": "ABCB1",
}


def _row(matrix, labels, name):
    e = labels.get(name)
    if e is None:
        return None
    idx = int(e["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _nearest_xml(parser, xmls, target_time):
    best_xml = None
    best_dt = float("inf")
    best_time = 0.0
    for xml in xmls:
        snap = parser._read_physicell_xml(xml)
        dt = abs(float(snap["time"]) - target_time)
        if dt < best_dt:
            best_dt = dt
            best_xml = xml
            best_time = float(snap["time"])
    return best_xml, best_time


def _sample_nearest(points, voxel_coords, values):
    if points.size == 0 or voxel_coords.size == 0:
        return np.array([], dtype=float)
    out = np.empty(points.shape[0], dtype=float)
    chunk = 500
    for i in range(0, points.shape[0], chunk):
        p = points[i : i + chunk]
        d2 = np.sum((p[:, None, :] - voxel_coords[None, :, :]) ** 2, axis=2)
        out[i : i + chunk] = values[np.argmin(d2, axis=1)]
    return out


def analyze_snapshot(parser, xml_path):
    snap = parser._read_physicell_xml(xml_path)
    t = snap["time"]
    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]
    micro_coords = snap["micro_coords"]
    micro_values = snap["micro_values"]

    cell_type = _row(matrix, labels, "cell_type")
    dead = _row(matrix, labels, "dead")
    death_model = _row(matrix, labels, "current_death_model")

    n_cells = matrix.shape[1]
    ctype_int = np.rint(cell_type).astype(int) if cell_type is not None else np.zeros(n_cells, int)

    tumor_mask = ctype_int == 0
    stroma_mask = ctype_int == 1
    dead_mask = (dead > 0.5) if dead is not None else np.zeros(n_cells, bool)
    apoptotic_mask = (np.rint(death_model).astype(int) == 100) if death_model is not None else np.zeros(n_cells, bool)
    live_tumor = tumor_mask & ~dead_mask & ~apoptotic_mask
    dead_tumor = tumor_mask & dead_mask

    n_live = int(np.sum(live_tumor))
    n_dead = int(np.sum(dead_tumor))
    n_stroma = int(np.sum(stroma_mask & ~dead_mask))

    # Positions
    pos_row = labels.get("position")
    if pos_row:
        pidx = int(pos_row["index"])
        positions = matrix[pidx : pidx + 3, :].T
    else:
        positions = np.empty((0, 3))

    tumor_pos = positions[live_tumor] if live_tumor.any() else np.empty((0, 3))
    centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.shape[0] > 0 else np.zeros(3)
    tumor_radius = float(np.max(np.linalg.norm(tumor_pos - centroid, axis=1))) if tumor_pos.shape[0] > 1 else 0.0

    m = {
        "time": float(t),
        "n_live": n_live,
        "n_dead": n_dead,
        "n_stroma": n_stroma,
        "tumor_radius": tumor_radius,
    }

    # Peritumoral ECM
    ecm = micro_values.get("ecm_density")
    if ecm is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        shell = (vd >= tumor_radius) & (vd <= tumor_radius + 100.0)
        m["peri_ecm"] = float(np.nanmean(ecm[shell])) if shell.any() else math.nan
    else:
        m["peri_ecm"] = math.nan

    # ECM at survivor locations
    if ecm is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        sampled_ecm = _sample_nearest(tumor_pos, micro_coords, ecm)
        m["ecm_at_survivors"] = float(np.nanmean(sampled_ecm))
    else:
        m["ecm_at_survivors"] = math.nan

    # ECM at all tumor (for pre-drug baseline)
    all_tumor_mask = tumor_mask & ~dead_mask
    all_tumor_pos = positions[all_tumor_mask] if all_tumor_mask.any() else np.empty((0, 3))
    if ecm is not None and micro_coords.size > 0 and all_tumor_pos.shape[0] > 0:
        sampled_ecm_all = _sample_nearest(all_tumor_pos, micro_coords, ecm)
        m["ecm_at_all_tumor"] = float(np.nanmean(sampled_ecm_all))
    else:
        m["ecm_at_all_tumor"] = math.nan

    # Death rate
    death_rates = _row(matrix, labels, "death_rates")
    m["mean_death_rate"] = float(np.nanmean(death_rates[live_tumor])) if (death_rates is not None and n_live > 0) else math.nan

    # ZEB1
    zeb1 = _row(matrix, labels, "zeb1_active")
    m["frac_zeb1"] = float(np.mean(zeb1[live_tumor] > 0.5)) if (zeb1 is not None and n_live > 0) else math.nan

    # HIF1A
    hif1a = _row(matrix, labels, "hif1a_active")
    m["frac_hif1a"] = float(np.mean(hif1a[live_tumor] > 0.5)) if (hif1a is not None and n_live > 0) else math.nan

    # ABCB1
    abcb1 = _row(matrix, labels, "abcb1_active")
    m["frac_abcb1"] = float(np.mean(abcb1[live_tumor] > 0.5)) if (abcb1 is not None and n_live > 0) else math.nan

    # NRF2
    nrf2 = _row(matrix, labels, "nrf2_active")
    m["frac_nrf2"] = float(np.mean(nrf2[live_tumor] > 0.5)) if (nrf2 is not None and n_live > 0) else math.nan

    # Drug fields
    drug_field = micro_values.get("drug")
    if drug_field is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        sampled_drug = _sample_nearest(tumor_pos, micro_coords, drug_field)
        m["drug_at_tumor"] = float(np.nanmean(sampled_drug))
        m["drug_at_tumor_max"] = float(np.nanmax(sampled_drug))
    else:
        m["drug_at_tumor"] = math.nan
        m["drug_at_tumor_max"] = math.nan

    if drug_field is not None and micro_coords.size > 0 and tumor_radius > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        peri = (vd >= tumor_radius - 50) & (vd <= tumor_radius + 100)
        core = vd <= tumor_radius * 0.3
        m["drug_periphery"] = float(np.nanmean(drug_field[peri])) if peri.any() else math.nan
        m["drug_core"] = float(np.nanmean(drug_field[core])) if core.any() else math.nan
    else:
        m["drug_periphery"] = math.nan
        m["drug_core"] = math.nan

    # Intracellular drug
    ic_drug = _row(matrix, labels, "intracellular_drug")
    if ic_drug is not None and n_live > 0:
        vals = ic_drug[live_tumor]
        m["ic_drug_mean"] = float(np.nanmean(vals))
        m["ic_drug_max"] = float(np.nanmax(vals))
        m["frac_ic_drug_gt0"] = float(np.mean(vals > 0.001))
    else:
        m["ic_drug_mean"] = math.nan
        m["ic_drug_max"] = math.nan
        m["frac_ic_drug_gt0"] = math.nan

    return m


def main():
    if not OUT_DIR.exists():
        print(f"ERROR: output not found: {OUT_DIR}")
        return 1

    parser = OutputParser(OUT_DIR)
    xmls = sorted(OUT_DIR.glob("output*.xml"))
    if not xmls:
        print("ERROR: no snapshots found")
        return 1

    print("=" * 72)
    print("  RC2 FULL DIAGNOSTIC — seed 42, 42-day schedule")
    print("  emt_death_increase=0.0001, dt_correct_drug_impedance_only=1")
    print("=" * 72)
    print(f"\n  Found {len(xmls)} snapshots\n")

    # Collect metrics at each checkpoint
    metrics = {}
    for label, target_t in CHECKPOINTS.items():
        xml, actual_t = _nearest_xml(parser, xmls, target_t)
        m = analyze_snapshot(parser, xml)
        metrics[label] = m

        print("  === %s (t=%.0f min, file=%s) ===" % (label, actual_t, xml.name))
        print("  Live tumor:        %d" % m["n_live"])
        print("  Dead tumor:        %d" % m["n_dead"])
        print("  Stroma:            %d" % m["n_stroma"])
        print("  Tumor radius:      %.1f um" % m["tumor_radius"])
        print("  Peri-ECM:          %.4f" % m["peri_ecm"])
        print("  ECM@survivors:     %.4f" % m["ecm_at_survivors"])
        print("  ECM@all_tumor:     %.4f" % m["ecm_at_all_tumor"])
        print("  Mean death rate:   %.6f" % m["mean_death_rate"])
        print("  ZEB1+ fraction:    %.4f" % m["frac_zeb1"])
        print("  HIF1A+ fraction:   %.4f" % m["frac_hif1a"])
        print("  ABCB1+ fraction:   %.4f" % m["frac_abcb1"])
        print("  NRF2+ fraction:    %.4f" % m["frac_nrf2"])
        print("  Drug periphery:    %.6f" % m.get("drug_periphery", float("nan")))
        print("  Drug core:         %.6f" % m.get("drug_core", float("nan")))
        print("  Drug@tumor mean:   %.6f" % m.get("drug_at_tumor", float("nan")))
        print("  Drug@tumor max:    %.6f" % m.get("drug_at_tumor_max", float("nan")))
        print("  IC drug mean:      %.6f" % m.get("ic_drug_mean", float("nan")))
        print("  IC drug max:       %.6f" % m.get("ic_drug_max", float("nan")))
        print("  Frac IC drug>0:    %.4f" % m.get("frac_ic_drug_gt0", float("nan")))
        print()

    # Time-series of tumor count (every 3 days)
    print("  --- Tumor Count Time-Series ---")
    print("  %8s  %8s  %8s  %s" % ("Day", "Live", "Dead", "Phase"))
    print("  " + "-" * 50)
    for day in [0, 3, 7, 10, 14, 17, 21, 24, 28, 31, 35, 38, 42]:
        target_t = day * 1440.0
        xml, actual_t = _nearest_xml(parser, xmls, target_t)
        snap = parser._read_physicell_xml(xml)
        mat = snap["cell_matrix"]
        lab = snap["label_name_map"]
        ct = _row(mat, lab, "cell_type")
        dd = _row(mat, lab, "dead")
        dm = _row(mat, lab, "current_death_model")
        n = mat.shape[1]
        cti = np.rint(ct).astype(int) if ct is not None else np.zeros(n, int)
        tumor = cti == 0
        dead_m = (dd > 0.5) if dd is not None else np.zeros(n, bool)
        apop = (np.rint(dm).astype(int) == 100) if dm is not None else np.zeros(n, bool)
        live = int(np.sum(tumor & ~dead_m & ~apop))
        dead_count = int(np.sum(tumor & dead_m))
        phase = "BARRIER" if day <= 14 else ("DRUG" if day <= 28 else "REGROWTH")
        print("  %8d  %8d  %8d  %s" % (day, live, dead_count, phase))
    print()

    # RC2 Criteria Evaluation
    pre = metrics["day 14 (pre-drug)"]
    treat = metrics["day 28 (treat-end)"]
    post = metrics["day 42 (post-withdrawal)"]

    print("  " + "=" * 60)
    print("  RC2 CRITERIA EVALUATION")
    print("  " + "=" * 60)

    # RC2-1: Partial response (>=10% reduction)
    if pre["n_live"] > 0:
        reduction = 1.0 - treat["n_live"] / pre["n_live"]
        rc2_1 = treat["n_live"] <= 0.90 * pre["n_live"]
        print("  [%s] RC2-1 Partial response: tumor %d->%d (%.1f%% reduction, need >=10%%)" %
              ("PASS" if rc2_1 else "FAIL", pre["n_live"], treat["n_live"], reduction * 100))
    else:
        rc2_1 = False
        print("  [FAIL] RC2-1 Partial response: no tumor at pre-drug")

    # RC2-2: Not eradicated
    rc2_2 = treat["n_live"] > 0
    print("  [%s] RC2-2 Not eradicated: %d live tumor at treat-end" %
          ("PASS" if rc2_2 else "FAIL", treat["n_live"]))

    # RC2-3: Barrier persists (ECM >= 80% of pre)
    if not math.isnan(pre["peri_ecm"]) and pre["peri_ecm"] > 0:
        barrier_frac = treat["peri_ecm"] / pre["peri_ecm"]
        rc2_3 = treat["peri_ecm"] >= 0.80 * pre["peri_ecm"]
        print("  [%s] RC2-3 Barrier persists: ECM %.4f -> %.4f (%.1f%% of pre, need >=80%%)" %
              ("PASS" if rc2_3 else "FAIL", pre["peri_ecm"], treat["peri_ecm"], barrier_frac * 100))
    else:
        rc2_3 = False
        print("  [FAIL] RC2-3 Barrier persists: peri-ECM unavailable")

    # RC2-4: Spatial sanctuary (ECM@survivors > ECM@all_tumor_pre)
    rc2_4 = treat["ecm_at_survivors"] > pre["ecm_at_all_tumor"]
    print("  [%s] RC2-4 Spatial sanctuary: ECM@survivors=%.4f > ECM@all_pre=%.4f" %
          ("PASS" if rc2_4 else "FAIL", treat["ecm_at_survivors"], pre["ecm_at_all_tumor"]))

    # RC2-5: ABCB1 emerges (SOFT — diagnostic only)
    abcb1_pre = pre["frac_abcb1"] if not math.isnan(pre["frac_abcb1"]) else 0
    abcb1_treat = treat["frac_abcb1"] if not math.isnan(treat["frac_abcb1"]) else 0
    rc2_5 = abcb1_treat > abcb1_pre
    print("  [%s] RC2-5 ABCB1 emerges [SOFT]: %.4f -> %.4f" %
          ("PASS" if rc2_5 else "FAIL", abcb1_pre, abcb1_treat))

    # RC2-6: Regrowth after withdrawal
    if treat["n_live"] > 0:
        rc2_6 = post["n_live"] > treat["n_live"]
        print("  [%s] RC2-6 Regrowth: tumor %d -> %d (treat-end -> post)" %
              ("PASS" if rc2_6 else "FAIL", treat["n_live"], post["n_live"]))
    else:
        rc2_6 = False
        print("  [FAIL] RC2-6 Regrowth: no survivors to regrow")

    print()
    hard_pass = sum([rc2_1, rc2_2, rc2_3, rc2_4, rc2_6])
    print("  HARD criteria: %d/5 PASS" % hard_pass)
    print("  SOFT (ABCB1):  %s" % ("PASS" if rc2_5 else "FAIL"))
    print()

    # Summary
    print("  " + "=" * 60)
    print("  SUMMARY")
    print("  " + "=" * 60)
    print("  partial_response:  %s" % rc2_1)
    print("  not_eradicated:    %s" % rc2_2)
    print("  barrier_persists:  %s" % rc2_3)
    print("  sanctuary:         %s" % rc2_4)
    print("  regrowth:          %s" % rc2_6)
    print("  abcb1_emerges:     %s (soft)" % rc2_5)
    overall = "GREEN" if hard_pass >= 4 else ("YELLOW" if hard_pass >= 3 else "RED")
    print("  OVERALL: %s (%d/5 hard)" % (overall, hard_pass))
    print("  " + "=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
