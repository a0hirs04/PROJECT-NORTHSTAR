#!/usr/bin/env python3
"""
Step 1 Diagnostic — Reads output from step1_validation runs and reports metrics.

Run A (RC1 baseline): tumor count at days 3, 7, 10, 14 + detailed metrics at days 10, 14
Run B (RC2 probe): drug penetration at day 16

Also runs the RC1 8-criteria evaluation on Run A output.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import OutputParser

RC1_DIR = PROJECT_ROOT / "build" / "step1_validation" / "rc1_baseline_seed42" / "output"
RC2_DIR = PROJECT_ROOT / "build" / "step1_validation" / "rc2_probe_seed42" / "output"

# Target times in minutes
RC1_CHECKPOINTS = {
    "day 3":  4320.0,
    "day 7":  10080.0,
    "day 10": 14400.0,
    "day 14": 20160.0,
}
RC2_CHECKPOINT = {"day 16": 23040.0}

# Stromal custom-data offset mapping (same as RC1/RC2 runners)
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
    """Find XML snapshot closest to target_time."""
    best_xml = None
    best_dt = float("inf")
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


def analyze_snapshot(parser, xml_path, detailed=False):
    """Analyze a single snapshot. Returns dict of metrics."""
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
    n_total_tumor = int(np.sum(tumor_mask))
    n_stroma = int(np.sum(stroma_mask & ~dead_mask))

    result = {
        "time": float(t),
        "n_live": n_live,
        "n_dead": n_dead,
        "n_total_tumor": n_total_tumor,
        "n_stroma": n_stroma,
    }

    if not detailed:
        return result

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

    # Peritumoral ECM
    ecm = micro_values.get("ecm_density")
    peri_ecm = math.nan
    if ecm is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        shell = (vd >= tumor_radius) & (vd <= tumor_radius + 100.0)
        if shell.any():
            peri_ecm = float(np.nanmean(ecm[shell]))
    result["peri_ecm"] = peri_ecm
    result["tumor_radius"] = tumor_radius

    # Cycle rate & death rate from custom_data
    cycle_rate = _row(matrix, labels, "cycle_rate")
    if cycle_rate is not None and n_live > 0:
        result["mean_cycle_rate"] = float(np.nanmean(cycle_rate[live_tumor]))
    else:
        result["mean_cycle_rate"] = math.nan

    death_rates = _row(matrix, labels, "death_rates")
    if death_rates is not None and n_live > 0:
        result["mean_death_rate"] = float(np.nanmean(death_rates[live_tumor]))
    else:
        result["mean_death_rate"] = math.nan

    # ZEB1 (EMT marker)
    zeb1 = _row(matrix, labels, "zeb1_active")
    if zeb1 is not None and n_live > 0:
        result["frac_zeb1"] = float(np.mean(zeb1[live_tumor] > 0.5))
    else:
        result["frac_zeb1"] = math.nan

    # HIF1A
    hif1a = _row(matrix, labels, "hif1a_active")
    if hif1a is not None and n_live > 0:
        result["frac_hif1a"] = float(np.mean(hif1a[live_tumor] > 0.5))
    else:
        result["frac_hif1a"] = math.nan

    # Contact inhibition (mechanical_pressure > 0.5)
    mech_pressure = _row(matrix, labels, "mechanical_pressure")
    if mech_pressure is not None and n_live > 0:
        result["frac_contact_inhibited"] = float(np.mean(mech_pressure[live_tumor] > 0.5))
    else:
        result["frac_contact_inhibited"] = math.nan

    # ABCB1 / NRF2
    abcb1 = _row(matrix, labels, "abcb1_active")
    nrf2 = _row(matrix, labels, "nrf2_active")
    result["frac_abcb1"] = float(np.mean(abcb1[live_tumor] > 0.5)) if (abcb1 is not None and n_live > 0) else math.nan
    result["frac_nrf2"] = float(np.mean(nrf2[live_tumor] > 0.5)) if (nrf2 is not None and n_live > 0) else math.nan

    # Drug fields
    drug_field = micro_values.get("drug")
    if drug_field is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        sampled = _sample_nearest(tumor_pos, micro_coords, drug_field)
        result["drug_at_tumor"] = float(np.nanmean(sampled))
        result["drug_at_tumor_max"] = float(np.nanmax(sampled))
    else:
        result["drug_at_tumor"] = math.nan
        result["drug_at_tumor_max"] = math.nan

    # Drug at periphery vs core
    if drug_field is not None and micro_coords.size > 0 and tumor_radius > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        peri_mask = (vd >= tumor_radius - 50) & (vd <= tumor_radius + 100)
        core_mask = vd <= tumor_radius * 0.3
        if peri_mask.any():
            result["drug_periphery"] = float(np.nanmean(drug_field[peri_mask]))
        else:
            result["drug_periphery"] = math.nan
        if core_mask.any():
            result["drug_core"] = float(np.nanmean(drug_field[core_mask]))
        else:
            result["drug_core"] = math.nan
    else:
        result["drug_periphery"] = math.nan
        result["drug_core"] = math.nan

    # Intracellular drug
    ic_drug = _row(matrix, labels, "intracellular_drug")
    if ic_drug is not None and n_live > 0:
        vals = ic_drug[live_tumor]
        result["ic_drug_mean"] = float(np.nanmean(vals))
        result["ic_drug_max"] = float(np.nanmax(vals))
    else:
        result["ic_drug_mean"] = math.nan
        result["ic_drug_max"] = math.nan

    # CAF activation (stromal ACTA2)
    acta2_label = STROMAL_TO_TUMOR_LABEL.get("acta2_active")
    if acta2_label:
        acta2_row = _row(matrix, labels, acta2_label)
        if acta2_row is not None:
            stroma_live = stroma_mask & ~dead_mask
            if stroma_live.any():
                result["n_caf_active"] = int(np.sum(acta2_row[stroma_live] > 0.5))

    return result


def rc1_criteria(metrics_d14):
    """Evaluate RC1 8-criteria against day-14 (or day-21) metrics."""
    criteria = {}
    # RC1-1: Barrier self-assembles (peri ECM > 0.5)
    criteria["RC1-1 Barrier"] = metrics_d14.get("peri_ecm", 0) > 0.5
    # RC1-5: Tumor-stroma ratio (tumor < 50%)
    n_t = metrics_d14.get("n_live", 0)
    n_s = metrics_d14.get("n_stroma", 0)
    total = n_t + n_s
    criteria["RC1-5 TSR"] = (n_t / total < 0.5) if total > 0 else False
    # RC1-7: CAFs present
    criteria["RC1-7 CAFs"] = metrics_d14.get("n_caf_active", 0) > 0
    # RC1-6: No instability (ECM <= 1.0, live cells > 0)
    criteria["RC1-6 Stability"] = n_t > 0 and metrics_d14.get("peri_ecm", 0) <= 1.0
    return criteria


def main():
    print("=" * 72)
    print("  STEP 1 DIAGNOSTIC: emt_death_increase = 0.0001")
    print("=" * 72)

    # ---- Run A: RC1 baseline ----
    print("\n--- RUN A: RC1 Baseline (14 days, no drug, seed 42) ---\n")

    if not RC1_DIR.exists():
        print(f"  ERROR: RC1 output not found: {RC1_DIR}")
        print("  Have the SLURM jobs finished? Check with sacct.")
        return 1

    parser_a = OutputParser(RC1_DIR)
    xmls_a = sorted(RC1_DIR.glob("output*.xml"))
    if not xmls_a:
        print("  ERROR: No snapshot XMLs found in RC1 output")
        return 1

    print(f"  Found {len(xmls_a)} snapshots")

    # Time-series summary (days 3, 7, 10, 14)
    print("\n  Time-series (live tumor count):")
    print("  %8s  %8s  %8s  %8s" % ("Day", "Live", "Dead", "Stroma"))
    print("  " + "-" * 40)

    for label, target_t in RC1_CHECKPOINTS.items():
        xml, actual_t = _nearest_xml(parser_a, xmls_a, target_t)
        m = analyze_snapshot(parser_a, xml, detailed=False)
        print("  %8s  %8d  %8d  %8d" % (label, m["n_live"], m["n_dead"], m["n_stroma"]))

    # Detailed metrics at days 10 and 14
    for label in ["day 10", "day 14"]:
        target_t = RC1_CHECKPOINTS[label]
        xml, actual_t = _nearest_xml(parser_a, xmls_a, target_t)
        m = analyze_snapshot(parser_a, xml, detailed=True)
        print("\n  === %s (t=%.0f min, file=%s) ===" % (label, actual_t, xml.name))
        print("  Live tumor:        %d" % m["n_live"])
        print("  Dead tumor:        %d" % m["n_dead"])
        print("  Stroma:            %d" % m["n_stroma"])
        print("  Tumor radius:      %.1f um" % m.get("tumor_radius", 0))
        print("  Peri-ECM (mean):   %.4f" % m.get("peri_ecm", float("nan")))
        print("  Mean cycle rate:   %.6f 1/min" % m.get("mean_cycle_rate", float("nan")))
        print("  Mean death rate:   %.6f 1/min" % m.get("mean_death_rate", float("nan")))
        print("  ZEB1+ fraction:    %.4f" % m.get("frac_zeb1", float("nan")))
        print("  HIF1A+ fraction:   %.4f" % m.get("frac_hif1a", float("nan")))
        print("  Contact inhibited: %.4f" % m.get("frac_contact_inhibited", float("nan")))
        print("  ABCB1+ fraction:   %.4f" % m.get("frac_abcb1", float("nan")))
        print("  NRF2+ fraction:    %.4f" % m.get("frac_nrf2", float("nan")))
        print("  CAFs active:       %d" % m.get("n_caf_active", 0))

    # RC1 criteria check
    xml_14, _ = _nearest_xml(parser_a, xmls_a, 20160.0)
    m14 = analyze_snapshot(parser_a, xml_14, detailed=True)
    criteria = rc1_criteria(m14)
    print("\n  --- RC1 Criteria (subset, from day 14) ---")
    for name, passed in criteria.items():
        tag = "PASS" if passed else "FAIL"
        print("    [%s] %s" % (tag, name))

    # ---- Run B: RC2 probe ----
    print("\n\n--- RUN B: RC2 Probe (16 days, drug ON at day 14, seed 42) ---\n")

    if not RC2_DIR.exists():
        print(f"  ERROR: RC2 output not found: {RC2_DIR}")
        print("  Have the SLURM jobs finished? Check with sacct.")
        return 1

    parser_b = OutputParser(RC2_DIR)
    xmls_b = sorted(RC2_DIR.glob("output*.xml"))
    if not xmls_b:
        print("  ERROR: No snapshot XMLs found in RC2 output")
        return 1

    print(f"  Found {len(xmls_b)} snapshots")

    # Day 14 (pre-drug) and day 16 (post-drug)
    for label, target_t in [("day 14 (pre-drug)", 20160.0), ("day 16 (post-drug)", 23040.0)]:
        xml, actual_t = _nearest_xml(parser_b, xmls_b, target_t)
        m = analyze_snapshot(parser_b, xml, detailed=True)
        print("\n  === %s (t=%.0f min, file=%s) ===" % (label, actual_t, xml.name))
        print("  Live tumor:        %d" % m["n_live"])
        print("  Peri-ECM (mean):   %.4f" % m.get("peri_ecm", float("nan")))
        print("  Drug at tumor:     %.6f (mean) / %.6f (max)" % (
            m.get("drug_at_tumor", float("nan")), m.get("drug_at_tumor_max", float("nan"))))
        print("  Drug periphery:    %.6f" % m.get("drug_periphery", float("nan")))
        print("  Drug core:         %.6f" % m.get("drug_core", float("nan")))
        print("  IC drug (mean):    %.6f" % m.get("ic_drug_mean", float("nan")))
        print("  IC drug (max):     %.6f" % m.get("ic_drug_max", float("nan")))

    print("\n" + "=" * 72)
    print("  STEP 1 DIAGNOSTIC COMPLETE")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
