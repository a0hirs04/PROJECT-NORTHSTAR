#!/usr/bin/env python3
"""
evaluate_resistance_sweep.py — Analyze results from sweep_resistance.py.

Reads the 3x3 sweep grid, evaluates Stage 1 gates, prints a truth table,
and recommends the best parameter combination.

Usage:
    python3 evaluate_resistance_sweep.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))
from python.wrapper.output_parser import OutputParser

SWEEP_DIR = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/resistance_sweep")
MANIFEST = SWEEP_DIR / "manifest.json"

DAY_MIN = 1440.0
CHECK_DAYS = [28, 29, 30, 31]

# Stage 1 gates
ABCB1_MEMORY_RATIO = 0.50   # d29-31 ABCB1 >= 50% of d28
SURVIVAL_RATIO = 0.70        # d31 n_live >= 70% of d28


def _row(matrix, labels, name: str):
    entry = labels.get(name)
    if entry is None:
        return None
    idx = int(entry["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _nearest_snapshot(parser: OutputParser, xmls: list[Path], target_time: float):
    best_xml = xmls[0]
    best_t = 0.0
    best_dt = float("inf")
    for xml in xmls:
        snap = parser._read_physicell_xml(xml)
        t = float(snap["time"])
        dt = abs(t - target_time)
        if dt < best_dt:
            best_dt = dt
            best_xml = xml
            best_t = t
    return best_xml, best_t


def _metrics_for_snapshot(parser: OutputParser, xml_path: Path) -> dict[str, float]:
    snap = parser._read_physicell_xml(xml_path)
    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]

    cell_type = _row(matrix, labels, "cell_type")
    dead = _row(matrix, labels, "dead")
    death_model = _row(matrix, labels, "current_death_model")
    abcb1 = _row(matrix, labels, "abcb1_active")
    ic_drug = _row(matrix, labels, "intracellular_drug")

    n_cells = matrix.shape[1]
    ctype = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1)
    tumor_mask = ctype == 0
    dead_mask = (dead > 0.5) if dead is not None else np.zeros(n_cells, dtype=bool)
    apoptotic_mask = (
        (np.rint(death_model).astype(int) == 100)
        if death_model is not None
        else np.zeros(n_cells, dtype=bool)
    )
    live_tumor = tumor_mask & ~dead_mask & ~apoptotic_mask

    n_live = int(np.sum(live_tumor))
    if n_live == 0:
        return {"n_live": 0.0, "abcb1_frac": float("nan"), "ic_mean": float("nan")}

    return {
        "n_live": float(n_live),
        "abcb1_frac": float(np.mean(abcb1[live_tumor] > 0.5)) if abcb1 is not None else float("nan"),
        "ic_mean": float(np.mean(ic_drug[live_tumor])) if ic_drug is not None else float("nan"),
    }


def _fmt(x: float) -> str:
    if isinstance(x, float) and np.isnan(x):
        return "  nan"
    return f"{x:.4f}"


def evaluate_variant(work_dir: Path) -> dict:
    """Evaluate a single variant. Returns dict with metrics and pass/fail."""
    out_dir = work_dir / "output"
    if not out_dir.exists():
        return {"status": "NO_OUTPUT", "rows": {}, "pass": False}

    xmls = sorted(out_dir.glob("output*.xml"))
    if len(xmls) < 5:
        return {"status": "INCOMPLETE", "rows": {}, "pass": False}

    parser = OutputParser(out_dir)
    rows = {}
    for day in CHECK_DAYS:
        t_target = day * DAY_MIN
        xml, t_actual = _nearest_snapshot(parser, xmls, t_target)
        m = _metrics_for_snapshot(parser, xml)
        m["t_actual"] = t_actual
        rows[day] = m

    # Gate 1: ABCB1 memory
    abcb1_28 = rows[28]["abcb1_frac"]
    pass_memory = True
    if np.isnan(abcb1_28) or abcb1_28 == 0:
        pass_memory = False
    else:
        for d in [29, 30, 31]:
            if np.isnan(rows[d]["abcb1_frac"]) or rows[d]["abcb1_frac"] < ABCB1_MEMORY_RATIO * abcb1_28:
                pass_memory = False
                break

    # Gate 2: Survival retention
    n28 = rows[28]["n_live"]
    n31 = rows[31]["n_live"]
    pass_survival = (n28 > 0) and (n31 >= SURVIVAL_RATIO * n28)

    overall = pass_memory and pass_survival

    return {
        "status": "PASS" if overall else "FAIL",
        "rows": rows,
        "pass": overall,
        "pass_memory": pass_memory,
        "pass_survival": pass_survival,
        "abcb1_retention": (rows[31]["abcb1_frac"] / abcb1_28) if (not np.isnan(abcb1_28) and abcb1_28 > 0) else 0.0,
        "survival_retention": (n31 / n28) if n28 > 0 else 0.0,
        "n_live_28": n28,
    }


def main() -> int:
    if not MANIFEST.exists():
        print(f"ERROR: manifest not found: {MANIFEST}")
        print("Run sweep_resistance.py first.")
        return 1

    manifest = json.loads(MANIFEST.read_text())
    variants = manifest["variants"]

    # Collect unique param values for the grid
    km_vals = sorted(set(v["drug_kill_multiplier"] for v in variants.values()))
    ap_vals = sorted(set(v.get("abcb1_production_rate", 0.0005) for v in variants.values()))

    results = {}
    for name, info in variants.items():
        work_dir = Path(info["work_dir"])
        results[name] = evaluate_variant(work_dir)
        results[name]["params"] = info

    # === Print 3x3 Truth Table ===
    print("=" * 72)
    print("  RESISTANCE SWEEP RESULTS — Stage 1 Micro-Sim")
    print("=" * 72)
    print()
    print(f"  drug_kill_multiplier \\ abcb1_production_rate")
    print()

    # Header row
    header = "           |"
    for ap in ap_vals:
        header += f"  ap={ap:<6}|"
    print(header)
    print("  ---------+" + "-----------+" * len(ap_vals))

    for km in km_vals:
        row = f"  km={km:<5} |"
        for ap in ap_vals:
            name = f"km{km}_ap{ap}"
            r = results.get(name, {"status": "MISSING"})
            status = r["status"]
            if status == "PASS":
                cell = "  PASS   "
            elif status == "FAIL":
                cell = "  FAIL   "
            else:
                cell = f"  {status[:7]:7s} "
            row += f"{cell} |"
        print(row)
    print("  ---------+" + "-----------+" * len(ap_vals))

    # === Detailed per-variant output ===
    print()
    print("-" * 72)
    print("  DETAILED RESULTS")
    print("-" * 72)

    for name in sorted(results.keys()):
        r = results[name]
        p = r["params"]
        print(f"\n  [{name}] km={p['drug_kill_multiplier']} ap={p.get('abcb1_production_rate', '?')}  "
              f"job={p['job_id']}  => {r['status']}")

        if r["status"] in ("NO_OUTPUT", "INCOMPLETE", "MISSING"):
            print(f"    {r['status']} — skipped")
            continue

        print(f"    {'day':>4}  {'t(min)':>7}  {'n_live':>7}  {'ABCB1+%':>8}  {'ic_mean':>8}")
        for day in CHECK_DAYS:
            m = r["rows"][day]
            abcb1_pct = 100.0 * m["abcb1_frac"] if not np.isnan(m["abcb1_frac"]) else float("nan")
            print(f"    {day:>4}  {m['t_actual']:>7.0f}  {int(m['n_live']):>7}  "
                  f"{_fmt(abcb1_pct):>8}  {_fmt(m['ic_mean']):>8}")

        print(f"    Memory gate:   {'PASS' if r.get('pass_memory') else 'FAIL'}")
        print(f"    Survival gate: {'PASS' if r.get('pass_survival') else 'FAIL'}")

    # === Best combo recommendation ===
    print()
    print("=" * 72)
    passing = [(name, r) for name, r in results.items() if r["pass"]]

    if not passing:
        print("  NO PASSING COMBINATIONS.")
        print("  All 9 combos failed. The issue may be elsewhere (not these parameters).")
        print("=" * 72)
        return 2

    # Rank by: (1) n_live at d28 closest to 20-30 (sweet spot), (2) highest ABCB1 retention
    def score(item):
        name, r = item
        n28 = r["n_live_28"]
        # Prefer ~25 survivors at d28 (sweet spot for RC2)
        n_score = -abs(n28 - 25)
        abcb1_score = r.get("abcb1_retention", 0)
        return (n_score, abcb1_score)

    passing.sort(key=score, reverse=True)
    best_name, best_r = passing[0]
    best_p = best_r["params"]

    print(f"  RECOMMENDED: {best_name}")
    print(f"    drug_kill_multiplier    = {best_p['drug_kill_multiplier']}")
    print(f"    abcb1_production_rate   = {best_p.get('abcb1_production_rate', '?')}")
    print(f"    drug_stress_threshold   = {best_p.get('drug_stress_threshold', '?')}")
    print(f"    Day 28 survivors:         {int(best_r['n_live_28'])}")
    print(f"    ABCB1 retention d31/d28:  {best_r.get('abcb1_retention', 0):.2%}")
    print(f"    Survival retention:       {best_r.get('survival_retention', 0):.2%}")
    print()
    print(f"  {len(passing)}/{len(results)} combos passed.")
    if len(passing) > 1:
        print(f"  All passing combos:")
        for name, r in passing:
            p = r["params"]
            print(f"    {name}: n28={int(r['n_live_28'])} "
                  f"abcb1_ret={r.get('abcb1_retention', 0):.2%} "
                  f"surv_ret={r.get('survival_retention', 0):.2%}")
    print("=" * 72)

    return 0


if __name__ == "__main__":
    sys.exit(main())
