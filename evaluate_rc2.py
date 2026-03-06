#!/usr/bin/env python3
"""Evaluate completed RC2 replicates (no re-submission)."""
from __future__ import annotations
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from run_reality_check_2 import (
    SEEDS, NUM_REPLICATES, QUORUM, CRITERIA_NAMES,
    T_PRE, T_TREAT_END, T_POST,
    PARTIAL_RESPONSE_MIN_REDUCTION, BOUNDARY_BARRIER_FRACTION,
    ReplicateResult, _evaluate,
)

WORK_DIR = PROJECT_ROOT / "build" / "reality_check_2"


def main() -> int:
    results = []
    for i, seed in enumerate(SEEDS):
        rep_dir = WORK_DIR / f"replicate_{i+1:02d}_seed{seed}"
        output_dir = rep_dir / "output"
        has_output = bool(list(output_dir.glob("output*.xml"))) if output_dir.exists() else False
        r = ReplicateResult(seed=seed, run_dir=rep_dir, success=has_output)
        if has_output:
            print(f"\n  Evaluating rep {i+1} (seed={seed}) ...")
            try:
                r = _evaluate(r)
            except Exception as exc:
                import traceback
                print(f"    ERROR: {exc}")
                traceback.print_exc()
                r.success = False
        else:
            print(f"\n  rep {i+1} (seed={seed}): NO OUTPUT")
        results.append(r)

    # Per-replicate
    print()
    print("=" * 72)
    print("  PER-REPLICATE RESULTS")
    print("=" * 72)
    for i, r in enumerate(results):
        hdr = f"Replicate {i+1}  (seed={r.seed})"
        if not r.success:
            print(f"\n  {hdr}  — SIMULATION FAILED")
            continue
        all_pass = all(r.criteria.get(k, False) for k, _ in CRITERIA_NAMES)
        tag = "ALL PASS" if all_pass else "SOME FAIL"
        print(f"\n  {hdr}  — {tag}")
        for key, label in CRITERIA_NAMES:
            ok = r.criteria.get(key, False)
            mark = "PASS" if ok else "FAIL"
            print(f"    [{mark}] {label}")
            print(f"           {r.details.get(key, '')}")
        if r.snap_pre and r.snap_treat_end and r.snap_post:
            sp, st, sw = r.snap_pre, r.snap_treat_end, r.snap_post
            print(f"    Timecourse:")
            print(f"      t={sp.time:.0f}  tumor={sp.n_tumor:4d}  caf={sp.n_caf:4d}  "
                  f"peri_ecm={sp.peri_ecm:.4f}  abcb1={sp.frac_abcb1:.4f}  "
                  f"ecm@tumor={sp.ecm_at_tumor:.4f}")
            print(f"      t={st.time:.0f}  tumor={st.n_tumor:4d}  caf={st.n_caf:4d}  "
                  f"peri_ecm={st.peri_ecm:.4f}  abcb1={st.frac_abcb1:.4f}  "
                  f"ecm@tumor={st.ecm_at_tumor:.4f}")
            print(f"      t={sw.time:.0f}  tumor={sw.n_tumor:4d}  caf={sw.n_caf:4d}  "
                  f"peri_ecm={sw.peri_ecm:.4f}  abcb1={sw.frac_abcb1:.4f}  "
                  f"ecm@tumor={sw.ecm_at_tumor:.4f}")

    # Aggregate
    print()
    print("=" * 72)
    print(f"  AGGREGATE  (>={QUORUM}/{NUM_REPLICATES} replicates must pass)")
    print("=" * 72)

    SOFT_CRITERIA = {"abcb1_emerges"}
    any_failure = False
    for key, label in CRITERIA_NAMES:
        passing = sum(1 for r in results if r.success and r.criteria.get(key, False))
        is_soft = key in SOFT_CRITERIA
        ok = passing >= QUORUM
        if is_soft:
            mark = "INFO" if not ok else "PASS"
        else:
            mark = "PASS" if ok else "FAIL"
            if not ok:
                any_failure = True
        print(f"  [{mark}] {label}  ({passing}/{NUM_REPLICATES})")

    print()
    if any_failure:
        print("  *** REALITY CHECK 2: FAIL ***")
        return 1
    else:
        print("  *** REALITY CHECK 2: PASS ***")
        return 0


if __name__ == "__main__":
    sys.exit(main())
