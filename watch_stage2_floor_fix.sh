#!/bin/bash
# Usage: bash watch_stage2_floor_fix.sh <SLURM_JOB_ID>
# Watches Stage 2 RC2 "floor fix" run — tracks cell counts, drug kill, ABCB1

set -euo pipefail

JOB_ID="${1:?Usage: bash watch_stage2_floor_fix.sh <JOB_ID>}"
ROOT="/home/a0hirs04/PROJECT-NORTHSTAR"
WORK_DIR="$ROOT/build/rc2_full_seed42/replicate_01_seed42"
OUT_DIR="$WORK_DIR/output"
OUT_LOG="$ROOT/logs/stage2_rc2_${JOB_ID}.out"
ERR_LOG="$ROOT/logs/stage2_rc2_${JOB_ID}.err"

job_state() {
  local jid="$1" st
  st=$(squeue -j "$jid" --format="%T" --noheader 2>/dev/null | head -1 | tr -d ' ')
  if [ -n "$st" ]; then printf "%s" "$st"; return; fi
  st=$(sacct -j "$jid" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ')
  printf "%s" "${st:-UNKNOWN}"
}

analyze_snapshots() {
  python3 - "$OUT_DIR" <<'PYEOF'
import sys, os, glob
import scipy.io, numpy as np

out_dir = sys.argv[1]
mats = sorted(glob.glob(os.path.join(out_dir, "output*_cells.mat")))
if not mats:
    print("  (no snapshots yet)")
    return

print(f"  {'Day':>5}  {'N_cells':>7}  {'Apo_rate':>10}  {'IntraDrug':>10}  {'NRF2':>8}  {'ABCB1':>8}  {'UptakeD':>8}")
print(f"  {'---':>5}  {'---':>7}  {'---':>10}  {'---':>10}  {'---':>8}  {'---':>8}  {'---':>8}")

# Sample every 4th snapshot (~1 day)
step = max(1, len(mats) // 20)
for mat_path in mats[::step]:
    snap_num = int(os.path.basename(mat_path).replace("output","").replace("_cells.mat",""))
    t_min = snap_num * 360
    t_day = t_min / 1440.0
    try:
        mat = scipy.io.loadmat(mat_path)
        c = mat['cells']
        n = c.shape[1]
        nr = c.shape[0]
        # Key rows from label map:
        # 28: death_rates[0] (apoptosis), 76: uptake_rates[3] (drug)
        # 172: nrf2_active, 173: abcb1_active, 175: intracellular_drug
        apo = c[28,:].mean() if nr > 28 else 0
        uptk = c[76,:].mean() if nr > 76 else 0
        nrf2 = c[172,:].mean() if nr > 172 else 0
        abcb1 = c[173,:].mean() if nr > 173 else 0
        idrug = c[175,:].mean() if nr > 175 else 0
        print(f"  {t_day:>5.1f}  {n:>7d}  {apo:>10.6f}  {idrug:>10.6f}  {nrf2:>8.4f}  {abcb1:>8.4f}  {uptk:>8.4f}")
    except Exception as e:
        print(f"  {t_day:>5.1f}  ERROR: {e}")

# Show last snapshot details
try:
    mat = scipy.io.loadmat(mats[-1])
    c = mat['cells']
    snap_num = int(os.path.basename(mats[-1]).replace("output","").replace("_cells.mat",""))
    t_day = snap_num * 360 / 1440.0
    n = c.shape[1]
    nr = c.shape[0]
    abcb1_vals = c[173,:] if nr > 173 else np.zeros(n)
    idrug_vals = c[175,:] if nr > 175 else np.zeros(n)
    apo_vals = c[28,:] if nr > 28 else np.zeros(n)
    print(f"\n  Latest snapshot (day {t_day:.1f}, {n} cells):")
    print(f"    intracellular_drug: mean={idrug_vals.mean():.6f}  max={idrug_vals.max():.6f}  nonzero={np.count_nonzero(idrug_vals)}/{n}")
    print(f"    abcb1_active:       mean={abcb1_vals.mean():.4f}  max={abcb1_vals.max():.4f}  >0.5={np.sum(abcb1_vals>0.5)}/{n}")
    print(f"    apoptosis_rate:     mean={apo_vals.mean():.6f}  max={apo_vals.max():.6f}  >0.001={np.sum(apo_vals>0.001)}/{n}")
except Exception:
    pass

PYEOF
}

while true; do
  clear 2>/dev/null || true
  STATE=$(job_state "$JOB_ID")
  N_SNAPS=$(ls "$OUT_DIR"/output*_cells.mat 2>/dev/null | wc -l)
  ELAPSED=$(squeue -j "$JOB_ID" --format="%M" --noheader 2>/dev/null | tr -d ' ')

  echo "╔══════════════════════════════════════════════════════════════╗"
  echo "║  Stage 2 RC2 — Floor Fix Watch        $(date '+%H:%M:%S')          ║"
  echo "╠══════════════════════════════════════════════════════════════╣"
  echo "║  Job: $JOB_ID   State: $STATE   Elapsed: ${ELAPSED:-n/a}"
  echo "║  Snapshots: $N_SNAPS   Output: $OUT_DIR"
  echo "╚══════════════════════════════════════════════════════════════╝"
  echo

  # Drug kill signature check
  if [ -f "$ERR_LOG" ]; then
    DRUG_ON=$(grep -c "Dirichlet.*drug\|drug_start\|DRUG_BC" "$ERR_LOG" 2>/dev/null || echo 0)
    CUSTOM_DBG=$(grep -c "CUSTOM_DATA_DEBUG" "$ERR_LOG" 2>/dev/null || echo 0)
    echo "  Drug BC events: $DRUG_ON   Custom data debug: $CUSTOM_DBG"
  fi

  echo
  echo "── Snapshot Timeline ──"
  analyze_snapshots 2>/dev/null || echo "  (waiting for data...)"

  echo
  if [ -f "$ERR_LOG" ]; then
    echo "── Recent stderr ──"
    tail -n 10 "$ERR_LOG" 2>/dev/null || true
  fi

  case "$STATE" in
    COMPLETED|FAILED|CANCELLED*|TIMEOUT)
      echo
      echo "════════════════════════════════════════"
      echo "  FINAL STATE: $STATE"
      echo "  Stdout: $OUT_LOG"
      echo "  Stderr: $ERR_LOG"
      echo "════════════════════════════════════════"
      echo
      echo "── Full Timeline ──"
      analyze_snapshots 2>/dev/null || true
      break
      ;;
  esac

  sleep 30
done
