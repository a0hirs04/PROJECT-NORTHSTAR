#!/bin/bash
# Watch Stage 2 RC2 with sweep-winning parameters
# Usage: bash watch_rc2_sweep_winner.sh [JOB_ID]

JOB_ID="${1:-10620}"
OUT_DIR="/home/a0hirs04/PROJECT-NORTHSTAR/build/rc2_full_seed42/replicate_01_seed42/output"

echo "Watching Stage 2 RC2 — job $JOB_ID (km=1.0, ap=0.0005, threshold=0.03)"
echo "================================================================"

while squeue -j "$JOB_ID" 2>/dev/null | grep -q "$JOB_ID"; do
    n_snap=$(ls "$OUT_DIR"/output*_cells.mat 2>/dev/null | wc -l)
    if [ "$n_snap" -gt 0 ]; then
        last=$(ls -t "$OUT_DIR"/output*_cells.mat | head -1)
        snap_num=$(basename "$last" _cells.mat | sed 's/output0*//')
        t_day=$(echo "scale=1; $snap_num * 360 / 1440" | bc 2>/dev/null)

        python3 -c "
import scipy.io as sio, numpy as np
mat = sio.loadmat('$last')
cells = mat['cells']
phase = cells[1,:]
live = phase < 100
n_live = int(np.sum(live))
n_dead = int(np.sum(~live))
if n_live > 0:
    abcb1 = cells[173, live]
    abcb1_pct = 100.0 * np.sum(abcb1 > 0.1) / n_live
    ic_drug = np.mean(cells[175, live])
    death_rate = np.mean(cells[28, live])
    print(f'  day ~$t_day | snap {$snap_num:>4} | live={n_live:>5} dead={n_dead:>5} | ABCB1+={abcb1_pct:5.1f}% | ic_drug={ic_drug:.4f} | death_rate={death_rate:.6f}')
else:
    print(f'  day ~$t_day | snap {$snap_num:>4} | ALL DEAD (n_dead={n_dead})')
" 2>/dev/null || echo "  snap $n_snap — parse error"
    else
        echo "  No snapshots yet..."
    fi
    sleep 60
done

echo ""
echo "Job $JOB_ID finished. Final snapshot count: $(ls "$OUT_DIR"/output*_cells.mat 2>/dev/null | wc -l)"
