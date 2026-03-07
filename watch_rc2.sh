#!/bin/bash
# Live monitor for RC2 full run.
# Usage: bash watch_rc2.sh <SLURM_JOB_ID>

set -euo pipefail

JOB_ID="${1:?Usage: bash watch_rc2.sh <JOB_ID>}"
OUT_DIR="/home/a0hirs04/PROJECT-NORTHSTAR/build/rc2_full_seed42/replicate_01_seed42/output"
SAVE_INTERVAL=360
TOTAL_SNAPS=168

SNAP_PRE=56
SNAP_TREAT=112

POLL=15
LAST_TUMOR_XML=""
TUMOR="?"
START_TIME=$(date +%s)

phase_of() {
    local s=$1
    if   [ "$s" -le "$SNAP_PRE" ];   then echo "BARRIER"
    elif [ "$s" -le "$SNAP_TREAT" ]; then echo "DRUG-ON"
    else                                   echo "REGROWTH"
    fi
}

job_state() {
    local jid=$1 st
    st=$(squeue -j "$jid" --format="%T" --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [ -n "$st" ]; then printf "%s" "$st"; return; fi
    st=$(sacct -j "$jid" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ')
    printf "%s" "${st:-UNKNOWN}"
}

tumor_count() {
    python3 -c "
import sys, numpy as np
sys.path.insert(0, '/home/a0hirs04/PROJECT-NORTHSTAR')
from python.wrapper.output_parser import OutputParser
p = OutputParser('$OUT_DIR')
s = p._read_physicell_xml('$1')
m = s['cell_matrix']; l = s['label_name_map']
ct = m[int(l['cell_type']['index']), :]
dead = m[int(l['dead']['index']), :] if 'dead' in l else np.zeros(m.shape[1])
dm = m[int(l['current_death_model']['index']), :] if 'current_death_model' in l else np.zeros(m.shape[1])
t = (np.rint(ct).astype(int)==0) & (dead<0.5) & (np.rint(dm).astype(int)!=100)
print(int(np.sum(t)))
" 2>/dev/null || echo "?"
}

clear
cat <<'BANNER'

  ╔══════════════════════════════════════════════════════════════════╗
  ║              RC2 FULL MONITOR — seed 42                        ║
  ╠══════════════════════════════════════════════════════════════════╣
  ║  BARRIER  d0-14   ░░░░░░░░░░░░░░  snap  0-56                  ║
  ║  DRUG-ON  d14-28  ░░░░░░░░░░░░░░  snap 56-112                 ║
  ║  REGROWTH d28-42  ░░░░░░░░░░░░░░  snap 112-168                ║
  ╚══════════════════════════════════════════════════════════════════╝

BANNER
printf "  Job: %s\n\n" "$JOB_ID"

while true; do
    STATE=$(job_state "$JOB_ID")
    N_SNAPS=$(ls "$OUT_DIR"/output*.xml 2>/dev/null | wc -l)
    SIM_DAY=$(echo "scale=1; $N_SNAPS * $SAVE_INTERVAL / 1440" | bc 2>/dev/null || echo "?")
    PHASE=$(phase_of "$N_SNAPS")
    NOW=$(date +%s)
    ELAPSED=$(( NOW - START_TIME ))
    ELAPSED_M=$(( ELAPSED / 60 ))
    ELAPSED_S=$(( ELAPSED % 60 ))

    # ETA
    if [ "$N_SNAPS" -gt 2 ]; then
        ETA_SECS=$(( ELAPSED * (TOTAL_SNAPS - N_SNAPS) / N_SNAPS ))
        ETA_M=$(( ETA_SECS / 60 ))
        ETA_STR="${ETA_M}m"
    else
        ETA_STR="..."
    fi

    # Tumor count (only re-read on new snapshot)
    LATEST_XML=$(ls "$OUT_DIR"/output*.xml 2>/dev/null | sort | tail -1)
    if [ -n "${LATEST_XML:-}" ] && [ "${LATEST_XML:-}" != "${LAST_TUMOR_XML:-}" ]; then
        TUMOR=$(tumor_count "$LATEST_XML")
        LAST_TUMOR_XML="$LATEST_XML"
    fi

    # Build 3-segment progress bar (width 60)
    # BARRIER=20 chars, DRUG=20 chars, REGROWTH=20 chars
    BAR_W=60
    SEG1=$(( SNAP_PRE * BAR_W / TOTAL_SNAPS ))       # 20
    SEG2=$(( (SNAP_TREAT - SNAP_PRE) * BAR_W / TOTAL_SNAPS ))  # 20
    SEG3=$(( BAR_W - SEG1 - SEG2 ))                   # 20
    FILLED=$(( N_SNAPS * BAR_W / TOTAL_SNAPS ))

    BAR=""
    for ((i=0; i<BAR_W; i++)); do
        if [ $i -eq $SEG1 ] || [ $i -eq $(( SEG1 + SEG2 )) ]; then
            BAR+="│"
        elif [ $i -lt $FILLED ]; then
            if   [ $i -lt $SEG1 ];                    then BAR+="▓"   # barrier (blue-ish)
            elif [ $i -lt $(( SEG1 + SEG2 )) ];       then BAR+="█"   # drug (solid)
            else                                            BAR+="▒"   # regrowth
            fi
        else
            BAR+="░"
        fi
    done

    # Color the state
    case "$STATE" in
        RUNNING)   STATE_FMT="\033[1;32m${STATE}\033[0m" ;;
        COMPLETED) STATE_FMT="\033[1;36m${STATE}\033[0m" ;;
        PENDING)   STATE_FMT="\033[1;33m${STATE}\033[0m" ;;
        *)         STATE_FMT="\033[1;31m${STATE}\033[0m" ;;
    esac

    # Move cursor up 6 lines and redraw
    printf "\033[6A"
    printf "  State:    %b\n" "$STATE_FMT"
    printf "  Elapsed:  %dm%02ds   ETA: %s\n" "$ELAPSED_M" "$ELAPSED_S" "$ETA_STR"
    printf "  \n"
    printf "  [%s]\n" "$BAR"
    printf "   BARRIER     │   DRUG-ON    │  REGROWTH\n"
    printf "  \n"
    printf "  Day: %-6s  Snaps: %3d/%-3d  Phase: %-10s  Tumor: %s   \n" \
        "$SIM_DAY" "$N_SNAPS" "$TOTAL_SNAPS" "$PHASE" "$TUMOR"

    if [ "$STATE" = "COMPLETED" ] || [ "$STATE" = "FAILED" ] || \
       [ "$STATE" = "CANCELLED" ] || [ "$STATE" = "TIMEOUT" ]; then
        echo ""
        if [ "$STATE" = "COMPLETED" ]; then
            echo "  Done! Run:  python diagnose_rc2_full.py"
        else
            echo "  Job ended: $STATE"
            echo "  Check: cat $(dirname "$OUT_DIR")/slurm_${JOB_ID}.err"
        fi
        break
    fi

    sleep "$POLL"
done
