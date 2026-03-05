#!/bin/bash
# RC2 live progress monitor — run with: bash tests/watch_rc2.sh
# Press Ctrl+C to stop.

RC2_DIR="/home/a0hirs04/PROJECT-NORTHSTAR/build/reality_check_2"
TOTAL_SNAPS=168   # 60480 / 360

while true; do
    clear
    echo "═══════════════════════════════════════════════════════════════"
    echo "  RC2 LIVE MONITOR  $(date '+%H:%M:%S')"
    echo "═══════════════════════════════════════════════════════════════"
    echo ""

    # SLURM job status
    echo "  SLURM JOBS:"
    squeue -u "$USER" --format="    %-10i %-20j %-8T %-10M %-15R" --noheader 2>/dev/null | grep -i rc2
    if [ $? -ne 0 ]; then
        echo "    (no rc2 jobs in queue — all finished or not started)"
    fi
    echo ""

    # Per-replicate snapshot progress
    echo "  SNAPSHOT PROGRESS (target: ${TOTAL_SNAPS}):"
    for d in "$RC2_DIR"/replicate_*/; do
        [ -d "$d" ] || continue
        name=$(basename "$d")
        n=$(ls "$d"/output/output*.xml 2>/dev/null | wc -l)
        pct=$((n * 100 / TOTAL_SNAPS))
        bar=""
        filled=$((pct / 5))
        for ((i=0; i<filled; i++)); do bar+="█"; done
        for ((i=filled; i<20; i++)); do bar+="░"; done
        printf "    %-28s  %3d/%d  [%s] %3d%%\n" "$name" "$n" "$TOTAL_SNAPS" "$bar" "$pct"
    done
    echo ""

    # Check for DRUG_SCHEDULE messages in stderr
    echo "  DRUG SCHEDULE EVENTS:"
    for d in "$RC2_DIR"/replicate_*/; do
        [ -d "$d" ] || continue
        name=$(basename "$d")
        err=$(ls "$d"/slurm_*.err 2>/dev/null | head -1)
        if [ -n "$err" ]; then
            msgs=$(grep "DRUG_SCHEDULE" "$err" 2>/dev/null)
            if [ -n "$msgs" ]; then
                echo "    $name:"
                echo "$msgs" | sed 's/^/      /'
            fi
        fi
    done
    echo ""

    # Check if all done
    all_done=true
    for d in "$RC2_DIR"/replicate_*/; do
        [ -d "$d" ] || continue
        n=$(ls "$d"/output/output*.xml 2>/dev/null | wc -l)
        if [ "$n" -lt "$TOTAL_SNAPS" ]; then
            all_done=false
            break
        fi
    done
    if $all_done; then
        echo "  *** ALL REPLICATES COMPLETE — run evaluation ***"
    fi

    sleep 30
done
