#!/bin/bash
# Watch fixH sweep (4 fixes x RC1+RC2)
WORK="/work/a0hirs04/PROJECT-NORTHSTAR/build/fixH"

while true; do
    clear
    echo "=== fixH Monitor  $(date) ==="
    echo ""

    echo "── SLURM Jobs ──"
    squeue -u "$USER" -o "%.10i %.22j %.2t %.10M" 2>/dev/null | grep -E "JOBID|fH_" || echo "  (no jobs in queue)"
    echo ""

    printf "%-6s  %11s %10s  %7s %8s\n" "Fix" "caf_contact" "drug_kill" "RC1" "RC2"
    echo "------------------------------------------------------"

    all_done=true
    for f in fixA fixB fixC fixD; do
        rc1_snaps=$(ls "$WORK/$f/rc1/replicate_01_seed42/output"/output*.xml 2>/dev/null | wc -l)
        rc2_snaps=$(ls "$WORK/$f/rc2/output"/output*.xml 2>/dev/null | wc -l)
        [[ $rc1_snaps -ge 85 ]] && rc1_st="DONE" || { rc1_st="${rc1_snaps}/85"; all_done=false; }
        [[ $rc2_snaps -ge 169 ]] && rc2_st="DONE" || { rc2_st="${rc2_snaps}/169"; all_done=false; }

        cfg="$WORK/$f/rc1/replicate_01_seed42/PhysiCell_settings.xml"
        caf=$(grep "ecm_emt_require_caf_contact" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
        dk=$(grep "drug_kill_coefficient" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')

        printf "%-6s  %11s %10s  %7s %8s\n" "$f" "$caf" "$dk" "$rc1_st" "$rc2_st"
    done

    echo ""
    if $all_done; then
        echo "══ ALL 4 FIXES COMPLETE ══"
        break
    fi

    sleep 60
done
