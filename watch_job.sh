#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════
#  PROJECT NORTHSTAR — Live Job Monitor
#  Usage:  bash watch_job.sh [SLURM_JOB_ID]
#          (auto-detects your running jobs if no ID given)
# ═══════════════════════════════════════════════════════════════════════
set -uo pipefail

# ──────────────────────── Configuration ────────────────────────────────
PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
DEFAULT_OUT_DIR="$PROJECT_ROOT/build/rc2_full_seed42/replicate_01_seed42/output"

SAVE_INTERVAL=360       # minutes per snapshot
TOTAL_TIME=60480        # total sim minutes (42 days)
TOTAL_SNAPS=168         # TOTAL_TIME / SAVE_INTERVAL

SNAP_PRE=56             # end of BARRIER  (day 14)
SNAP_TREAT=112          # end of DRUG-ON  (day 28)

POLL=10                 # seconds between refreshes
BAR_WIDTH=62            # characters for progress bar

# ──────────────────────── Colors / Symbols ─────────────────────────────
RST="\033[0m"
BOLD="\033[1m"
DIM="\033[2m"
RED="\033[1;31m"
GRN="\033[1;32m"
YEL="\033[1;33m"
CYN="\033[1;36m"
BLU="\033[1;34m"
MAG="\033[1;35m"
WHT="\033[1;37m"
BG_GRN="\033[42m"
BG_YEL="\033[43m"
BG_BLU="\033[44m"
BG_RED="\033[41m"

SPIN_CHARS=('⠋' '⠙' '⠹' '⠸' '⠼' '⠴' '⠦' '⠧' '⠇' '⠏')
SPARK_CHARS=('▁' '▂' '▃' '▄' '▅' '▆' '▇' '█')

# ──────────────────────── Auto-detect Job ──────────────────────────────
if [ -n "${1:-}" ]; then
    JOB_ID="$1"
else
    JOB_ID=$(squeue -u "$USER" --format="%i" --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [ -z "$JOB_ID" ]; then
        echo -e "${RED}ERROR:${RST} No running jobs found. Pass a job ID: bash watch_job.sh <JOB_ID>"
        exit 1
    fi
    echo -e "${DIM}Auto-detected job: ${JOB_ID}${RST}"
fi

# ──────────────────────── Detect output dir ────────────────────────────
OUT_DIR="$DEFAULT_OUT_DIR"
# Try to find from scontrol
SLURM_CMD=$(scontrol show job "$JOB_ID" 2>/dev/null | grep -oP 'Command=\K\S+' || true)
if [ -n "$SLURM_CMD" ] && [ -f "$SLURM_CMD" ]; then
    # Try to extract output dir from the slurm script directory
    SLURM_DIR=$(dirname "$SLURM_CMD")
    if [ -d "$SLURM_DIR/output" ]; then
        OUT_DIR="$SLURM_DIR/output"
    fi
fi

# ──────────────────────── Helper Functions ─────────────────────────────

job_info() {
    # Returns: STATE|ELAPSED|TIMELIMIT|PARTITION|NODES|NODELIST|JOBNAME|NCPUS|MEM
    squeue -j "$JOB_ID" --format="%T|%M|%l|%P|%D|%N|%j|%C|%m" --noheader 2>/dev/null | head -1
}

job_state_sacct() {
    sacct -j "$JOB_ID" --format=State --noheader 2>/dev/null | head -1 | tr -d ' '
}

phase_of() {
    local s=$1
    if   (( s <= SNAP_PRE ));   then echo "BARRIER"
    elif (( s <= SNAP_TREAT )); then echo "DRUG-ON"
    else                              echo "REGROWTH"
    fi
}

phase_color() {
    case "$1" in
        BARRIER)  echo "$BLU" ;;
        DRUG-ON)  echo "$MAG" ;;
        REGROWTH) echo "$CYN" ;;
        *)        echo "$WHT" ;;
    esac
}

state_color() {
    case "$1" in
        RUNNING)    echo "$GRN" ;;
        COMPLETED)  echo "$CYN" ;;
        PENDING)    echo "$YEL" ;;
        FAILED)     echo "$RED" ;;
        CANCELLED*) echo "$RED" ;;
        TIMEOUT)    echo "$RED" ;;
        *)          echo "$WHT" ;;
    esac
}

# Tumor cell count via Python
tumor_count() {
    python3 -c "
import sys, numpy as np
sys.path.insert(0, '$PROJECT_ROOT')
from python.wrapper.output_parser import OutputParser
p = OutputParser('$OUT_DIR')
s = p._read_physicell_xml('$1')
m = s['cell_matrix']; l = s['label_name_map']
ct = m[int(l['cell_type']['index']), :]
dead = m[int(l['dead']['index']), :] if 'dead' in l else np.zeros(m.shape[1])
dm = m[int(l['current_death_model']['index']), :] if 'current_death_model' in l else np.zeros(m.shape[1])
tumor = (np.rint(ct).astype(int)==0) & (dead<0.5) & (np.rint(dm).astype(int)!=100)
stroma = (np.rint(ct).astype(int)==1) & (dead<0.5) & (np.rint(dm).astype(int)!=100)
total = int(m.shape[1])
print(f'{int(np.sum(tumor))}|{int(np.sum(stroma))}|{total}')
" 2>/dev/null || echo "?|?|?"
}

# Last lines of stderr log
last_errors() {
    local slurm_err
    slurm_err=$(find "$(dirname "$OUT_DIR")" -name "slurm_${JOB_ID}.err" 2>/dev/null | head -1)
    if [ -n "$slurm_err" ] && [ -f "$slurm_err" ]; then
        tail -3 "$slurm_err" 2>/dev/null
    else
        echo ""
    fi
}

# Disk usage of output directory
disk_usage() {
    du -sh "$OUT_DIR" 2>/dev/null | cut -f1
}

# Build sparkline from tumor count history
# Globals: TUMOR_HISTORY (array)
sparkline() {
    local arr=("$@")
    local n=${#arr[@]}
    if (( n < 2 )); then echo ""; return; fi

    local min=${arr[0]} max=${arr[0]}
    for v in "${arr[@]}"; do
        (( v < min )) && min=$v
        (( v > max )) && max=$v
    done

    local range=$(( max - min ))
    (( range == 0 )) && range=1

    local line=""
    # Show last 30 values max
    local start=0
    (( n > 30 )) && start=$(( n - 30 ))
    for (( i=start; i<n; i++ )); do
        local idx=$(( (arr[i] - min) * 7 / range ))
        (( idx > 7 )) && idx=7
        (( idx < 0 )) && idx=0
        line+="${SPARK_CHARS[$idx]}"
    done
    echo "$line"
}

# Build the multi-segment progress bar
build_progress_bar() {
    local n_snaps=$1
    local seg1=$(( SNAP_PRE * BAR_WIDTH / TOTAL_SNAPS ))
    local seg2=$(( (SNAP_TREAT - SNAP_PRE) * BAR_WIDTH / TOTAL_SNAPS ))
    local filled=$(( n_snaps * BAR_WIDTH / TOTAL_SNAPS ))

    local bar=""
    for (( i=0; i<BAR_WIDTH; i++ )); do
        if (( i == seg1 )) || (( i == seg1 + seg2 )); then
            bar+="${WHT}│${RST}"
        elif (( i < filled )); then
            if   (( i < seg1 ));          then bar+="${BLU}█${RST}"
            elif (( i < seg1 + seg2 ));   then bar+="${MAG}█${RST}"
            else                               bar+="${CYN}█${RST}"
            fi
        else
            bar+="${DIM}░${RST}"
        fi
    done
    echo -e "$bar"
}

# Format seconds as HH:MM:SS
fmt_time() {
    local s=$1
    printf "%02d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60))
}

# ──────────────────────── State Variables ──────────────────────────────
LAST_TUMOR_XML=""
TUMOR_LIVE="?"
STROMA_LIVE="?"
TOTAL_CELLS="?"
declare -a TUMOR_HISTORY=()
declare -a SNAP_TIMESTAMPS=()
START_TIME=$(date +%s)
SPIN_IDX=0
TICK=0

# ──────────────────────── Initial Banner ───────────────────────────────
clear
cat <<'EOF'

    ╔══════════════════════════════════════════════════════════════════════╗
    ║                                                                    ║
    ║    ███╗   ██╗ ██████╗ ██████╗ ████████╗██╗  ██╗███████╗████████╗   ║
    ║    ████╗  ██║██╔═══██╗██╔══██╗╚══██╔══╝██║  ██║██╔════╝╚══██╔══╝   ║
    ║    ██╔██╗ ██║██║   ██║██████╔╝   ██║   ███████║███████╗   ██║      ║
    ║    ██║╚██╗██║██║   ██║██╔══██╗   ██║   ██╔══██║╚════██║   ██║      ║
    ║    ██║ ╚████║╚██████╔╝██║  ██║   ██║   ██║  ██║███████║   ██║      ║
    ║    ╚═╝  ╚═══╝ ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚══════╝   ╚═╝      ║
    ║                                                                    ║
    ║               A R  ── Live Job Monitor ──                          ║
    ╚══════════════════════════════════════════════════════════════════════╝

EOF

# Print placeholder lines that we'll overwrite
DISPLAY_LINES=20
for (( i=0; i<DISPLAY_LINES; i++ )); do echo ""; done

# ──────────────────────── Main Loop ────────────────────────────────────
while true; do
    NOW=$(date +%s)
    WALL_ELAPSED=$(( NOW - START_TIME ))
    SPIN_IDX=$(( (SPIN_IDX + 1) % ${#SPIN_CHARS[@]} ))
    SPINNER="${SPIN_CHARS[$SPIN_IDX]}"
    TICK=$(( TICK + 1 ))

    # ── Fetch job info ─────────────────────────────────────────────
    INFO=$(job_info)
    if [ -n "$INFO" ]; then
        IFS='|' read -r STATE JOB_ELAPSED TIMELIMIT PARTITION NODES NODELIST JOBNAME NCPUS MEM <<< "$INFO"
        STATE=$(echo "$STATE" | tr -d ' ')
    else
        STATE=$(job_state_sacct)
        JOB_ELAPSED="N/A"
        TIMELIMIT="N/A"
        PARTITION="?"
        NODES="?"
        NODELIST="?"
        JOBNAME="?"
        NCPUS="?"
        MEM="?"
    fi

    S_COLOR=$(state_color "$STATE")

    # ── Snapshot count ─────────────────────────────────────────────
    N_SNAPS=$(find "$OUT_DIR" -maxdepth 1 -name "output*.xml" 2>/dev/null | wc -l)
    N_SVG=$(find "$OUT_DIR" -maxdepth 1 -name "snapshot*.svg" 2>/dev/null | wc -l)
    SIM_MIN=$(( N_SNAPS * SAVE_INTERVAL ))
    SIM_DAY=$(echo "scale=1; $SIM_MIN / 1440" | bc 2>/dev/null || echo "?")
    PHASE=$(phase_of "$N_SNAPS")
    P_COLOR=$(phase_color "$PHASE")
    PCT=$(( N_SNAPS * 100 / TOTAL_SNAPS ))

    # ── ETA calculation ────────────────────────────────────────────
    if (( N_SNAPS > 3 )); then
        SNAP_RATE=$(echo "scale=4; $N_SNAPS / $WALL_ELAPSED" | bc 2>/dev/null || echo "0")
        REMAINING=$(( TOTAL_SNAPS - N_SNAPS ))
        if [ "$SNAP_RATE" != "0" ] && [ -n "$SNAP_RATE" ]; then
            ETA_SECS=$(echo "scale=0; $REMAINING / $SNAP_RATE" | bc 2>/dev/null || echo "0")
            ETA_STR=$(fmt_time "${ETA_SECS:-0}")
        else
            ETA_STR="calculating..."
        fi
        # Speed: snapshots per minute
        SPEED=$(echo "scale=2; $SNAP_RATE * 60" | bc 2>/dev/null || echo "?")
    else
        ETA_STR="calculating..."
        SPEED="?"
    fi

    # ── Tumor count (only on new snapshot) ─────────────────────────
    LATEST_XML=$(find "$OUT_DIR" -maxdepth 1 -name "output*.xml" 2>/dev/null | sort | tail -1)
    if [ -n "${LATEST_XML:-}" ] && [ "${LATEST_XML:-}" != "${LAST_TUMOR_XML:-}" ]; then
        COUNTS=$(tumor_count "$LATEST_XML")
        IFS='|' read -r TUMOR_LIVE STROMA_LIVE TOTAL_CELLS <<< "$COUNTS"
        LAST_TUMOR_XML="$LATEST_XML"
        # Record history for sparkline
        if [ "$TUMOR_LIVE" != "?" ]; then
            TUMOR_HISTORY+=("$TUMOR_LIVE")
        fi
    fi

    # ── Disk usage (every 6 ticks) ─────────────────────────────────
    if (( TICK % 6 == 1 )); then
        DISK=$(disk_usage)
    fi

    # ── Build sparkline ────────────────────────────────────────────
    SPARK=""
    if (( ${#TUMOR_HISTORY[@]} > 1 )); then
        SPARK=$(sparkline "${TUMOR_HISTORY[@]}")
    fi

    # ── Build progress bar ─────────────────────────────────────────
    PBAR=$(build_progress_bar "$N_SNAPS")

    # ── Time limit parsing for time-used bar ───────────────────────
    TLIMIT_BAR=""
    if [[ "$TIMELIMIT" =~ ^([0-9]+):([0-9]+):([0-9]+)$ ]]; then
        TL_SECS=$(( ${BASH_REMATCH[1]} * 3600 + ${BASH_REMATCH[2]} * 60 + ${BASH_REMATCH[3]} ))
        if [[ "$JOB_ELAPSED" =~ ^([0-9]+):([0-9]+):([0-9]+)$ ]]; then
            JE_SECS=$(( ${BASH_REMATCH[1]} * 3600 + ${BASH_REMATCH[2]} * 60 + ${BASH_REMATCH[3]} ))
        elif [[ "$JOB_ELAPSED" =~ ^([0-9]+):([0-9]+)$ ]]; then
            JE_SECS=$(( ${BASH_REMATCH[1]} * 60 + ${BASH_REMATCH[2]} ))
        else
            JE_SECS=0
        fi
        if (( TL_SECS > 0 )); then
            TIME_PCT=$(( JE_SECS * 100 / TL_SECS ))
            TIME_FILLED=$(( JE_SECS * 20 / TL_SECS ))
            TLIMIT_BAR=""
            for (( i=0; i<20; i++ )); do
                if (( i < TIME_FILLED )); then
                    if (( TIME_PCT > 80 )); then
                        TLIMIT_BAR+="${RED}█${RST}"
                    elif (( TIME_PCT > 60 )); then
                        TLIMIT_BAR+="${YEL}█${RST}"
                    else
                        TLIMIT_BAR+="${GRN}█${RST}"
                    fi
                else
                    TLIMIT_BAR+="${DIM}░${RST}"
                fi
            done
            TLIMIT_BAR="[$(echo -e "$TLIMIT_BAR")] ${TIME_PCT}%"
        fi
    fi

    # ── Draw ───────────────────────────────────────────────────────
    printf "\033[${DISPLAY_LINES}A"   # move cursor up

    printf "  ${DIM}─────────────────────────── Job Info ──────────────────────────────${RST}\n"
    printf "  ${BOLD}Job ID:${RST}    %-12s  ${BOLD}Name:${RST}  %-20s  ${BOLD}State:${RST} %b %s\n" \
        "$JOB_ID" "$JOBNAME" "${S_COLOR}● ${STATE}${RST}" "$SPINNER"
    printf "  ${BOLD}Partition:${RST} %-12s  ${BOLD}Node:${RST}  %-20s  ${BOLD}CPUs:${RST}  %s  ${BOLD}Mem:${RST} %s\n" \
        "$PARTITION" "$NODELIST" "$NCPUS" "$MEM"
    printf "  ${BOLD}Elapsed:${RST}   %-12s  ${BOLD}Limit:${RST} %-20s  %b\n" \
        "$JOB_ELAPSED" "$TIMELIMIT" "$TLIMIT_BAR"
    echo ""

    printf "  ${DIM}─────────────────────── Simulation Progress ────────────────────────${RST}\n"
    printf "  %b\n" "$PBAR"
    printf "  ${BLU}■${RST} BARRIER (d0-14)  ${WHT}│${RST}  ${MAG}■${RST} DRUG-ON (d14-28)  ${WHT}│${RST}  ${CYN}■${RST} REGROWTH (d28-42)\n"
    echo ""
    printf "  ${BOLD}Day:${RST} %-8s ${BOLD}Snap:${RST} %3d / %-3d  [${BOLD}%3d%%${RST}]   ${BOLD}Phase:${RST} %b%-10s${RST}  ${BOLD}ETA:${RST} %s\n" \
        "$SIM_DAY" "$N_SNAPS" "$TOTAL_SNAPS" "$PCT" "$P_COLOR" "$PHASE" "$ETA_STR"
    printf "  ${BOLD}Speed:${RST} %s snap/min                        ${BOLD}SVGs:${RST} %d   ${BOLD}Disk:${RST} %s\n" \
        "$SPEED" "$N_SVG" "${DISK:-...}"
    echo ""

    printf "  ${DIM}────────────────────────── Cell Counts ─────────────────────────────${RST}\n"
    printf "  ${BOLD}🧬 Tumor:${RST}  %-8s  ${BOLD}🔬 Stroma:${RST} %-8s  ${BOLD}📊 Total:${RST} %-8s\n" \
        "$TUMOR_LIVE" "$STROMA_LIVE" "$TOTAL_CELLS"

    if [ -n "$SPARK" ]; then
        printf "  ${BOLD}Tumor trend:${RST} ${GRN}%s${RST}\n" "$SPARK"
    else
        printf "  ${DIM}Tumor trend: (collecting data...)${RST}\n"
    fi
    echo ""

    printf "  ${DIM}Monitor: %-8s │ Refresh: ${POLL}s │ %s${RST}\n" \
        "$(fmt_time $WALL_ELAPSED)" "$(date '+%H:%M:%S')"
    printf "  ${DIM}Ctrl+C to detach (job continues running)${RST}\n"

    # ── Check terminal state ───────────────────────────────────────
    if [ "$STATE" = "COMPLETED" ] || [ "$STATE" = "FAILED" ] || \
       [[ "$STATE" == CANCELLED* ]] || [ "$STATE" = "TIMEOUT" ]; then
        echo ""
        if [ "$STATE" = "COMPLETED" ]; then
            printf "  ${GRN}✓ Job completed successfully!${RST}\n"
            printf "  ${BOLD}Next:${RST} python diagnose_rc2_full.py\n"
        else
            printf "  ${RED}✗ Job ended: %s${RST}\n" "$STATE"
            ERR_FILE=$(find "$(dirname "$OUT_DIR")" -name "slurm_${JOB_ID}.err" 2>/dev/null | head -1)
            if [ -n "$ERR_FILE" ] && [ -f "$ERR_FILE" ]; then
                printf "  ${BOLD}Last errors:${RST}\n"
                tail -5 "$ERR_FILE" 2>/dev/null | sed 's/^/    /'
            fi
        fi
        break
    fi

    sleep "$POLL"
done
