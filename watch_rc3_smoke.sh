#!/bin/bash
# Live monitor for the current RC3 smoke test.
# Usage:
#   bash watch_rc3_smoke.sh
#   bash watch_rc3_smoke.sh --once
#   bash watch_rc3_smoke.sh LABEL:JOB_ID:RUN_DIR [...]

set -euo pipefail

ONCE=0
if [[ "${1:-}" == "--once" ]]; then
    ONCE=1
    shift
fi

DEFAULT_DIRS=(
    "ArmA|/home/a0hirs04/PROJECT-NORTHSTAR/build/rc3_vismodegib_smoke/Arm_A_Control/seed_42"
    "ArmB|/home/a0hirs04/PROJECT-NORTHSTAR/build/rc3_vismodegib_smoke/Arm_B_SHH_Only/seed_42"
)

latest_job_id_for_run_dir() {
    local run_dir="$1"
    local latest
    latest=$(find "$run_dir" -maxdepth 1 -type f -name 'slurm_*.out' -printf '%T@ %f\n' 2>/dev/null | \
        sort -nr | head -1 | awk '{print $2}')
    if [[ -z "${latest:-}" ]]; then
        return 1
    fi
    latest="${latest#slurm_}"
    latest="${latest%.out}"
    printf "%s" "$latest"
}

SPECS=()
if [[ "$#" -gt 0 ]]; then
    for arg in "$@"; do
        IFS=':' read -r label job_id run_dir <<< "$arg"
        if [[ -z "${label:-}" || -z "${job_id:-}" || -z "${run_dir:-}" ]]; then
            printf "Invalid spec: %s\n" "$arg" >&2
            printf "Expected: LABEL:JOB_ID:RUN_DIR\n" >&2
            exit 1
        fi
        SPECS+=("${label}|${job_id}|${run_dir}")
    done
else
    for entry in "${DEFAULT_DIRS[@]}"; do
        IFS='|' read -r label run_dir <<< "$entry"
        job_id=$(latest_job_id_for_run_dir "$run_dir") || {
            printf "Could not find a smoke-test slurm log in %s\n" "$run_dir" >&2
            exit 1
        }
        SPECS+=("${label}|${job_id}|${run_dir}")
    done
fi

POLL_SECONDS=15

job_state() {
    local job_id="$1"
    local state
    state=$(squeue -j "$job_id" --format="%T" --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [[ -n "$state" ]]; then
        printf "%s" "$state"
        return
    fi

    state=$(sacct -j "$job_id" --format=State --noheader --parsable2 2>/dev/null | \
        head -1 | cut -d'|' -f1 | tr -d ' ')
    printf "%s" "${state:-UNKNOWN}"
}

format_seconds() {
    python3 - "$1" <<'PY'
import sys

secs = float(sys.argv[1])
if secs < 0:
    secs = 0.0

days = int(secs // 86400)
secs -= days * 86400
hours = int(secs // 3600)
secs -= hours * 3600
minutes = int(secs // 60)
secs -= minutes * 60
seconds = int(round(secs))

if seconds == 60:
    seconds = 0
    minutes += 1
if minutes == 60:
    minutes = 0
    hours += 1
if hours == 24:
    hours = 0
    days += 1

parts = []
if days:
    parts.append(f"{days}d")
if hours:
    parts.append(f"{hours}h")
if minutes:
    parts.append(f"{minutes}m")
if not parts or seconds:
    parts.append(f"{seconds}s")
print(" ".join(parts))
PY
}

parse_progress() {
    python3 - "$1" <<'PY'
import os
import re
import sys

path = sys.argv[1]
if not os.path.exists(path):
    print("0\t80640\t0\t0\tSTARTUP")
    raise SystemExit(0)

sim_re = re.compile(r"current simulated time:\s*([0-9.eE+-]+)\s*min\s*\(max:\s*([0-9.eE+-]+)\s*min\)")
agents_re = re.compile(r"total agents:\s*(\d+)")
wall_re = re.compile(
    r"total wall time:\s*(\d+)\s*days,\s*(\d+)\s*hours,\s*(\d+)\s*minutes,\s*and\s*([0-9.eE+-]+)\s*seconds"
)

sim_time = 0.0
max_time = 80640.0
agents = 0
wall_seconds = 0.0

with open(path, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        m = sim_re.search(line)
        if m:
            sim_time = float(m.group(1))
            max_time = float(m.group(2))
            continue

        m = agents_re.search(line)
        if m:
            agents = int(m.group(1))
            continue

        m = wall_re.search(line)
        if m:
            days = int(m.group(1))
            hours = int(m.group(2))
            minutes = int(m.group(3))
            seconds = float(m.group(4))
            wall_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds

if sim_time < 20160.0:
    phase = "Wk0-4"
elif sim_time < max_time:
    phase = "Wk4-8"
else:
    phase = "DONE"

print(f"{sim_time}\t{max_time}\t{wall_seconds}\t{agents}\t{phase}")
PY
}

terminal_state() {
    case "$1" in
        COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

all_terminal() {
    local state
    for state in "$@"; do
        if ! terminal_state "$state"; then
            return 1
        fi
    done
    return 0
}

while true; do
    if [[ "$ONCE" -eq 0 ]]; then
        if [[ -n "${TERM:-}" ]]; then
            clear 2>/dev/null || printf '\033[2J\033[H'
        else
            printf '\033[2J\033[H'
        fi
    fi

    printf "RC3 Smoke Monitor\n\n"
    printf "%-6s %-6s %-10s %-8s %-13s %-8s %-8s %-11s %-20s\n" \
        "Label" "Job" "State" "Phase" "SimWeek" "Agents" "Snaps" "ETA" "Finish (UTC)"
    printf "%-6s %-6s %-10s %-8s %-13s %-8s %-8s %-11s %-20s\n" \
        "------" "------" "----------" "--------" "-------------" "--------" "--------" "-----------" "--------------------"

    overall_eta=-1
    overall_finish="unknown"
    states=()

    for spec in "${SPECS[@]}"; do
        IFS='|' read -r label job_id run_dir <<< "$spec"
        out_dir="${run_dir}/output"
        out_log="${run_dir}/slurm_${job_id}.out"

        state=$(job_state "$job_id")
        states+=("$state")
        IFS=$'\t' read -r sim_time max_time wall_seconds agents phase <<< "$(parse_progress "$out_log")"

        snaps=$(find "$out_dir" -maxdepth 1 -name 'output*.xml' 2>/dev/null | wc -l | tr -d ' ')
        sim_week=$(python3 - "$sim_time" "$max_time" <<'PY'
import sys
sim = float(sys.argv[1]) / 10080.0
max_week = float(sys.argv[2]) / 10080.0
print(f"{sim:4.2f}/{max_week:4.2f}")
PY
)

        eta_seconds=$(python3 - "$sim_time" "$max_time" "$wall_seconds" <<'PY'
import sys
sim_time = float(sys.argv[1])
max_time = float(sys.argv[2])
wall_seconds = float(sys.argv[3])
if sim_time <= 0 or wall_seconds <= 0 or max_time <= 0:
    print(-1)
else:
    total_est = wall_seconds * (max_time / sim_time)
    print(max(0.0, total_est - wall_seconds))
PY
)

        if python3 - "$eta_seconds" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= 0 else 1)
PY
        then
            eta_fmt=$(format_seconds "$eta_seconds")
            finish_ts=$(python3 - "$eta_seconds" <<'PY'
import sys
import time
print(int(time.time() + float(sys.argv[1])))
PY
)
            finish_at=$(date -u -d "@${finish_ts}" '+%Y-%m-%d %H:%M:%S')
            if python3 - "$eta_seconds" "$overall_eta" <<'PY'
import sys
curr = float(sys.argv[1])
best = float(sys.argv[2])
if best < 0 or curr > best:
    raise SystemExit(0)
raise SystemExit(1)
PY
            then
                overall_eta="$eta_seconds"
                overall_finish="${finish_at} UTC"
            fi
        else
            eta_fmt="estimating"
            finish_at="unknown"
        fi

        printf "%-6s %-6s %-10s %-8s %-13s %-8s %-8s %-11s %-20s\n" \
            "$label" "$job_id" "$state" "$phase" "$sim_week" "$agents" "$snaps" "$eta_fmt" "$finish_at"
        printf "  run: %s\n" "$run_dir"
        printf "  log: %s\n" "$out_log"
    done

    printf "\n"
    if python3 - "$overall_eta" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= 0 else 1)
PY
    then
        printf "Batch ETA:      %s\n" "$(format_seconds "$overall_eta")"
        printf "Batch Finish:   %s\n" "$overall_finish"
    else
        printf "Batch ETA:      estimating\n"
        printf "Batch Finish:   unknown\n"
    fi

    if all_terminal "${states[@]}"; then
        printf "Smoke Status:   all jobs reached terminal states\n"
        exit 0
    fi

    if [[ "$ONCE" -eq 1 ]]; then
        exit 0
    fi

    sleep "$POLL_SECONDS"
done
