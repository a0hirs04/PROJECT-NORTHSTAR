#!/bin/bash
# Live monitor for the 5-seed RC2 gold benchmark batch.
# Usage:
#   bash watch_rc2_gold_benchmark.sh
#   bash watch_rc2_gold_benchmark.sh --once

set -euo pipefail

ROOT_DIR="${RC2_GOLD_ROOT:-/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_gold_benchmark_5seed}"
POLL_SECONDS=15
ONCE=0

if [[ "${1:-}" == "--once" ]]; then
    ONCE=1
fi

format_seconds() {
    python3 - "$1" <<'PY'
import sys
secs = max(0.0, float(sys.argv[1]))
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

job_state() {
    local job_id="$1"
    local state
    state=$(squeue -j "$job_id" --format="%T" --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [[ -n "$state" ]]; then
        printf "%s" "$state"
        return
    fi
    state=$(sacct -j "$job_id" --format=State --noheader --parsable2 2>/dev/null | head -1 | cut -d'|' -f1 | tr -d ' ')
    printf "%s" "${state:-UNKNOWN}"
}

job_elapsed() {
    local job_id="$1"
    local elapsed
    elapsed=$(squeue -j "$job_id" --format="%M" --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [[ -n "$elapsed" ]]; then
        printf "%s" "$elapsed"
        return
    fi
    elapsed=$(sacct -j "$job_id" --format=Elapsed --noheader --parsable2 2>/dev/null | head -1 | cut -d'|' -f1 | tr -d ' ')
    printf "%s" "${elapsed:-unknown}"
}

latest_job_id() {
    local run_dir="$1"
    local latest
    latest=$(find "$run_dir" -maxdepth 1 -type f -name 'slurm_*.out' -printf '%T@ %f\n' 2>/dev/null | sort -nr | head -1 | awk '{print $2}')
    if [[ -z "${latest:-}" ]]; then
        return 1
    fi
    latest="${latest#slurm_}"
    latest="${latest%.out}"
    printf "%s" "$latest"
}

job_name_for_run_dir() {
    local run_dir="$1"
    awk -F= '/^#SBATCH --job-name=/{print $2; exit}' "$run_dir/run.slurm.sh" 2>/dev/null
}

job_id_from_name() {
    local job_name="$1"
    sacct -S 2026-03-29T00:00:00 -P --format=JobID,JobName,State 2>/dev/null | \
        awk -F'|' -v target="$job_name" '
            $2 == target && $1 !~ /\./ {
                id = $1 + 0
                if (id > best) {
                    best = id
                }
            }
            END {
                if (best > 0) {
                    print best
                }
            }
        '
}

parse_progress() {
    python3 - "$1" <<'PY'
import os
import re
import sys

path = sys.argv[1]
if not os.path.exists(path):
    print("0\t60480\t0\t0\tSTARTUP")
    raise SystemExit(0)

sim_re = re.compile(r"current simulated time:\s*([0-9.eE+-]+)\s*min\s*\(max:\s*([0-9.eE+-]+)\s*min\)")
agents_re = re.compile(r"total agents:\s*(\d+)")
wall_re = re.compile(r"total wall time:\s*(\d+)\s*days,\s*(\d+)\s*hours,\s*(\d+)\s*minutes,\s*and\s*([0-9.eE+-]+)\s*seconds")

sim_time = 0.0
max_time = 60480.0
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
    phase = "BARRIER"
elif sim_time < 40320.0:
    phase = "DRUG-ON"
elif sim_time < max_time:
    phase = "REGROWTH"
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

discover_specs() {
    [[ -d "$ROOT_DIR" ]] || return 0
    find "$ROOT_DIR" -mindepth 1 -maxdepth 1 -type d -name 'replicate_*_seed*' | sort | while read -r run_dir; do
        local seed label job_id job_name
        seed=$(basename "$run_dir")
        label="${seed#replicate_}"
        job_id=$(latest_job_id "$run_dir" || true)
        if [[ -z "${job_id:-}" ]]; then
            job_name=$(job_name_for_run_dir "$run_dir" || true)
            [[ -n "${job_name:-}" ]] && job_id=$(job_id_from_name "$job_name" || true)
        fi
        [[ -z "${job_id:-}" ]] && continue
        printf "%s|%s|%s\n" "$label" "$job_id" "$run_dir"
    done
}

while true; do
    mapfile -t SPECS < <(discover_specs)

    if [[ "$ONCE" -eq 0 ]]; then
        if [[ -n "${TERM:-}" ]]; then
            clear 2>/dev/null || printf '\033[2J\033[H'
        else
            printf '\033[2J\033[H'
        fi
    fi

    printf "RC2 Gold Benchmark Monitor\n"
    printf "Root: %s\n\n" "$ROOT_DIR"
    printf "%-13s %-6s %-10s %-10s %-13s %-8s %-8s %-9s %-11s %-20s\n" \
        "Replicate" "Job" "State" "Phase" "SimDay" "Agents" "Snaps" "%Done" "ETA" "Finish (UTC)"
    printf "%-13s %-6s %-10s %-10s %-13s %-8s %-8s %-9s %-11s %-20s\n" \
        "-------------" "------" "----------" "----------" "-------------" "--------" "--------" "---------" "-----------" "--------------------"

    running_etas=()
    completed_count=0
    states=()
    progress_values=()

    for spec in "${SPECS[@]}"; do
        IFS='|' read -r label job_id run_dir <<< "$spec"
        out_dir="${run_dir}/output"
        out_log="${run_dir}/slurm_${job_id}.out"
        state=$(job_state "$job_id")
        states+=("$state")
        IFS=$'\t' read -r sim_time max_time wall_seconds agents phase <<< "$(parse_progress "$out_log")"

        snaps=$(find "$out_dir" -maxdepth 1 -name 'output*.xml' 2>/dev/null | wc -l | tr -d ' ')
        sim_day=$(python3 - "$sim_time" "$max_time" <<'PY'
import sys
sim = float(sys.argv[1]) / 1440.0
max_day = float(sys.argv[2]) / 1440.0
print(f"{sim:5.2f}/{max_day:5.2f}")
PY
)

        pct_done=$(python3 - "$sim_time" "$max_time" "$state" <<'PY'
import sys
sim_time = float(sys.argv[1])
max_time = float(sys.argv[2])
state = sys.argv[3]
terminal = {"COMPLETED","FAILED","CANCELLED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE"}
if state in terminal:
    pct = 100.0
elif max_time <= 0:
    pct = 0.0
else:
    pct = max(0.0, min(100.0, 100.0 * sim_time / max_time))
print(f"{pct:.1f}%")
print(pct)
PY
)
        pct_fmt=$(printf "%s" "$pct_done" | sed -n '1p')
        pct_value=$(printf "%s" "$pct_done" | sed -n '2p')
        progress_values+=("$pct_value")

        eta_seconds=$(python3 - "$sim_time" "$max_time" "$wall_seconds" "$state" <<'PY'
import sys
sim_time = float(sys.argv[1])
max_time = float(sys.argv[2])
wall_seconds = float(sys.argv[3])
state = sys.argv[4]
terminal = {"COMPLETED","FAILED","CANCELLED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE"}
if state in terminal:
    print(0.0)
elif sim_time <= 0 or wall_seconds <= 0 or max_time <= 0:
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
        else
            eta_fmt="estimating"
            finish_at="unknown"
        fi

        if terminal_state "$state"; then
            completed_count=$((completed_count + 1))
        else
            running_etas+=("$eta_seconds")
        fi

        printf "%-13s %-6s %-10s %-10s %-13s %-8s %-8s %-9s %-11s %-20s\n" \
            "$label" "$job_id" "$state" "$phase" "$sim_day" "$agents" "$snaps" "$pct_fmt" "$eta_fmt" "${finish_at} UTC"
    done

    batch_progress=$(python3 - "${progress_values[@]:-0}" <<'PY'
import sys
vals = [float(v) for v in sys.argv[1:] if v.strip()]
if not vals:
    print("0.0%")
else:
    print(f"{sum(vals)/len(vals):.1f}%")
PY
)

    overall_eta=$(python3 - "${running_etas[@]:--1}" <<'PY'
import sys
vals = [float(v) for v in sys.argv[1:] if float(v) >= 0]
if not vals:
    print(-1)
else:
    print(max(vals))
PY
)

    if python3 - "$overall_eta" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= 0 else 1)
PY
    then
        overall_eta_fmt=$(format_seconds "$overall_eta")
        overall_finish_ts=$(python3 - "$overall_eta" <<'PY'
import sys
import time
print(int(time.time() + float(sys.argv[1])))
PY
)
        overall_finish=$(date -u -d "@${overall_finish_ts}" '+%Y-%m-%d %H:%M:%S UTC')
    else
        overall_eta_fmt="estimating"
        overall_finish="unknown"
    fi

    printf "\nCompleted: %d / %d\n" "$completed_count" "${#SPECS[@]}"
    printf "Batch Progress: %s\n" "$batch_progress"
    printf "Batch ETA: %s\n" "$overall_eta_fmt"
    printf "Projected Finish: %s\n" "$overall_finish"

    all_terminal=1
    for state in "${states[@]}"; do
        if ! terminal_state "$state"; then
            all_terminal=0
            break
        fi
    done

    if [[ "$ONCE" -eq 1 || "$all_terminal" -eq 1 ]]; then
        break
    fi

    sleep "$POLL_SECONDS"
done
