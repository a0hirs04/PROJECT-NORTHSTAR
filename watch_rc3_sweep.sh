#!/bin/bash
# Live monitor for the current RC3 full sweep.
# Usage:
#   bash watch_rc3_sweep.sh
#   bash watch_rc3_sweep.sh --once

set -euo pipefail

ROOT_DIR_OVERRIDE="${RC3_SWEEP_ROOT:-}"
ROOT_CANDIDATES=(
    "/home/a0hirs04/PROJECT-NORTHSTAR/build/rc3_vismodegib"
    "/work/a0hirs04/PROJECT-NORTHSTAR/build/rc3_vismodegib"
)
POLL_SECONDS=15
ONCE=0

if [[ "${1:-}" == "--once" ]]; then
    ONCE=1
fi

resolve_root_dir() {
    if [[ -n "$ROOT_DIR_OVERRIDE" ]]; then
        printf "%s" "$ROOT_DIR_OVERRIDE"
        return 0
    fi

    local best_root=""
    local best_score=-1
    local root score logs outputs
    for root in "${ROOT_CANDIDATES[@]}"; do
        [[ -d "$root" ]] || continue
        logs=$(find "$root" -maxdepth 3 -type f -name 'slurm_*.out' 2>/dev/null | wc -l | tr -d ' ')
        outputs=$(find "$root" -maxdepth 4 -type f -path '*/output/output*.xml' 2>/dev/null | wc -l | tr -d ' ')
        score=$(( logs * 1000 + outputs ))
        if (( score > best_score )); then
            best_root="$root"
            best_score=$score
        fi
    done

    if [[ -n "$best_root" ]]; then
        printf "%s" "$best_root"
        return 0
    fi

    printf "%s" "${ROOT_CANDIDATES[0]}"
}

ROOT_DIR="$(resolve_root_dir)"

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
    sacct -S 2026-03-28T00:00:00 -P --format=JobID,JobName,State 2>/dev/null | \
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
    print("0\t80640\t0\t0\tSTARTUP")
    raise SystemExit(0)

sim_re = re.compile(r"current simulated time:\s*([0-9.eE+-]+)\s*min\s*\(max:\s*([0-9.eE+-]+)\s*min\)")
agents_re = re.compile(r"total agents:\s*(\d+)")
wall_re = re.compile(r"total wall time:\s*(\d+)\s*days,\s*(\d+)\s*hours,\s*(\d+)\s*minutes,\s*and\s*([0-9.eE+-]+)\s*seconds")

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
elif sim_time < 40320.0:
    phase = "Wk4-6"
elif sim_time < max_time:
    phase = "Wk6-8"
else:
    phase = "DONE"

print(f"{sim_time}\t{max_time}\t{wall_seconds}\t{agents}\t{phase}")
PY
}

discover_specs() {
    [[ -d "$ROOT_DIR" ]] || return 0
    find "$ROOT_DIR" -mindepth 2 -maxdepth 2 -type d | sort | while read -r run_dir; do
        local arm seed label job_id job_name
        arm=$(basename "$(dirname "$run_dir")")
        seed=$(basename "$run_dir")
        job_id=$(latest_job_id "$run_dir" || true)
        if [[ -z "${job_id:-}" ]]; then
            job_name=$(job_name_for_run_dir "$run_dir" || true)
            [[ -n "${job_name:-}" ]] && job_id=$(job_id_from_name "$job_name" || true)
        fi
        [[ -z "${job_id:-}" ]] && continue
        label="${arm#Arm_}/${seed#seed_}"
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

    printf "RC3 Full Sweep Monitor\n"
    printf "Root: %s\n\n" "$ROOT_DIR"
    printf "%-18s %-6s %-10s %-8s %-13s %-8s %-8s %-9s %-11s %-20s\n" \
        "Label" "Job" "State" "Phase" "SimWeek" "Agents" "Snaps" "%Done" "ETA" "Finish (UTC)"
    printf "%-18s %-6s %-10s %-8s %-13s %-8s %-8s %-9s %-11s %-20s\n" \
        "------------------" "------" "----------" "--------" "-------------" "--------" "--------" "---------" "-----------" "--------------------"

    running_etas=()
    completed_count=0
    pending_count=0
    states=()
    mean_runtime_guess=0
    batch_progress_values=()

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
print(f"{pct:6.1f}%")
PY
)
        eta_seconds=$(python3 - "$sim_time" "$max_time" "$wall_seconds" "$state" <<'PY'
import sys
sim_time = float(sys.argv[1])
max_time = float(sys.argv[2])
wall_seconds = float(sys.argv[3])
state = sys.argv[4]
terminal = {"COMPLETED","FAILED","CANCELLED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE"}
if state in terminal:
    print(0)
elif sim_time <= 0 or wall_seconds <= 0 or max_time <= 0:
    print(-1)
else:
    total_est = wall_seconds * (max_time / sim_time)
    print(max(0.0, total_est - wall_seconds))
PY
)
        total_runtime_est=$(python3 - "$sim_time" "$max_time" "$wall_seconds" <<'PY'
import sys
sim_time = float(sys.argv[1])
max_time = float(sys.argv[2])
wall_seconds = float(sys.argv[3])
if sim_time <= 0 or wall_seconds <= 0 or max_time <= 0:
    print(-1)
else:
    print(wall_seconds * (max_time / sim_time))
PY
)

        if [[ "$state" == "RUNNING" ]] && python3 - "$eta_seconds" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= 0 else 1)
PY
        then
            running_etas+=("$eta_seconds")
            if python3 - "$total_runtime_est" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) > 0 else 1)
PY
            then
                mean_runtime_guess=$(python3 - "$mean_runtime_guess" "$total_runtime_est" <<'PY'
import sys
curr = float(sys.argv[1]); new = float(sys.argv[2])
print(new if curr <= 0 else (curr + new) / 2.0)
PY
)
            fi
        fi

        if [[ "$state" == "PENDING" ]]; then
            pending_count=$((pending_count + 1))
        fi
        if [[ "$state" == "COMPLETED" ]]; then
            completed_count=$((completed_count + 1))
        fi

        progress_value=$(python3 - "$sim_time" "$max_time" "$state" <<'PY'
import sys
sim_time = float(sys.argv[1])
max_time = float(sys.argv[2])
state = sys.argv[3]
terminal = {"COMPLETED","FAILED","CANCELLED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE"}
if state in terminal:
    print(100.0)
elif max_time <= 0:
    print(0.0)
else:
    print(max(0.0, min(100.0, 100.0 * sim_time / max_time)))
PY
)
        batch_progress_values+=("$progress_value")

        if python3 - "$eta_seconds" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= 0 else 1)
PY
        then
            eta_fmt=$(format_seconds "$eta_seconds")
            finish_ts=$(python3 - "$eta_seconds" <<'PY'
import sys, time
print(int(time.time() + float(sys.argv[1])))
PY
)
            finish_at=$(date -u -d "@${finish_ts}" '+%Y-%m-%d %H:%M:%S')
        else
            eta_fmt="unknown"
            finish_at="unknown"
        fi

        printf "%-18s %-6s %-10s %-8s %-13s %-8s %-8s %-9s %-11s %-20s\n" \
            "$label" "$job_id" "$state" "$phase" "$sim_week" "$agents" "$snaps" "$pct_done" "$eta_fmt" "$finish_at"
    done

    batch_eta=$(python3 - "$pending_count" "$mean_runtime_guess" "${running_etas[*]:-}" <<'PY'
import math
import sys
pending = int(float(sys.argv[1]))
mean_runtime = float(sys.argv[2])
etas = [float(x) for x in sys.argv[3].split() if x.strip()]
if not etas:
    print(-1)
    raise SystemExit(0)
running_max = max(etas)
if pending <= 0 or mean_runtime <= 0:
    print(running_max)
    raise SystemExit(0)
waves = math.ceil(pending / 8.0)
print(running_max + waves * mean_runtime)
PY
)

    if python3 - "$batch_eta" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= 0 else 1)
PY
    then
        batch_finish_ts=$(python3 - "$batch_eta" <<'PY'
import sys, time
print(int(time.time() + float(sys.argv[1])))
PY
)
        batch_finish=$(date -u -d "@${batch_finish_ts}" '+%Y-%m-%d %H:%M:%S UTC')
        printf "\nBatch ETA:    %s\n" "$(format_seconds "$batch_eta")"
        printf "Batch Finish: %s\n" "$batch_finish"
    else
        printf "\nBatch ETA:    unknown\n"
        printf "Batch Finish: unknown\n"
    fi

    batch_pct=$(python3 - "${batch_progress_values[*]:-}" <<'PY'
import sys
vals = [float(x) for x in sys.argv[1].split() if x.strip()]
if not vals:
    print("0.0%")
else:
    print(f"{sum(vals)/len(vals):.1f}%")
PY
)

    printf "Completed:    %d / %d\n" "$completed_count" "${#SPECS[@]}"
    printf "Pending:      %d\n" "$pending_count"
    printf "Progress:     %s\n" "$batch_pct"

    if [[ "$ONCE" -eq 1 ]]; then
        break
    fi

    sleep "$POLL_SECONDS"
done
