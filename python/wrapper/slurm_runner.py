from __future__ import annotations

import logging
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    from .physicell_runner import PhysiCellRunner, SimulationResult
except Exception:  # pragma: no cover
    from python.wrapper.physicell_runner import PhysiCellRunner, SimulationResult  # type: ignore


LOGGER = logging.getLogger(__name__)


@dataclass
class SlurmResources:
    partition: str = "gpu"
    nodes: int = 1
    ntasks_per_node: int = 1
    cpus_per_task: int = 32
    mem: str = "32G"
    time_limit: str = "02:00:00"
    account: Optional[str] = None
    qos: Optional[str] = None


@dataclass
class _ActiveJob:
    index: int
    run_dir: Path
    output_path: Path
    start_time: float


class SlurmRunner(PhysiCellRunner):
    """
    SLURM-native simulation runner for multi-node distribution on the
    University of Louisville Zurada system.

    This runner submits each PhysiCell simulation as a separate SLURM job:
    1) write per-simulation run script
    2) submit via sbatch
    3) poll sacct for completion
    4) collect SimulationResult objects
    """

    TERMINAL_STATES = {
        "COMPLETED",
        "FAILED",
        "CANCELLED",
        "TIMEOUT",
        "OUT_OF_MEMORY",
        "NODE_FAIL",
        "PREEMPTED",
        "BOOT_FAIL",
        "DEADLINE",
    }

    DEFAULT_ZURADA_MODULES = [
        "module purge || true",
        "module load gcc/12",
        "module load python/3.10",
    ]

    def __init__(
        self,
        binary_path: Path | str,
        config_path: Path | str,
        output_dir: Path | str,
        timeout_seconds: int = 3600,
        resources: Optional[SlurmResources] = None,
        module_commands: Optional[List[str]] = None,
        project_root: Optional[Path | str] = None,
        poll_interval_seconds: int = 15,
        sbatch_args: Optional[List[str]] = None,
    ):
        super().__init__(
            binary_path=binary_path,
            config_path=config_path,
            output_dir=output_dir,
            timeout_seconds=timeout_seconds,
        )
        self.resources = resources or SlurmResources()
        self.module_commands = list(module_commands or self.DEFAULT_ZURADA_MODULES)
        self.project_root = (
            Path(project_root).expanduser().resolve() if project_root is not None else Path.cwd().resolve()
        )
        self.slurm_poll_interval_seconds = max(1, int(poll_interval_seconds))
        self.default_sbatch_args = list(sbatch_args or [])

    def run_batch(
        self,
        intervention_jsons: list,
        max_parallel: int = 50,
        sbatch_args: Optional[List[str]] = None,
    ) -> List[SimulationResult]:
        return self.run_batch_slurm(
            intervention_jsons=intervention_jsons,
            max_parallel=max_parallel,
            sbatch_args=sbatch_args,
            poll_interval_seconds=self.slurm_poll_interval_seconds,
        )

    def run_batch_slurm(
        self,
        intervention_jsons: list,
        max_parallel: int = 50,
        sbatch_args: Optional[List[str]] = None,
        poll_interval_seconds: Optional[int] = None,
    ) -> List[SimulationResult]:
        if max_parallel < 1:
            raise ValueError("max_parallel must be >= 1")

        poll_seconds = max(1, int(poll_interval_seconds or self.slurm_poll_interval_seconds))
        merged_sbatch_args = [*self.default_sbatch_args, *(sbatch_args or [])]

        prepared: List[Dict[str, Any]] = []
        for idx, intervention_json in enumerate(intervention_jsons):
            item = self._prepare_slurm_run(intervention_json)
            item["index"] = idx
            prepared.append(item)

        results: List[Optional[SimulationResult]] = [None] * len(prepared)
        pending = list(prepared)
        active: Dict[str, _ActiveJob] = {}

        while pending or active:
            while pending and len(active) < max_parallel:
                item = pending.pop(0)
                idx = int(item["index"])
                run_dir = Path(item["run_dir"])
                script_path = self._write_zurada_job_script(
                    run_dir=run_dir,
                    config_path=Path(item["config_path"]),
                    intervention_path=Path(item["intervention_path"]),
                    output_path=Path(item["output_path"]),
                )
                job_id = self._submit_sbatch(script_path, merged_sbatch_args)
                if job_id is None:
                    result = SimulationResult(
                        run_dir=run_dir,
                        exit_code=-1,
                        wall_time=0.0,
                        success=False,
                        output_files={},
                    )
                    if self.cleanup_failed_runs:
                        self._cleanup_failed_run(result, reason="sbatch_submit_failed")
                    results[idx] = result
                    continue

                active[job_id] = _ActiveJob(
                    index=idx,
                    run_dir=run_dir,
                    output_path=Path(item["output_path"]),
                    start_time=time.perf_counter(),
                )
                LOGGER.info("Submitted simulation job %s for index=%d (%s)", job_id, idx, run_dir)

            if not active:
                continue

            time.sleep(poll_seconds)

            for job_id in list(active.keys()):
                job = active[job_id]
                elapsed = time.perf_counter() - job.start_time

                if elapsed > self.timeout_seconds:
                    subprocess.run(["scancel", job_id], capture_output=True, text=True, check=False)
                    result = SimulationResult(
                        run_dir=job.run_dir,
                        exit_code=-1,
                        wall_time=elapsed,
                        success=False,
                        output_files={},
                    )
                    if self.cleanup_failed_runs:
                        self._cleanup_failed_run(result, reason="slurm_timeout")
                    results[job.index] = result
                    del active[job_id]
                    continue

                state, exit_code, terminal = self._query_sacct_state(job_id)
                if not terminal:
                    continue

                output_files = self._collect_output_files(job.output_path)
                success = (state == "COMPLETED") and (exit_code == 0) and self._output_has_files(job.output_path)
                result = SimulationResult(
                    run_dir=job.run_dir,
                    exit_code=exit_code,
                    wall_time=elapsed,
                    success=success,
                    output_files=output_files,
                )
                if not success and self.cleanup_failed_runs:
                    self._cleanup_failed_run(result, reason=f"slurm_{state.lower()}")
                results[job.index] = result
                del active[job_id]

        return [r for r in results if r is not None]

    def _submit_sbatch(self, script_path: Path, sbatch_args: List[str]) -> Optional[str]:
        cmd = ["sbatch", "--parsable", *sbatch_args, str(script_path)]
        submit = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if submit.returncode != 0:
            LOGGER.error("sbatch failed for %s: %s", script_path, submit.stderr.strip())
            return None
        job_id = submit.stdout.strip().split(";")[0]
        if not job_id:
            LOGGER.error("Unable to parse SLURM job id from sbatch output: %s", submit.stdout.strip())
            return None
        return job_id

    def _write_zurada_job_script(
        self,
        run_dir: Path,
        config_path: Path,
        intervention_path: Path,
        output_path: Path,
    ) -> Path:
        script_path = run_dir / "run.slurm.sh"
        log_out = run_dir / "slurm_%j.out"
        log_err = run_dir / "slurm_%j.err"

        sbatch_lines = [
            "#!/bin/bash",
            "#SBATCH --job-name=stroma_world_sim",
            f"#SBATCH --partition={self.resources.partition}",
            f"#SBATCH --nodes={self.resources.nodes}",
            f"#SBATCH --ntasks-per-node={self.resources.ntasks_per_node}",
            f"#SBATCH --cpus-per-task={self.resources.cpus_per_task}",
            f"#SBATCH --mem={self.resources.mem}",
            f"#SBATCH --time={self.resources.time_limit}",
            f"#SBATCH --output={log_out}",
            f"#SBATCH --error={log_err}",
        ]
        if self.resources.account:
            sbatch_lines.append(f"#SBATCH --account={self.resources.account}")
        if self.resources.qos:
            sbatch_lines.append(f"#SBATCH --qos={self.resources.qos}")

        primary_cmd = " ".join(
            [
                shlex.quote(str(self.binary_path)),
                "--config",
                shlex.quote(str(config_path)),
                "--intervention",
                shlex.quote(str(intervention_path)),
                "--output",
                shlex.quote(str(output_path)),
            ]
        )
        legacy_cmd = " ".join(
            [
                shlex.quote(str(self.binary_path)),
                shlex.quote(str(config_path)),
                shlex.quote(str(intervention_path)),
            ]
        )

        body_lines = [
            "set -euo pipefail",
            "source /etc/profile >/dev/null 2>&1 || true",
            f'export ZURADA_PROJECT_ROOT="{self.project_root}"',
            'export PYTHONPATH="$ZURADA_PROJECT_ROOT:${PYTHONPATH:-}"',
            'export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"',
            "if command -v module >/dev/null 2>&1; then",
        ]
        for line in self.module_commands:
            body_lines.append(f"  {line}")
        body_lines.extend(
            [
                "fi",
                'cd "$ZURADA_PROJECT_ROOT"',
                f"{primary_cmd} || {legacy_cmd}",
                "",
            ]
        )

        script_path.write_text("\n".join([*sbatch_lines, *body_lines]), encoding="utf-8")
        script_path.chmod(0o755)
        return script_path

    def _query_sacct_state(self, job_id: str) -> Tuple[str, int, bool]:
        cmd = ["sacct", "-j", job_id, "--format=JobIDRaw,State,ExitCode", "-n", "-P"]
        probe = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if probe.returncode == 0 and probe.stdout.strip():
            state, exit_code = self._parse_sacct_output(job_id, probe.stdout)
            if state:
                return state, exit_code, state in self.TERMINAL_STATES

        # Fallback: scheduler still knows about the job but sacct is lagging.
        sq = subprocess.run(["squeue", "-j", job_id, "-h", "-o", "%T"], capture_output=True, text=True, check=False)
        if sq.returncode == 0 and sq.stdout.strip():
            return sq.stdout.strip().splitlines()[0], -1, False
        return "UNKNOWN", -1, True

    def _parse_sacct_output(self, job_id: str, output: str) -> Tuple[str, int]:
        # Prefer exact JobIDRaw match, then fallback to job steps.
        exact_matches: List[Tuple[str, int]] = []
        step_matches: List[Tuple[str, int]] = []
        for line in output.strip().splitlines():
            fields = line.split("|")
            if len(fields) < 3:
                continue
            jid = fields[0].strip()
            state_raw = fields[1].strip()
            exit_raw = fields[2].strip()
            if not jid:
                continue

            state = state_raw.split()[0].split("+")[0]
            code = self._parse_exit_code(exit_raw)
            if jid == job_id:
                exact_matches.append((state, code))
            elif jid.startswith(f"{job_id}."):
                step_matches.append((state, code))

        if exact_matches:
            return exact_matches[0]
        if step_matches:
            return step_matches[0]
        return "", -1

    @staticmethod
    def _parse_exit_code(raw: str) -> int:
        if ":" in raw:
            raw = raw.split(":", 1)[0]
        try:
            return int(raw)
        except ValueError:
            return -1
