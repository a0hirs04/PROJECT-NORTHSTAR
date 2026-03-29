#!/usr/bin/env python3
"""Summarize a clean 5-seed RC2 gold benchmark set.

This leaves `evaluate_rc2.py` untouched and reuses its parsing helpers so the
benchmark summary is consistent with the main RC2 evaluator.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import re
import statistics
import xml.etree.ElementTree as ET
from pathlib import Path

from evaluate_rc2 import OutputParser, _criteria, _load_snapshot, _parse_snapshot, _report_one


EXPECTED_SEEDS = [42, 43, 44, 45, 46]
DEFAULT_ROOT = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_gold_benchmark_5seed")


def _seed_from_name(path: Path) -> int:
    match = re.search(r"seed(\d+)", path.name)
    if not match:
        raise ValueError(f"Could not infer seed from directory name: {path}")
    return int(match.group(1))


def _read_float(root: ET.Element, xpath: str) -> float:
    node = root.find(xpath)
    if node is None or node.text is None:
        raise ValueError(f"Missing XML node for xpath: {xpath}")
    return float(node.text)


def _timing_from_config(config_path: Path) -> dict[str, float]:
    root = ET.parse(config_path).getroot()
    return {
        "drug_start_time": _read_float(root, ".//user_parameters/drug_start_time"),
        "drug_end_time": _read_float(root, ".//user_parameters/drug_end_time"),
        "max_time": _read_float(root, ".//overall/max_time"),
    }


def _fmt_float(value: float, ndigits: int = 4) -> str:
    if value is None or not math.isfinite(value):
        return "nan"
    return f"{value:.{ndigits}f}"


def _discover_runs(root_dir: Path) -> list[tuple[int, Path, Path]]:
    runs = []
    for rep_dir in sorted(root_dir.glob("replicate_*_seed*")):
        config_path = rep_dir / "config.xml"
        out_dir = rep_dir / "output"
        if config_path.exists() and out_dir.is_dir():
            runs.append((_seed_from_name(rep_dir), rep_dir, out_dir))
    if not runs:
        raise FileNotFoundError(f"No RC2 benchmark runs found under {root_dir}")
    found = sorted(seed for seed, _, _ in runs)
    if found != EXPECTED_SEEDS:
        raise ValueError(f"Expected seeds {EXPECTED_SEEDS}, found {found}")
    return runs


def summarize_run(seed: int, rep_dir: Path, out_dir: Path) -> dict[str, object]:
    timing = _timing_from_config(rep_dir / "config.xml")
    parser = OutputParser(out_dir)

    pre_snap, _ = _load_snapshot(parser, out_dir, timing["drug_start_time"])
    treat_snap, _ = _load_snapshot(parser, out_dir, timing["drug_end_time"])
    post_snap, _ = _load_snapshot(parser, out_dir, timing["max_time"])

    pre = _parse_snapshot(parser, pre_snap)
    treat = _parse_snapshot(parser, treat_snap)
    post = _parse_snapshot(parser, post_snap)
    criteria = _criteria(pre, treat, post)

    with contextlib.redirect_stdout(io.StringIO()):
        hard_score, verdict = _report_one(f"seed_{seed}", out_dir, timing, fast_mode=False)

    return {
        "seed": seed,
        "day14_live": pre["n_live_tumor"],
        "day28_live": treat["n_live_tumor"],
        "day42_live": post["n_live_tumor"],
        "peri_d14": pre["peri_ecm"],
        "peri_d28": treat["peri_ecm"],
        "ecm_surv_d14": pre["ecm_at_survivors"],
        "ecm_surv_d28": treat["ecm_at_survivors"],
        "criteria": criteria,
        "hard_score": hard_score,
        "verdict": verdict,
    }


def render_summary(root_dir: Path, rows: list[dict[str, object]]) -> str:
    day28_values = [int(row["day28_live"]) for row in rows]
    median_day28 = statistics.median(day28_values)

    lines: list[str] = []
    lines.append("RC2 GOLD BENCHMARK SUMMARY")
    lines.append(f"Root: {root_dir}")
    lines.append("")
    lines.append(
        "Seed  d14_live  d28_live  d42_live  peri_d14  peri_d28  ecm_surv_d14  "
        "ecm_surv_d28  RC2-1  RC2-2  RC2-3  RC2-4  RC2-5  RC2-6  Verdict"
    )

    for row in rows:
        criteria = row["criteria"]
        assert isinstance(criteria, dict)
        lines.append(
            f"{int(row['seed']):<4}  "
            f"{int(row['day14_live']):>8}  "
            f"{int(row['day28_live']):>8}  "
            f"{int(row['day42_live']):>8}  "
            f"{_fmt_float(float(row['peri_d14'])):>8}  "
            f"{_fmt_float(float(row['peri_d28'])):>8}  "
            f"{_fmt_float(float(row['ecm_surv_d14'])):>12}  "
            f"{_fmt_float(float(row['ecm_surv_d28'])):>12}  "
            f"{criteria['RC2-1']['status']:>5}  "
            f"{criteria['RC2-2']['status']:>5}  "
            f"{criteria['RC2-3']['status']:>5}  "
            f"{criteria['RC2-4']['status']:>5}  "
            f"{criteria['RC2-5']['status']:>5}  "
            f"{criteria['RC2-6']['status']:>5}  "
            f"{str(row['verdict'])}"
        )

    lines.append("")
    lines.append(f"Benchmark median day-28 live tumor count: {median_day28}")
    lines.append("")
    lines.append("Benchmark recommendation for final RC3 C3.5:")
    lines.append(f"--rc2-benchmark {root_dir}")
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize a clean 5-seed RC2 gold benchmark set.")
    parser.add_argument("--root-dir", type=Path, default=DEFAULT_ROOT, help="Benchmark root directory.")
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional summary output path. Defaults to <root-dir>/benchmark_summary.txt",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root_dir = args.root_dir.expanduser().resolve()
    summary_out = (
        args.summary_out.expanduser().resolve()
        if args.summary_out is not None
        else (root_dir / "benchmark_summary.txt").resolve()
    )

    try:
        runs = _discover_runs(root_dir)
        rows = [summarize_run(seed, rep_dir, out_dir) for seed, rep_dir, out_dir in runs]
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1

    rows.sort(key=lambda row: int(row["seed"]))
    text = render_summary(root_dir, rows)
    summary_out.write_text(text, encoding="utf-8")
    print(text, end="")
    print(f"Saved summary: {summary_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
