"""
PhysiCell wrapper package for Stroma World.

Handles launching PhysiCell simulations, parsing MultiCellDS output,
and generating per-run configuration files from EA individuals.

Modules:
    physicell_runner  — Subprocess management and batch execution (local + SLURM)
    output_parser     — MultiCellDS XML/MAT output parser and metrics extraction
    config_generator  — Per-run XML config and intervention JSON generation
    slurm_runner      — SLURM HPC batch submission (Zurada/LARCC)
"""
