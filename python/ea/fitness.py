from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict

try:
    from ..wrapper.output_parser import SimulationMetrics
except Exception:  # pragma: no cover
    from python.wrapper.output_parser import SimulationMetrics  # type: ignore


@dataclass
class FitnessConfig:
    """
    Tunable fitness weights and normalization constants for Cycle 1.

    Biological rationale:
    - Tumor kill and stromal disruption are the top priorities for this cycle.
    - Drug penetration is rewarded because barrier removal is only useful if drug
      reaches tumor interior.
    - Resistance and invasion are penalized to avoid strategies that appear to
      "open" stroma while selecting for NRF2/ABCB1-mediated drug resistance or
      promoting spatial spread.
    """

    tumor_kill_weight: float = 0.35
    stroma_disruption_weight: float = 0.25
    drug_penetration_weight: float = 0.20
    resistance_penalty_weight: float = -0.10
    invasion_penalty_weight: float = -0.10

    # Normalization constants
    initial_tumor_cells: int = 50
    max_expected_drug: float = 1.0
    domain_size: float = 2000.0

    @classmethod
    def from_dict(cls, payload: Dict[str, Any]) -> "FitnessConfig":
        valid_fields = set(cls.__dataclass_fields__.keys())  # type: ignore[attr-defined]
        filtered = {k: v for k, v in payload.items() if k in valid_fields}
        return cls(**filtered)

    @classmethod
    def from_json(cls, filepath: Path | str) -> "FitnessConfig":
        path = Path(filepath).expanduser().resolve()
        with path.open("r", encoding="utf-8") as f:
            payload = json.load(f)
        if not isinstance(payload, dict):
            raise ValueError(f"Fitness config JSON must be an object: {path}")
        return cls.from_dict(payload)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


DEFAULT_FITNESS_CONFIG = FitnessConfig()
_ACTIVE_FITNESS_CONFIG = FitnessConfig()


def set_fitness_config(config: FitnessConfig) -> None:
    global _ACTIVE_FITNESS_CONFIG
    _ACTIVE_FITNESS_CONFIG = config


def load_fitness_config(filepath: Path | str) -> FitnessConfig:
    config = FitnessConfig.from_json(filepath)
    set_fitness_config(config)
    return config


def get_fitness_config() -> FitnessConfig:
    return _ACTIVE_FITNESS_CONFIG


def validate_metrics(metrics: SimulationMetrics) -> bool:
    """
    Basic sanity checks to reject corrupted / crashed simulation outputs.

    Checks include:
    - Non-negative cell counts
    - Live tumor cells <= total tumor cells
    - Fractions and ECM metrics constrained to biologically valid ranges
    - Non-negative spread/penetration values
    - Finite numeric values
    """

    int_fields = {
        "total_tumor_cells": metrics.total_tumor_cells,
        "total_stromal_cells": metrics.total_stromal_cells,
        "live_tumor_cells": metrics.live_tumor_cells,
        "activated_cafs": metrics.activated_cafs,
        "mesenchymal_tumor_cells": metrics.mesenchymal_tumor_cells,
    }
    if any(v < 0 for v in int_fields.values()):
        return False
    if metrics.live_tumor_cells > metrics.total_tumor_cells:
        return False
    if metrics.activated_cafs > metrics.total_stromal_cells:
        return False
    if metrics.mesenchymal_tumor_cells > metrics.total_tumor_cells:
        return False

    float_ranges = [
        (metrics.mean_ecm_density, 0.0, 1.0),
        (metrics.max_ecm_density, 0.0, 1.0),
        (metrics.stroma_barrier_score, 0.0, 1.0),
        (metrics.mean_tumor_drug_sensitivity, 0.0, 1.0),
        (metrics.hypoxic_fraction, 0.0, 1.0),
    ]
    for value, lo, hi in float_ranges:
        if not math.isfinite(value) or value < lo or value > hi:
            return False

    nonnegative = [metrics.tumor_extent, metrics.drug_penetration]
    if any((not math.isfinite(v)) or v < 0.0 for v in nonnegative):
        return False

    return True


def compute_fitness_detailed(metrics: SimulationMetrics) -> dict:
    """
    Compute decomposed Cycle 1 fitness components for logging/analysis.

    Components:
    1. Tumor kill: direct anti-tumor efficacy.
    2. Stroma disruption: ECM barrier removal.
    3. Drug penetration: delivery into tumor interior.
    4. Resistance penalty: discourages NRF2/ABCB1-associated low sensitivity.
    5. Invasion penalty: discourages spread caused by aggressive stroma removal.
    """

    cfg = get_fitness_config()
    valid = validate_metrics(metrics)
    if not valid:
        return {
            "valid": False,
            "fitness": 0.0,
            "tumor_kill_score": 0.0,
            "stroma_score": 0.0,
            "penetration_score": 0.0,
            "resistance_penalty": 1.0,
            "invasion_penalty": 1.0,
            "weighted_tumor_kill": 0.0,
            "weighted_stroma": 0.0,
            "weighted_penetration": 0.0,
            "weighted_resistance_penalty": cfg.resistance_penalty_weight,
            "weighted_invasion_penalty": cfg.invasion_penalty_weight,
            "config": cfg.to_dict(),
        }

    tumor_kill_score = _clamp01(
        1.0 - (_safe_div(float(metrics.live_tumor_cells), float(max(1, cfg.initial_tumor_cells))))
    )
    stroma_score = _clamp01(1.0 - metrics.stroma_barrier_score)
    penetration_score = _clamp01(_safe_div(metrics.drug_penetration, max(cfg.max_expected_drug, 1e-12)))
    resistance_penalty = _clamp01(1.0 - metrics.mean_tumor_drug_sensitivity)
    invasion_penalty = _clamp01(_safe_div(metrics.tumor_extent, max(cfg.domain_size, 1e-12)))

    weighted_tumor_kill = cfg.tumor_kill_weight * tumor_kill_score
    weighted_stroma = cfg.stroma_disruption_weight * stroma_score
    weighted_penetration = cfg.drug_penetration_weight * penetration_score
    weighted_resistance_penalty = cfg.resistance_penalty_weight * resistance_penalty
    weighted_invasion_penalty = cfg.invasion_penalty_weight * invasion_penalty

    raw_fitness = (
        weighted_tumor_kill
        + weighted_stroma
        + weighted_penetration
        + weighted_resistance_penalty
        + weighted_invasion_penalty
    )
    fitness = _clamp01(raw_fitness)

    return {
        "valid": True,
        "fitness": fitness,
        "raw_fitness": raw_fitness,
        "tumor_kill_score": tumor_kill_score,
        "stroma_score": stroma_score,
        "penetration_score": penetration_score,
        "resistance_penalty": resistance_penalty,
        "invasion_penalty": invasion_penalty,
        "weighted_tumor_kill": weighted_tumor_kill,
        "weighted_stroma": weighted_stroma,
        "weighted_penetration": weighted_penetration,
        "weighted_resistance_penalty": weighted_resistance_penalty,
        "weighted_invasion_penalty": weighted_invasion_penalty,
        "config": cfg.to_dict(),
    }


def compute_fitness(metrics: SimulationMetrics) -> float:
    """
    Scalar Cycle 1 fitness objective (higher is better).

    Formula:
        0.35*tumor_kill + 0.25*stroma_disruption + 0.20*penetration
      - 0.10*resistance_penalty - 0.10*invasion_penalty

    Final value is clamped to [0, 1] for stable EA selection pressure.
    """

    detailed = compute_fitness_detailed(metrics)
    return float(_clamp01(detailed.get("fitness", 0.0)))


def _safe_div(num: float, den: float) -> float:
    if den == 0.0:
        return 0.0
    return num / den


def _clamp01(value: float) -> float:
    if not math.isfinite(value):
        return 0.0
    if value < 0.0:
        return 0.0
    if value > 1.0:
        return 1.0
    return float(value)
