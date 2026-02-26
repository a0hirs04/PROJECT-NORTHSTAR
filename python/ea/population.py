"""
Population management utilities for the Stroma World evolutionary algorithm.

This module provides helpers for creating and validating individuals in the
EA population. Each individual is a list of intervention dicts:

    [
        {"gene": "EGFR",   "effect": "INHIBIT",  "strength": 0.73},
        {"gene": "HAS2",   "effect": "INHIBIT",  "strength": 0.45},
        ...
    ]

The list is variable-length (1 to max_interventions). Each intervention
targets one druggable gene with a specified effect type and continuous
strength in [0.1, 1.0].
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional

# Druggable gene pool — matches EAConfig.druggable_genes default.
DEFAULT_DRUGGABLE_GENES: List[str] = [
    "EGFR",
    "BCL_XL",
    "MMP2",
    "TGFB1",
    "SHH",
    "GLI1",
    "HAS2",
    "ABCB1",
]

# Effect types available to the EA.  INHIBIT is weighted more heavily (80%)
# because the majority of clinically-relevant interventions are inhibitory.
ALLOWED_EFFECTS: tuple = ("INHIBIT", "ACTIVATE")
INHIBIT_WEIGHT: float = 0.8  # probability of choosing INHIBIT over ACTIVATE


def random_intervention(
    druggable_genes: Optional[List[str]] = None,
    inhibit_weight: float = INHIBIT_WEIGHT,
) -> Dict:
    """Return one randomly generated intervention dict.

    Args:
        druggable_genes: Pool of gene names the EA can target.
            Defaults to DEFAULT_DRUGGABLE_GENES.
        inhibit_weight: Probability of choosing INHIBIT effect (vs ACTIVATE).

    Returns:
        dict with keys ``gene``, ``effect``, ``strength``.
    """
    genes = druggable_genes if druggable_genes is not None else DEFAULT_DRUGGABLE_GENES
    gene = random.choice(genes)
    effect = "INHIBIT" if random.random() < inhibit_weight else "ACTIVATE"
    strength = round(random.uniform(0.1, 1.0), 4)
    return {"gene": gene, "effect": effect, "strength": strength}


def random_individual(
    max_interventions: int = 5,
    min_interventions: int = 1,
    druggable_genes: Optional[List[str]] = None,
    inhibit_weight: float = INHIBIT_WEIGHT,
) -> List[Dict]:
    """Create one random individual (list of intervention dicts).

    The number of interventions is uniformly drawn from
    [min_interventions, max_interventions]. Gene targets are deduplicated
    so no gene appears more than once.

    Args:
        max_interventions: Maximum interventions per individual.
        min_interventions: Minimum interventions per individual.
        druggable_genes: Pool of gene names available for targeting.
        inhibit_weight: Probability of choosing INHIBIT effect.

    Returns:
        List of intervention dicts (length in [min_interventions, max_interventions]).
    """
    genes = druggable_genes if druggable_genes is not None else DEFAULT_DRUGGABLE_GENES
    n = random.randint(min_interventions, min(max_interventions, len(genes)))

    # Sample distinct gene targets, then assign effect + strength.
    selected_genes = random.sample(genes, n)
    individual = []
    for gene in selected_genes:
        effect = "INHIBIT" if random.random() < inhibit_weight else "ACTIVATE"
        strength = round(random.uniform(0.1, 1.0), 4)
        individual.append({"gene": gene, "effect": effect, "strength": strength})
    return individual


def initialize_population(
    population_size: int,
    max_interventions: int = 5,
    min_interventions: int = 1,
    druggable_genes: Optional[List[str]] = None,
    inhibit_weight: float = INHIBIT_WEIGHT,
) -> List[List[Dict]]:
    """Create a full initial population of random individuals.

    Args:
        population_size: Number of individuals to generate.
        max_interventions: Maximum interventions per individual.
        min_interventions: Minimum interventions per individual.
        druggable_genes: Gene pool for targeting.
        inhibit_weight: Probability of INHIBIT over ACTIVATE.

    Returns:
        List of ``population_size`` random individuals.
    """
    return [
        random_individual(
            max_interventions=max_interventions,
            min_interventions=min_interventions,
            druggable_genes=druggable_genes,
            inhibit_weight=inhibit_weight,
        )
        for _ in range(population_size)
    ]


def validate_individual(
    individual: List[Dict],
    druggable_genes: Optional[List[str]] = None,
    max_interventions: int = 5,
) -> bool:
    """Check that an individual is structurally valid.

    Validation rules:
    - Must be a non-empty list with length <= max_interventions.
    - Each entry must have ``gene``, ``effect``, and ``strength`` keys.
    - ``gene`` must be in the druggable_genes pool.
    - ``effect`` must be one of ALLOWED_EFFECTS.
    - ``strength`` must be in [0.0, 1.0].
    - No two entries may target the same gene (deduplication constraint).

    Args:
        individual: The individual to validate.
        druggable_genes: Allowed gene targets.
        max_interventions: Maximum allowed interventions.

    Returns:
        True if valid, False otherwise.
    """
    genes = set(druggable_genes) if druggable_genes is not None else set(DEFAULT_DRUGGABLE_GENES)
    if not individual or len(individual) > max_interventions:
        return False
    seen_genes: set = set()
    for entry in individual:
        if not isinstance(entry, dict):
            return False
        if not {"gene", "effect", "strength"}.issubset(entry):
            return False
        if entry["gene"] not in genes:
            return False
        if entry["effect"] not in ALLOWED_EFFECTS:
            return False
        s = entry["strength"]
        if not (isinstance(s, (int, float)) and 0.0 <= float(s) <= 1.0):
            return False
        if entry["gene"] in seen_genes:
            return False  # duplicate gene target
        seen_genes.add(entry["gene"])
    return True


def individual_to_json_payload(individual: List[Dict]) -> dict:
    """Convert an individual to the intervention JSON format consumed by PhysiCell.

    Returns:
        dict matching the schema expected by BooleanNetwork::load_from_json(),
        e.g.::

            {
                "interventions": [
                    {"gene": "EGFR", "effect": "INHIBIT", "strength": 0.8, "name": "EA_EGFR_0"},
                    ...
                ]
            }
    """
    entries = []
    for i, iv in enumerate(individual):
        entries.append(
            {
                "gene": iv["gene"],
                "effect": iv["effect"],
                "strength": float(iv["strength"]),
                "name": f"EA_{iv['gene']}_{i}",
            }
        )
    return {"interventions": entries}
