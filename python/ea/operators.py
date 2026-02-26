"""
Genetic operators for the Stroma World evolutionary algorithm.

Implements the crossover, mutation, and selection operators used by
StromaWorldEA.  All operators work on the variable-length individual
representation: a list of intervention dicts.

Operator specifications (from implementation.pdf, Prompt 4.4):

MUTATION (mutate_individual):
  For each intervention in the individual, with probability mutation_rate:
    - 40%: mutate strength by ±0.1 (Gaussian, clamped to [0.1, 1.0])
    - 20%: change target gene to a random (different) druggable gene
    - 10%: change effect type (INHIBIT ↔ ACTIVATE)
    - 15%: add a new random intervention (if below max_interventions)
    - 15%: remove a random intervention (if above 1)

CROSSOVER (crossover_individuals):
  Uniform crossover at the intervention level:
    - Pool all interventions from both parents.
    - Each child randomly selects from the pool.
    - Deduplicate: no two interventions on the same gene in one child.

SELECTION (tournament_select):
  Tournament selection — draw ``tournament_size`` individuals, return the
  best fitness.
"""

from __future__ import annotations

import copy
import random
from typing import Dict, List, Optional, Tuple

from .population import (
    ALLOWED_EFFECTS,
    DEFAULT_DRUGGABLE_GENES,
    random_intervention,
    validate_individual,
)


# ---------------------------------------------------------------------------
# Mutation
# ---------------------------------------------------------------------------

def mutate_individual(
    individual: List[Dict],
    mutation_rate: float = 0.2,
    max_interventions: int = 5,
    druggable_genes: Optional[List[str]] = None,
) -> List[Dict]:
    """Apply mutation operators to an individual in-place (returns a copy).

    With probability ``mutation_rate``, one of five sub-operators is applied
    to each intervention.  The individual is returned as a new list so the
    original is not modified (DEAP expects a fresh object).

    Args:
        individual: List of intervention dicts.
        mutation_rate: Per-intervention probability of mutation.
        max_interventions: Upper bound on individual length.
        druggable_genes: Available gene targets.

    Returns:
        Mutated copy of the individual.
    """
    genes = druggable_genes if druggable_genes is not None else DEFAULT_DRUGGABLE_GENES
    ind = copy.deepcopy(individual)

    for i in range(len(ind)):
        if random.random() >= mutation_rate:
            continue

        # Choose sub-operator proportionally.
        r = random.random()

        if r < 0.40:
            # Strength perturbation (Gaussian ±0.1, clamped to [0.1, 1.0]).
            delta = random.gauss(0.0, 0.1)
            new_strength = float(ind[i]["strength"]) + delta
            ind[i]["strength"] = round(max(0.1, min(1.0, new_strength)), 4)

        elif r < 0.60:  # 20% bucket
            # Change target gene to a different druggable gene.
            current_genes = {iv["gene"] for iv in ind}
            candidates = [g for g in genes if g != ind[i]["gene"] and g not in current_genes]
            if candidates:
                ind[i]["gene"] = random.choice(candidates)

        elif r < 0.70:  # 10% bucket
            # Toggle effect type.
            ind[i]["effect"] = "ACTIVATE" if ind[i]["effect"] == "INHIBIT" else "INHIBIT"

        elif r < 0.85:  # 15% bucket — add a new intervention
            if len(ind) < max_interventions:
                existing_genes = {iv["gene"] for iv in ind}
                candidates = [g for g in genes if g not in existing_genes]
                if candidates:
                    new_gene = random.choice(candidates)
                    new_effect = "INHIBIT" if random.random() < 0.8 else "ACTIVATE"
                    new_strength = round(random.uniform(0.1, 1.0), 4)
                    ind.append({"gene": new_gene, "effect": new_effect, "strength": new_strength})

        else:  # 15% bucket — remove a random intervention
            if len(ind) > 1:
                del ind[random.randrange(len(ind))]
                # i may now be out of bounds; loop will terminate naturally
                break  # restart loop safely

    return ind


# ---------------------------------------------------------------------------
# Crossover
# ---------------------------------------------------------------------------

def crossover_individuals(
    ind1: List[Dict],
    ind2: List[Dict],
    max_interventions: int = 5,
) -> Tuple[List[Dict], List[Dict]]:
    """Uniform crossover at the intervention level.

    Algorithm:
    1. Pool all interventions from both parents (deduplicate by gene name,
       keeping the last occurrence to avoid bias).
    2. Each child randomly samples from the combined pool.
    3. Ensure no two interventions in one child target the same gene.
    4. Each child has at most ``max_interventions`` interventions.

    Args:
        ind1: First parent individual.
        ind2: Second parent individual.
        max_interventions: Maximum interventions per child.

    Returns:
        Tuple of two child individuals (new lists, parents unchanged).
    """
    # Build a combined pool keyed by gene (last-writer-wins for duplicates).
    pool: Dict[str, Dict] = {}
    for iv in ind1:
        pool[iv["gene"]] = copy.deepcopy(iv)
    for iv in ind2:
        pool[iv["gene"]] = copy.deepcopy(iv)

    pool_list = list(pool.values())
    random.shuffle(pool_list)

    def _sample_child() -> List[Dict]:
        seen: set = set()
        child: List[Dict] = []
        for iv in pool_list:
            if iv["gene"] not in seen and len(child) < max_interventions:
                child.append(copy.deepcopy(iv))
                seen.add(iv["gene"])
        # Ensure at least 1 intervention.
        if not child and pool_list:
            child = [copy.deepcopy(random.choice(pool_list))]
        return child

    child1 = _sample_child()
    # Re-shuffle so child2 gets a different ordering.
    random.shuffle(pool_list)
    child2 = _sample_child()

    return child1, child2


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------

def tournament_select(
    population: List,
    fitnesses: List[float],
    tournament_size: int = 3,
) -> object:
    """Select one individual via tournament selection.

    Randomly samples ``tournament_size`` indices, returns the individual
    with the highest fitness among them.

    Args:
        population: List of individuals.
        fitnesses: Corresponding fitness values (same length as population).
        tournament_size: Number of contestants per tournament.

    Returns:
        The winning individual (not a copy — caller should deepcopy if needed).
    """
    contestants = random.sample(range(len(population)), min(tournament_size, len(population)))
    winner_idx = max(contestants, key=lambda i: fitnesses[i])
    return population[winner_idx]


def elitist_survivors(
    population: List,
    fitnesses: List[float],
    elite_fraction: float = 0.05,
) -> Tuple[List, List[float]]:
    """Extract the elite subset of the population for carry-forward.

    Args:
        population: Current population.
        fitnesses: Corresponding fitness values.
        elite_fraction: Fraction of population to keep (e.g. 0.05 = top 5%).

    Returns:
        Tuple (elite_individuals, elite_fitnesses) sorted by fitness descending.
    """
    n_elite = max(1, int(len(population) * elite_fraction))
    ranked = sorted(zip(population, fitnesses), key=lambda x: x[1], reverse=True)
    elites = ranked[:n_elite]
    elite_inds = [copy.deepcopy(e[0]) for e in elites]
    elite_fits = [e[1] for e in elites]
    return elite_inds, elite_fits
