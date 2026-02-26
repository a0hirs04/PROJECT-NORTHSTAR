from __future__ import annotations

from python.ea.fitness import compute_fitness, validate_metrics


def test_compute_fitness_known_good_metrics_high(metrics_factory):
    metrics = metrics_factory(
        live_tumor_cells=5,
        stroma_barrier_score=0.1,
        drug_penetration=0.9,
        mean_tumor_drug_sensitivity=0.9,
        tumor_extent=200.0,
        mean_ecm_density=0.3,
        max_ecm_density=0.5,
    )
    score = compute_fitness(metrics)
    assert 0.6 <= score <= 1.0


def test_compute_fitness_known_bad_metrics_low(metrics_factory):
    metrics = metrics_factory(
        live_tumor_cells=50,
        stroma_barrier_score=0.95,
        drug_penetration=0.05,
        mean_tumor_drug_sensitivity=0.1,
        tumor_extent=1800.0,
        mean_ecm_density=0.95,
        max_ecm_density=1.0,
    )
    score = compute_fitness(metrics)
    assert 0.0 <= score <= 0.15


def test_resistance_penalty_lowers_fitness(metrics_factory):
    sensitive = metrics_factory(mean_tumor_drug_sensitivity=0.9)
    resistant = metrics_factory(mean_tumor_drug_sensitivity=0.1)

    sensitive_score = compute_fitness(sensitive)
    resistant_score = compute_fitness(resistant)
    assert resistant_score < sensitive_score


def test_invasion_penalty_lowers_fitness(metrics_factory):
    low_invasion = metrics_factory(tumor_extent=100.0)
    high_invasion = metrics_factory(tumor_extent=1900.0)

    low_score = compute_fitness(low_invasion)
    high_score = compute_fitness(high_invasion)
    assert high_score < low_score


def test_boundary_conditions_zero_cells_and_invalid_ecm(metrics_factory):
    zero_cells = metrics_factory(
        total_tumor_cells=0,
        live_tumor_cells=0,
        total_stromal_cells=0,
        activated_cafs=0,
        mesenchymal_tumor_cells=0,
        mean_ecm_density=0.0,
        max_ecm_density=0.0,
        stroma_barrier_score=0.0,
        drug_penetration=0.0,
        mean_tumor_drug_sensitivity=0.5,
        hypoxic_fraction=0.0,
        tumor_extent=0.0,
    )
    assert validate_metrics(zero_cells)
    assert 0.0 <= compute_fitness(zero_cells) <= 1.0

    invalid_ecm = metrics_factory(mean_ecm_density=1.2, max_ecm_density=1.2)
    assert not validate_metrics(invalid_ecm)
    assert compute_fitness(invalid_ecm) == 0.0
