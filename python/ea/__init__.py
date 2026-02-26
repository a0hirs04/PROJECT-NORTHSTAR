"""
Evolutionary Algorithm package for Stroma World.

Provides the DEAP-based EA that evolves intervention strategies
(drug combinations and strengths) to optimize the fitness function
against the PhysiCell PDAC stroma simulation.

Modules:
    evolutionary_algorithm — Main EA loop (StromaWorldEA, EAConfig, EAResult)
    fitness                — Fitness function (compute_fitness, FitnessConfig)
    population             — Individual representation and initialization helpers
    operators              — Crossover, mutation, and selection operators
"""
