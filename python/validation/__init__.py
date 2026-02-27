"""
Biological validation package for Stroma World.

Runs a dependency-ordered set of known-outcome simulations to verify the
PhysiCell model reproduces PDAC "Reality Anchors". Serves as a regression
test suite for emergent biology wiring.

Modules:
    validate_biology — Ten anchor scenarios with dependency gating:
                       1) self-assembly
                       2) penetration vs maturity
                       3) SHH paradox
                       4) SMAD4 asymmetry
                       5) central hypoxia
                       6) peripheral EMT
                       7) ECM-degradation limits
                       8) two-component barrier proxy
                       9) density prognostic behavior
                       10) spatial sanctuary integration
                       Includes a 34-test directional pass/fail matrix and
                       replicate-median decision protocol.
"""
