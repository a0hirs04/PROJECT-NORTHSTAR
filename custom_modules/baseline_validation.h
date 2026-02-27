#ifndef BASELINE_VALIDATION_H
#define BASELINE_VALIDATION_H

// ============================================================================
// baseline_validation.h — Stroma World / PROJECT-NORTHSTAR
//
// §8.1 Baseline Behavioral Validation Suite
//       "Tumor-Intrinsic Baseline Traits v1.0"
//
// PURPOSE:
//   Verifies that the BooleanNetwork gene-network logic faithfully implements
//   all eight tumor-intrinsic baseline traits before any EA evolution begins.
//
//   Each trait maps directly to one or more BooleanNetwork behaviors that
//   can be checked analytically (no full PhysiCell simulation required).
//   Tests operate on BooleanNetwork instances with default ThresholdConfig
//   values — no XML parsing, no Cell/Phenotype objects needed.
//
// CALL SITE:
//   // In setup_tissue() immediately after cell types are registered:
//   if (!run_baseline_validation()) {
//       std::cerr << "[FATAL] Baseline traits failed — do not proceed.\n";
//       exit(EXIT_FAILURE);
//   }
//
// TRAITS TESTED:
//   1. Constitutive Paracrine Secretion (TGFB1, SHH driven by KRAS=GOF)
//   2. TGF-β Insensitivity (SMAD4=LOF disconnects arrest arm; invasion arm intact)
//   3. Rapid Proliferation / Broken Checkpoints (KRAS+/MYC+, braking=NONE)
//   4. Apoptotic Resistance at Baseline (BCL-XL+/TP53-; apoptosis << proliferation)
//   5. Inducible EMT (hypoxic/TGF-β periphery mesenchymal; normoxic core epithelial)
//   6. Hypoxia-Responsive Phenotype Switching (HIF1A activates; secretion amplified)
//   7. Drug-Inducible Efflux — not constitutive (NRF2→ABCB1 only when drug present)
//   8. ECM Compaction — mechanical not signaling (solid_stress formula verified)
// ============================================================================

/// Run all 8 baseline behavioral tests.
///
/// Writes a PASS/FAIL summary table to std::cerr.
/// @returns  true if all tests pass; false if any test fails.
bool run_baseline_validation();

#endif // BASELINE_VALIDATION_H
