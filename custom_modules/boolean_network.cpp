// ============================================================================
// boolean_network.cpp — Stroma World / PROJECT-NORTHSTAR
//
// Implementation of the continuous Boolean gene network for PDAC tumor-stroma
// simulation. See boolean_network.h for design overview.
// ============================================================================

#include "boolean_network.h"
#include "../modules/PhysiCell_settings.h"   // PhysiCell::parameters
#include <algorithm>       // std::fill
#include <cstring>         // std::memset
#include <fstream>         // std::ifstream
#include <stdexcept>       // std::runtime_error
#include <unordered_map>   // gene name → index lookup table
#include <iostream>        // std::cerr (intervention logging)

// nlohmann/json — bundled single-header in BioFVM/
#include "../BioFVM/json.hpp"
using json = nlohmann::json;

// ============================================================================
// ThresholdConfig::load_from_xml
//
// Reads authoritative threshold values from PhysiCell's parsed XML parameters.
// Call once at simulation setup AFTER PhysiCell has loaded the XML (i.e.,
// after setup_microenvironment() or at the top of setup_tissue()).
//
// Pattern:
//   ThresholdConfig cfg = ThresholdConfig::load_from_xml();
//   // pass cfg to BooleanNetwork::update() for every cell every timestep
//
// If a parameter is absent from the XML, the default value (matching the
// compile-time constant in gene_definitions.h) is silently retained.
// ============================================================================

ThresholdConfig ThresholdConfig::load_from_xml()
{
    ThresholdConfig cfg;  // starts with compile-time fallback defaults

    // Helper: read a double parameter by name, warn if absent, keep default.
    auto read_double = [&](const char* name, double& target)
    {
        try
        {
            target = PhysiCell::parameters.doubles(name);
        }
        catch (...)
        {
            std::cerr << "[ThresholdConfig] WARNING: parameter '" << name
                      << "' not found in XML; using default = " << target << "\n";
        }
    };

    read_double("oxygen_hypoxia_threshold",   cfg.hypoxia_hif1a_threshold);
    read_double("oxygen_necrosis_threshold",  cfg.hypoxia_necrosis_threshold);
    read_double("tgfb_activation_threshold",  cfg.tgfb_activation_threshold);
    read_double("shh_activation_threshold",   cfg.shh_activation_threshold);

    // ecm_high_threshold, contact_inhibition_threshold, emt_zeb1_threshold,
    // stress_nrf2_threshold have no XML user_parameters yet — keep defaults.

    std::cerr << "[ThresholdConfig] Loaded from XML:"
              << "  HIF1A O2=" << cfg.hypoxia_hif1a_threshold
              << " mmHg  Necrosis O2=" << cfg.hypoxia_necrosis_threshold
              << " mmHg  TGFB=" << cfg.tgfb_activation_threshold
              << "  SHH=" << cfg.shh_activation_threshold << "\n";

    return cfg;
}

// ============================================================================
// TAU ASSIGNMENT RATIONALE
//
// Time constants reflect biological turnover rates of each gene product:
//
//   6.0  min — signaling molecules: rapid post-translational modification,
//              protein phosphorylation/dephosphorylation (e.g., HIF1A protein
//              stabilization under hypoxia occurs within minutes).
//
//  60.0  min — transcription factors: nuclear import, transcription, mRNA
//              half-lives (~30-90 min for most TF mRNAs). Covers KRAS effectors
//              (MYC, ZEB1), cell-death regulators (BCL_XL, BAX), and pathway TFs.
//
// 360.0  min — structural / secreted proteins: synthesis, folding, secretion,
//              and ECM incorporation are slow processes. Collagen fibers have
//              half-lives of hours to days in vivo; hyaluronic acid ~2-5 hours.
// ============================================================================

static constexpr double TAU_SIGNALING  =   6.0;  // rapid signaling proteins
static constexpr double TAU_TF         =  60.0;  // transcription factors / regulators
static constexpr double TAU_STRUCTURAL = 360.0;  // secreted / structural proteins

// ============================================================================
// GeneParams
// ============================================================================

GeneParams::GeneParams()
{
    std::fill(state_offset,   state_offset   + GENE_COUNT, 0.0);
    std::fill(tau_scale,      tau_scale      + GENE_COUNT, 1.0);
    std::fill(drug_inhibition, drug_inhibition + GENE_COUNT, 0.0);
}

// ============================================================================
// BooleanNetwork constructor
// ============================================================================

BooleanNetwork::BooleanNetwork()
    : cell_type_(CellType::TUMOR), params_()
{
    std::fill(gene_states, gene_states + GENE_COUNT, 0.0);
    std::fill(is_mutant,   is_mutant   + GENE_COUNT, false);
    assign_default_tau();
}

// ============================================================================
// assign_default_tau
//
// Classifies each gene into a biological turnover category and assigns the
// corresponding time constant. Applied at construction; scaled by params_.tau_scale
// during update().
// ============================================================================

void BooleanNetwork::assign_default_tau()
{
    // Signaling proteins — rapid kinetics
    tau[HIF1A]  = TAU_SIGNALING;   // HIF1A: O2-regulated prolyl hydroxylation, ~minutes
    tau[NRF2]   = TAU_SIGNALING;   // NRF2:  Keap1-mediated ubiquitination, rapid

    // Transcription factors and their direct targets
    tau[KRAS]   = TAU_TF;          // KRAS:  GTPase cycling (locked mutant, tau irrelevant)
    tau[MYC]    = TAU_TF;          // MYC:   short-lived TF (~30 min mRNA half-life)
    tau[EGFR]   = TAU_TF;          // EGFR:  RTK activation/internalization
    tau[TP53]   = TAU_TF;          // TP53:  MDM2-regulated stability (locked mutant)
    tau[BCL_XL] = TAU_TF;          // BCL-XL: anti-apoptotic; KRAS-regulated
    tau[BAX]    = TAU_TF;          // BAX:   pro-apoptotic; TP53-regulated
    tau[CDKN2A] = TAU_TF;          // CDKN2A: epigenetically silenced (locked mutant)
    tau[SMAD4]  = TAU_TF;          // SMAD4: nuclear TGFb effector (locked mutant)
    tau[RB1]    = TAU_TF;          // RB1:   phosphorylation state changes over hours
    tau[ZEB1]   = TAU_TF;          // ZEB1:  EMT TF, hours to switch state
    tau[CDH1]   = TAU_TF;          // CDH1:  E-cadherin mRNA turnover ~hours
    tau[ACTA2]  = TAU_TF;          // ACTA2: alpha-SMA induction takes hours in CAFs
    tau[TGFB1]  = TAU_TF;          // TGFB1: mRNA expression proxy
    tau[SHH]    = TAU_TF;          // SHH:   ligand expression regulated transcriptionally
    tau[GLI1]   = TAU_TF;          // GLI1:  Hedgehog target TF, nuclear translocation
    tau[ABCB1]  = TAU_TF;          // ABCB1: drug-induced upregulation, hours
    tau[MMP2]   = TAU_TF;          // MMP2:  secreted protease, moderate turnover

    // Structural / secreted ECM proteins — slow turnover
    tau[HAS2]   = TAU_STRUCTURAL;  // HAS2:  HA synthase; HA t1/2 ~2-5 hours
    tau[COL1A1] = TAU_STRUCTURAL;  // COL1A1: collagen fiber assembly, hours to days
}

// ============================================================================
// initialize
//
// Sets gene states from GENE_INFO defaults for the given cell_type, applies
// EA parameter offsets, locks PDAC mutant genes, and scales tau values.
// ============================================================================

void BooleanNetwork::initialize(CellType cell_type, const GeneParams& params)
{
    cell_type_ = cell_type;
    params_    = params;

    assign_default_tau();

    // ---- Set initial gene states from metadata defaults --------------------
    for (int i = 0; i < GENE_COUNT; ++i)
    {
        double base = (cell_type == CellType::TUMOR)
                      ? GENE_INFO[i].default_tumor_state
                      : GENE_INFO[i].default_stroma_state;

        // Apply EA-evolvable offset (clamped to [0,1])
        gene_states[i] = clamp(base + params.state_offset[i], 0.0, 1.0);

        // Apply tau scaling from params (EA can tune dynamics)
        tau[i] *= params.tau_scale[i];

        is_mutant[i] = false;
    }

    // ---- Lock PDAC-defining mutations (tumor cells only) -------------------
    // These represent irreversible genomic alterations that cannot be reversed
    // by drug interventions in the current model. The EA may not evolve these.
    if (cell_type == CellType::TUMOR)
    {
        // KRAS G12D/V: constitutively active GTPase — always ON
        gene_states[KRAS]   = 1.0;
        is_mutant[KRAS]     = true;

        // TP53 loss-of-function mutation — always OFF (apoptosis gatekeeper lost)
        gene_states[TP53]   = 0.0;
        is_mutant[TP53]     = true;

        // CDKN2A homozygous deletion — always OFF (G1 brake lost)
        gene_states[CDKN2A] = 0.0;
        is_mutant[CDKN2A]   = true;

        // SMAD4 loss — always OFF (TGFb anti-proliferative arm lost)
        gene_states[SMAD4]  = 0.0;
        is_mutant[SMAD4]    = true;
    }
}

// ============================================================================
// sync_from_cell
//
// Reads current gene states from pCell->custom_data[] into the local
// gene_states[] array. Called at the start of each phenotype timestep so
// that the network reflects any external modifications made by PhysiCell
// (e.g., internalized drug affecting ABCB1, contact-inhibition signals).
// ============================================================================

void BooleanNetwork::sync_from_cell(PhysiCell::Cell* pCell)
{
    for (int i = 0; i < GENE_COUNT; ++i)
    {
        gene_states[i] = pCell->custom_data[i];
    }
}

// ============================================================================
// sync_to_cell
//
// Writes the updated gene_states[] back into pCell->custom_data[] so that
// PhysiCell's phenotype mapping functions (tumor_cell.cpp, stromal_cell.cpp)
// can read the current network state when computing cycle rates, apoptosis
// rates, secretion rates, and motility parameters.
// ============================================================================

void BooleanNetwork::sync_to_cell(PhysiCell::Cell* pCell) const
{
    for (int i = 0; i < GENE_COUNT; ++i)
    {
        pCell->custom_data[i] = gene_states[i];
    }
}

// ============================================================================
// update
//
// Main network integration step. Called every dt_phenotype (6 min default).
//
// Algorithm:
//   1. Compute target[] from Boolean rules (compute_*_targets)
//   2. Apply drug inhibition: target[i] *= (1 - drug_inhibition[i])
//   3. Apply first-order update: gene[i] += dt * (target[i] - gene[i]) / tau[i]
//   4. Clamp all states to [0, 1]
//   5. Re-lock mutant genes to their fixed values
// ============================================================================

void BooleanNetwork::update(double dt,
                            double oxygen,
                            double tgfb_local,
                            double shh_local,
                            double drug_local,
                            const ThresholdConfig& cfg)
{
    double target[GENE_COUNT];
    std::fill(target, target + GENE_COUNT, 0.0);

    // ---- Step 1: Compute Boolean rule targets ------------------------------
    if (cell_type_ == CellType::TUMOR)
        compute_tumor_targets(target, oxygen, tgfb_local, shh_local, drug_local, cfg);
    else
        compute_stroma_targets(target, oxygen, tgfb_local, shh_local, drug_local, cfg);

    // ---- Step 2: Apply drug inhibition (competitive / graded) -------------
    // drug_inhibition[i] in [0, 1]: fraction by which the target is reduced.
    // 0.0 = no effect; 1.0 = complete blockade (drives target toward 0).
    // This represents the EA's ability to intervene on specific gene nodes.
    for (int i = 0; i < GENE_COUNT; ++i)
    {
        target[i] = clamp(target[i] * (1.0 - params_.drug_inhibition[i]), 0.0, 1.0);
    }

    // ---- Step 3 & 4: First-order update + clamp ----------------------------
    for (int i = 0; i < GENE_COUNT; ++i)
    {
        if (is_mutant[i]) continue;  // mutant genes are permanently locked

        const double effective_tau = tau[i] > 0.0 ? tau[i] : 1.0;  // guard /0
        gene_states[i] += dt * (target[i] - gene_states[i]) / effective_tau;
        gene_states[i]  = clamp(gene_states[i], 0.0, 1.0);
    }

    // ---- Step 5: Re-lock mutant genes (safety: prevents drift from sync) ---
    if (cell_type_ == CellType::TUMOR)
    {
        gene_states[KRAS]   = 1.0;
        gene_states[TP53]   = 0.0;
        gene_states[CDKN2A] = 0.0;
        gene_states[SMAD4]  = 0.0;
    }
}

// ============================================================================
// compute_tumor_targets
//
// Evaluates the Boolean-logic rules for tumor cells. Each rule computes a
// target value in [0, 1] representing the gene's "desired" steady state given
// current inputs. The actual state relaxes toward this target with time constant
// tau (see update()).
//
// Shorthand aliases (const references) are used for readability.
// ============================================================================

void BooleanNetwork::compute_tumor_targets(double* target,
                                           double oxygen,
                                           double tgfb_local,
                                           double shh_local,
                                           double drug_local,
                                           const ThresholdConfig& cfg) const
{
    // Shorthand aliases for current gene states — read-only
    const double KRAS_   = gene_states[KRAS];
    const double EGFR_   = gene_states[EGFR];
    const double TP53_   = gene_states[TP53];
    const double NRF2_   = gene_states[NRF2];
    const double HIF1A_  = gene_states[HIF1A];
    const double ZEB1_   = gene_states[ZEB1];
    const double CDH1_   = gene_states[CDH1];
    const double CDKN2A_ = gene_states[CDKN2A];
    const double SMAD4_  = gene_states[SMAD4];

    // =========================================================================
    // GROWTH AXIS
    // =========================================================================

    // MYC — driven by KRAS, modulated by EGFR availability.
    //
    // Biological basis:
    //   KRAS activates ERK and AKT, both of which stabilize MYC protein and
    //   promote MYC transcription. EGFR reinforces RAS/MAPK signaling; when
    //   EGFR is pharmacologically blocked (erlotinib), MYC is partially
    //   suppressed even in the presence of mutant KRAS (~50% reduction).
    //   The factor (1 - 0.5*(1-EGFR)) = 0.5 + 0.5*EGFR: ranges from 0.5
    //   (EGFR fully blocked) to 1.0 (EGFR fully active).
    target[MYC] = clamp(KRAS_ * (0.5 + 0.5 * EGFR_), 0.0, 1.0);

    // =========================================================================
    // DEATH AXIS
    // =========================================================================

    // BCL_XL — anti-apoptotic; basally expressed, boosted by survival signals.
    //
    // Biological basis:
    //   BCL-XL is transcriptionally activated by NF-kB (downstream of KRAS)
    //   and by NRF2 (oxidative stress survival). The 0.3 basal term reflects
    //   constitutive expression in many cell types. KRAS contributes 0.5
    //   via PI3K/AKT/NF-kB. NRF2 contributes 0.2 (stress-induced survival).
    //   Sum can reach 1.0 in worst-case (KRAS ON + NRF2 ON).
    target[BCL_XL] = clamp(0.3 + 0.5 * KRAS_ + 0.2 * NRF2_, 0.0, 1.0);

    // BAX — pro-apoptotic; TP53-dependent, basally present at low levels.
    //
    // Biological basis:
    //   BAX is a direct transcriptional target of p53. In PDAC, TP53 is lost,
    //   so BAX is only basally expressed (0.3 floor from p53-independent
    //   pathways). When TP53 is restored (therapeutic scenario), BAX rises
    //   strongly, tipping the BCL_XL/BAX balance toward apoptosis.
    target[BAX] = clamp(0.3 + 0.5 * TP53_, 0.0, 1.0);

    // =========================================================================
    // BRAKING AXIS
    // =========================================================================

    // RB1 — functional G1/S checkpoint; depends on upstream brakes.
    //
    // Biological basis:
    //   pRB (RB1) is maintained in its active (hypophosphorylated) form by
    //   CDK4/6 inhibition via p16/CDKN2A. When CDKN2A is deleted, CDK4/6
    //   constitutively phosphorylates and inactivates RB1 (even if RB1 gene
    //   is intact). SMAD4 also promotes RB1-mediated quiescence via
    //   p15/CDKN2B induction. Both lost in PDAC → RB1 functionally ~0.
    //   Weights: CDKN2A contributes 70% of RB1 functional support; SMAD4 30%.
    target[RB1] = clamp(CDKN2A_ * 0.7 + SMAD4_ * 0.3, 0.0, 1.0);

    // =========================================================================
    // INVASION / EMT AXIS
    // =========================================================================

    // ZEB1 — EMT master regulator; bistable switch with CDH1.
    //
    // Biological basis:
    //   ZEB1 is activated by TGF-beta (via non-canonical SMAD-independent
    //   pathways, since SMAD4 is lost in PDAC), by hypoxia (HIF1A binds ZEB1
    //   promoter), and by KRAS (via ERK). Crucially, CDH1/E-cadherin suppresses
    //   ZEB1 through miR-200 family induction (CDH1 → miR-200 → represses ZEB1).
    //   CDH1_feedback is the current CDH1 value used as a repressive input.
    //   The ZEB1/CDH1 double-negative feedback creates a bistable switch:
    //   epithelial (ZEB1↓, CDH1↑) or mesenchymal (ZEB1↑, CDH1↓) states.
    {
        const double CDH1_feedback = CDH1_;
        target[ZEB1] = clamp(0.3 * tgfb_local
                           + 0.3 * HIF1A_
                           + 0.2 * KRAS_
                           - 0.3 * CDH1_feedback,
                             0.0, 1.0);
    }

    // CDH1 — E-cadherin; epithelial identity marker; ZEB1-repressed.
    //
    // Biological basis:
    //   ZEB1 directly represses CDH1 transcription. The mutual repression
    //   (ZEB1 → CDH1 off; CDH1 → ZEB1 off) creates the EMT bistable switch.
    //   CDH1 = 1 - ZEB1 implements a continuous toggle: as ZEB1 rises,
    //   E-cadherin falls, weakening cell-cell adhesion and enabling invasion.
    target[CDH1] = clamp(1.0 - ZEB1_, 0.0, 1.0);

    // MMP2 — matrix metalloproteinase-2; ECM degradation / invasion.
    //
    // Biological basis:
    //   MT1-MMP on tumor cells activates pro-MMP2 in the pericellular space.
    //   MMP2 expression is driven by ZEB1 (EMT-associated invasion program),
    //   HIF1A (hypoxia upregulates MMP2 via AP-1 and HIF binding sites),
    //   and KRAS (ERK-dependent transcriptional activation).
    target[MMP2] = clamp(0.2 * ZEB1_
                       + 0.3 * HIF1A_
                       + 0.2 * KRAS_,
                         0.0, 1.0);

    // =========================================================================
    // SECRETION AXIS
    // =========================================================================

    // TGFB1 — tumor-secreted TGF-beta; drives CAF activation in stroma.
    //
    // Biological basis:
    //   PDAC tumor cells constitutively secrete TGF-beta (0.1 basal term)
    //   through KRAS-driven transcription (0.3 * KRAS). Hypoxia further
    //   amplifies TGF-beta secretion via HIF1A binding to the TGFB1 promoter
    //   (0.2 * HIF1A). This creates a feedforward loop: hypoxia → TGF-beta
    //   → stromal activation → desmoplasia → further hypoxia.
    target[TGFB1] = clamp(0.1 + 0.3 * KRAS_ + 0.2 * HIF1A_, 0.0, 1.0);

    // SHH — Sonic Hedgehog ligand; activates GLI1 in stromal cells.
    //
    // Biological basis:
    //   SHH expression in PDAC is directly driven by oncogenic KRAS via
    //   transcription factors downstream of RAS/MAPK. SHH secreted by tumor
    //   cells binds Patched (PTCH1) on stromal fibroblasts, de-represses
    //   Smoothened (SMO), and activates GLI1 transcription — the primary
    //   route to CAF activation and desmoplasia in PDAC.
    //   Basal 0.1 reflects leaky expression even without strong KRAS signal.
    target[SHH] = clamp(0.1 + 0.3 * KRAS_, 0.0, 1.0);

    // =========================================================================
    // STRESS AXIS
    // =========================================================================

    // HIF1A — oxygen sensor; activated under severe hypoxia.
    //
    // Biological basis:
    //   Under normoxia, HIF1A protein is hydroxylated by PHD enzymes (require
    //   O2) and degraded via VHL/26S proteasome. Below ~5 mmHg, PHD activity
    //   collapses and HIF1A accumulates rapidly (within minutes). The target
    //   scales linearly from 0 at threshold to 1 at anoxia (0 mmHg), modeling
    //   the graded PHD oxygen dependence (Km ~5 mmHg for PHD2).
    if (oxygen < cfg.hypoxia_hif1a_threshold)
        target[HIF1A] = clamp(1.0 - oxygen / cfg.hypoxia_hif1a_threshold, 0.0, 1.0);
    else
        target[HIF1A] = 0.0;

    // NRF2 — oxidative/drug stress sensor; multi-input activation.
    //
    // Biological basis:
    //   NRF2 is normally sequestered by KEAP1 for proteasomal degradation.
    //   Reactive oxygen species (ROS) — elevated in hypoxia and by drug
    //   metabolism — modify KEAP1 cysteines, releasing NRF2 to translocate
    //   to the nucleus. KRAS also directly activates NRF2 transcription via
    //   an AP-1 site. Drug exposure generates ROS and electrophiles that
    //   further stabilize NRF2, making it a key node for drug resistance.
    //   drug_local > 0 term uses a step function: any non-zero drug activates.
    target[NRF2] = clamp(0.1 * KRAS_
                       + 0.2 * HIF1A_
                       + 0.3 * (drug_local > 0.0 ? 1.0 : 0.0),
                         0.0, 1.0);

    // ABCB1 — MDR1/P-glycoprotein drug efflux pump.
    //
    // Biological basis:
    //   The MDR1 gene promoter contains NRF2-binding antioxidant response
    //   elements (AREs) and HIF1A hypoxia response elements (HREs). Both
    //   stress signals transcriptionally upregulate ABCB1. The pump actively
    //   exports cytotoxic drugs, reducing intracellular drug concentration.
    //   High ABCB1 → reduced drug_sensitivity (handled in tumor_cell.cpp).
    target[ABCB1] = clamp(0.5 * NRF2_ + 0.3 * HIF1A_, 0.0, 1.0);

    // =========================================================================
    // GENES WITH NO TUMOR RULE (states persist from previous step)
    // =========================================================================
    // ACTA2, GLI1, HAS2, COL1A1 — stroma-specific; not updated in tumor cells.
    // Their target defaults to 0.0 (already zero-initialized above), so if any
    // residual state exists, it will decay toward 0 with tau_TF.
    // KRAS, TP53, CDKN2A, SMAD4 — locked mutants, handled in update().
    // EGFR — not updated here; its state is set externally (drug intervention)
    //         or held at its initial value. A rule could be added if needed.
}

// ============================================================================
// compute_stroma_targets
//
// Boolean-logic rules for stromal cells (CAFs / pancreatic stellate cells).
// Stromal cells are activated from a quiescent PSC state to an active CAF
// state by paracrine TGF-beta and SHH signals from tumor cells.
// Activated CAFs produce ECM (collagen, HA), secrete autocrine TGF-beta,
// and contribute to the desmoplastic barrier.
// ============================================================================

void BooleanNetwork::compute_stroma_targets(double* target,
                                            double oxygen,
                                            double tgfb_local,
                                            double shh_local,
                                            double /*drug_local*/,
                                            const ThresholdConfig& cfg) const
{
    // Shorthand aliases for current stromal gene states
    const double ACTA2_ = gene_states[ACTA2];
    const double HIF1A_ = gene_states[HIF1A];
    const double GLI1_  = gene_states[GLI1];

    // ---- Apply activation thresholds to paracrine signals -------------------
    // Below the threshold, the signal is too weak to trigger receptor activation
    // (noise floor). Above the threshold, the effective signal ramps linearly
    // from 0 to its maximum. This mirrors the oxygen/HIF1A threshold pattern
    // and ensures that trace concentrations of TGF-beta or SHH do not
    // prematurely activate quiescent stellate cells.
    //
    // Thresholds are read from XML user_parameters via ThresholdConfig:
    //   cfg.tgfb_activation_threshold  (XML default: 0.1)
    //   cfg.shh_activation_threshold   (XML default: 0.05)

    const double effective_tgfb = (tgfb_local > cfg.tgfb_activation_threshold)
        ? (tgfb_local - cfg.tgfb_activation_threshold)
          / (1.0 - cfg.tgfb_activation_threshold)   // ramp [threshold,1] → [0,1]
        : 0.0;

    const double effective_shh = (shh_local > cfg.shh_activation_threshold)
        ? (shh_local - cfg.shh_activation_threshold)
          / (1.0 - cfg.shh_activation_threshold)     // ramp [threshold,1] → [0,1]
        : 0.0;

    // =========================================================================
    // ACTIVATION AXIS — quiescent PSC → activated CAF transition
    // =========================================================================

    // GLI1 — Hedgehog pathway transcription factor; SHH signal transducer.
    //
    // Biological basis:
    //   SHH secreted by tumor cells binds PTCH1 on stromal fibroblasts.
    //   This de-represses Smoothened (SMO), which activates the GLI1/2/3
    //   transcription factor cascade. GLI1 is both a direct target and an
    //   amplifier of the Hedgehog pathway (positive feedback via GLI1 binding
    //   to its own promoter). The linear relationship (0.8 * effective_shh)
    //   models the near-linear dose-response of GLI1 mRNA to SHH concentration
    //   observed in in vitro stellate cell assays.
    //
    // NOTE: effective_shh is zero below cfg.shh_activation_threshold (XML: 0.05),
    //   then ramps linearly to 1.0 at shh_local = 1.0. This prevents trace SHH
    //   from activating GLI1 in distant stromal cells.
    target[GLI1] = clamp(0.8 * effective_shh, 0.0, 1.0);

    // ACTA2 — alpha-smooth muscle actin; canonical CAF activation marker.
    //
    // Biological basis:
    //   TGF-beta (0.6 weight) is the dominant driver of PSC-to-CAF
    //   transdifferentiation, acting through both SMAD-dependent (pSMAD2/3)
    //   and SMAD-independent (RhoA, p38-MAPK) pathways to induce alpha-SMA.
    //   GLI1 (0.3 weight) reinforces activation: activated Hedgehog signaling
    //   synergizes with TGF-beta to promote fibroblast activation and
    //   collagen production in the pancreatic stroma (Hwang et al. 2012).
    //   Coefficient sum can reach 0.9 under maximal stimulation, leaving
    //   room below 1.0 to avoid saturation artefacts.
    //
    // NOTE: effective_tgfb is zero below cfg.tgfb_activation_threshold (XML: 0.1),
    //   then ramps linearly. This ensures quiescent PSCs remain quiescent until
    //   TGF-beta reaches a biologically meaningful concentration.
    target[ACTA2] = clamp(0.6 * effective_tgfb + 0.3 * GLI1_, 0.0, 1.0);

    // =========================================================================
    // ECM PRODUCTION — only operative when ACTA2 (CAF activation) is high
    // =========================================================================

    // HAS2 — hyaluronan synthase-2; HA ECM production.
    //
    // Biological basis:
    //   HAS2 expression is driven by activated CAF state (ACTA2) and amplified
    //   by hypoxia (HIF1A binds HAS2 promoter). Hyaluronic acid creates a
    //   hydrophilic, space-filling ECM component that increases interstitial
    //   pressure and impairs drug delivery. PEGPH20 targets this mechanism.
    target[HAS2] = clamp(0.7 * ACTA2_ + 0.2 * HIF1A_, 0.0, 1.0);

    // COL1A1 — collagen I alpha-1; fibrous desmoplastic ECM.
    //
    // Biological basis:
    //   Type I collagen is the dominant ECM component of the PDAC desmoplastic
    //   stroma (~80% of the fibrous content). Its deposition is a direct output
    //   of activated CAF state, driven primarily by TGF-beta/ACTA2 signaling.
    //   No independent HIF1A contribution here (collagen prolyl hydroxylation
    //   actually requires O2, so hypoxia may reduce collagen production — this
    //   simplification uses ACTA2 as the single gating input).
    target[COL1A1] = clamp(0.8 * ACTA2_, 0.0, 1.0);

    // MMP2 — matrix metalloproteinase-2; ECM remodeling.
    //
    // Biological basis:
    //   Activated CAFs secrete MMP2 to remodel the fibrous ECM and create
    //   invasion conduits for tumor cell migration. CAF-derived MMP2 is lower
    //   than tumor-cell MMP2 in the PDAC context (coefficient 0.3 vs 0.4+ in
    //   tumor rules), reflecting the primarily structural role of stromal MMP2.
    target[MMP2] = clamp(0.3 * ACTA2_, 0.0, 1.0);

    // =========================================================================
    // PARACRINE FEEDBACK — activated CAFs amplify the TGF-beta signal
    // =========================================================================

    // TGFB1 — autocrine/paracrine TGF-beta secretion from activated stroma.
    //
    // Biological basis:
    //   Activated CAFs are themselves a major source of TGF-beta in the tumor
    //   microenvironment, creating a positive feedback loop: tumor TGF-beta
    //   activates CAFs → CAF TGF-beta further activates other PSCs and
    //   suppresses anti-tumor immune responses. The 0.4 coefficient on ACTA2
    //   (lower than tumor's 0.3+0.1 baseline) reflects the secondary nature of
    //   stromal TGF-beta secretion relative to tumor-derived TGF-beta.
    target[TGFB1] = clamp(0.4 * ACTA2_, 0.0, 1.0);

    // =========================================================================
    // STRESS AXIS
    // =========================================================================

    // HIF1A — identical O2-sensing rule as tumor cells.
    //
    // Biological basis:
    //   Stromal cells in the PDAC core experience the same hypoxic environment
    //   as tumor cells. Stromal HIF1A activation further drives HAS2 (see above),
    //   creating a hypoxia → desmoplasia positive feedback distinct from the
    //   tumor cell hypoxia → TGF-beta → desmoplasia route.
    if (oxygen < cfg.hypoxia_hif1a_threshold)
        target[HIF1A] = clamp(1.0 - oxygen / cfg.hypoxia_hif1a_threshold, 0.0, 1.0);
    else
        target[HIF1A] = 0.0;

    // =========================================================================
    // STROMA-ABSENT GENES (tumor-specific)
    // =========================================================================
    // KRAS, MYC, EGFR, TP53, BCL_XL, BAX, CDKN2A, SMAD4, RB1,
    // ZEB1, CDH1, SHH, NRF2, ABCB1 — not expressed in stromal cells.
    // Their targets default to 0.0 (already zero-initialized), so residual
    // state (if any) decays to 0 over tau_TF ~ 60 min.
}

// ============================================================================
// apply_interventions
//
// Applies a vector of Intervention objects to gene_states[] after update().
//
// Effect semantics:
//
//   INHIBIT:
//     gene *= (1 - strength)
//     Multiplicative suppression. Models competitive inhibitors, blocking
//     antibodies, and partial knockdowns where some residual activity remains.
//     Strength 0.0 = no effect; 1.0 = complete elimination of current state.
//     Reversible: removing the intervention from the list restores dynamics.
//
//   ACTIVATE:
//     gene += strength * (1 - gene)
//     Multiplicative activation toward 1.0. Models agonists, activating
//     mutations introduced therapeutically, or forced overexpression where
//     the maximum achievable state is 1.0. Respects the ceiling naturally.
//     Reversible: dynamics resume from the elevated state on next update().
//
//   KNOCKDOWN:
//     gene = clamp(gene - strength, 0, 1)
//     is_mutant[i] = true  (permanent lock)
//     Models RNAi, CRISPR loss-of-function, or irreversible silencing.
//     Pins the gene at a reduced level and prevents further update() changes.
//     The EA can evolve strength to find the minimum effective knockdown.
//
//   OVEREXPRESS:
//     gene = clamp(gene + strength, 0, 1)
//     is_mutant[i] = true  (permanent lock)
//     Models stable transgene insertion, viral vector delivery, or CRISPR
//     activation (CRISPRa). Pins the gene at an elevated state permanently.
//
// IMPORTANT: Interventions on genomically locked mutant genes (KRAS, TP53,
// CDKN2A, SMAD4 in tumor cells) are silently skipped for INHIBIT/ACTIVATE
// since update() re-locks them anyway. KNOCKDOWN/OVEREXPRESS on mutant genes
// ARE applied (allowing therapeutic scenarios like partial KRAS suppression
// via direct KRAS inhibitors, which the EA may explore).
// ============================================================================

void BooleanNetwork::apply_interventions(const std::vector<Intervention>& interventions)
{
    for (const auto& iv : interventions)
    {
        const int i = static_cast<int>(iv.gene_index);

        // Bounds check (defensive)
        if (i < 0 || i >= GENE_COUNT) continue;

        // Clamp strength to [0, 1]
        const double s = clamp(iv.strength, 0.0, 1.0);

        switch (iv.effect_type)
        {
            // ---- Reversible: drug-like ----------------------------------------

            case EffectType::INHIBIT:
                // Skip if already a genomic mutant — update() would re-lock anyway.
                if (is_mutant[i]) break;
                gene_states[i] = clamp(gene_states[i] * (1.0 - s), 0.0, 1.0);
                break;

            case EffectType::ACTIVATE:
                // Skip if already a genomic mutant.
                if (is_mutant[i]) break;
                gene_states[i] = clamp(gene_states[i] + s * (1.0 - gene_states[i]),
                                       0.0, 1.0);
                break;

            // ---- Irreversible: gene-therapy-like --------------------------------

            case EffectType::KNOCKDOWN:
                // Lock the gene and reduce its state permanently.
                // Allows KRAS partial suppression (e.g., MRTX1133 direct KRAS G12D
                // inhibitor): even if KRAS is marked mutant, KNOCKDOWN overrides.
                gene_states[i] = clamp(gene_states[i] - s, 0.0, 1.0);
                is_mutant[i]   = true;
                break;

            case EffectType::OVEREXPRESS:
                gene_states[i] = clamp(gene_states[i] + s, 0.0, 1.0);
                is_mutant[i]   = true;
                break;
        }
    }
}

// ============================================================================
// load_from_json
//
// Parses an intervention JSON file generated by the Python EA and returns a
// vector of resolved Intervention structs ready for apply_interventions().
//
// JSON schema (all fields required per entry):
//   "gene"     : string matching GENE_INFO[i].name exactly (case-sensitive)
//   "effect"   : "INHIBIT" | "ACTIVATE" | "KNOCKDOWN" | "OVEREXPRESS"
//   "strength" : float in [0.0, 1.0]
//   "name"     : string label (any; used for logging only)
//
// The gene name → GeneIndex map is built lazily from GENE_INFO[] on first
// call and reused (static local). Thread-safety: safe if called before any
// parallel PhysiCell threads start (i.e., during setup, not during simulate).
//
// Throws std::runtime_error on:
//   - File not found / unreadable
//   - Malformed JSON (nlohmann exception propagated)
//   - Unknown gene name
//   - Unknown effect string
//   - Missing required JSON fields
// ============================================================================

std::vector<Intervention> BooleanNetwork::load_from_json(const std::string& filepath)
{
    // ---- Build gene-name → GeneIndex lookup table (built once) --------------
    // Constructed from the compile-time GENE_INFO array so it is always in
    // sync with the enum without manual maintenance.
    static const std::unordered_map<std::string, GeneIndex> NAME_TO_INDEX = []()
    {
        std::unordered_map<std::string, GeneIndex> m;
        m.reserve(GENE_COUNT);
        for (int i = 0; i < GENE_COUNT; ++i)
        {
            m.emplace(std::string(GENE_INFO[i].name), static_cast<GeneIndex>(i));
        }
        return m;
    }();

    // ---- Build effect-string → EffectType lookup table ---------------------
    static const std::unordered_map<std::string, EffectType> NAME_TO_EFFECT = {
        { "INHIBIT",     EffectType::INHIBIT     },
        { "ACTIVATE",    EffectType::ACTIVATE    },
        { "KNOCKDOWN",   EffectType::KNOCKDOWN   },
        { "OVEREXPRESS", EffectType::OVEREXPRESS },
    };

    // ---- Open and parse the JSON file --------------------------------------
    std::ifstream ifs(filepath);
    if (!ifs.is_open())
    {
        throw std::runtime_error(
            "BooleanNetwork::load_from_json: cannot open file: " + filepath);
    }

    json root;
    try
    {
        ifs >> root;
    }
    catch (const json::parse_error& e)
    {
        throw std::runtime_error(
            std::string("BooleanNetwork::load_from_json: JSON parse error in '")
            + filepath + "': " + e.what());
    }

    // ---- Validate top-level structure --------------------------------------
    if (!root.contains("interventions") || !root["interventions"].is_array())
    {
        throw std::runtime_error(
            "BooleanNetwork::load_from_json: missing or non-array 'interventions' "
            "key in: " + filepath);
    }

    // ---- Parse each intervention entry ------------------------------------
    std::vector<Intervention> result;
    result.reserve(root["interventions"].size());

    for (std::size_t idx = 0; idx < root["interventions"].size(); ++idx)
    {
        const json& entry = root["interventions"][idx];

        // Helper lambda: throw with context on missing field
        auto require = [&](const char* field)
        {
            if (!entry.contains(field))
                throw std::runtime_error(
                    std::string("BooleanNetwork::load_from_json: entry [")
                    + std::to_string(idx) + "] missing field '" + field
                    + "' in: " + filepath);
        };

        require("gene");
        require("effect");
        require("strength");
        require("name");

        // ---- Resolve gene name → GeneIndex ----------------------------------
        const std::string gene_name = entry["gene"].get<std::string>();
        auto git = NAME_TO_INDEX.find(gene_name);
        if (git == NAME_TO_INDEX.end())
        {
            throw std::runtime_error(
                std::string("BooleanNetwork::load_from_json: unknown gene '")
                + gene_name + "' at entry [" + std::to_string(idx)
                + "] in: " + filepath);
        }

        // ---- Resolve effect string → EffectType ----------------------------
        const std::string effect_str = entry["effect"].get<std::string>();
        auto eit = NAME_TO_EFFECT.find(effect_str);
        if (eit == NAME_TO_EFFECT.end())
        {
            throw std::runtime_error(
                std::string("BooleanNetwork::load_from_json: unknown effect '")
                + effect_str + "' at entry [" + std::to_string(idx)
                + "] — valid values: INHIBIT, ACTIVATE, KNOCKDOWN, OVEREXPRESS"
                + " in: " + filepath);
        }

        // ---- Parse and validate strength ------------------------------------
        const double strength = entry["strength"].get<double>();
        if (strength < 0.0 || strength > 1.0)
        {
            std::cerr << "[BooleanNetwork] WARNING: intervention ["
                      << idx << "] strength " << strength
                      << " out of [0,1]; clamping.\n";
        }

        // ---- Warn if targeting a non-druggable gene -------------------------
        if (!GENE_INFO[git->second].is_druggable)
        {
            std::cerr << "[BooleanNetwork] WARNING: intervention ["
                      << idx << "] targets non-druggable gene '"
                      << gene_name << "' — allowed but biologically unusual.\n";
        }

        result.push_back({
            git->second,
            eit->second,
            strength,
            entry["name"].get<std::string>()
        });
    }

    std::cerr << "[BooleanNetwork] Loaded " << result.size()
              << " intervention(s) from: " << filepath << "\n";

    return result;
}
