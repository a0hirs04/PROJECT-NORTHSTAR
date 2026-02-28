#include <cassert>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

void set_left_boundary_drug_source(double value)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        if (c[0] <= -900.0)
        {
            microenvironment.density_vector(n)[drug_index] = value;
        }
    }
}

void set_cohort_abcb1(const std::unordered_set<int>& ids, double value)
{
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        if (ids.find(pCell->ID) == ids.end()) continue;
        const int abcb1_idx = custom_index(pCell, "abcb1_active");
        if (abcb1_idx >= 0) pCell->custom_data[abcb1_idx] = value;
    }
}

struct Stats
{
    int live_n = 0;
    double intracellular_mean = 0.0;
    double apoptosis_mean = 0.0;
    int n = 0;
};

Stats summarize(const std::unordered_set<int>& ids)
{
    Stats s;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        if (ids.find(pCell->ID) == ids.end()) continue;
        if (pCell->phenotype.death.dead) continue;

        ++s.live_n;
        const int intra_idx = custom_index(pCell, "intracellular_drug");
        const double intra = (intra_idx >= 0) ? pCell->custom_data[intra_idx] : 0.0;
        const int apop_idx = pCell->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
        const double apop = (apop_idx >= 0) ? pCell->phenotype.death.rates[apop_idx] : 0.0;

        s.intracellular_mean += intra;
        s.apoptosis_mean += apop;
        ++s.n;
    }

    if (s.n > 0)
    {
        s.intracellular_mean /= static_cast<double>(s.n);
        s.apoptosis_mean /= static_cast<double>(s.n);
    }
    return s;
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);
    assert(drug_index >= 0);

    // Isolate E25: manually force ABCB1 states.
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.05;
    parameters.doubles("efflux_drug_reduction") = 0.7;

    parameters.doubles("drug_uptake_rate") = 0.03;
    parameters.doubles("drug_stress_threshold") = 2.0;
    parameters.doubles("efflux_induction_delay") = 1e9;
    parameters.doubles("efflux_strength") = 0.1;
    parameters.doubles("nrf2_decay_rate") = 0.0;
    parameters.doubles("drug_natural_decay_rate") = 0.0;
    parameters.doubles("hif1a_nrf2_priming_bonus") = 0.0;

    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 1e9;
    parameters.doubles("emt_induction_threshold") = 10.0;

    microenvironment.diffusion_coefficients[drug_index] = 200.0;
    microenvironment.decay_rates[drug_index] = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    // A: drug + ABCB1=0, B: drug + ABCB1=1, C: no drug + ABCB1=1, D: no drug + ABCB1=0.
    place_tumor_cluster(pTumor, 20, -880.0, -90.0, 25.0); // A
    place_tumor_cluster(pTumor, 20, -880.0, 90.0, 25.0);  // B
    place_tumor_cluster(pTumor, 20, 700.0, -90.0, 35.0);  // C
    place_tumor_cluster(pTumor, 20, 700.0, 90.0, 35.0);   // D

    std::unordered_set<int> A_ids;
    std::unordered_set<int> B_ids;
    std::unordered_set<int> C_ids;
    std::unordered_set<int> D_ids;

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->type != pTumor->type) continue;
        pCell->phenotype.motility.is_motile = false;

        if (pCell->position[0] < 0.0)
        {
            if (pCell->position[1] < 0.0) A_ids.insert(pCell->ID);
            else B_ids.insert(pCell->ID);
        }
        else
        {
            if (pCell->position[1] < 0.0) C_ids.insert(pCell->ID);
            else D_ids.insert(pCell->ID);
        }
    }

    set_cohort_abcb1(A_ids, 0.0);
    set_cohort_abcb1(B_ids, 1.0);
    set_cohort_abcb1(C_ids, 1.0);
    set_cohort_abcb1(D_ids, 0.0);

    double t = 0.0;
    const double dt = 1.0;
    Cell_Container* cc = cell_container();

    Stats A50, B50;
    Stats A100, B100, C100, D100;

    for (int step = 1; step <= 100; ++step)
    {
        set_left_boundary_drug_source(1.0);
        microenvironment.simulate_diffusion_decay(dt);
        cc->update_all_cells(t, dt, dt, dt);
        t += dt;

        // Preserve forced ABCB1 states throughout run.
        set_cohort_abcb1(A_ids, 0.0);
        set_cohort_abcb1(B_ids, 1.0);
        set_cohort_abcb1(C_ids, 1.0);
        set_cohort_abcb1(D_ids, 0.0);

        if (step == 50)
        {
            A50 = summarize(A_ids);
            B50 = summarize(B_ids);
        }
    }

    A100 = summarize(A_ids);
    B100 = summarize(B_ids);
    C100 = summarize(C_ids);
    D100 = summarize(D_ids);

    const bool pass_survival = (B100.live_n > A100.live_n);
    const bool pass_intra = (B50.intracellular_mean < A50.intracellular_mean);
    const bool pass_apop = (B50.apoptosis_mean < A50.apoptosis_mean);
    const bool pass_no_drug_inert = (C100.live_n == D100.live_n);

    std::cout << "E25 measurements"
              << " A50_intra=" << A50.intracellular_mean
              << " B50_intra=" << B50.intracellular_mean
              << " A50_apop=" << A50.apoptosis_mean
              << " B50_apop=" << B50.apoptosis_mean
              << " A100_live=" << A100.live_n
              << " B100_live=" << B100.live_n
              << " C100_live_noDrug_abcb1On=" << C100.live_n
              << " D100_live_noDrug_abcb1Off=" << D100.live_n
              << std::endl;

    std::cout << "E25 checks"
              << " B_survives_gt_A=" << (pass_survival ? 1 : 0)
              << " B_intra_lt_A=" << (pass_intra ? 1 : 0)
              << " B_apop_lt_A=" << (pass_apop ? 1 : 0)
              << " abcb1_inert_without_drug=" << (pass_no_drug_inert ? 1 : 0)
              << std::endl;

    if (!(pass_survival && pass_intra && pass_apop && pass_no_drug_inert))
    {
        std::cout << "FAIL E25 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E25 integration" << std::endl;
    return 0;
}
