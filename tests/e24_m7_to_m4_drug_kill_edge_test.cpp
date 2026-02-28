#include <algorithm>
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

struct CohortSummary
{
    int live_count = 0;
    double apoptosis_mean = 0.0;
    double intracellular_mean = 0.0;
    double apoptosis_high_intra_mean = 0.0;
    double apoptosis_low_intra_mean = 0.0;
    int n = 0;
};

CohortSummary summarize(const std::unordered_set<int>& ids)
{
    CohortSummary s;
    std::vector<double> intra;
    std::vector<double> apop;

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        if (ids.find(pCell->ID) == ids.end()) continue;
        if (pCell->phenotype.death.dead) continue;

        ++s.live_count;
        const int intra_idx = custom_index(pCell, "intracellular_drug");
        const double intra_val = (intra_idx >= 0) ? pCell->custom_data[intra_idx] : 0.0;
        const int apop_idx = pCell->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
        const double apop_val = (apop_idx >= 0) ? pCell->phenotype.death.rates[apop_idx] : 0.0;

        intra.push_back(intra_val);
        apop.push_back(apop_val);
        s.intracellular_mean += intra_val;
        s.apoptosis_mean += apop_val;
    }

    s.n = static_cast<int>(intra.size());
    if (s.n > 0)
    {
        s.intracellular_mean /= static_cast<double>(s.n);
        s.apoptosis_mean /= static_cast<double>(s.n);

        std::vector<double> sorted_intra = intra;
        std::sort(sorted_intra.begin(), sorted_intra.end());
        const double median = sorted_intra[s.n / 2];

        double high_sum = 0.0;
        double low_sum = 0.0;
        int high_n = 0;
        int low_n = 0;
        for (int i = 0; i < s.n; ++i)
        {
            if (intra[i] >= median)
            {
                high_sum += apop[i];
                ++high_n;
            }
            else
            {
                low_sum += apop[i];
                ++low_n;
            }
        }
        s.apoptosis_high_intra_mean = (high_n > 0) ? (high_sum / static_cast<double>(high_n)) : 0.0;
        s.apoptosis_low_intra_mean = (low_n > 0) ? (low_sum / static_cast<double>(low_n)) : 0.0;
    }

    return s;
}

void set_abcb1_for_cohort(const std::unordered_set<int>& ids, double value)
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

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);
    assert(drug_index >= 0);

    // Isolate E24 (drug kill path M7 -> M4).
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.2;
    parameters.doubles("efflux_drug_reduction") = 0.7;

    parameters.doubles("drug_uptake_rate") = 0.08;
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

    // Left cohort: drug-exposed. Right cohort: no-drug control (far from source).
    place_tumor_cluster(pTumor, 20, -880.0, 0.0, 25.0);
    place_tumor_cluster(pTumor, 20, 700.0, 0.0, 40.0);

    std::unordered_set<int> drug_ids;
    std::unordered_set<int> ctrl_ids;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->type != pTumor->type) continue;
        pCell->phenotype.motility.is_motile = false;

        if (pCell->position[0] < 0.0) drug_ids.insert(pCell->ID);
        else ctrl_ids.insert(pCell->ID);
    }

    set_abcb1_for_cohort(drug_ids, 0.0);
    set_abcb1_for_cohort(ctrl_ids, 0.0);

    const int init_drug = static_cast<int>(drug_ids.size());
    const int init_ctrl = static_cast<int>(ctrl_ids.size());

    double t = 0.0;
    const double dt = 1.0;
    Cell_Container* cc = cell_container();

    CohortSummary step10_drug;
    CohortSummary step10_ctrl;
    CohortSummary step25_drug;
    CohortSummary step25_ctrl;
    CohortSummary step100_drug;
    CohortSummary step100_ctrl;

    for (int step = 1; step <= 100; ++step)
    {
        set_left_boundary_drug_source(2.0);
        microenvironment.simulate_diffusion_decay(dt);
        cc->update_all_cells(t, dt, dt, dt);
        t += dt;

        if (step == 10)
        {
            step10_drug = summarize(drug_ids);
            step10_ctrl = summarize(ctrl_ids);
        }
        if (step == 25)
        {
            step25_drug = summarize(drug_ids);
            step25_ctrl = summarize(ctrl_ids);
        }
    }

    step100_drug = summarize(drug_ids);
    step100_ctrl = summarize(ctrl_ids);

    const bool pass_drug_decrease = (step100_drug.live_count < init_drug);
    const bool pass_ctrl_stable = (step100_ctrl.live_count >= init_ctrl);
    const bool pass_apop_higher =
        (step10_drug.n > 0) &&
        (step10_ctrl.n > 0) &&
        (step10_drug.apoptosis_mean > step10_ctrl.apoptosis_mean);
    const bool pass_proportional =
        (step10_drug.n > 2) &&
        (step10_drug.apoptosis_high_intra_mean > step10_drug.apoptosis_low_intra_mean);

    std::cout << "E24 measurements"
              << " init_drug_n=" << init_drug
              << " step100_drug_n=" << step100_drug.live_count
              << " init_ctrl_n=" << init_ctrl
              << " step100_ctrl_n=" << step100_ctrl.live_count
              << " step10_drug_n=" << step10_drug.live_count
              << " step10_ctrl_n=" << step10_ctrl.live_count
              << " apop_mean_drug_step10=" << step10_drug.apoptosis_mean
              << " apop_mean_ctrl_step10=" << step10_ctrl.apoptosis_mean
              << " intra_mean_drug_step10=" << step10_drug.intracellular_mean
              << " intra_mean_ctrl_step10=" << step10_ctrl.intracellular_mean
              << " apop_high_intra_mean_step10=" << step10_drug.apoptosis_high_intra_mean
              << " apop_low_intra_mean_step10=" << step10_drug.apoptosis_low_intra_mean
              << " step25_drug_n=" << step25_drug.live_count
              << " step25_ctrl_n=" << step25_ctrl.live_count
              << " apop_mean_drug_step25=" << step25_drug.apoptosis_mean
              << " apop_mean_ctrl_step25=" << step25_ctrl.apoptosis_mean
              << " intra_mean_drug_step25=" << step25_drug.intracellular_mean
              << " intra_mean_ctrl_step25=" << step25_ctrl.intracellular_mean
              << " apop_high_intra_mean_step25=" << step25_drug.apoptosis_high_intra_mean
              << " apop_low_intra_mean_step25=" << step25_drug.apoptosis_low_intra_mean
              << std::endl;

    std::cout << "E24 checks"
              << " drug_cells_decrease=" << (pass_drug_decrease ? 1 : 0)
              << " ctrl_stable_or_growth=" << (pass_ctrl_stable ? 1 : 0)
              << " apop_drug_gt_ctrl=" << (pass_apop_higher ? 1 : 0)
              << " apop_proportional_to_intracellular=" << (pass_proportional ? 1 : 0)
              << std::endl;

    if (!(pass_drug_decrease && pass_ctrl_stable && pass_apop_higher && pass_proportional))
    {
        std::cout << "FAIL E24 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E24 integration" << std::endl;
    return 0;
}
