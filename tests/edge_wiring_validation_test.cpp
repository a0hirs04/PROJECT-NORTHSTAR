#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

bool nearly_equal(double a, double b, double eps = 1e-10)
{
    return std::fabs(a - b) <= eps;
}

int idx(Cell* pCell, const char* name)
{
    const int i = pCell->custom_data.find_variable_index(name);
    assert(i >= 0);
    return i;
}

double local_density(Cell* pCell, int substrate_index)
{
    std::vector<double> p = pCell->position;
    const int voxel = microenvironment.nearest_voxel_index(p);
    return microenvironment.density_vector(voxel)[substrate_index];
}

void set_local_density(Cell* pCell, int substrate_index, double value)
{
    std::vector<double> p = pCell->position;
    const int voxel = microenvironment.nearest_voxel_index(p);
    microenvironment.density_vector(voxel)[substrate_index] = value;
}

int voxel_for_cell(Cell* pCell)
{
    std::vector<double> p = pCell->position;
    return microenvironment.nearest_voxel_index(p);
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    assert(oxygen_index >= 0);
    assert(tgfb_index >= 0);
    assert(shh_index >= 0);
    assert(ecm_index >= 0);
    assert(drug_index >= 0);

    const double dt = 6.0;
    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.6);

    // ---------------------------------------------------------------------
    // E10 (load-bearing wall #1): CAF activation -> ECM production
    // ---------------------------------------------------------------------
    parameters.doubles("tgfb_activation_threshold") = 0.6;
    parameters.doubles("shh_activation_threshold") = 0.05;
    parameters.doubles("ecm_production_rate_base") = 0.01;
    parameters.doubles("ecm_production_rate_boosted") = 0.015;
    parameters.doubles("ecm_ha_fraction_default") = 0.6;
    parameters.doubles("mmp2_degradation_rate") = 0.0;

    Cell* stroma_e10 = create_cell(*pStroma);
    stroma_e10->assign_position(std::vector<double>{100.0, 100.0, 0.0});
    const int stroma_voxel = voxel_for_cell(stroma_e10);
    set_local_density(stroma_e10, tgfb_index, 0.5);
    set_local_density(stroma_e10, shh_index, 0.3);
    microenvironment.density_vector(stroma_voxel)[ecm_index] = 0.0;

    module3_stromal_activation(stroma_e10, stroma_e10->phenotype, dt, ModulePhase::SENSING);
    const int acta2_idx = idx(stroma_e10, "acta2_active");
    const int gli1_idx = idx(stroma_e10, "gli1_active");
    assert(nearly_equal(stroma_e10->custom_data[acta2_idx], 1.0));
    assert(nearly_equal(stroma_e10->custom_data[gli1_idx], 1.0));

    module6_ecm_production(stroma_e10, stroma_e10->phenotype, dt, ModulePhase::WRITE);
    const double ecm_after_acta2_on = microenvironment.density_vector(stroma_voxel)[ecm_index];
    assert(ecm_after_acta2_on > 0.0); // direction: activation increases ECM

    // Gating enforcement: no ACTA2 -> no stromal ECM production.
    stroma_e10->custom_data[acta2_idx] = 0.0;
    microenvironment.density_vector(stroma_voxel)[ecm_index] = 0.0;
    module6_ecm_production(stroma_e10, stroma_e10->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(microenvironment.density_vector(stroma_voxel)[ecm_index], 0.0));

    std::cout << "PASS E10 wiring/gating/direction" << std::endl;

    // ---------------------------------------------------------------------
    // E06 (load-bearing wall #2): TGF-beta field -> stromal activation
    // Source: Module2 secretion rates -> diffusion/field -> Module3 reads tgfb
    // ---------------------------------------------------------------------
    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.6);
    parameters.doubles("tgfb_secretion_rate") = 1.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.0;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.0;
    parameters.doubles("tgfb_activation_threshold") = 0.005;

    Cell* tumor_e06 = create_cell(*pTumor);
    tumor_e06->assign_position(std::vector<double>{200.0, 200.0, 0.0});
    Cell* stroma_near_e06 = create_cell(*pStroma);
    stroma_near_e06->assign_position(std::vector<double>{220.0, 200.0, 0.0});
    Cell* stroma_far_e06 = create_cell(*pStroma);
    stroma_far_e06->assign_position(std::vector<double>{800.0, 800.0, 0.0});

    const int near_acta2_idx = idx(stroma_near_e06, "acta2_active");
    const int far_acta2_idx = idx(stroma_far_e06, "acta2_active");
    stroma_near_e06->custom_data[near_acta2_idx] = 0.0;
    stroma_far_e06->custom_data[far_acta2_idx] = 0.0;
    tumor_e06->custom_data[idx(tumor_e06, "hif1a_active")] = 0.0;

    double t = 0.0;
    bool near_activated = false;
    for (int step = 0; step < 30; ++step)
    {
        advance_steps(1, dt, t);
        if (stroma_near_e06->custom_data[near_acta2_idx] == 1.0)
        {
            near_activated = true;
            break;
        }
    }

    const double tgfb_near = local_density(stroma_near_e06, tgfb_index);
    const double tgfb_far = local_density(stroma_far_e06, tgfb_index);
    assert(tgfb_near > 0.0);
    assert(tgfb_near > tgfb_far);                    // direction in space
    assert(near_activated); // field drives activation

    std::cout << "PASS E06 wiring/gating/direction" << std::endl;

    // ---------------------------------------------------------------------
    // E21 (load-bearing wall #3): ECM -> diffusion block
    // Source: Module6 writes ECM -> coupling reads ECM -> effective D decreases
    // ---------------------------------------------------------------------
    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.8);
    parameters.doubles("ecm_production_rate_base") = 0.02;
    parameters.doubles("ecm_production_rate_boosted") = 0.02;

    Cell* stroma_e21 = create_cell(*pStroma);
    stroma_e21->assign_position(std::vector<double>{300.0, 300.0, 0.0});
    const int e21_acta2 = idx(stroma_e21, "acta2_active");
    const int e21_gli1 = idx(stroma_e21, "gli1_active");
    stroma_e21->custom_data[e21_acta2] = 1.0;
    stroma_e21->custom_data[e21_gli1] = 0.0;
    const int e21_voxel = voxel_for_cell(stroma_e21);

    update_ecm_effective_diffusion_coefficients(microenvironment);
    const double drug_base = microenvironment.diffusion_coefficients[drug_index];
    const double d0 = get_effective_diffusion_coefficient(drug_index, e21_voxel);
    assert(nearly_equal(d0, drug_base));

    module6_ecm_production(stroma_e21, stroma_e21->phenotype, dt, ModulePhase::WRITE);
    const double ecm_after_first_write = microenvironment.density_vector(e21_voxel)[ecm_index];
    update_ecm_effective_diffusion_coefficients(microenvironment);
    const double d1 = get_effective_diffusion_coefficient(drug_index, e21_voxel);
    assert(ecm_after_first_write > 0.0);
    assert(d1 < d0); // direction: higher ECM lowers drug diffusion

    module6_ecm_production(stroma_e21, stroma_e21->phenotype, dt, ModulePhase::WRITE);
    update_ecm_effective_diffusion_coefficients(microenvironment);
    const double d2 = get_effective_diffusion_coefficient(drug_index, e21_voxel);
    assert(d2 <= d1); // monotone with added ECM

    std::cout << "PASS E21 wiring/gating/direction" << std::endl;

    // ---------------------------------------------------------------------
    // E08 / E09 split validation
    // E08: TGF-beta -> EMT is SMAD4-independent (never gated by SMAD4)
    // E09: TGF-beta -> growth arrest is SMAD4-dependent (only WT)
    // ---------------------------------------------------------------------
    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.6);
    parameters.doubles("emt_induction_threshold") = 0.4;
    parameters.ints("emt_phenotype_extent") = 2;
    parameters.doubles("hif1a_emt_boost") = 0.0;
    parameters.doubles("motility_epithelial") = 0.1;
    parameters.doubles("motility_mesenchymal_med") = 1.0;
    parameters.doubles("adhesion_epithelial") = 5.0;
    parameters.doubles("adhesion_mesenchymal_med") = 1.5;

    Cell* tumor_lof = create_cell(*pTumor);
    tumor_lof->assign_position(std::vector<double>{400.0, 400.0, 0.0});
    Cell* tumor_wt = create_cell(*pTumor);
    tumor_wt->assign_position(std::vector<double>{420.0, 400.0, 0.0});

    const int smad4_lof = idx(tumor_lof, "SMAD4");
    const int smad4_wt = idx(tumor_wt, "SMAD4");
    const int zeb1_lof = idx(tumor_lof, "zeb1_active");
    const int zeb1_wt = idx(tumor_wt, "zeb1_active");
    tumor_lof->custom_data[smad4_lof] = 0.0;
    tumor_wt->custom_data[smad4_wt] = 1.0;
    tumor_lof->custom_data[idx(tumor_lof, "hif1a_active")] = 0.0;
    tumor_wt->custom_data[idx(tumor_wt, "hif1a_active")] = 0.0;

    set_local_density(tumor_lof, tgfb_index, 0.6);
    set_local_density(tumor_wt, tgfb_index, 0.6);
    module5_emt_engine(tumor_lof, tumor_lof->phenotype, dt, ModulePhase::DECISION);
    module5_emt_engine(tumor_wt, tumor_wt->phenotype, dt, ModulePhase::DECISION);
    assert(nearly_equal(tumor_lof->custom_data[zeb1_lof], 1.0));
    assert(nearly_equal(tumor_wt->custom_data[zeb1_wt], 1.0));

    // E09 gating + direction in Module 4
    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("go_grow_penalty") = 0.0;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 1e9;
    parameters.doubles("tgfb_brake_sensitivity") = 0.5;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;
    parameters.doubles("hypoxia_death_resistance_bonus") = 0.0;

    tumor_lof->custom_data[idx(tumor_lof, "zeb1_active")] = 0.0;
    tumor_wt->custom_data[idx(tumor_wt, "zeb1_active")] = 0.0;
    tumor_lof->custom_data[idx(tumor_lof, "mechanical_pressure")] = 0.0;
    tumor_wt->custom_data[idx(tumor_wt, "mechanical_pressure")] = 0.0;
    tumor_lof->custom_data[idx(tumor_lof, "intracellular_drug")] = 0.0;
    tumor_wt->custom_data[idx(tumor_wt, "intracellular_drug")] = 0.0;
    tumor_lof->custom_data[idx(tumor_lof, "abcb1_active")] = 0.0;
    tumor_wt->custom_data[idx(tumor_wt, "abcb1_active")] = 0.0;

    set_local_density(tumor_lof, tgfb_index, 1.0);
    set_local_density(tumor_wt, tgfb_index, 1.0);
    module4_proliferation_death(tumor_lof, tumor_lof->phenotype, dt, ModulePhase::DECISION);
    module4_proliferation_death(tumor_wt, tumor_wt->phenotype, dt, ModulePhase::DECISION);

    const double prolif_lof_high_tgfb = tumor_lof->phenotype.cycle.data.transition_rate(0, 0);
    const double prolif_wt_high_tgfb = tumor_wt->phenotype.cycle.data.transition_rate(0, 0);
    assert(nearly_equal(prolif_lof_high_tgfb, 0.01)); // no E09 when SMAD4=LOF
    assert(prolif_wt_high_tgfb < prolif_lof_high_tgfb); // E09 active when SMAD4=WT

    set_local_density(tumor_wt, tgfb_index, 0.0);
    module4_proliferation_death(tumor_wt, tumor_wt->phenotype, dt, ModulePhase::DECISION);
    const double prolif_wt_zero_tgfb = tumor_wt->phenotype.cycle.data.transition_rate(0, 0);
    assert(prolif_wt_zero_tgfb > prolif_wt_high_tgfb); // direction: TGFB decreases proliferation

    std::cout << "PASS E08/E09 split wiring" << std::endl;
    std::cout << "PASS edge_wiring_validation_test" << std::endl;
    return 0;
}
