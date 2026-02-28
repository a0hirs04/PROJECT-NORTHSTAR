#include <cassert>
#include <iostream>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // E17/E18 setup.
    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("go_grow_penalty") = 0.3;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 1e9;
    parameters.doubles("apoptosis_resistance") = 0.85;
    parameters.doubles("hypoxia_death_resistance_bonus") = 0.0;
    parameters.doubles("drug_kill_coefficient") = 0.01;
    parameters.doubles("efflux_drug_reduction") = 0.0;

    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    Cell* zeb_off = create_cell(*pTumor);
    zeb_off->assign_position(std::vector<double>{-100.0, 0.0, 0.0});
    zeb_off->phenotype.motility.is_motile = false;

    Cell* zeb_on = create_cell(*pTumor);
    zeb_on->assign_position(std::vector<double>{100.0, 0.0, 0.0});
    zeb_on->phenotype.motility.is_motile = false;

    zeb_off->custom_data[custom_index(zeb_off, "zeb1_active")] = 0.0;
    zeb_on->custom_data[custom_index(zeb_on, "zeb1_active")] = 1.0;

    zeb_off->custom_data[custom_index(zeb_off, "hif1a_active")] = 0.0;
    zeb_on->custom_data[custom_index(zeb_on, "hif1a_active")] = 0.0;

    zeb_off->custom_data[custom_index(zeb_off, "abcb1_active")] = 0.0;
    zeb_on->custom_data[custom_index(zeb_on, "abcb1_active")] = 0.0;

    zeb_off->custom_data[custom_index(zeb_off, "intracellular_drug")] = 0.5;
    zeb_on->custom_data[custom_index(zeb_on, "intracellular_drug")] = 0.5;

    zeb_off->custom_data[custom_index(zeb_off, "mechanical_pressure")] = 0.0;
    zeb_on->custom_data[custom_index(zeb_on, "mechanical_pressure")] = 0.0;

    module4_proliferation_death(zeb_off, zeb_off->phenotype, 1.0, ModulePhase::DECISION);
    module4_proliferation_death(zeb_on, zeb_on->phenotype, 1.0, ModulePhase::DECISION);

    const double prolif_off = zeb_off->phenotype.cycle.data.transition_rate(0, 0);
    const double prolif_on = zeb_on->phenotype.cycle.data.transition_rate(0, 0);

    const int apop_idx_off = zeb_off->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    const int apop_idx_on = zeb_on->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    const double apop_off = zeb_off->phenotype.death.rates[apop_idx_off];
    const double apop_on = zeb_on->phenotype.death.rates[apop_idx_on];

    const bool pass_e17 = (prolif_on < prolif_off);
    const bool pass_e18 = (apop_on < apop_off);
    const bool pass_tradeoff = pass_e17 && pass_e18;

    std::cout << "P5D measurements"
              << " prolif_zeb_off=" << prolif_off
              << " prolif_zeb_on=" << prolif_on
              << " apop_zeb_off=" << apop_off
              << " apop_zeb_on=" << apop_on
              << std::endl;

    std::cout << "P5D checks"
              << " E17_prolif_zeb_on_lt_off=" << (pass_e17 ? 1 : 0)
              << " E18_apop_zeb_on_lt_off=" << (pass_e18 ? 1 : 0)
              << " E17E18_simultaneous=" << (pass_tradeoff ? 1 : 0)
              << std::endl;

    if (!pass_tradeoff)
    {
        std::cout << "FAIL P5 Cluster D" << std::endl;
        return 1;
    }

    std::cout << "PASS P5 Cluster D" << std::endl;
    return 0;
}
