#include <cassert>
#include <iostream>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.0;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.0;

    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("emt_induction_threshold") = 10.0;

    microenvironment.decay_rates[tgfb_index] = 0.001;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double left_cx = -350.0;
    const double left_cy = 0.0;
    const double right_cx = 350.0;
    const double right_cy = 0.0;

    place_tumor_cluster(pTumor, 20, left_cx, left_cy, 40.0);   // with ECM ring
    place_tumor_cluster(pTumor, 20, right_cx, right_cy, 40.0);  // no ECM ring

    set_annulus_ecm(left_cx, left_cy, 60.0, 140.0, 0.8, 0.6);

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        pCell->phenotype.motility.is_motile = false;
    }

    double t = 0.0;
    const double dt = 1.0;
    advance_steps(100, dt, t);

    const double interface_with_ecm = mean_field_in_annulus(tgfb_index, left_cx, left_cy, 60.0, 140.0);
    const double interface_no_ecm = mean_field_in_annulus(tgfb_index, right_cx, right_cy, 60.0, 140.0);

    const double beyond_with_ecm = mean_field_in_annulus(tgfb_index, left_cx, left_cy, 170.0, 260.0);
    const double beyond_no_ecm = mean_field_in_annulus(tgfb_index, right_cx, right_cy, 170.0, 260.0);

    const bool pass_interface = (interface_with_ecm > interface_no_ecm);
    const bool pass_beyond = (beyond_with_ecm < beyond_no_ecm);

    std::cout << "P5E measurements"
              << " interface_with_ecm=" << interface_with_ecm
              << " interface_no_ecm=" << interface_no_ecm
              << " beyond_with_ecm=" << beyond_with_ecm
              << " beyond_no_ecm=" << beyond_no_ecm
              << std::endl;

    std::cout << "P5E checks"
              << " E23_interface_trapping=" << (pass_interface ? 1 : 0)
              << " E23_beyond_escape_reduced=" << (pass_beyond ? 1 : 0)
              << std::endl;

    if (!(pass_interface && pass_beyond))
    {
        std::cout << "FAIL P5 Cluster E" << std::endl;
        return 1;
    }

    std::cout << "PASS P5 Cluster E" << std::endl;
    return 0;
}
