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

bool nearly_equal(double a, double b, double eps = 1e-12)
{
    return std::fabs(a - b) <= eps;
}

} // namespace

int main()
{
    initialize_world();

    assert(tgfb_index >= 0);
    assert(shh_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.3;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.5;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    Cell* pCell = create_cell(*pTumor);
    pCell->assign_position(std::vector<double>{0.0, 0.0, 0.0});
    pCell->phenotype.motility.is_motile = false;

    const int hif_idx = custom_index(pCell, "hif1a_active");
    assert(hif_idx >= 0);

    // Condition A: HIF1A OFF.
    pCell->custom_data[hif_idx] = 0.0;
    module2_paracrine_secretion(pCell, pCell->phenotype, 6.0, ModulePhase::WRITE);
    const double tgfb_a = pCell->phenotype.secretion.secretion_rates[tgfb_index];
    const double shh_a = pCell->phenotype.secretion.secretion_rates[shh_index];

    // Condition B: HIF1A ON.
    pCell->custom_data[hif_idx] = 1.0;
    module2_paracrine_secretion(pCell, pCell->phenotype, 6.0, ModulePhase::WRITE);
    const double tgfb_b = pCell->phenotype.secretion.secretion_rates[tgfb_index];
    const double shh_b = pCell->phenotype.secretion.secretion_rates[shh_index];

    const double ratio = (tgfb_a > 0.0) ? (tgfb_b / tgfb_a) : 0.0;

    const bool pass_tgfb_increase = (tgfb_b > tgfb_a);
    const bool pass_ratio = nearly_equal(ratio, 1.5);
    const bool pass_shh_unchanged = nearly_equal(shh_a, shh_b) && nearly_equal(shh_a, 0.3);

    std::cout << "E02 measurements"
              << " tgfb_hif0=" << tgfb_a
              << " tgfb_hif1=" << tgfb_b
              << " tgfb_ratio=" << ratio
              << " shh_hif0=" << shh_a
              << " shh_hif1=" << shh_b
              << std::endl;

    std::cout << "E02 checks"
              << " tgfb_increases=" << (pass_tgfb_increase ? 1 : 0)
              << " ratio_matches_amp_factor=" << (pass_ratio ? 1 : 0)
              << " shh_unchanged=" << (pass_shh_unchanged ? 1 : 0)
              << std::endl;

    if (!(pass_tgfb_increase && pass_ratio && pass_shh_unchanged))
    {
        std::cout << "FAIL E02 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E02 integration" << std::endl;
    return 0;
}
