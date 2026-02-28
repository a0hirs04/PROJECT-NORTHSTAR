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

    assert(oxygen_index >= 0);
    assert(tgfb_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // E01 setup.
    parameters.doubles("hypoxia_response_threshold") = 0.5; // mapped threshold = 0.05
    parameters.doubles("emt_induction_threshold") = 0.4;
    parameters.doubles("hif1a_emt_boost") = 0.2;
    parameters.ints("emt_phenotype_extent") = 2;
    parameters.doubles("motility_epithelial") = 0.1;
    parameters.doubles("motility_mesenchymal_med") = 1.0;
    parameters.doubles("adhesion_epithelial") = 5.0;
    parameters.doubles("adhesion_mesenchymal_med") = 1.5;
    parameters.doubles("go_grow_penalty") = 0.0;

    reset_all_fields(0.08, 0.0, 0.0, 0.0, 0.0, 0.6);

    Cell* cell_a = create_cell(*pTumor);
    cell_a->assign_position(std::vector<double>{-200.0, 0.0, 0.0});
    cell_a->phenotype.motility.is_motile = false;

    Cell* cell_b = create_cell(*pTumor);
    cell_b->assign_position(std::vector<double>{200.0, 0.0, 0.0});
    cell_b->phenotype.motility.is_motile = false;

    // Same TGF-beta input; only oxygen differs.
    cell_a->nearest_density_vector()[tgfb_index] = 0.3;
    cell_b->nearest_density_vector()[tgfb_index] = 0.3;

    cell_a->nearest_density_vector()[oxygen_index] = 0.08; // Condition A (normoxic)
    cell_b->nearest_density_vector()[oxygen_index] = 0.01; // Condition B (hypoxic)

    module1_oxygen_sensing(cell_a, cell_a->phenotype, 6.0, ModulePhase::SENSING);
    module1_oxygen_sensing(cell_b, cell_b->phenotype, 6.0, ModulePhase::SENSING);

    module5_emt_engine(cell_a, cell_a->phenotype, 6.0, ModulePhase::DECISION);
    module5_emt_engine(cell_b, cell_b->phenotype, 6.0, ModulePhase::DECISION);

    const int hif_a = custom_index(cell_a, "hif1a_active");
    const int hif_b = custom_index(cell_b, "hif1a_active");
    const int zeb_a = custom_index(cell_a, "zeb1_active");
    const int zeb_b = custom_index(cell_b, "zeb1_active");
    assert(hif_a >= 0 && hif_b >= 0 && zeb_a >= 0 && zeb_b >= 0);

    const double induction_a = 0.3;
    const double induction_b = 0.3 + 0.2;

    const bool pass_a = nearly_equal(cell_a->custom_data[hif_a], 0.0) && nearly_equal(cell_a->custom_data[zeb_a], 0.0);
    const bool pass_b = nearly_equal(cell_b->custom_data[hif_b], 1.0) && nearly_equal(cell_b->custom_data[zeb_b], 1.0);

    std::cout << "E01 measurements"
              << " A_o2=0.08"
              << " A_hif1a=" << cell_a->custom_data[hif_a]
              << " A_induction_signal=" << induction_a
              << " A_zeb1=" << cell_a->custom_data[zeb_a]
              << " B_o2=0.01"
              << " B_hif1a=" << cell_b->custom_data[hif_b]
              << " B_induction_signal=" << induction_b
              << " B_zeb1=" << cell_b->custom_data[zeb_b]
              << std::endl;

    std::cout << "E01 checks"
              << " condA_no_emt=" << (pass_a ? 1 : 0)
              << " condB_emt_via_boost=" << (pass_b ? 1 : 0)
              << std::endl;

    if (!(pass_a && pass_b))
    {
        std::cout << "FAIL E01 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E01 integration" << std::endl;
    return 0;
}
