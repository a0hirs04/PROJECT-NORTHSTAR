#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../Stroma_world/PhysiCell/core/PhysiCell.h"
#include "../Stroma_world/PhysiCell/modules/PhysiCell_standard_modules.h"
#include "../custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

namespace
{

bool nearly_equal(double a, double b, double eps = 1e-12)
{
    return std::fabs(a - b) <= eps;
}

} // namespace

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    assert(tgfb_index >= 0);
    assert(shh_index >= 0);

    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.3;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.5;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.3;

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    // Test A — Tumor cell, normoxic.
    Cell* tumor = create_cell(*pTumor);
    tumor->assign_position(std::vector<double>{100.0, 100.0, 0.0});
    const int hif1a_idx = tumor->custom_data.find_variable_index("hif1a_active");
    assert(hif1a_idx >= 0);
    tumor->custom_data[hif1a_idx] = 0.0;

    module2_paracrine_secretion(tumor, tumor->phenotype, 1.0, ModulePhase::WRITE);
    assert(nearly_equal(tumor->phenotype.secretion.secretion_rates[tgfb_index], 0.5));
    assert(nearly_equal(tumor->phenotype.secretion.secretion_rates[shh_index], 0.3));
    std::cout << "PASS Test A" << std::endl;

    // Test B — Tumor cell, hypoxic.
    tumor->custom_data[hif1a_idx] = 1.0;
    module2_paracrine_secretion(tumor, tumor->phenotype, 1.0, ModulePhase::WRITE);
    assert(nearly_equal(tumor->phenotype.secretion.secretion_rates[tgfb_index], 0.75));
    assert(nearly_equal(tumor->phenotype.secretion.secretion_rates[shh_index], 0.3));
    std::cout << "PASS Test B" << std::endl;

    // Test C — Activated CAF.
    Cell* caf = create_cell(*pStroma);
    caf->assign_position(std::vector<double>{140.0, 100.0, 0.0});
    const int acta2_idx_c = caf->custom_data.find_variable_index("acta2_active");
    assert(acta2_idx_c >= 0);
    caf->custom_data[acta2_idx_c] = 1.0;

    module2_paracrine_secretion(caf, caf->phenotype, 1.0, ModulePhase::WRITE);
    assert(nearly_equal(caf->phenotype.secretion.secretion_rates[tgfb_index], 0.3));
    assert(nearly_equal(caf->phenotype.secretion.secretion_rates[shh_index], 0.0));
    std::cout << "PASS Test C" << std::endl;

    // Test D — Quiescent PSC.
    Cell* psc = create_cell(*pStroma);
    psc->assign_position(std::vector<double>{180.0, 100.0, 0.0});
    const int acta2_idx_d = psc->custom_data.find_variable_index("acta2_active");
    assert(acta2_idx_d >= 0);
    psc->custom_data[acta2_idx_d] = 0.0;

    module2_paracrine_secretion(psc, psc->phenotype, 1.0, ModulePhase::WRITE);
    assert(nearly_equal(psc->phenotype.secretion.secretion_rates[tgfb_index], 0.0));
    assert(nearly_equal(psc->phenotype.secretion.secretion_rates[shh_index], 0.0));
    std::cout << "PASS Test D" << std::endl;

    std::cout << "PASS module2_paracrine_secretion_test" << std::endl;
    return 0;
}
