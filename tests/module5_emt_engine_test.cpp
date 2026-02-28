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

    parameters.doubles("emt_induction_threshold") = 0.4;
    parameters.ints("emt_phenotype_extent") = 2;
    parameters.doubles("hif1a_emt_boost") = 0.2;
    parameters.doubles("motility_epithelial") = 0.1;
    parameters.doubles("motility_mesenchymal_low") = 0.5;
    parameters.doubles("motility_mesenchymal_med") = 1.0;
    parameters.doubles("motility_mesenchymal_high") = 2.0;
    parameters.doubles("adhesion_epithelial") = 5.0;
    parameters.doubles("adhesion_mesenchymal_low") = 3.0;
    parameters.doubles("adhesion_mesenchymal_med") = 1.5;
    parameters.doubles("adhesion_mesenchymal_high") = 0.5;

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    // Test A — Tumor cell, signal above threshold, MEDIUM extent.
    Cell* tumor_a = create_cell(*pTumor);
    tumor_a->assign_position(std::vector<double>{100.0, 100.0, 0.0});

    const int hif1a_idx_a = tumor_a->custom_data.find_variable_index("hif1a_active");
    const int zeb1_idx_a = tumor_a->custom_data.find_variable_index("zeb1_active");
    const int cdh1_idx_a = tumor_a->custom_data.find_variable_index("cdh1_expressed");
    const int mmp2_idx_a = tumor_a->custom_data.find_variable_index("mmp2_active");
    assert(hif1a_idx_a >= 0);
    assert(zeb1_idx_a >= 0);
    assert(cdh1_idx_a >= 0);
    assert(mmp2_idx_a >= 0);

    std::vector<double>& rho_a = tumor_a->nearest_density_vector();
    rho_a[tgfb_index] = 0.5;
    tumor_a->custom_data[hif1a_idx_a] = 0.0;

    module5_emt_engine(tumor_a, tumor_a->phenotype, 1.0, ModulePhase::DECISION);

    assert(nearly_equal(tumor_a->custom_data[zeb1_idx_a], 1.0));
    assert(nearly_equal(tumor_a->custom_data[cdh1_idx_a], 0.0));
    assert(nearly_equal(tumor_a->custom_data[mmp2_idx_a], 1.0));
    assert(nearly_equal(tumor_a->phenotype.motility.migration_speed, 1.0));
    assert(nearly_equal(tumor_a->phenotype.mechanics.cell_cell_adhesion_strength, 1.5));
    std::cout << "PASS Test A" << std::endl;

    // Test B — Same cell, signal drops below threshold (reversibility).
    rho_a[tgfb_index] = 0.1;
    tumor_a->custom_data[hif1a_idx_a] = 0.0;
    module5_emt_engine(tumor_a, tumor_a->phenotype, 1.0, ModulePhase::DECISION);

    assert(nearly_equal(tumor_a->custom_data[zeb1_idx_a], 0.0));
    assert(nearly_equal(tumor_a->custom_data[cdh1_idx_a], 1.0));
    assert(nearly_equal(tumor_a->custom_data[mmp2_idx_a], 0.0));
    assert(nearly_equal(tumor_a->phenotype.motility.migration_speed, 0.1));
    std::cout << "PASS Test B" << std::endl;

    // Test C — Hypoxia boosts induction past threshold.
    Cell* tumor_c = create_cell(*pTumor);
    tumor_c->assign_position(std::vector<double>{140.0, 100.0, 0.0});
    const int hif1a_idx_c = tumor_c->custom_data.find_variable_index("hif1a_active");
    const int zeb1_idx_c = tumor_c->custom_data.find_variable_index("zeb1_active");
    assert(hif1a_idx_c >= 0);
    assert(zeb1_idx_c >= 0);

    std::vector<double>& rho_c = tumor_c->nearest_density_vector();
    rho_c[tgfb_index] = 0.3;
    tumor_c->custom_data[hif1a_idx_c] = 1.0;
    module5_emt_engine(tumor_c, tumor_c->phenotype, 1.0, ModulePhase::DECISION);

    assert(nearly_equal(tumor_c->custom_data[zeb1_idx_c], 1.0));
    std::cout << "PASS Test C" << std::endl;

    // Test D — Stromal cell is ignored.
    Cell* stroma_d = create_cell(*pStroma);
    stroma_d->assign_position(std::vector<double>{180.0, 100.0, 0.0});
    const int zeb1_idx_d = stroma_d->custom_data.find_variable_index("zeb1_active");
    const int mmp2_idx_d = stroma_d->custom_data.find_variable_index("mmp2_active");
    assert(zeb1_idx_d >= 0);
    assert(mmp2_idx_d >= 0);

    const double zeb1_before = stroma_d->custom_data[zeb1_idx_d];
    const double mmp2_before = stroma_d->custom_data[mmp2_idx_d];
    std::vector<double>& rho_d = stroma_d->nearest_density_vector();
    rho_d[tgfb_index] = 0.9;

    module5_emt_engine(stroma_d, stroma_d->phenotype, 1.0, ModulePhase::DECISION);

    assert(nearly_equal(stroma_d->custom_data[zeb1_idx_d], zeb1_before));
    assert(nearly_equal(stroma_d->custom_data[mmp2_idx_d], mmp2_before));
    std::cout << "PASS Test D" << std::endl;

    // Test E — LOW extent does NOT activate MMP2.
    parameters.ints("emt_phenotype_extent") = 1;
    Cell* tumor_e = create_cell(*pTumor);
    tumor_e->assign_position(std::vector<double>{220.0, 100.0, 0.0});
    const int hif1a_idx_e = tumor_e->custom_data.find_variable_index("hif1a_active");
    const int zeb1_idx_e = tumor_e->custom_data.find_variable_index("zeb1_active");
    const int mmp2_idx_e = tumor_e->custom_data.find_variable_index("mmp2_active");
    assert(hif1a_idx_e >= 0);
    assert(zeb1_idx_e >= 0);
    assert(mmp2_idx_e >= 0);

    std::vector<double>& rho_e = tumor_e->nearest_density_vector();
    rho_e[tgfb_index] = 0.6;
    tumor_e->custom_data[hif1a_idx_e] = 0.0;
    module5_emt_engine(tumor_e, tumor_e->phenotype, 1.0, ModulePhase::DECISION);

    assert(nearly_equal(tumor_e->custom_data[zeb1_idx_e], 1.0));
    assert(nearly_equal(tumor_e->custom_data[mmp2_idx_e], 0.0));
    std::cout << "PASS Test E" << std::endl;

    std::cout << "PASS module5_emt_engine_test" << std::endl;
    return 0;
}
