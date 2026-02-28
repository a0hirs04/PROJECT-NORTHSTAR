#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
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

struct ConditionResult
{
    std::string label;
    std::string smad4;
    double tgfb = 0.0;
    double proliferation = 0.0;
    double zeb1 = 0.0;
};

ConditionResult run_condition(Cell* pCell,
                              const std::string& label,
                              const std::string& smad4_label,
                              double tgfb)
{
    pCell->nearest_density_vector()[tgfb_index] = tgfb;

    pCell->custom_data[custom_index(pCell, "hif1a_active")] = 0.0;
    pCell->custom_data[custom_index(pCell, "mechanical_pressure")] = 0.0;
    pCell->custom_data[custom_index(pCell, "abcb1_active")] = 0.0;
    pCell->custom_data[custom_index(pCell, "intracellular_drug")] = 0.0;

    module5_emt_engine(pCell, pCell->phenotype, 6.0, ModulePhase::DECISION);
    module4_proliferation_death(pCell, pCell->phenotype, 6.0, ModulePhase::DECISION);

    ConditionResult out;
    out.label = label;
    out.smad4 = smad4_label;
    out.tgfb = tgfb;
    out.proliferation = pCell->phenotype.cycle.data.transition_rate(0, 0);
    out.zeb1 = pCell->custom_data[custom_index(pCell, "zeb1_active")];
    return out;
}

} // namespace

int main()
{
    initialize_world();

    assert(tgfb_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // E08 + E09 isolation.
    parameters.doubles("emt_induction_threshold") = 0.4;
    parameters.ints("emt_phenotype_extent") = 2;
    parameters.doubles("hif1a_emt_boost") = 0.0;
    parameters.doubles("motility_epithelial") = 0.1;
    parameters.doubles("motility_mesenchymal_med") = 1.0;
    parameters.doubles("adhesion_epithelial") = 5.0;
    parameters.doubles("adhesion_mesenchymal_med") = 1.5;

    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("go_grow_penalty") = 0.0;              // isolate E09 from ZEB1 go-grow
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 1e9; // disable contact arrest
    parameters.doubles("tgfb_brake_sensitivity") = 0.5;

    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;
    parameters.doubles("hypoxia_death_resistance_bonus") = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    // Four primary conditions.
    Cell* A = create_cell(*pTumor);
    A->assign_position(std::vector<double>{-300.0, 0.0, 0.0});
    Cell* B = create_cell(*pTumor);
    B->assign_position(std::vector<double>{-100.0, 0.0, 0.0});
    Cell* C = create_cell(*pTumor);
    C->assign_position(std::vector<double>{100.0, 0.0, 0.0});
    Cell* D = create_cell(*pTumor);
    D->assign_position(std::vector<double>{300.0, 0.0, 0.0});

    A->custom_data[custom_index(A, "SMAD4")] = 1.0; // WT
    B->custom_data[custom_index(B, "SMAD4")] = 1.0; // WT
    C->custom_data[custom_index(C, "SMAD4")] = 0.0; // LOF
    D->custom_data[custom_index(D, "SMAD4")] = 0.0; // LOF

    const ConditionResult rA = run_condition(A, "A", "WT", 0.0);
    const ConditionResult rB = run_condition(B, "B", "WT", 0.6);
    const ConditionResult rC = run_condition(C, "C", "LOF", 0.0);
    const ConditionResult rD = run_condition(D, "D", "LOF", 0.6);

    const double base_rate = parameters.doubles("base_proliferation_rate");

    const bool pass_A = nearly_equal(rA.proliferation, base_rate) && nearly_equal(rA.zeb1, 0.0);
    const bool pass_B = (rB.proliferation < base_rate) && nearly_equal(rB.zeb1, 1.0);
    const bool pass_C = nearly_equal(rC.proliferation, base_rate) && nearly_equal(rC.zeb1, 0.0);
    const bool pass_D = nearly_equal(rD.proliferation, base_rate) && nearly_equal(rD.zeb1, 1.0);

    // Mandatory negative sweep.
    std::vector<double> lof_tgfb = {0.3, 0.6, 0.9};
    std::vector<double> lof_prolif;
    Cell* lof_sweep = create_cell(*pTumor);
    lof_sweep->assign_position(std::vector<double>{500.0, 0.0, 0.0});
    lof_sweep->custom_data[custom_index(lof_sweep, "SMAD4")] = 0.0;

    for (size_t i = 0; i < lof_tgfb.size(); ++i)
    {
        const ConditionResult rr = run_condition(lof_sweep, "LOF_sweep", "LOF", lof_tgfb[i]);
        lof_prolif.push_back(rr.proliferation);
    }

    std::vector<double> wt_tgfb = {0.3, 0.6, 0.9};
    std::vector<double> wt_prolif;
    Cell* wt_sweep = create_cell(*pTumor);
    wt_sweep->assign_position(std::vector<double>{700.0, 0.0, 0.0});
    wt_sweep->custom_data[custom_index(wt_sweep, "SMAD4")] = 1.0;

    for (size_t i = 0; i < wt_tgfb.size(); ++i)
    {
        const ConditionResult rr = run_condition(wt_sweep, "WT_sweep", "WT", wt_tgfb[i]);
        wt_prolif.push_back(rr.proliferation);
    }

    bool pass_lof_no_brake = true;
    for (size_t i = 0; i < lof_prolif.size(); ++i)
    {
        if (!nearly_equal(lof_prolif[i], base_rate)) pass_lof_no_brake = false;
    }

    bool pass_wt_has_brake = true;
    for (size_t i = 0; i < wt_prolif.size(); ++i)
    {
        if (!(wt_prolif[i] < base_rate)) pass_wt_has_brake = false;
    }

    const bool pass =
        pass_A &&
        pass_B &&
        pass_C &&
        pass_D &&
        pass_lof_no_brake &&
        pass_wt_has_brake;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "E08_E09 primary_conditions" << std::endl;
    std::cout << "Condition\tSMAD4\tTGFb\tProliferation\tZEB1" << std::endl;
    std::cout << rA.label << "\t" << rA.smad4 << "\t" << rA.tgfb << "\t" << rA.proliferation << "\t" << rA.zeb1 << std::endl;
    std::cout << rB.label << "\t" << rB.smad4 << "\t" << rB.tgfb << "\t" << rB.proliferation << "\t" << rB.zeb1 << std::endl;
    std::cout << rC.label << "\t" << rC.smad4 << "\t" << rC.tgfb << "\t" << rC.proliferation << "\t" << rC.zeb1 << std::endl;
    std::cout << rD.label << "\t" << rD.smad4 << "\t" << rD.tgfb << "\t" << rD.proliferation << "\t" << rD.zeb1 << std::endl;

    std::cout << "E08_E09 negative_sweep"
              << " lof_tgfb_0.3=" << lof_prolif[0]
              << " lof_tgfb_0.6=" << lof_prolif[1]
              << " lof_tgfb_0.9=" << lof_prolif[2]
              << " wt_tgfb_0.3=" << wt_prolif[0]
              << " wt_tgfb_0.6=" << wt_prolif[1]
              << " wt_tgfb_0.9=" << wt_prolif[2]
              << std::endl;

    std::cout << "E08_E09 checks"
              << " A=" << (pass_A ? 1 : 0)
              << " B=" << (pass_B ? 1 : 0)
              << " C=" << (pass_C ? 1 : 0)
              << " D=" << (pass_D ? 1 : 0)
              << " LOF_no_brake=" << (pass_lof_no_brake ? 1 : 0)
              << " WT_has_brake=" << (pass_wt_has_brake ? 1 : 0)
              << std::endl;

    if (!pass)
    {
        std::cout << "FAIL E08/E09 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E08/E09 integration" << std::endl;
    return 0;
}
