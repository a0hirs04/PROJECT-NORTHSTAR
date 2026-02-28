#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

int voxel_for_xy(double x, double y)
{
    std::vector<double> p(3, 0.0);
    p[0] = x;
    p[1] = y;
    p[2] = 0.0;
    return microenvironment.nearest_voxel_index(p);
}

} // namespace

int main()
{
    std::ofstream null_out("/dev/null");
    std::streambuf* cout_backup = std::cout.rdbuf(null_out.rdbuf());
    initialize_world();
    std::cout.rdbuf(cout_backup);
    null_out.close();

    assert(ecm_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // Isolate E19 (M5 -> M6 via mmp2_active flag): no production/compaction.
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0015;
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.0;
    parameters.doubles("hif1a_emt_boost") = 0.0;
    parameters.doubles("emt_induction_threshold") = 0.4;
    parameters.ints("emt_phenotype_extent") = 2;
    parameters.doubles("hypoxia_response_threshold") = 0.0;

    // Keep TGF-beta localized so one cohort is EMT/MMP2+ and the control is EMT/MMP2-.
    microenvironment.diffusion_coefficients[tgfb_index] = 0.0;
    microenvironment.decay_rates[tgfb_index] = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double x_on = -220.0;
    const double y_on = 0.0;
    const double x_off = 220.0;
    const double y_off = 0.0;

    const int voxel_on = voxel_for_xy(x_on, y_on);
    const int voxel_off = voxel_for_xy(x_off, y_off);
    const int voxel_on_adj = voxel_for_xy(x_on + 80.0, y_on);
    const int voxel_off_adj = voxel_for_xy(x_off + 80.0, y_off);

    microenvironment.density_vector(voxel_on)[ecm_index] = 0.5;
    microenvironment.density_vector(voxel_off)[ecm_index] = 0.5;
    microenvironment.density_vector(voxel_on_adj)[ecm_index] = 0.5;
    microenvironment.density_vector(voxel_off_adj)[ecm_index] = 0.5;
    microenvironment.density_vector(voxel_on)[tgfb_index] = 0.6;
    microenvironment.density_vector(voxel_off)[tgfb_index] = 0.0;

    Cell* mmp2_on = create_cell(*pTumor);
    mmp2_on->assign_position(std::vector<double>{x_on, y_on, 0.0});
    mmp2_on->phenotype.motility.is_motile = false;

    Cell* mmp2_off = create_cell(*pTumor);
    mmp2_off->assign_position(std::vector<double>{x_off, y_off, 0.0});
    mmp2_off->phenotype.motility.is_motile = false;

    const int mmp_on_idx = mmp2_on->custom_data.find_variable_index("mmp2_active");
    const int mmp_off_idx = mmp2_off->custom_data.find_variable_index("mmp2_active");
    assert(mmp_on_idx >= 0);
    assert(mmp_off_idx >= 0);

    mmp2_on->custom_data[mmp_on_idx] = 1.0;
    mmp2_off->custom_data[mmp_off_idx] = 0.0;

    double t = 0.0;
    const double dt = 6.0;

    double on_step10 = 0.0;
    double on_step25 = 0.0;
    double on_step50 = 0.0;

    double off_step10 = 0.0;
    double off_step25 = 0.0;
    double off_step50 = 0.0;

    for (int step = 1; step <= 50; ++step)
    {
        advance_steps(1, dt, t);

        if (step == 10)
        {
            on_step10 = microenvironment.density_vector(voxel_on)[ecm_index];
            off_step10 = microenvironment.density_vector(voxel_off)[ecm_index];
        }
        if (step == 25)
        {
            on_step25 = microenvironment.density_vector(voxel_on)[ecm_index];
            off_step25 = microenvironment.density_vector(voxel_off)[ecm_index];
        }
        if (step == 50)
        {
            on_step50 = microenvironment.density_vector(voxel_on)[ecm_index];
            off_step50 = microenvironment.density_vector(voxel_off)[ecm_index];
        }
    }

    const double on_adj_final = microenvironment.density_vector(voxel_on_adj)[ecm_index];
    const double off_adj_final = microenvironment.density_vector(voxel_off_adj)[ecm_index];

    const bool pass_on_decreases = (on_step50 < on_step25) && (on_step25 < on_step10);
    const bool pass_off_constant =
        (std::fabs(off_step10 - 0.5) < 1e-12) &&
        (std::fabs(off_step25 - 0.5) < 1e-12) &&
        (std::fabs(off_step50 - 0.5) < 1e-12);
    const bool pass_locality =
        (std::fabs(on_adj_final - 0.5) < 1e-12) &&
        (std::fabs(off_adj_final - 0.5) < 1e-12);
    const bool pass_floor = (on_step50 >= 0.0);

    std::cout << "E19 measurements"
              << " on_step10=" << on_step10
              << " on_step25=" << on_step25
              << " on_step50=" << on_step50
              << " off_step10=" << off_step10
              << " off_step25=" << off_step25
              << " off_step50=" << off_step50
              << " on_adj_final=" << on_adj_final
              << " off_adj_final=" << off_adj_final
              << std::endl;

    std::cout << "E19 checks"
              << " mmp2_on_decreases=" << (pass_on_decreases ? 1 : 0)
              << " mmp2_off_constant=" << (pass_off_constant ? 1 : 0)
              << " local_only=" << (pass_locality ? 1 : 0)
              << " nonnegative_floor=" << (pass_floor ? 1 : 0)
              << std::endl;

    if (!(pass_on_decreases && pass_off_constant && pass_locality && pass_floor))
    {
        std::cout << "FAIL E19 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E19 integration" << std::endl;
    return 0;
}
