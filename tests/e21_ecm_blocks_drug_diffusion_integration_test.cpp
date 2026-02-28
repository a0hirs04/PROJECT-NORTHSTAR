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

void set_drug_source_left_boundary(double value)
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

double mean_drug_in_band(double x_min, double x_max, double y_abs_max = 60.0)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    double sum = 0.0;
    int count = 0;
    for (int n = 0; n < n_vox; ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        if (c[0] >= x_min && c[0] <= x_max && std::fabs(c[1]) <= y_abs_max)
        {
            sum += microenvironment.density_vector(n)[drug_index];
            ++count;
        }
    }
    return (count > 0) ? (sum / static_cast<double>(count)) : 0.0;
}

double run_condition(bool with_ecm_band)
{
    // Reset fields and tumor internal state.
    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.8);

    if (with_ecm_band)
    {
        // Dense barrier between source (left) and tumor (center).
        const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
        for (int n = 0; n < n_vox; ++n)
        {
            const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
            if (c[0] >= -600.0 && c[0] <= -200.0)
            {
                microenvironment.density_vector(n)[ecm_index] = 0.8;
                set_ecm_ha_fraction(n, 0.8);
            }
        }
    }

    Cell_Container* cc = (Cell_Container*)microenvironment.agent_container;
    double t = 0.0;
    const double dt = 1.0;
    for (int step = 0; step < 200; ++step)
    {
        set_drug_source_left_boundary(1.0);
        microenvironment.simulate_diffusion_decay(dt);
        cc->update_all_cells(t, dt, dt, dt);
        t += dt;
    }

    // Tumor location is fixed at center.
    const std::vector<double> center{0.0, 0.0, 0.0};
    const double drug_at_tumor = substrate_at_position(drug_index, center);
    return drug_at_tumor;
}

} // namespace

int main()
{
    std::ofstream null_out("/dev/null");
    std::streambuf* cout_backup = std::cout.rdbuf(null_out.rdbuf());
    initialize_world();
    std::cout.rdbuf(cout_backup);
    null_out.close();

    assert(drug_index >= 0);
    assert(ecm_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // Isolation for E21: only ECM differs between conditions.
    parameters.doubles("drug_uptake_rate") = 0.0;
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;

    // Single tumor cell at center.
    Cell* tumor = create_cell(*pTumor);
    tumor->assign_position(std::vector<double>{0.0, 0.0, 0.0});
    tumor->phenotype.motility.is_motile = false;
    const int intra_idx = tumor->custom_data.find_variable_index("intracellular_drug");
    const int texp_idx = tumor->custom_data.find_variable_index("time_since_drug_exposure");
    if (intra_idx >= 0) tumor->custom_data[intra_idx] = 0.0;
    if (texp_idx >= 0) tumor->custom_data[texp_idx] = -1.0;

    const double drug_A = run_condition(false);
    const double drug_B = run_condition(true);

    // Gradient in the ECM band (condition B) should be steep.
    const double band_left = mean_drug_in_band(-590.0, -560.0);
    const double band_right = mean_drug_in_band(-240.0, -210.0);

    const bool pass_less_drug_with_ecm = (drug_B < drug_A);
    const bool pass_30pct_reduction = (drug_B <= 0.7 * drug_A);
    const bool pass_steep_gradient = (band_left > band_right);

    std::cout << "E21 measurements:"
              << " drug_at_tumor_no_ecm=" << drug_A
              << " drug_at_tumor_with_ecm=" << drug_B
              << " relative_ratio_B_over_A=" << (drug_A > 0.0 ? (drug_B / drug_A) : 0.0)
              << " band_left_drug=" << band_left
              << " band_right_drug=" << band_right
              << std::endl;
    std::cout << "E21 checks:"
              << " less_drug_with_ecm=" << (pass_less_drug_with_ecm ? 1 : 0)
              << " reduction_ge_30pct=" << (pass_30pct_reduction ? 1 : 0)
              << " steep_band_gradient=" << (pass_steep_gradient ? 1 : 0)
              << std::endl;

    if (!(pass_less_drug_with_ecm && pass_30pct_reduction && pass_steep_gradient))
    {
        std::cout << "FAIL E21 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E21 integration" << std::endl;
    return 0;
}
