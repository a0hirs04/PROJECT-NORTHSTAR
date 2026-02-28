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

double mean_ecm_at_voxels(const std::vector<int>& voxels)
{
    double sum = 0.0;
    for (size_t i = 0; i < voxels.size(); ++i)
    {
        sum += microenvironment.density_vector(voxels[i])[ecm_index];
    }
    return voxels.empty() ? 0.0 : sum / static_cast<double>(voxels.size());
}

int voxel_for_pos(double x, double y)
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

    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pStroma != NULL);

    // Isolate E10 with manual CAF/PSC states.
    parameters.doubles("ecm_production_rate_base") = 0.001;
    parameters.doubles("ecm_production_rate_boosted") = 0.001;
    parameters.doubles("ecm_ha_fraction_default") = 0.6;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("caf_proliferation_rate") = 0.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;

    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.6);

    std::vector<int> caf_voxels;
    std::vector<int> psc_voxels;
    caf_voxels.reserve(5);
    psc_voxels.reserve(5);

    const double two_pi = 6.283185307179586;
    for (int i = 0; i < 5; ++i)
    {
        const double th = two_pi * (static_cast<double>(i) / 5.0);
        const double x = 120.0 * std::cos(th);
        const double y = 120.0 * std::sin(th);
        Cell* caf = create_cell(*pStroma);
        caf->assign_position(std::vector<double>{x, y, 0.0});
        caf->phenotype.motility.is_motile = false;
        const int acta2_idx = caf->custom_data.find_variable_index("acta2_active");
        const int gli1_idx = caf->custom_data.find_variable_index("gli1_active");
        assert(acta2_idx >= 0);
        assert(gli1_idx >= 0);
        caf->custom_data[acta2_idx] = 1.0;
        caf->custom_data[gli1_idx] = 0.0;
        caf_voxels.push_back(voxel_for_pos(x, y));
    }

    for (int i = 0; i < 5; ++i)
    {
        const double th = two_pi * (static_cast<double>(i) / 5.0);
        const double x = 420.0 * std::cos(th);
        const double y = 420.0 * std::sin(th);
        Cell* psc = create_cell(*pStroma);
        psc->assign_position(std::vector<double>{x, y, 0.0});
        psc->phenotype.motility.is_motile = false;
        const int acta2_idx = psc->custom_data.find_variable_index("acta2_active");
        const int gli1_idx = psc->custom_data.find_variable_index("gli1_active");
        assert(acta2_idx >= 0);
        assert(gli1_idx >= 0);
        psc->custom_data[acta2_idx] = 0.0;
        psc->custom_data[gli1_idx] = 0.0;
        psc_voxels.push_back(voxel_for_pos(x, y));
    }

    // Ensure zero start.
    for (size_t i = 0; i < caf_voxels.size(); ++i) microenvironment.density_vector(caf_voxels[i])[ecm_index] = 0.0;
    for (size_t i = 0; i < psc_voxels.size(); ++i) microenvironment.density_vector(psc_voxels[i])[ecm_index] = 0.0;

    const int empty_voxel = voxel_for_pos(850.0, 850.0);
    microenvironment.density_vector(empty_voxel)[ecm_index] = 0.0;

    const double dt = 6.0;
    double t = 0.0;
    double caf25 = 0.0;
    double caf50 = 0.0;
    double psc50 = 0.0;
    double empty50 = 0.0;

    for (int step = 1; step <= 50; ++step)
    {
        advance_steps(1, dt, t);
        if (step == 25) caf25 = mean_ecm_at_voxels(caf_voxels);
        if (step == 50)
        {
            caf50 = mean_ecm_at_voxels(caf_voxels);
            psc50 = mean_ecm_at_voxels(psc_voxels);
            empty50 = microenvironment.density_vector(empty_voxel)[ecm_index];
        }
    }

    const bool pass_caf_increase = (caf50 > 0.0);
    const bool pass_psc_zero = std::fabs(psc50) < 1e-12;
    const bool pass_empty_zero = std::fabs(empty50) < 1e-12;
    const bool pass_ongoing = (caf50 > caf25);

    std::cout << "E10 measurements:"
              << " ecm_caf_step25=" << caf25
              << " ecm_caf_step50=" << caf50
              << " ecm_psc_step50=" << psc50
              << " ecm_empty_step50=" << empty50
              << std::endl;
    std::cout << "E10 checks:"
              << " caf_increase=" << (pass_caf_increase ? 1 : 0)
              << " psc_zero=" << (pass_psc_zero ? 1 : 0)
              << " empty_zero=" << (pass_empty_zero ? 1 : 0)
              << " ongoing=" << (pass_ongoing ? 1 : 0)
              << std::endl;

    if (!(pass_caf_increase && pass_psc_zero && pass_empty_zero && pass_ongoing))
    {
        std::cout << "FAIL E10 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E10 integration" << std::endl;
    return 0;
}
