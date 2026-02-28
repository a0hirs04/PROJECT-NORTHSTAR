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

void set_oxygen_boundary_source(double oxygen_value)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        if (std::fabs(c[0]) >= 900.0 || std::fabs(c[1]) >= 900.0)
        {
            microenvironment.density_vector(n)[oxygen_index] = oxygen_value;
        }
    }
}

std::vector<double> centroid_for_half(bool left_half)
{
    std::vector<double> c(3, 0.0);
    int n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (left_half && pCell->position[0] >= 0.0) continue;
        if (!left_half && pCell->position[0] <= 0.0) continue;

        c[0] += pCell->position[0];
        c[1] += pCell->position[1];
        c[2] += pCell->position[2];
        ++n;
    }
    if (n > 0)
    {
        c[0] /= static_cast<double>(n);
        c[1] /= static_cast<double>(n);
        c[2] /= static_cast<double>(n);
    }
    return c;
}

double hif_fraction_for_half(bool left_half)
{
    int total = 0;
    int active = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (left_half && pCell->position[0] >= 0.0) continue;
        if (!left_half && pCell->position[0] <= 0.0) continue;

        const int hif_idx = custom_index(pCell, "hif1a_active");
        if (hif_idx < 0) continue;

        ++total;
        if (pCell->custom_data[hif_idx] == 1.0) ++active;
    }

    if (total == 0) return 0.0;
    return static_cast<double>(active) / static_cast<double>(total);
}

} // namespace

int main()
{
    std::ofstream null_out("/dev/null");
    std::streambuf* cout_backup = std::cout.rdbuf(null_out.rdbuf());
    initialize_world();
    std::cout.rdbuf(cout_backup);
    null_out.close();

    assert(oxygen_index >= 0);
    assert(ecm_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // Isolate E20: only ECM barrier differs between cohorts.
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.0;
    parameters.doubles("drug_uptake_rate") = 0.0;

    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;

    parameters.doubles("hypoxia_response_threshold") = 0.6;

    reset_all_fields(0.06, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double left_cx = -400.0;
    const double left_cy = 0.0;
    const double right_cx = 400.0;
    const double right_cy = 0.0;

    place_tumor_cluster(pTumor, 45, left_cx, left_cy, 55.0);
    place_tumor_cluster(pTumor, 45, right_cx, right_cy, 55.0);

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        pCell->phenotype.motility.is_motile = false;
    }

    // Condition B (left): dense ECM ring around tumor; Condition A (right): no ECM.
    set_annulus_ecm(left_cx, left_cy, 60.0, 220.0, 0.9, 0.8);

    double t = 0.0;
    const double dt = 6.0;
    Cell_Container* cc = cell_container();

    for (int step = 0; step < 200; ++step)
    {
        set_oxygen_boundary_source(0.06);
        microenvironment.simulate_diffusion_decay(dt);
        cc->update_all_cells(t, dt, dt, dt);
        t += dt;
    }

    const std::vector<double> left_centroid = centroid_for_half(true);
    const std::vector<double> right_centroid = centroid_for_half(false);

    const double o2_with_ecm = substrate_at_position(oxygen_index, left_centroid);
    const double o2_no_ecm = substrate_at_position(oxygen_index, right_centroid);

    const double hif_with_ecm = hif_fraction_for_half(true);
    const double hif_no_ecm = hif_fraction_for_half(false);

    const bool pass_o2_drop = (o2_with_ecm < o2_no_ecm);
    const bool pass_hif_rise = (hif_with_ecm > hif_no_ecm);
    const bool pass_majority_hypoxic = (hif_with_ecm >= 0.5);

    std::cout << "E20 measurements"
              << " o2_no_ecm=" << o2_no_ecm
              << " o2_with_ecm=" << o2_with_ecm
              << " hif_fraction_no_ecm=" << hif_no_ecm
              << " hif_fraction_with_ecm=" << hif_with_ecm
              << std::endl;

    std::cout << "E20 checks"
              << " o2_lower_with_ecm=" << (pass_o2_drop ? 1 : 0)
              << " hif_higher_with_ecm=" << (pass_hif_rise ? 1 : 0)
              << " majority_hypoxic_with_ecm=" << (pass_majority_hypoxic ? 1 : 0)
              << std::endl;

    if (!(pass_o2_drop && pass_hif_rise && pass_majority_hypoxic))
    {
        std::cout << "FAIL E20 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E20 integration" << std::endl;
    return 0;
}
