#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "feedback_loop_common.h"

using namespace feedback_loop_common;

namespace
{

struct FieldStats
{
    double min_value = std::numeric_limits<double>::infinity();
    double max_value = -std::numeric_limits<double>::infinity();
    int nan_count = 0;
    int negative_count = 0;
};

void update_stats(FieldStats& stats, double value)
{
    if (std::isnan(value))
    {
        ++stats.nan_count;
        return;
    }
    stats.min_value = std::min(stats.min_value, value);
    stats.max_value = std::max(stats.max_value, value);
    if (value < -1e-12) ++stats.negative_count;
}

FieldStats scan_substrate_stats(int substrate_index)
{
    FieldStats stats;
    const int n_voxels = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_voxels; ++n)
    {
        const std::vector<double>& rho = microenvironment.density_vector(n);
        if (substrate_index < 0 || substrate_index >= static_cast<int>(rho.size())) continue;
        update_stats(stats, rho[substrate_index]);
    }
    return stats;
}

int count_live_cells_of_type(int type_id)
{
    int count = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        if (pCell->phenotype.death.dead) continue;
        if (pCell->type == type_id) ++count;
    }
    return count;
}

int count_tumor_zeb1_on(int tumor_type_id)
{
    int count = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != tumor_type_id) continue;
        const int zeb1_idx = custom_index(pCell, "zeb1_active");
        if (zeb1_idx >= 0 && pCell->custom_data[zeb1_idx] == 1.0) ++count;
    }
    return count;
}

int count_cells_in_annulus(double cx, double cy, double r_inner, double r_outer)
{
    int count = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        const double dx = pCell->position[0] - cx;
        const double dy = pCell->position[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r >= r_inner && r < r_outer) ++count;
    }
    return count;
}

void print_snapshot_row(double cx, double cy, double r0, double r1)
{
    const double o2 = mean_field_in_annulus(oxygen_index, cx, cy, r0, r1);
    const double tgfb = mean_field_in_annulus(tgfb_index, cx, cy, r0, r1);
    const double shh = mean_field_in_annulus(shh_index, cx, cy, r0, r1);
    const double ecm = mean_field_in_annulus(ecm_index, cx, cy, r0, r1);
    const int cells = count_cells_in_annulus(cx, cy, r0, r1);
    std::cout << "SNAPSHOT annulus[" << r0 << "," << r1 << ")"
              << " cells=" << cells
              << " o2=" << o2
              << " tgfb=" << tgfb
              << " shh=" << shh
              << " ecm=" << ecm << std::endl;
}

} // namespace

int main()
{
    // PhysiCell core currently emits very verbose per-cell debug output on stdout.
    // Mute stdout during setup/simulation and restore it for the smoke summary.
    std::ofstream null_out("/dev/null");
    std::streambuf* cout_backup = std::cout.rdbuf(null_out.rdbuf());

    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    if (pTumor == NULL || pStroma == NULL)
    {
        std::cerr << "FAIL: missing tumor_cell or stromal_cell definition." << std::endl;
        return 1;
    }

    const int tumor_type_id = pTumor->type;
    const int stromal_type_id = pStroma->type;

    // Smoke-test setup requested by user.
    const double cx = 0.0;
    const double cy = 0.0;
    const int initial_tumor_cells = 20;
    const int initial_psc_cells = 50;
    const double dt = 1.0;
    const int n_steps = 500;

    // Stabilize long smoke run while preserving growth/activation dynamics.
    parameters.doubles("base_proliferation_rate") = 0.003;
    parameters.doubles("caf_proliferation_rate") = 0.001;
    parameters.doubles("psc_proliferation_rate") = 0.0;

    // No drug.
    parameters.doubles("drug_uptake_rate") = 0.0;
    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.6);

    place_tumor_cluster(pTumor, initial_tumor_cells, cx, cy, 45.0);
    place_stroma_ring(pStroma, initial_psc_cells, cx, cy, 95.0, 170.0);

    const int tumor_start = count_live_cells_of_type(tumor_type_id);
    const int stromal_start = count_live_cells_of_type(stromal_type_id);
    const int caf_start = count_activated_cafs(stromal_type_id);
    const double ecm_near_tumor_start = mean_ecm_in_annulus(cx, cy, 0.0, 200.0);

    bool any_nan = false;
    bool any_negative = false;

    double t = 0.0;
    for (int step = 1; step <= n_steps; ++step)
    {
        advance_steps(1, dt, t);

        const FieldStats o2_stats = scan_substrate_stats(oxygen_index);
        const FieldStats tgfb_stats = scan_substrate_stats(tgfb_index);
        const FieldStats shh_stats = scan_substrate_stats(shh_index);
        const FieldStats ecm_stats = scan_substrate_stats(ecm_index);
        const FieldStats drug_stats = scan_substrate_stats(drug_index);

        any_nan = any_nan ||
                  o2_stats.nan_count > 0 ||
                  tgfb_stats.nan_count > 0 ||
                  shh_stats.nan_count > 0 ||
                  ecm_stats.nan_count > 0 ||
                  drug_stats.nan_count > 0;
        any_negative = any_negative ||
                       o2_stats.negative_count > 0 ||
                       tgfb_stats.negative_count > 0 ||
                       shh_stats.negative_count > 0 ||
                       ecm_stats.negative_count > 0 ||
                       drug_stats.negative_count > 0;
    }

    std::cout.rdbuf(cout_backup);
    null_out.close();

    const int tumor_end = count_live_cells_of_type(tumor_type_id);
    const int stromal_end = count_live_cells_of_type(stromal_type_id);
    const int caf_end = count_activated_cafs(stromal_type_id);
    const int tumor_zeb1_on = count_tumor_zeb1_on(tumor_type_id);

    const double ecm_near_tumor_end = mean_ecm_in_annulus(cx, cy, 0.0, 200.0);
    const std::vector<double> center{cx, cy, 0.0};
    const double o2_center = substrate_at_position(oxygen_index, center);
    const double o2_boundary = mean_field_in_annulus(oxygen_index, cx, cy, 820.0, 980.0);
    const double tgfb_center = mean_field_in_annulus(tgfb_index, cx, cy, 0.0, 200.0);
    const double tgfb_boundary = mean_field_in_annulus(tgfb_index, cx, cy, 820.0, 980.0);
    const double shh_center = mean_field_in_annulus(shh_index, cx, cy, 0.0, 200.0);
    const double shh_boundary = mean_field_in_annulus(shh_index, cx, cy, 820.0, 980.0);

    const FieldStats o2_final = scan_substrate_stats(oxygen_index);
    const FieldStats tgfb_final = scan_substrate_stats(tgfb_index);
    const FieldStats shh_final = scan_substrate_stats(shh_index);
    const FieldStats ecm_final = scan_substrate_stats(ecm_index);
    const FieldStats drug_final = scan_substrate_stats(drug_index);

    const bool pass_no_nan_or_negative = (!any_nan && !any_negative);
    const bool pass_tumor_growth = (tumor_end > tumor_start);
    const bool pass_caf_activation = (caf_end > 0);
    const bool pass_ecm_growth = (ecm_near_tumor_end > ecm_near_tumor_start);
    const bool pass_o2_gradient = (o2_center < o2_boundary);
    const bool pass_emt = (tumor_zeb1_on > 0);
    const bool pass_tgfb_gradient = (tgfb_center > tgfb_boundary);
    const bool pass_shh_gradient = (shh_center > shh_boundary);

    std::cout << "SMOKE step_count=" << n_steps
              << " dt=" << dt
              << " elapsed_time_min=" << t << std::endl;
    std::cout << "COUNTS tumor_start=" << tumor_start
              << " tumor_end=" << tumor_end
              << " stromal_start=" << stromal_start
              << " stromal_end=" << stromal_end
              << " caf_start=" << caf_start
              << " caf_end=" << caf_end
              << " tumor_zeb1_on=" << tumor_zeb1_on << std::endl;
    std::cout << "RANGES oxygen=[" << o2_final.min_value << "," << o2_final.max_value << "]"
              << " tgfb=[" << tgfb_final.min_value << "," << tgfb_final.max_value << "]"
              << " shh=[" << shh_final.min_value << "," << shh_final.max_value << "]"
              << " ecm=[" << ecm_final.min_value << "," << ecm_final.max_value << "]"
              << " drug=[" << drug_final.min_value << "," << drug_final.max_value << "]"
              << std::endl;
    std::cout << "SPATIAL center_vs_boundary"
              << " o2_center=" << o2_center << " o2_boundary=" << o2_boundary
              << " tgfb_center=" << tgfb_center << " tgfb_boundary=" << tgfb_boundary
              << " shh_center=" << shh_center << " shh_boundary=" << shh_boundary
              << " ecm_near_start=" << ecm_near_tumor_start
              << " ecm_near_end=" << ecm_near_tumor_end
              << std::endl;

    // One snapshot of spatial state (radial profile at final step).
    print_snapshot_row(cx, cy, 0.0, 120.0);
    print_snapshot_row(cx, cy, 120.0, 240.0);
    print_snapshot_row(cx, cy, 240.0, 400.0);
    print_snapshot_row(cx, cy, 400.0, 700.0);
    print_snapshot_row(cx, cy, 700.0, 980.0);

    std::cout << "CHECK no_nan_or_negative=" << (pass_no_nan_or_negative ? 1 : 0)
              << " tumor_growth=" << (pass_tumor_growth ? 1 : 0)
              << " caf_activation=" << (pass_caf_activation ? 1 : 0)
              << " ecm_growth=" << (pass_ecm_growth ? 1 : 0)
              << " o2_gradient=" << (pass_o2_gradient ? 1 : 0)
              << " emt_present=" << (pass_emt ? 1 : 0)
              << " tgfb_gradient=" << (pass_tgfb_gradient ? 1 : 0)
              << " shh_gradient=" << (pass_shh_gradient ? 1 : 0)
              << std::endl;

    const bool all_pass =
        pass_no_nan_or_negative &&
        pass_tumor_growth &&
        pass_caf_activation &&
        pass_ecm_growth &&
        pass_o2_gradient &&
        pass_emt &&
        pass_tgfb_gradient &&
        pass_shh_gradient;

    if (!all_pass)
    {
        std::cout << "FAIL smoke_integration_test" << std::endl;
        return 1;
    }

    std::cout << "PASS smoke_integration_test" << std::endl;
    return 0;
}
