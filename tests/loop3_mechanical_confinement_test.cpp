#include <cassert>
#include <cmath>
#include <iostream>
#include "feedback_loop_common.h"

using namespace feedback_loop_common;

namespace
{

struct CohortMetrics
{
    int live_count = 0;
    double boundary_pressure_mean = 0.0;
    double boundary_prolif_mean = 0.0;
    double interior_prolif_mean = 0.0;
    int boundary_n = 0;
    int interior_n = 0;
};

CohortMetrics compute_cohort_metrics(bool left_half, double cx, double cy, double threshold)
{
    CohortMetrics out;
    double boundary_high_prolif_sum = 0.0;
    int boundary_high_n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (left_half && pCell->position[0] >= 0.0) continue;
        if (!left_half && pCell->position[0] <= 0.0) continue;

        ++out.live_count;

        const int pressure_idx = custom_index(pCell, "mechanical_pressure");
        const double pressure = (pressure_idx >= 0) ? pCell->custom_data[pressure_idx] : 0.0;
        const double prolif = pCell->phenotype.cycle.data.transition_rate(0, 0);

        const double dx = pCell->position[0] - cx;
        const double dy = pCell->position[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);

        if (r >= 45.0 && r <= 100.0)
        {
            out.boundary_pressure_mean += pressure;
            out.boundary_prolif_mean += prolif;
            ++out.boundary_n;
            if (pressure > threshold)
            {
                boundary_high_prolif_sum += prolif;
                ++boundary_high_n;
            }
        }
        if (r < 25.0)
        {
            out.interior_prolif_mean += prolif;
            ++out.interior_n;
        }
    }

    if (out.boundary_n > 0)
    {
        out.boundary_pressure_mean /= static_cast<double>(out.boundary_n);
    }
    if (boundary_high_n > 0)
    {
        out.boundary_prolif_mean = boundary_high_prolif_sum / static_cast<double>(boundary_high_n);
    }
    if (out.interior_n > 0)
    {
        out.interior_prolif_mean /= static_cast<double>(out.interior_n);
    }
    return out;
}

} // namespace

int main()
{
    initialize_world();

    assert(ecm_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // Loop-3-specific parameterization.
    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("go_grow_penalty") = 0.0;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 0.2;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;

    parameters.doubles("mechanical_compaction_strength") = 1.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;

    parameters.doubles("emt_induction_threshold") = 10.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    // Two cohorts in one simulation: left has ECM ring, right has no ring.
    const double left_cx = -300.0;
    const double left_cy = 0.0;
    const double right_cx = 300.0;
    const double right_cy = 0.0;

    const int n_each = 50;
    place_tumor_cluster(pTumor, n_each, left_cx, left_cy, 60.0);
    place_tumor_cluster(pTumor, n_each, right_cx, right_cy, 60.0);

    // Collagen-dominant ring around left cohort.
    set_annulus_ecm(left_cx, left_cy, 45.0, 110.0, 0.8, 0.1);

    double t = 0.0;
    const double dt = 1.0;
    advance_steps(200, dt, t);

    const double threshold = parameters.doubles("contact_inhibition_threshold");
    const CohortMetrics ring = compute_cohort_metrics(true, left_cx, left_cy, threshold);
    const CohortMetrics control = compute_cohort_metrics(false, right_cx, right_cy, threshold);

    const bool pass =
        (ring.live_count < control.live_count) &&
        (ring.boundary_n > 0) &&
        (ring.interior_n > 0) &&
        (ring.boundary_pressure_mean > threshold) &&
        (ring.boundary_prolif_mean < 1e-4) &&
        (ring.interior_prolif_mean > 1e-4);

    std::cout << "LOOP3 ring_live_count=" << ring.live_count
              << " control_live_count=" << control.live_count << std::endl;
    std::cout << "LOOP3 boundary_pressure_mean=" << ring.boundary_pressure_mean
              << " contact_threshold=" << threshold << std::endl;
    std::cout << "LOOP3 boundary_prolif_mean=" << ring.boundary_prolif_mean
              << " interior_prolif_mean=" << ring.interior_prolif_mean
              << " boundary_n=" << ring.boundary_n
              << " interior_n=" << ring.interior_n << std::endl;

    std::cout << "LOOP3 " << (pass ? "PASS" : "FAIL") << std::endl;
    if (!pass) return 1;

    std::cout << "PASS loop3_mechanical_confinement_test" << std::endl;
    return 0;
}
