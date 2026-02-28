#include <cassert>
#include <cmath>
#include <iostream>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

struct CohortMetrics
{
    int live_count = 0;
    double boundary_pressure_mean = 0.0;
    double boundary_prolif_when_over_threshold_mean = 0.0;
    double interior_prolif_mean = 0.0;
    int boundary_n = 0;
    int boundary_over_threshold_n = 0;
    int interior_n = 0;
};

CohortMetrics compute_metrics(bool left_half, double cx, double cy, double threshold)
{
    CohortMetrics m;

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;

        if (left_half && pCell->position[0] >= 0.0) continue;
        if (!left_half && pCell->position[0] <= 0.0) continue;

        ++m.live_count;

        const int pressure_idx = custom_index(pCell, "mechanical_pressure");
        const double pressure = (pressure_idx >= 0) ? pCell->custom_data[pressure_idx] : 0.0;
        const double prolif = pCell->phenotype.cycle.data.transition_rate(0, 0);

        const double dx = pCell->position[0] - cx;
        const double dy = pCell->position[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);

        if (r >= 45.0 && r <= 105.0)
        {
            m.boundary_pressure_mean += pressure;
            ++m.boundary_n;

            if (pressure > threshold)
            {
                m.boundary_prolif_when_over_threshold_mean += prolif;
                ++m.boundary_over_threshold_n;
            }
        }

        if (r < 30.0)
        {
            m.interior_prolif_mean += prolif;
            ++m.interior_n;
        }
    }

    if (m.boundary_n > 0)
    {
        m.boundary_pressure_mean /= static_cast<double>(m.boundary_n);
    }
    if (m.boundary_over_threshold_n > 0)
    {
        m.boundary_prolif_when_over_threshold_mean /= static_cast<double>(m.boundary_over_threshold_n);
    }
    if (m.interior_n > 0)
    {
        m.interior_prolif_mean /= static_cast<double>(m.interior_n);
    }

    return m;
}

} // namespace

int main()
{
    initialize_world();

    assert(ecm_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // Isolate E26 (M8 -> M4) with matched control cohort.
    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("go_grow_penalty") = 0.0;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 0.5;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;
    parameters.doubles("tgfb_brake_sensitivity") = 0.0;

    parameters.doubles("mechanical_compaction_strength") = 1.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;

    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.2);

    const double left_cx = -320.0;
    const double left_cy = 0.0;
    const double right_cx = 320.0;
    const double right_cy = 0.0;

    place_tumor_cluster(pTumor, 30, left_cx, left_cy, 55.0);
    place_tumor_cluster(pTumor, 30, right_cx, right_cy, 55.0);

    // ECM-ring cohort (left): dense, collagen-dominant ring.
    set_annulus_ecm(left_cx, left_cy, 45.0, 110.0, 0.9, 0.2);

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        pCell->phenotype.motility.is_motile = false;
    }

    double t = 0.0;
    const double dt = 6.0;
    advance_steps(100, dt, t);

    const double threshold = parameters.doubles("contact_inhibition_threshold");
    const CohortMetrics ecm_cohort = compute_metrics(true, left_cx, left_cy, threshold);
    const CohortMetrics control_cohort = compute_metrics(false, right_cx, right_cy, threshold);

    const bool pass_pressure_over_threshold =
        (ecm_cohort.boundary_over_threshold_n > 0) &&
        (ecm_cohort.boundary_pressure_mean > threshold);

    // One-step pressure->proliferation delay means we check near-zero arrest.
    const bool pass_boundary_arrest =
        (ecm_cohort.boundary_over_threshold_n > 0) &&
        (ecm_cohort.boundary_prolif_when_over_threshold_mean <= 2e-3);

    const bool pass_interior_divides =
        (ecm_cohort.interior_n > 0) &&
        (ecm_cohort.interior_prolif_mean > 1e-4);

    const bool pass_population_slowed =
        (ecm_cohort.live_count < control_cohort.live_count);

    const bool pass =
        pass_pressure_over_threshold &&
        pass_boundary_arrest &&
        pass_interior_divides &&
        pass_population_slowed;

    std::cout << "E26 measurements"
              << " ecm_live_count=" << ecm_cohort.live_count
              << " control_live_count=" << control_cohort.live_count
              << " boundary_pressure_mean=" << ecm_cohort.boundary_pressure_mean
              << " contact_threshold=" << threshold
              << " boundary_prolif_over_threshold_mean=" << ecm_cohort.boundary_prolif_when_over_threshold_mean
              << " interior_prolif_mean=" << ecm_cohort.interior_prolif_mean
              << " boundary_over_threshold_n=" << ecm_cohort.boundary_over_threshold_n
              << " interior_n=" << ecm_cohort.interior_n
              << std::endl;

    std::cout << "E26 checks"
              << " pressure_over_threshold=" << (pass_pressure_over_threshold ? 1 : 0)
              << " boundary_arrest=" << (pass_boundary_arrest ? 1 : 0)
              << " interior_divides=" << (pass_interior_divides ? 1 : 0)
              << " slowed_vs_control=" << (pass_population_slowed ? 1 : 0)
              << std::endl;

    if (!pass)
    {
        std::cout << "FAIL E26 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E26 integration" << std::endl;
    return 0;
}
