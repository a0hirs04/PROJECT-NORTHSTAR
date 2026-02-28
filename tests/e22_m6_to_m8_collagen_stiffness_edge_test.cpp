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
    double boundary_pressure = 0.0;
    double boundary_ecm = 0.0;
    double tumor_radius = 0.0;
    int boundary_n = 0;
    int live_n = 0;
};

CohortMetrics measure_cohort(bool left_half, double cx, double cy)
{
    CohortMetrics m;
    m.boundary_ecm = mean_ecm_in_annulus(cx, cy, 30.0, 130.0);

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;

        if (left_half && pCell->position[0] >= 0.0) continue;
        if (!left_half && pCell->position[0] <= 0.0) continue;

        ++m.live_n;

        const double dx = pCell->position[0] - cx;
        const double dy = pCell->position[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);
        m.tumor_radius = std::max(m.tumor_radius, r);

        if (r >= 30.0 && r <= 130.0)
        {
            const int pressure_idx = custom_index(pCell, "mechanical_pressure");
            if (pressure_idx >= 0)
            {
                m.boundary_pressure += pCell->custom_data[pressure_idx];
                ++m.boundary_n;
            }
        }
    }

    if (m.boundary_n > 0)
    {
        m.boundary_pressure /= static_cast<double>(m.boundary_n);
    }
    return m;
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // Isolate E22 with same total ECM density, different composition.
    parameters.doubles("base_proliferation_rate") = 0.02;
    parameters.doubles("go_grow_penalty") = 0.0;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 0.2;
    parameters.doubles("apoptosis_resistance") = 1.0;

    parameters.doubles("mechanical_compaction_strength") = 1.5;
    parameters.doubles("compaction_ecm_increment") = 0.01;

    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("drug_uptake_rate") = 0.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("emt_induction_threshold") = 10.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.5);

    const double left_cx = -350.0;
    const double left_cy = 0.0;
    const double right_cx = 350.0;
    const double right_cy = 0.0;

    place_tumor_cluster(pTumor, 20, left_cx, left_cy, 80.0);   // A collagen-dominant
    place_tumor_cluster(pTumor, 20, right_cx, right_cy, 80.0);  // B HA-dominant

    // Same ECM density, different composition.
    set_annulus_ecm(left_cx, left_cy, 30.0, 130.0, 0.7, 0.2);   // collagen fraction 0.8
    set_annulus_ecm(right_cx, right_cy, 30.0, 130.0, 0.7, 0.8); // collagen fraction 0.2

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        pCell->phenotype.motility.is_motile = false;
    }

    double t = 0.0;
    const double dt = 1.0;

    advance_steps(10, dt, t);
    const CohortMetrics A10 = measure_cohort(true, left_cx, left_cy);
    const CohortMetrics B10 = measure_cohort(false, right_cx, right_cy);

    advance_steps(140, dt, t);
    const CohortMetrics A100 = measure_cohort(true, left_cx, left_cy);
    const CohortMetrics B100 = measure_cohort(false, right_cx, right_cy);

    const double A_compaction = A100.boundary_ecm - A10.boundary_ecm;
    const double B_compaction = B100.boundary_ecm - B10.boundary_ecm;

    const bool pass_pressure = (A100.boundary_pressure > B100.boundary_pressure);
    const bool pass_radius = (A100.tumor_radius < B100.tumor_radius);
    const bool pass_compaction = (A_compaction > B_compaction);

    std::cout << "E22 measurements"
              << " A_boundary_pressure=" << A100.boundary_pressure
              << " B_boundary_pressure=" << B100.boundary_pressure
              << " A_tumor_radius=" << A100.tumor_radius
              << " B_tumor_radius=" << B100.tumor_radius
              << " A_compaction=" << A_compaction
              << " B_compaction=" << B_compaction
              << " A_ecm_density_init=0.7"
              << " B_ecm_density_init=0.7"
              << " A_ha_fraction=0.2"
              << " B_ha_fraction=0.8"
              << std::endl;

    std::cout << "E22 checks"
              << " pressure_A_gt_B=" << (pass_pressure ? 1 : 0)
              << " radius_A_lt_B=" << (pass_radius ? 1 : 0)
              << " compaction_A_gt_B=" << (pass_compaction ? 1 : 0)
              << std::endl;

    if (!(pass_pressure && pass_radius && pass_compaction))
    {
        std::cout << "FAIL E22 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E22 integration" << std::endl;
    return 0;
}
