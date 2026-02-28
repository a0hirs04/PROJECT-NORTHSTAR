#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

double parse_base_prolif(int argc, char** argv)
{
    double val = 0.01;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg.rfind("--base-prolif=", 0) == 0)
        {
            val = std::stod(arg.substr(std::string("--base-prolif=").size()));
        }
    }
    return val;
}

struct Metrics
{
    double boundary_pressure = 0.0;
    double boundary_ecm = 0.0;
    int boundary_n = 0;
};

Metrics measure(double cx, double cy)
{
    Metrics m;
    m.boundary_ecm = mean_ecm_in_annulus(cx, cy, 30.0, 120.0);

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;

        const double dx = pCell->position[0] - cx;
        const double dy = pCell->position[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r < 30.0 || r > 120.0) continue;

        const int pressure_idx = custom_index(pCell, "mechanical_pressure");
        if (pressure_idx >= 0)
        {
            m.boundary_pressure += pCell->custom_data[pressure_idx];
            ++m.boundary_n;
        }
    }

    if (m.boundary_n > 0)
    {
        m.boundary_pressure /= static_cast<double>(m.boundary_n);
    }
    return m;
}

} // namespace

int main(int argc, char** argv)
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    const double base_prolif = parse_base_prolif(argc, argv);

    // Isolate E14: stress generation from proliferation in confinement.
    parameters.doubles("base_proliferation_rate") = base_prolif;
    parameters.doubles("go_grow_penalty") = 0.0;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 1e9;
    parameters.doubles("apoptosis_resistance") = 1.0;

    parameters.doubles("mechanical_compaction_strength") = 1.0;
    parameters.doubles("compaction_ecm_increment") = 0.01;

    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("drug_uptake_rate") = 0.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("emt_induction_threshold") = 10.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.2);

    const double cx = 0.0;
    const double cy = 0.0;

    place_tumor_cluster(pTumor, 20, cx, cy, 80.0);
    set_annulus_ecm(cx, cy, 30.0, 120.0, 0.8, 0.2);

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        pCell->phenotype.motility.is_motile = false;
    }

    double t = 0.0;
    const double dt = 1.0;

    advance_steps(10, dt, t);
    const Metrics m10 = measure(cx, cy);

    advance_steps(90, dt, t);
    const Metrics m100 = measure(cx, cy);

    std::cout << "E14 scenario"
              << " base_prolif=" << base_prolif
              << " pressure_step10=" << m10.boundary_pressure
              << " pressure_step100=" << m100.boundary_pressure
              << " ecm_step10=" << m10.boundary_ecm
              << " ecm_step100=" << m100.boundary_ecm
              << " boundary_n_step100=" << m100.boundary_n
              << std::endl;

    // Per-scenario sanity checks; cross-condition check is done in runner.
    const bool local_pass =
        (m100.boundary_n > 0) &&
        (m10.boundary_ecm >= 0.0) &&
        (m100.boundary_ecm >= 0.0);

    std::cout << "E14 scenario_check=" << (local_pass ? 1 : 0) << std::endl;
    if (!local_pass)
    {
        std::cout << "FAIL E14 scenario" << std::endl;
        return 1;
    }

    std::cout << "PASS E14 scenario" << std::endl;
    return 0;
}
