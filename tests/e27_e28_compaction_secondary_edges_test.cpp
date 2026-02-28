#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

bool parse_compaction_on(int argc, char** argv)
{
    bool enabled = true;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--compaction=on") enabled = true;
        if (arg == "--compaction=off") enabled = false;
    }
    return enabled;
}

struct Metrics
{
    double boundary_ecm = 0.0;
    double centroid_o2 = 0.0;
    int live_tumor_n = 0;
};

Metrics measure(double cx, double cy, int tumor_type)
{
    Metrics m;
    m.boundary_ecm = mean_ecm_in_annulus(cx, cy, 30.0, 120.0);

    const std::vector<double> center{cx, cy, 0.0};
    m.centroid_o2 = substrate_at_position(oxygen_index, center);

    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type == tumor_type) ++m.live_tumor_n;
    }

    return m;
}

} // namespace

int main(int argc, char** argv)
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);
    assert(oxygen_index >= 0);
    assert(ecm_index >= 0);

    const bool compaction_on = parse_compaction_on(argc, argv);

    // Isolate E27/E28: no CAFs, no ECM production, only passive compaction.
    parameters.doubles("base_proliferation_rate") = 0.02;
    parameters.doubles("go_grow_penalty") = 0.0;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 1e9;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;
    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("emt_induction_threshold") = 10.0;
    parameters.doubles("mechanical_compaction_strength") = compaction_on ? 2.0 : 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.03;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.2);

    // Keep ECM changes attributable only to compaction, not PDE turnover.
    microenvironment.diffusion_coefficients[ecm_index] = 0.0;
    microenvironment.decay_rates[ecm_index] = 0.0;

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
    const Metrics m10 = measure(cx, cy, pTumor->type);

    advance_steps(90, dt, t);
    const Metrics m100 = measure(cx, cy, pTumor->type);

    const double ecm_delta = m100.boundary_ecm - m10.boundary_ecm;

    const bool local_pass =
        (m10.live_tumor_n > 0) &&
        (m100.live_tumor_n > 0) &&
        (m10.boundary_ecm >= 0.0) &&
        (m100.boundary_ecm >= 0.0) &&
        (m100.boundary_ecm <= 1.0);

    std::cout << "E27E28 scenario"
              << " compaction_on=" << (compaction_on ? 1 : 0)
              << " boundary_ecm_step10=" << m10.boundary_ecm
              << " boundary_ecm_step100=" << m100.boundary_ecm
              << " boundary_ecm_delta=" << ecm_delta
              << " centroid_o2_step10=" << m10.centroid_o2
              << " centroid_o2_step100=" << m100.centroid_o2
              << " live_tumor_step100=" << m100.live_tumor_n
              << std::endl;

    std::cout << "E27E28 scenario_check=" << (local_pass ? 1 : 0) << std::endl;
    if (!local_pass)
    {
        std::cout << "FAIL E27/E28 scenario" << std::endl;
        return 1;
    }

    std::cout << "PASS E27/E28 scenario" << std::endl;
    return 0;
}
