#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

double live_tgfb_flux(int tumor_type)
{
    double sum = 0.0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != tumor_type) continue;
        sum += pCell->phenotype.secretion.secretion_rates[tgfb_index];
    }
    return sum;
}

int live_count(int tumor_type)
{
    int n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type == tumor_type) ++n;
    }
    return n;
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.0;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.0;

    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double two_pi = 6.283185307179586;
    std::vector<Cell*> tumors;
    for (int i = 0; i < 20; ++i)
    {
        const double th = two_pi * (static_cast<double>(i) / 20.0);
        Cell* p = create_cell(*pTumor);
        p->assign_position(std::vector<double>{60.0 * std::cos(th), 60.0 * std::sin(th), 0.0});
        p->phenotype.motility.is_motile = false;
        p->custom_data[custom_index(p, "hif1a_active")] = 0.0;
        tumors.push_back(p);
    }

    double t = 0.0;
    const double dt = 1.0;
    advance_steps(1, dt, t);

    const double flux_before = live_tgfb_flux(pTumor->type);
    const int live_before = live_count(pTumor->type);

    // Kill 10 cells.
    const int apop_model = tumors.front()->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    std::vector<Cell*> killed;
    for (int i = 0; i < 10; ++i)
    {
        tumors[i]->start_death(apop_model);
        killed.push_back(tumors[i]);
    }

    advance_steps(1, dt, t);

    const double flux_after = live_tgfb_flux(pTumor->type);
    const int live_after = live_count(pTumor->type);

    bool dead_zero = true;
    for (size_t i = 0; i < killed.size(); ++i)
    {
        Cell* pCell = killed[i];
        const double sec = pCell->phenotype.secretion.secretion_rates[tgfb_index];
        if (std::fabs(sec) > 1e-12) dead_zero = false;
    }

    const bool pass_e15 = (flux_after > 0.0);
    const bool pass_e16 = (flux_after < flux_before) && (live_after < live_before);
    const bool pass_dead_zero = dead_zero;

    const bool pass = pass_e15 && pass_e16 && pass_dead_zero;

    std::cout << "P5C measurements"
              << " live_before=" << live_before
              << " live_after=" << live_after
              << " flux_before=" << flux_before
              << " flux_after=" << flux_after
              << std::endl;

    std::cout << "P5C checks"
              << " E15_survivors_still_secrete=" << (pass_e15 ? 1 : 0)
              << " E16_total_flux_decreases=" << (pass_e16 ? 1 : 0)
              << " dead_cells_exact_zero_secretion=" << (pass_dead_zero ? 1 : 0)
              << std::endl;

    if (!pass)
    {
        std::cout << "FAIL P5 Cluster C" << std::endl;
        return 1;
    }

    std::cout << "PASS P5 Cluster C" << std::endl;
    return 0;
}
