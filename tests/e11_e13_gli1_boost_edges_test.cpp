#include <cassert>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

bool nearly_equal(double a, double b, double eps = 1e-6)
{
    return std::fabs(a - b) <= eps;
}

int count_live_ids(const std::unordered_set<int>& ids)
{
    int n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (ids.find(pCell->ID) != ids.end()) ++n;
    }
    return n;
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pStroma != NULL);
    assert(ecm_index >= 0);

    // Shared isolation.
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;
    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    // E11: ECM production gli1 ON > OFF in activated CAFs.
    parameters.doubles("ecm_production_rate_base") = 0.01;
    parameters.doubles("ecm_production_rate_boosted") = 0.015;
    parameters.doubles("ecm_ha_fraction_default") = 0.6;

    Cell* caf_off = create_cell(*pStroma);
    caf_off->assign_position(std::vector<double>{-150.0, 0.0, 0.0});
    caf_off->phenotype.motility.is_motile = false;
    Cell* caf_on = create_cell(*pStroma);
    caf_on->assign_position(std::vector<double>{150.0, 0.0, 0.0});
    caf_on->phenotype.motility.is_motile = false;

    caf_off->custom_data[custom_index(caf_off, "acta2_active")] = 1.0;
    caf_off->custom_data[custom_index(caf_off, "gli1_active")] = 0.0;
    caf_on->custom_data[custom_index(caf_on, "acta2_active")] = 1.0;
    caf_on->custom_data[custom_index(caf_on, "gli1_active")] = 1.0;

    const int voxel_off = voxel_index_for_cell(caf_off);
    const int voxel_on = voxel_index_for_cell(caf_on);
    microenvironment.density_vector(voxel_off)[ecm_index] = 0.0;
    microenvironment.density_vector(voxel_on)[ecm_index] = 0.0;

    double t = 0.0;
    const double dt = 1.0;
    advance_steps(50, dt, t);

    const double ecm_off = microenvironment.density_vector(voxel_off)[ecm_index];
    const double ecm_on = microenvironment.density_vector(voxel_on)[ecm_index];
    const bool pass_e11 = (ecm_on > ecm_off);

    // E13: CAF proliferation boost via GLI1.
    parameters.doubles("caf_proliferation_rate") = 0.03;
    parameters.doubles("gli1_proliferation_boost") = 2.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;

    std::unordered_set<int> off_ids;
    std::unordered_set<int> on_ids;

    const double two_pi = 6.283185307179586;
    for (int i = 0; i < 10; ++i)
    {
        const double th = two_pi * (static_cast<double>(i) / 10.0);
        Cell* c0 = create_cell(*pStroma);
        c0->assign_position(std::vector<double>{-500.0 + 20.0 * std::cos(th), 120.0 + 20.0 * std::sin(th), 0.0});
        c0->phenotype.motility.is_motile = false;
        c0->custom_data[custom_index(c0, "acta2_active")] = 1.0;
        c0->custom_data[custom_index(c0, "gli1_active")] = 0.0;
        off_ids.insert(c0->ID);

        Cell* c1 = create_cell(*pStroma);
        c1->assign_position(std::vector<double>{500.0 + 20.0 * std::cos(th), 120.0 + 20.0 * std::sin(th), 0.0});
        c1->phenotype.motility.is_motile = false;
        c1->custom_data[custom_index(c1, "acta2_active")] = 1.0;
        c1->custom_data[custom_index(c1, "gli1_active")] = 1.0;
        on_ids.insert(c1->ID);
    }

    const int init_off = count_live_ids(off_ids);
    const int init_on = count_live_ids(on_ids);

    advance_steps(1, dt, t);

    double rate_off = 0.0;
    double rate_on = 0.0;
    bool found_off = false;
    bool found_on = false;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (!found_off && off_ids.find(pCell->ID) != off_ids.end())
        {
            rate_off = pCell->phenotype.cycle.data.transition_rate(0, 0);
            found_off = true;
        }
        if (!found_on && on_ids.find(pCell->ID) != on_ids.end())
        {
            rate_on = pCell->phenotype.cycle.data.transition_rate(0, 0);
            found_on = true;
        }
        if (found_off && found_on) break;
    }

    advance_steps(99, dt, t);

    const int final_off = count_live_ids(off_ids);
    const int final_on = count_live_ids(on_ids);

    const double expected_boost = parameters.doubles("gli1_proliferation_boost");
    const double rate_ratio = (rate_off > 0.0) ? (rate_on / rate_off) : 0.0;

    const bool pass_e13_ratio = nearly_equal(rate_ratio, expected_boost);
    const bool pass_e13_count = (final_on > final_off) && (final_on > init_on) && (final_off >= init_off);

    const bool pass = pass_e11 && pass_e13_ratio && pass_e13_count;

    std::cout << "P5B measurements"
              << " ecm_gli1_off=" << ecm_off
              << " ecm_gli1_on=" << ecm_on
              << " rate_off=" << rate_off
              << " rate_on=" << rate_on
              << " rate_ratio=" << rate_ratio
              << " expected_boost=" << expected_boost
              << " init_off_n=" << init_off
              << " init_on_n=" << init_on
              << " final_off_n=" << final_off
              << " final_on_n=" << final_on
              << std::endl;

    std::cout << "P5B checks"
              << " E11_ecm_on_gt_off=" << (pass_e11 ? 1 : 0)
              << " E13_rate_ratio_matches_boost=" << (pass_e13_ratio ? 1 : 0)
              << " E13_count_on_gt_off=" << (pass_e13_count ? 1 : 0)
              << std::endl;

    if (!pass)
    {
        std::cout << "FAIL P5 Cluster B" << std::endl;
        return 1;
    }

    std::cout << "PASS P5 Cluster B" << std::endl;
    return 0;
}
