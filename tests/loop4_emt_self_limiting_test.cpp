#include <cassert>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "feedback_loop_common.h"

using namespace feedback_loop_common;

namespace
{

bool has_decrease(const std::vector<double>& values, double eps = 1e-9)
{
    for (size_t i = 1; i < values.size(); ++i)
    {
        if (values[i] + eps < values[i - 1]) return true;
    }
    return false;
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // Loop-4-specific parameterization.
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;

    parameters.doubles("hypoxia_response_threshold") = 1.0; // mapped threshold=0.1
    parameters.doubles("emt_induction_threshold") = 0.25;
    parameters.ints("emt_phenotype_extent") = 2;
    parameters.doubles("hif1a_emt_boost") = 0.2;

    parameters.doubles("mmp2_degradation_rate") = 0.015;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;

    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.0;

    reset_all_fields(1.0, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double cx = 0.0;
    const double cy = 0.0;
    place_tumor_cluster(pTumor, 60, cx, cy, 70.0);

    // Mature, collagen-leaning barrier and hypoxic/tgfb-rich tumor boundary.
    set_annulus_ecm(cx, cy, 20.0, 110.0, 0.8, 0.4);
    set_annulus_field(tgfb_index, cx, cy, 20.0, 110.0, 0.35);
    set_annulus_field(oxygen_index, cx, cy, 0.0, 110.0, 0.0);

    double t = 0.0;
    const double dt = 6.0;

    // Warm-up to establish boundary MMP2+ state.
    advance_steps(1, dt, t);

    std::unordered_set<int> tracked_ids;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != pTumor->type) continue;
        const int hif_idx = custom_index(pCell, "hif1a_active");
        const int zeb_idx = custom_index(pCell, "zeb1_active");
        const int mmp_idx = custom_index(pCell, "mmp2_active");
        if (hif_idx < 0 || zeb_idx < 0 || mmp_idx < 0) continue;

        const double dx = pCell->position[0] - cx;
        const double dy = pCell->position[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r < 20.0 || r > 70.0) continue;

        if (pCell->custom_data[hif_idx] == 1.0 &&
            pCell->custom_data[zeb_idx] == 1.0 &&
            pCell->custom_data[mmp_idx] == 1.0)
        {
            tracked_ids.insert(pCell->ID);
        }
    }

    if (tracked_ids.empty())
    {
        int seeded = 0;
        for (size_t i = 0; i < all_cells->size(); ++i)
        {
            Cell* pCell = (*all_cells)[i];
            if (pCell == NULL || pCell->phenotype.death.dead) continue;
            if (pCell->type != pTumor->type) continue;

            const double dx = pCell->position[0] - cx;
            const double dy = pCell->position[1] - cy;
            const double r = std::sqrt(dx * dx + dy * dy);
            if (r < 20.0 || r > 70.0) continue;

            const int hif_idx = custom_index(pCell, "hif1a_active");
            const int zeb_idx = custom_index(pCell, "zeb1_active");
            const int mmp_idx = custom_index(pCell, "mmp2_active");
            const int cdh_idx = custom_index(pCell, "cdh1_expressed");
            if (hif_idx < 0 || zeb_idx < 0 || mmp_idx < 0 || cdh_idx < 0) continue;

            pCell->custom_data[hif_idx] = 1.0;
            pCell->custom_data[zeb_idx] = 1.0;
            pCell->custom_data[mmp_idx] = 1.0;
            pCell->custom_data[cdh_idx] = 0.0;
            tracked_ids.insert(pCell->ID);
            ++seeded;
            if (seeded >= 12) break;
        }
    }

    assert(!tracked_ids.empty());

    std::vector<double> tracked_ecm_series;
    std::vector<double> tracked_o2_series;
    std::vector<double> boundary_mes_fraction_series;
    std::unordered_set<int> reverted_ids;

    for (int step = 0; step < 100; ++step)
    {
        double ecm_sum = 0.0;
        double o2_sum = 0.0;
        int tracked_live_n = 0;

        int boundary_n = 0;
        int boundary_mes_n = 0;

        for (size_t i = 0; i < all_cells->size(); ++i)
        {
            Cell* pCell = (*all_cells)[i];
            if (pCell == NULL || pCell->phenotype.death.dead) continue;
            if (pCell->type != pTumor->type) continue;

            const int zeb_idx = custom_index(pCell, "zeb1_active");
            const int mmp_idx = custom_index(pCell, "mmp2_active");
            if (zeb_idx < 0 || mmp_idx < 0) continue;

            const double dx = pCell->position[0] - cx;
            const double dy = pCell->position[1] - cy;
            const double r = std::sqrt(dx * dx + dy * dy);
            if (r >= 20.0 && r <= 110.0)
            {
                ++boundary_n;
                if (pCell->custom_data[zeb_idx] == 1.0) ++boundary_mes_n;
            }

            if (tracked_ids.find(pCell->ID) != tracked_ids.end())
            {
                const int voxel = voxel_index_for_cell(pCell);
                const std::vector<double>& rho = microenvironment.density_vector(voxel);
                ecm_sum += rho[ecm_index];
                o2_sum += rho[oxygen_index];
                ++tracked_live_n;

                if (pCell->custom_data[zeb_idx] == 0.0 || pCell->custom_data[mmp_idx] == 0.0)
                {
                    reverted_ids.insert(pCell->ID);
                }
            }
        }

        if (tracked_live_n > 0)
        {
            tracked_ecm_series.push_back(ecm_sum / static_cast<double>(tracked_live_n));
            tracked_o2_series.push_back(o2_sum / static_cast<double>(tracked_live_n));
        }
        else
        {
            tracked_ecm_series.push_back(0.0);
            tracked_o2_series.push_back(0.0);
        }

        const double boundary_mes_frac =
            (boundary_n > 0) ? (static_cast<double>(boundary_mes_n) / static_cast<double>(boundary_n)) : 0.0;
        boundary_mes_fraction_series.push_back(boundary_mes_frac);

        advance_steps(1, dt, t);
    }

    const double ecm_start = tracked_ecm_series.front();
    const double ecm_end = tracked_ecm_series.back();
    const double o2_start = tracked_o2_series.front();
    const double o2_end = tracked_o2_series.back();

    double mes_max = 0.0;
    for (size_t i = 0; i < boundary_mes_fraction_series.size(); ++i)
    {
        mes_max = std::max(mes_max, boundary_mes_fraction_series[i]);
    }
    const double mes_end = boundary_mes_fraction_series.back();

    const int reverted_count = static_cast<int>(reverted_ids.size());

    const bool pass =
        (ecm_end < ecm_start) &&
        (o2_end > o2_start) &&
        (reverted_count > 0) &&
        (has_decrease(boundary_mes_fraction_series) || mes_end < mes_max);

    std::cout << "LOOP4 tracked_ecm_start=" << ecm_start
              << " tracked_ecm_end=" << ecm_end << std::endl;
    std::cout << "LOOP4 tracked_o2_start=" << o2_start
              << " tracked_o2_end=" << o2_end << std::endl;
    std::cout << "LOOP4 reverted_count=" << reverted_count
              << " boundary_mes_end=" << mes_end
              << " boundary_mes_max=" << mes_max << std::endl;

    std::cout << "LOOP4 " << (pass ? "PASS" : "FAIL") << std::endl;
    if (!pass) return 1;

    std::cout << "PASS loop4_emt_self_limiting_test" << std::endl;
    return 0;
}
