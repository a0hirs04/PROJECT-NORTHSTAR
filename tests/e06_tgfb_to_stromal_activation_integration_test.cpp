#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

struct PscRecord
{
    Cell* cell = NULL;
    double radius = 0.0;
    int activation_step = -1;
    double tgfb_at_activation = 0.0;
};

double local_tgfb(Cell* pCell)
{
    std::vector<double> p = pCell->position;
    const int voxel = microenvironment.nearest_voxel_index(p);
    return microenvironment.density_vector(voxel)[tgfb_index];
}

double mean_tgfb(const std::vector<PscRecord>& ps)
{
    double sum = 0.0;
    for (size_t i = 0; i < ps.size(); ++i) sum += local_tgfb(ps[i].cell);
    return ps.empty() ? 0.0 : sum / static_cast<double>(ps.size());
}

int count_activated(const std::vector<PscRecord>& ps, int acta2_idx)
{
    int count = 0;
    for (size_t i = 0; i < ps.size(); ++i)
    {
        if (ps[i].cell->custom_data[acta2_idx] == 1.0) ++count;
    }
    return count;
}

} // namespace

int main()
{
    // Mute verbose PhysiCell setup prints; keep test metrics readable.
    std::ofstream null_out("/dev/null");
    std::streambuf* cout_backup = std::cout.rdbuf(null_out.rdbuf());
    initialize_world();
    std::cout.rdbuf(cout_backup);
    null_out.close();

    assert(tgfb_index >= 0);
    assert(shh_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    // Isolation of E06 (TGF-beta arm only).
    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.0;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.0;
    parameters.doubles("tgfb_activation_threshold") = 0.3;
    parameters.doubles("shh_activation_threshold") = 1e9;
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("caf_proliferation_rate") = 0.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;
    parameters.doubles("drug_uptake_rate") = 0.0;
    microenvironment.diffusion_coefficients[tgfb_index] = 20.0;
    microenvironment.decay_rates[tgfb_index] = 0.0;

    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double cx = 0.0;
    const double cy = 0.0;

    // 5 tumor cells at center.
    const double two_pi = 6.283185307179586;
    for (int i = 0; i < 5; ++i)
    {
        const double frac = (static_cast<double>(i) + 0.5) / 5.0;
        const double r = 20.0 * std::sqrt(frac);
        const double theta = two_pi * frac;
        Cell* tcell = create_cell(*pTumor);
        tcell->assign_position(std::vector<double>{cx + r * std::cos(theta), cy + r * std::sin(theta), 0.0});
        tcell->phenotype.motility.is_motile = false;
        if (tgfb_index >= 0 && tgfb_index < static_cast<int>(tcell->phenotype.secretion.secretion_rates.size()))
        {
            tcell->phenotype.secretion.secretion_rates[tgfb_index] = 0.5;
        }
        if (shh_index >= 0 && shh_index < static_cast<int>(tcell->phenotype.secretion.secretion_rates.size()))
        {
            tcell->phenotype.secretion.secretion_rates[shh_index] = 0.0;
        }
    }

    // 10 PSCs around ~100um; near ring (80um) and far ring (120um).
    std::vector<PscRecord> pscs;
    pscs.reserve(10);
    for (int i = 0; i < 10; ++i)
    {
        const double theta = two_pi * (static_cast<double>(i) / 10.0);
        const double radius = (i < 5) ? 80.0 : 120.0;
        Cell* p = create_cell(*pStroma);
        p->assign_position(std::vector<double>{cx + radius * std::cos(theta), cy + radius * std::sin(theta), 0.0});
        p->phenotype.motility.is_motile = false;
        const int acta2_idx_local = p->custom_data.find_variable_index("acta2_active");
        assert(acta2_idx_local >= 0);
        p->custom_data[acta2_idx_local] = 0.0;
        pscs.push_back(PscRecord{p, radius, -1});
    }

    const int acta2_idx = pscs.front().cell->custom_data.find_variable_index("acta2_active");
    assert(acta2_idx >= 0);

    const double dt = 10.0;
    double t = 0.0;
    Cell_Container* cc = cell_container();

    double tgfb10 = 0.0;
    double tgfb50 = 0.0;
    double tgfb100 = 0.0;
    int act10 = 0;
    int act50 = 0;
    int act100 = 0;

    for (int step = 1; step <= 100; ++step)
    {
        microenvironment.simulate_diffusion_decay(dt);
        // Enforce SHH=0 to isolate TGF-beta-only E06 arm.
        const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
        for (int n = 0; n < n_vox; ++n)
        {
            microenvironment.density_vector(n)[shh_index] = 0.0;
        }
        cc->update_all_cells(t, dt, dt, dt);
        t += dt;

        for (size_t i = 0; i < pscs.size(); ++i)
        {
            if (pscs[i].activation_step < 0 && pscs[i].cell->custom_data[acta2_idx] == 1.0)
            {
                pscs[i].activation_step = step;
                pscs[i].tgfb_at_activation = local_tgfb(pscs[i].cell);
            }
        }

        if (step == 10)
        {
            tgfb10 = mean_tgfb(pscs);
            act10 = count_activated(pscs, acta2_idx);
        }
        if (step == 50)
        {
            tgfb50 = mean_tgfb(pscs);
            act50 = count_activated(pscs, acta2_idx);
        }
        if (step == 100)
        {
            tgfb100 = mean_tgfb(pscs);
            act100 = count_activated(pscs, acta2_idx);
        }
    }

    // Near-vs-far activation timing.
    int near_first = std::numeric_limits<int>::max();
    int far_first = std::numeric_limits<int>::max();
    for (size_t i = 0; i < pscs.size(); ++i)
    {
        if (pscs[i].activation_step < 0) continue;
        if (pscs[i].radius <= 100.0) near_first = std::min(near_first, pscs[i].activation_step);
        else far_first = std::min(far_first, pscs[i].activation_step);
    }

    // Threshold-consistency check at final time.
    const double activation_threshold = parameters.doubles("tgfb_activation_threshold");
    bool threshold_rule_ok = true;
    for (size_t i = 0; i < pscs.size(); ++i)
    {
        const double tg = local_tgfb(pscs[i].cell);
        const bool active = (pscs[i].cell->custom_data[acta2_idx] == 1.0);
        if (tg > activation_threshold && !active) threshold_rule_ok = false;
        if (!active && tg < activation_threshold) continue;
        if (!active && tg >= activation_threshold) threshold_rule_ok = false;
        if (active && tg < activation_threshold &&
            !(pscs[i].activation_step >= 0 && pscs[i].tgfb_at_activation > activation_threshold))
        {
            threshold_rule_ok = false;
        }
    }

    if (!threshold_rule_ok)
    {
        for (size_t i = 0; i < pscs.size(); ++i)
        {
            const double tg = local_tgfb(pscs[i].cell);
            const bool active = (pscs[i].cell->custom_data[acta2_idx] == 1.0);
            std::cout << "E06 debug psc=" << i
                      << " radius=" << pscs[i].radius
                      << " tgfb_final=" << tg
                      << " active=" << (active ? 1 : 0)
                      << " activation_step=" << pscs[i].activation_step
                      << " tgfb_at_activation=" << pscs[i].tgfb_at_activation
                      << std::endl;
        }
    }

    const bool pass_tgfb_increases = (tgfb10 < tgfb50 && tgfb50 < tgfb100);
    const bool pass_some_activate = (act100 > 0);
    const bool near_has_activation = (near_first != std::numeric_limits<int>::max());
    const bool far_has_activation = (far_first != std::numeric_limits<int>::max());
    const bool pass_near_before_far =
        near_has_activation && (!far_has_activation || near_first < far_first);
    const bool pass_threshold_rule = threshold_rule_ok;

    std::cout << "E06 measurements:"
              << " tgfb_step10=" << tgfb10
              << " tgfb_step50=" << tgfb50
              << " tgfb_step100=" << tgfb100
              << " act_step10=" << act10
              << " act_step50=" << act50
              << " act_step100=" << act100
              << " near_first_activation_step=" << (near_first == std::numeric_limits<int>::max() ? -1 : near_first)
              << " far_first_activation_step=" << (far_first == std::numeric_limits<int>::max() ? -1 : far_first)
              << std::endl;
    std::cout << "E06 checks:"
              << " tgfb_increases=" << (pass_tgfb_increases ? 1 : 0)
              << " some_activate=" << (pass_some_activate ? 1 : 0)
              << " near_before_far=" << (pass_near_before_far ? 1 : 0)
              << " threshold_rule=" << (pass_threshold_rule ? 1 : 0)
              << std::endl;

    if (!(pass_tgfb_increases && pass_some_activate && pass_near_before_far && pass_threshold_rule))
    {
        std::cout << "FAIL E06 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E06 integration" << std::endl;
    return 0;
}
