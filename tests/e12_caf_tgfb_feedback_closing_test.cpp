#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

bool nondecreasing(const std::vector<double>& values, double eps = 1e-12)
{
    for (size_t i = 1; i < values.size(); ++i)
    {
        if (values[i] + eps < values[i - 1]) return false;
    }
    return true;
}

double parse_caf_rate(int argc, char** argv)
{
    double caf_rate = 0.3;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg.rfind("--caf-rate=", 0) == 0)
        {
            caf_rate = std::stod(arg.substr(std::string("--caf-rate=").size()));
        }
    }
    return caf_rate;
}

bool parse_control_mode(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--control") return true;
    }
    return false;
}

} // namespace

int main(int argc, char** argv)
{
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

    const double caf_tgfb_rate = parse_caf_rate(argc, argv);
    const bool control_mode = parse_control_mode(argc, argv);

    // Isolate E12 (M3 -> M2 CAF TGF-beta feedback).
    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.0;
    parameters.doubles("caf_tgfb_secretion_rate") = caf_tgfb_rate;

    parameters.doubles("tgfb_activation_threshold") = 0.25;
    parameters.doubles("shh_activation_threshold") = 1e9;

    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("caf_proliferation_rate") = 0.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;

    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;

    parameters.doubles("drug_uptake_rate") = 0.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;

    microenvironment.diffusion_coefficients[tgfb_index] = 20.0;
    microenvironment.decay_rates[tgfb_index] = 0.0005;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double cx = 0.0;
    const double cy = 0.0;
    place_tumor_cluster(pTumor, 5, cx, cy, 25.0);

    const double two_pi = 6.283185307179586;
    for (int i = 0; i < 20; ++i)
    {
        const double theta = two_pi * (static_cast<double>(i) / 20.0);
        const double ring_blend = (i % 4 == 0) ? 0.1 : ((i % 4 == 1) ? 0.35 : ((i % 4 == 2) ? 0.65 : 0.9));
        const double r = 70.0 + 60.0 * ring_blend;
        Cell* p = create_cell(*pStroma);
        p->assign_position(std::vector<double>{cx + r * std::cos(theta), cy + r * std::sin(theta), 0.0});
        p->phenotype.motility.is_motile = false;
        const int acta2_idx = p->custom_data.find_variable_index("acta2_active");
        assert(acta2_idx >= 0);
        p->custom_data[acta2_idx] = 0.0;
    }

    double t = 0.0;
    const double dt = 6.0;

    double tgfb_50 = 0.0;
    double tgfb_100 = 0.0;
    double tgfb_200 = 0.0;
    int caf_50 = 0;
    int caf_100 = 0;
    int caf_200 = 0;

    std::vector<double> tgfb_series;
    std::vector<double> caf_series;

    for (int step = 1; step <= 200; ++step)
    {
        advance_steps(1, dt, t);
        if (step == 50 || step == 100 || step == 200)
        {
            const double tgfb_boundary = mean_field_in_annulus(tgfb_index, cx, cy, 60.0, 150.0);
            const int caf_count = count_activated_cafs(pStroma->type);

            tgfb_series.push_back(tgfb_boundary);
            caf_series.push_back(static_cast<double>(caf_count));

            if (step == 50)
            {
                tgfb_50 = tgfb_boundary;
                caf_50 = caf_count;
            }
            else if (step == 100)
            {
                tgfb_100 = tgfb_boundary;
                caf_100 = caf_count;
            }
            else
            {
                tgfb_200 = tgfb_boundary;
                caf_200 = caf_count;
            }
        }
    }

    const double estimated_caf_feedback_200 = static_cast<double>(caf_200) * caf_tgfb_rate;

    const bool pass_tgfb_amp = (tgfb_200 > tgfb_100) && (tgfb_100 > tgfb_50);
    const bool pass_caf_increase = (caf_200 > caf_50) && (caf_100 >= caf_50);
    const bool pass_tgfb_mono = nondecreasing(tgfb_series);
    const bool pass_caf_mono = nondecreasing(caf_series);
    const bool pass_any_caf = (caf_200 > 0);

    bool scenario_pass = true;
    if (!control_mode)
    {
        scenario_pass =
            pass_tgfb_amp &&
            pass_caf_increase &&
            pass_tgfb_mono &&
            pass_caf_mono &&
            pass_any_caf;
    }

    std::cout << "E12 scenario measurements"
              << " caf_tgfb_rate=" << caf_tgfb_rate
              << " tgfb_step50=" << tgfb_50
              << " tgfb_step100=" << tgfb_100
              << " tgfb_step200=" << tgfb_200
              << " caf_step50=" << caf_50
              << " caf_step100=" << caf_100
              << " caf_step200=" << caf_200
              << " est_caf_feedback_step200=" << estimated_caf_feedback_200
              << std::endl;

    std::cout << "E12 scenario checks"
              << " tgfb_amp=" << (pass_tgfb_amp ? 1 : 0)
              << " caf_increase=" << (pass_caf_increase ? 1 : 0)
              << " tgfb_mono=" << (pass_tgfb_mono ? 1 : 0)
              << " caf_mono=" << (pass_caf_mono ? 1 : 0)
              << " any_caf=" << (pass_any_caf ? 1 : 0)
              << std::endl;

    std::cout << "E12 scenario tgfb200=" << tgfb_200 << " caf200=" << caf_200 << std::endl;

    if (!scenario_pass)
    {
        std::cout << "FAIL E12 scenario" << std::endl;
        return 1;
    }

    if (control_mode)
    {
        std::cout << "PASS E12 control scenario" << std::endl;
    }
    else
    {
        std::cout << "PASS E12 scenario" << std::endl;
    }
    return 0;
}
