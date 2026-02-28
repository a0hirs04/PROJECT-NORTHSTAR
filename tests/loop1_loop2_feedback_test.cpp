#include <cassert>
#include <iostream>
#include <vector>

#include "feedback_loop_common.h"

using namespace feedback_loop_common;

namespace
{

struct Loop1Metrics
{
    int step = 0;
    double tgfb_boundary = 0.0;
    int caf_count = 0;
    double ecm_peri = 0.0;
};

bool nondecreasing(const std::vector<double>& values, double eps = 1e-9)
{
    for (size_t i = 1; i < values.size(); ++i)
    {
        if (values[i] + eps < values[i - 1]) return false;
    }
    return true;
}

double tumor_hif1a_fraction(int tumor_type)
{
    int total = 0;
    int active = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != tumor_type) continue;
        const int hif_idx = custom_index(pCell, "hif1a_active");
        if (hif_idx < 0) continue;
        ++total;
        if (pCell->custom_data[hif_idx] == 1.0) ++active;
    }
    if (total == 0) return 0.0;
    return static_cast<double>(active) / static_cast<double>(total);
}

std::vector<double> tumor_centroid(int tumor_type)
{
    std::vector<double> c(3, 0.0);
    int total = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != tumor_type) continue;
        c[0] += pCell->position[0];
        c[1] += pCell->position[1];
        c[2] += pCell->position[2];
        ++total;
    }
    if (total == 0) return c;
    c[0] /= static_cast<double>(total);
    c[1] /= static_cast<double>(total);
    c[2] /= static_cast<double>(total);
    return c;
}

void tumor_tgfb_secretion_means(int tumor_type,
                                double& hypoxic_mean,
                                double& normoxic_mean,
                                int& hypoxic_n,
                                int& normoxic_n)
{
    hypoxic_mean = 0.0;
    normoxic_mean = 0.0;
    hypoxic_n = 0;
    normoxic_n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != tumor_type) continue;
        const int hif_idx = custom_index(pCell, "hif1a_active");
        if (hif_idx < 0) continue;
        const double sec = (tgfb_index >= 0 &&
                            tgfb_index < static_cast<int>(pCell->phenotype.secretion.secretion_rates.size()))
            ? pCell->phenotype.secretion.secretion_rates[tgfb_index]
            : 0.0;
        if (pCell->custom_data[hif_idx] == 1.0)
        {
            hypoxic_mean += sec;
            ++hypoxic_n;
        }
        else
        {
            normoxic_mean += sec;
            ++normoxic_n;
        }
    }
    if (hypoxic_n > 0) hypoxic_mean /= static_cast<double>(hypoxic_n);
    if (normoxic_n > 0) normoxic_mean /= static_cast<double>(normoxic_n);
}

} // namespace

int main()
{
    initialize_world();

    assert(tgfb_index >= 0);
    assert(shh_index >= 0);
    assert(ecm_index >= 0);
    assert(oxygen_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);
    const int tumor_type = pTumor->type;
    const int stromal_type = pStroma->type;

    // Promote clear loop dynamics and suppress unrelated pathways.
    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.5;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.8;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.30;

    parameters.doubles("tgfb_activation_threshold") = 0.22;
    parameters.doubles("shh_activation_threshold") = 0.05;
    parameters.doubles("ecm_production_rate_base") = 0.004;
    parameters.doubles("ecm_production_rate_boosted") = 0.007;
    parameters.doubles("ecm_ha_fraction_default") = 0.6;

    parameters.doubles("hypoxia_response_threshold") = 0.9; // mapped threshold=0.09
    parameters.doubles("emt_induction_threshold") = 10.0;   // keep EMT/MMP2 out of Loops 1-2
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.ints("emt_phenotype_extent") = 2;
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("caf_proliferation_rate") = 0.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;

    parameters.doubles("drug_uptake_rate") = 0.0;
    parameters.doubles("drug_stress_threshold") = 10.0;
    parameters.doubles("efflux_induction_delay") = 1e9;
    parameters.doubles("efflux_strength") = 0.0;
    parameters.doubles("nrf2_decay_rate") = 0.0;
    parameters.doubles("drug_natural_decay_rate") = 0.0;
    parameters.doubles("hif1a_nrf2_priming_bonus") = 0.0;

    // Slow TGF-beta/SHH decay so positive feedback can amplify over 200 steps.
    microenvironment.decay_rates[tgfb_index] = 0.001;
    microenvironment.decay_rates[shh_index] = 0.001;
    microenvironment.decay_rates[oxygen_index] = 0.01;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double cx = 0.0;
    const double cy = 0.0;
    place_tumor_cluster(pTumor, 10, cx, cy, 35.0);
    place_stroma_ring(pStroma, 30, cx, cy, 60.0, 110.0);

    double t = 0.0;
    const double dt = 6.0;

    std::vector<Loop1Metrics> checkpoints;
    for (int step = 1; step <= 200; ++step)
    {
        advance_steps(1, dt, t);
        if (step == 50 || step == 100 || step == 150 || step == 200)
        {
            Loop1Metrics m;
            m.step = step;
            m.tgfb_boundary = mean_field_in_annulus(tgfb_index, cx, cy, 50.0, 130.0);
            m.caf_count = count_activated_cafs(stromal_type);
            m.ecm_peri = mean_ecm_in_annulus(cx, cy, 50.0, 130.0);
            checkpoints.push_back(m);
        }
    }

    assert(checkpoints.size() == 4);
    const Loop1Metrics m50 = checkpoints[0];
    const Loop1Metrics m200 = checkpoints[3];

    std::vector<double> tgfb_series;
    std::vector<double> caf_series;
    std::vector<double> ecm_series;
    for (size_t i = 0; i < checkpoints.size(); ++i)
    {
        tgfb_series.push_back(checkpoints[i].tgfb_boundary);
        caf_series.push_back(static_cast<double>(checkpoints[i].caf_count));
        ecm_series.push_back(checkpoints[i].ecm_peri);
    }

    const bool loop1_pass =
        (m200.tgfb_boundary > m50.tgfb_boundary) &&
        (m200.caf_count > m50.caf_count) &&
        (m200.ecm_peri > m50.ecm_peri) &&
        nondecreasing(tgfb_series) &&
        nondecreasing(caf_series) &&
        nondecreasing(ecm_series);

    std::cout << "LOOP1 step50 tgfb_boundary=" << m50.tgfb_boundary
              << " caf_count=" << m50.caf_count
              << " ecm_peri=" << m50.ecm_peri << std::endl;
    std::cout << "LOOP1 step200 tgfb_boundary=" << m200.tgfb_boundary
              << " caf_count=" << m200.caf_count
              << " ecm_peri=" << m200.ecm_peri << std::endl;
    std::cout << "LOOP1 monotonic tgfb=" << (nondecreasing(tgfb_series) ? 1 : 0)
              << " caf=" << (nondecreasing(caf_series) ? 1 : 0)
              << " ecm=" << (nondecreasing(ecm_series) ? 1 : 0)
              << std::endl;

    // Loop 2 starts from Loop 1 endpoint.
    const std::vector<double> center_pos = {cx, cy, 0.0};
    const double o2_200 = substrate_at_position(oxygen_index, center_pos);
    const double hif_frac_200 = tumor_hif1a_fraction(tumor_type);

    for (int step = 201; step <= 400; ++step)
    {
        advance_steps(1, dt, t);
    }

    const double o2_400 = substrate_at_position(oxygen_index, center_pos);
    const double hif_frac_400 = tumor_hif1a_fraction(tumor_type);

    double hypoxic_sec = 0.0;
    double normoxic_sec = 0.0;
    int hypoxic_n = 0;
    int normoxic_n = 0;
    tumor_tgfb_secretion_means(tumor_type, hypoxic_sec, normoxic_sec, hypoxic_n, normoxic_n);

    const bool loop2_pass =
        (o2_400 < o2_200) &&
        (hif_frac_400 > hif_frac_200) &&
        (hypoxic_n > 0) &&
        (normoxic_n > 0) &&
        (hypoxic_sec > normoxic_sec);

    std::cout << "LOOP2 step200 o2_centroid=" << o2_200
              << " hif1a_fraction=" << hif_frac_200 << std::endl;
    std::cout << "LOOP2 step400 o2_centroid=" << o2_400
              << " hif1a_fraction=" << hif_frac_400 << std::endl;
    std::cout << "LOOP2 secretion hypoxic_mean=" << hypoxic_sec
              << " normoxic_mean=" << normoxic_sec
              << " hypoxic_n=" << hypoxic_n
              << " normoxic_n=" << normoxic_n << std::endl;

    std::cout << "LOOP1 " << (loop1_pass ? "PASS" : "FAIL") << std::endl;
    std::cout << "LOOP2 " << (loop2_pass ? "PASS" : "FAIL") << std::endl;

    if (!(loop1_pass && loop2_pass))
    {
        return 1;
    }

    std::cout << "PASS loop1_loop2_feedback_test" << std::endl;
    return 0;
}
