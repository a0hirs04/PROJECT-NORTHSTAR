#ifndef FEEDBACK_LOOP_COMMON_H
#define FEEDBACK_LOOP_COMMON_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "../Stroma_world/PhysiCell/core/PhysiCell.h"
#include "../Stroma_world/PhysiCell/modules/PhysiCell_standard_modules.h"
#include "../custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

namespace feedback_loop_common
{

inline void initialize_world()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    intervention_state.ha_degrade_active = false;
    intervention_state.col_degrade_active = false;
    intervention_state.ha_degrade_strength = 0.0;
    intervention_state.col_degrade_strength = 0.0;
}

inline Cell_Container* cell_container()
{
    return (Cell_Container*)microenvironment.agent_container;
}

inline void advance_steps(int n_steps, double dt, double& t)
{
    Cell_Container* cc = cell_container();
    for (int n = 0; n < n_steps; ++n)
    {
        microenvironment.simulate_diffusion_decay(dt);
        cc->update_all_cells(t, dt, dt, dt);
        t += dt;
    }
}

inline int voxel_index_for_position(double x, double y, double z = 0.0)
{
    std::vector<double> p(3, 0.0);
    p[0] = x;
    p[1] = y;
    p[2] = z;
    return microenvironment.nearest_voxel_index(p);
}

inline int voxel_index_for_cell(Cell* pCell)
{
    std::vector<double> p = pCell->position;
    return microenvironment.nearest_voxel_index(p);
}

inline void reset_all_fields(double oxygen_value,
                             double tgfb_value,
                             double shh_value,
                             double ecm_value,
                             double drug_value,
                             double ha_fraction_value)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        std::vector<double>& rho = microenvironment.density_vector(n);
        if (oxygen_index >= 0 && oxygen_index < static_cast<int>(rho.size())) rho[oxygen_index] = oxygen_value;
        if (tgfb_index >= 0 && tgfb_index < static_cast<int>(rho.size())) rho[tgfb_index] = tgfb_value;
        if (shh_index >= 0 && shh_index < static_cast<int>(rho.size())) rho[shh_index] = shh_value;
        if (ecm_index >= 0 && ecm_index < static_cast<int>(rho.size())) rho[ecm_index] = ecm_value;
        if (drug_index >= 0 && drug_index < static_cast<int>(rho.size())) rho[drug_index] = drug_value;
        set_ecm_ha_fraction(n, ha_fraction_value);
    }
}

inline void set_annulus_field(int substrate_index,
                              double cx,
                              double cy,
                              double r_inner,
                              double r_outer,
                              double value)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        const double dx = c[0] - cx;
        const double dy = c[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r >= r_inner && r <= r_outer)
        {
            std::vector<double>& rho = microenvironment.density_vector(n);
            if (substrate_index >= 0 && substrate_index < static_cast<int>(rho.size()))
            {
                rho[substrate_index] = value;
            }
        }
    }
}

inline void set_annulus_ecm(double cx,
                            double cy,
                            double r_inner,
                            double r_outer,
                            double ecm_value,
                            double ha_fraction_value)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        const double dx = c[0] - cx;
        const double dy = c[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r >= r_inner && r <= r_outer)
        {
            std::vector<double>& rho = microenvironment.density_vector(n);
            if (ecm_index >= 0 && ecm_index < static_cast<int>(rho.size()))
            {
                rho[ecm_index] = ecm_value;
            }
            set_ecm_ha_fraction(n, ha_fraction_value);
        }
    }
}

inline double mean_field_in_annulus(int substrate_index,
                                    double cx,
                                    double cy,
                                    double r_inner,
                                    double r_outer)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    double sum = 0.0;
    int count = 0;
    for (int n = 0; n < n_vox; ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        const double dx = c[0] - cx;
        const double dy = c[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r >= r_inner && r <= r_outer)
        {
            const std::vector<double>& rho = microenvironment.density_vector(n);
            if (substrate_index >= 0 && substrate_index < static_cast<int>(rho.size()))
            {
                sum += rho[substrate_index];
                ++count;
            }
        }
    }
    if (count == 0) return 0.0;
    return sum / static_cast<double>(count);
}

inline double mean_ecm_in_annulus(double cx, double cy, double r_inner, double r_outer)
{
    return mean_field_in_annulus(ecm_index, cx, cy, r_inner, r_outer);
}

inline int custom_index(Cell* pCell, const char* name)
{
    return pCell->custom_data.find_variable_index(name);
}

inline std::vector<Cell*> live_cells_of_type(int type_id)
{
    std::vector<Cell*> out;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL) continue;
        if (pCell->phenotype.death.dead) continue;
        if (pCell->type == type_id) out.push_back(pCell);
    }
    return out;
}

inline double substrate_at_position(int substrate_index, const std::vector<double>& pos)
{
    std::vector<double> p = pos;
    const int voxel = microenvironment.nearest_voxel_index(p);
    const std::vector<double>& rho = microenvironment.density_vector(voxel);
    if (substrate_index < 0 || substrate_index >= static_cast<int>(rho.size())) return 0.0;
    return rho[substrate_index];
}

inline void place_tumor_cluster(Cell_Definition* pTumorDef,
                                int n_cells,
                                double cx,
                                double cy,
                                double radius)
{
    const double two_pi = 6.283185307179586;
    for (int i = 0; i < n_cells; ++i)
    {
        const double frac = (static_cast<double>(i) + 0.5) / static_cast<double>(n_cells);
        const double r = radius * std::sqrt(frac);
        const double theta = two_pi * frac * 2.0;
        Cell* pCell = create_cell(*pTumorDef);
        pCell->assign_position(std::vector<double>{cx + r * std::cos(theta), cy + r * std::sin(theta), 0.0});
    }
}

inline void place_stroma_ring(Cell_Definition* pStromaDef,
                              int n_cells,
                              double cx,
                              double cy,
                              double r_inner,
                              double r_outer)
{
    const double two_pi = 6.283185307179586;
    for (int i = 0; i < n_cells; ++i)
    {
        const double frac = (static_cast<double>(i) + 0.5) / static_cast<double>(n_cells);
        const double theta = two_pi * frac;
        const double blend = (i % 3 == 0) ? 0.15 : ((i % 3 == 1) ? 0.5 : 0.85);
        const double r = r_inner + (r_outer - r_inner) * blend;
        Cell* pCell = create_cell(*pStromaDef);
        pCell->assign_position(std::vector<double>{cx + r * std::cos(theta), cy + r * std::sin(theta), 0.0});
    }
}

inline int count_activated_cafs(int stromal_type_id)
{
    int count = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != stromal_type_id) continue;
        const int acta2_idx = custom_index(pCell, "acta2_active");
        if (acta2_idx >= 0 && pCell->custom_data[acta2_idx] == 1.0) ++count;
    }
    return count;
}

} // namespace feedback_loop_common

#endif
