#ifndef CUSTOM_H
#define CUSTOM_H

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

#include "stromal_cell.h"
#include "tumor_cell.h"

// ---------------------------------------------------------------------------
// Global substrate indices (set in setup_microenvironment()).
// Exposed so phenotype modules can reuse a single source of truth.
// ---------------------------------------------------------------------------
extern int oxygen_index;
extern int tgfb_index;
extern int shh_index;
extern int ecm_index;
extern int drug_index;

// ---------------------------------------------------------------------------
// PhysiCell custom module entry points
// ---------------------------------------------------------------------------
void setup_microenvironment(void);
void create_cell_types(void);
void setup_tissue(void);

// Dispatcher for phenotype update.
void custom_function(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt);

// Optional coloring callback used by standard PhysiCell main.cpp templates.
std::vector<std::string> my_coloring_function(PhysiCell::Cell* pCell);

// ---------------------------------------------------------------------------
// ECM barrier microenvironment modifier
// ---------------------------------------------------------------------------
void ecm_dependent_diffusion(double dt);
void ecm_dependent_diffusion_solver(BioFVM::Microenvironment& M, double dt);
void register_ecm_dependent_diffusion_solver(void);

#endif // CUSTOM_H
