#ifndef TUMOR_CELL_H
#define TUMOR_CELL_H

#include "../core/PhysiCell.h"
#include "boolean_network.h"

#include <vector>

// Tumor-cell phenotype callback (invoked by custom_function dispatcher through
// update_phenotype registration).
void tumor_phenotype_update(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt);

// Shared runtime helpers used by both tumor and stromal phenotype functions.
BooleanNetwork* get_boolean_network(PhysiCell::Cell* pCell, CellType cell_type);

const std::vector<Intervention>& get_current_interventions();
void set_current_interventions(const std::vector<Intervention>& interventions);

const ThresholdConfig& get_threshold_config();
void set_threshold_config(const ThresholdConfig& cfg);

#endif // TUMOR_CELL_H
