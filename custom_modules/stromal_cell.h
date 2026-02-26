#ifndef STROMAL_CELL_H
#define STROMAL_CELL_H

#include "../core/PhysiCell.h"

// Stromal-cell phenotype callback (invoked by custom_function dispatcher through
// update_phenotype registration).
void stromal_phenotype_update(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt);

#endif // STROMAL_CELL_H
