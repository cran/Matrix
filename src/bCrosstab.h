#ifndef MATRIX_BCROSSTAB_H
#define MATRIX_BCROSSTAB_H

#include "Mutils.h"
#include "dgCMatrix.h"
#include "Metis_utils.h"
#include "triplet_to_col.h"
#include "R_ldl.h"
#include "dsCMatrix.h"
#include "dtCMatrix.h"

void lmer_populate(SEXP bCtab);

#endif
