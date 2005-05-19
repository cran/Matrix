#ifndef MATRIX_LSCMATRIX_H
#define MATRIX_LSCMATRIX_H

#include <Rdefines.h>
#include "Mutils.h"
#include "triplet_to_col.h"
#include "Metis_utils.h"
#include "R_ldl.h"

SEXP lsCMatrix_validate(SEXP x);
SEXP lsCMatrix_trans(SEXP x);
SEXP lsCMatrix_chol(SEXP x, SEXP pivot);

#endif
