#ifndef MATRIX_TRT_H
#define MATRIX_TRT_H

#include "Mutils.h"
#include "chm_common.h"
#include "triplet_to_col.h"

SEXP dtTMatrix_validate(SEXP x);
SEXP dtTMatrix_as_dtrMatrix(SEXP x);
SEXP dtTMatrix_as_dgCMatrix(SEXP x);

#endif
