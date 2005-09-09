#ifndef MATRIX_TRIPLET_H
#define MATRIX_TRIPLET_H

#include "Mutils.h"
#include "triplet_to_col.h"

SEXP dgTMatrix_validate(SEXP x);
SEXP dgTMatrix_to_dgCMatrix(SEXP x);
SEXP dgTMatrix_to_dgeMatrix(SEXP x);
SEXP dgTMatrix_to_matrix(SEXP x);
SEXP graphNEL_as_dgTMatrix(SEXP x, SEXP symmetric);

#endif
