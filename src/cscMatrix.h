#ifndef MATRIX_CSC_H
#define MATRIX_CSC_H

#include <Rdefines.h>
#include "Mutils.h"
#include "R_ldl.h"

SEXP csc_crossprod(SEXP x);
SEXP csc_matrix_crossprod(SEXP x, SEXP y);
SEXP csc_validate(SEXP x);
SEXP csc_to_triplet(SEXP x);
SEXP csc_to_matrix(SEXP x);
SEXP csc_to_geMatrix(SEXP x);
SEXP matrix_to_csc(SEXP A);
SEXP triplet_to_csc(SEXP triplet);
SEXP csc_getDiag(SEXP x);
SEXP csc_transpose(SEXP x);
SEXP csc_matrix_mm(SEXP a, SEXP b);
SEXP csc_col_permute(SEXP x, SEXP perm);

#endif
