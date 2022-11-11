#ifndef MATRIX_TSC_H
#define MATRIX_TSC_H

#include "Mutils.h"
#include "dgCMatrix.h"

SEXP dtCMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed);
SEXP dtCMatrix_sparse_solve(SEXP a, SEXP b);

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0
SEXP tCMatrix_validate(SEXP x);
SEXP tRMatrix_validate(SEXP x);
#endif /* MJ */

#endif
