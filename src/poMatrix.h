#ifndef MATRIX_POMATRIX_H
#define MATRIX_POMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP poMatrix_rcond(SEXP obj, SEXP type);
SEXP poMatrix_solve(SEXP a);
SEXP poMatrix_matrix_solve(SEXP a, SEXP b);
SEXP poMatrix_geMatrix_solve(SEXP a, SEXP b);
SEXP poMatrix_chol(SEXP x);
double get_norm_sy(SEXP obj, char *typstr);

#endif
