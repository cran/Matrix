#ifndef MATRIX_SYMATRIX_H
#define MATRIX_SYMATRIX_H

#include "geMatrix.h"
#include "R_ext/Lapack.h"

SEXP syMatrix_validate(SEXP obj);
double get_norm_sy(SEXP obj, char *typstr);
SEXP syMatrix_norm(SEXP obj, SEXP type);
SEXP syMatrix_rcond(SEXP obj, SEXP type);
SEXP syMatrix_solve(SEXP a);
SEXP syMatrix_matrix_solve(SEXP a, SEXP b);
SEXP syMatrix_as_geMatrix(SEXP from);
SEXP syMatrix_as_matrix(SEXP from);
SEXP syMatrix_geMatrix_mm(SEXP a, SEXP b);
SEXP syMatrix_geMatrix_mm_R(SEXP a, SEXP b);

#endif
