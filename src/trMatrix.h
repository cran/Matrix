#ifndef MATRIX_TRMATRIX_H
#define MATRIX_TRMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP trMatrix_validate(SEXP obj);
SEXP trMatrix_norm(SEXP obj, SEXP type);
SEXP trMatrix_rcond(SEXP obj, SEXP type);
SEXP trMatrix_solve(SEXP a);
SEXP trMatrix_matrix_solve(SEXP a, SEXP b);
SEXP trMatrix_as_geMatrix(SEXP from);
SEXP trMatrix_as_matrix(SEXP from);
SEXP trMatrix_getDiag(SEXP x);
void make_array_triangular(double *x, SEXP from);
SEXP trMatrix_geMatrix_mm(SEXP a, SEXP b);
SEXP trMatrix_geMatrix_mm_R(SEXP a, SEXP b);

#endif
