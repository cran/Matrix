#ifndef MATRIX_TRMATRIX_H
#define MATRIX_TRMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP dtrMatrix_validate(SEXP obj);
SEXP dtrMatrix_norm(SEXP obj, SEXP type);
SEXP dtrMatrix_rcond(SEXP obj, SEXP type);
SEXP dtrMatrix_solve(SEXP a);
SEXP dtrMatrix_matrix_solve(SEXP a, SEXP b);
SEXP dtrMatrix_as_dgeMatrix(SEXP from);
SEXP dtrMatrix_as_matrix(SEXP from);
SEXP dtrMatrix_getDiag(SEXP x);
void make_array_triangular(double *x, SEXP from);
SEXP dtrMatrix_dgeMatrix_mm(SEXP a, SEXP b);
SEXP dtrMatrix_dgeMatrix_mm_R(SEXP a, SEXP b);

#endif
