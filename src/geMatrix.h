#ifndef MATRIX_GEMATRIX_H
#define MATRIX_GEMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP geMatrix_validate(SEXP obj);
SEXP geMatrix_norm(SEXP obj, SEXP norm);
SEXP geMatrix_crossprod(SEXP x);
SEXP geMatrix_geMatrix_crossprod(SEXP x, SEXP y);
SEXP geMatrix_matrix_crossprod(SEXP x, SEXP y);
SEXP geMatrix_getDiag(SEXP x);
SEXP geMatrix_LU(SEXP x);
SEXP geMatrix_determinant(SEXP x, SEXP logarithm);
SEXP geMatrix_solve(SEXP a);
SEXP geMatrix_geMatrix_mm(SEXP a, SEXP b);

#endif
