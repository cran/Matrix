#ifndef MATRIX_SYMATRIX_H
#define MATRIX_SYMATRIX_H

#include "dgeMatrix.h"
#include "R_ext/Lapack.h"

SEXP dsyMatrix_validate(SEXP obj);
double get_norm_sy(SEXP obj, char *typstr);
SEXP dsyMatrix_norm(SEXP obj, SEXP type);
SEXP dsyMatrix_rcond(SEXP obj, SEXP type);
SEXP dsyMatrix_solve(SEXP a);
SEXP dsyMatrix_matrix_solve(SEXP a, SEXP b);
SEXP dsyMatrix_dgeMatrix_solve(SEXP a, SEXP b);
SEXP dsyMatrix_as_dgeMatrix(SEXP from);
SEXP dsyMatrix_as_dspMatrix(SEXP from);
SEXP dsyMatrix_as_matrix(SEXP from);
SEXP dsyMatrix_dgeMatrix_mm(SEXP a, SEXP b);
SEXP dsyMatrix_dgeMatrix_mm_R(SEXP a, SEXP b);
SEXP dsyMatrix_trf(SEXP x);

#endif
