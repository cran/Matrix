#ifndef MATRIX_SYMATRIX_H
#define MATRIX_SYMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm); /* factorizations.c */

SEXP dsyMatrix_trf_(SEXP obj,  int warn);
SEXP dsyMatrix_trf (SEXP obj, SEXP warn);
SEXP    matrix_trf_(SEXP obj,  int warn, char uplo);
SEXP    matrix_trf (SEXP obj, SEXP warn, SEXP uplo);

double get_norm_dsy(SEXP obj, const char *typstr);
SEXP dsyMatrix_norm(SEXP obj, SEXP type);
SEXP dsyMatrix_rcond(SEXP obj);
SEXP dsyMatrix_determinant(SEXP obj, SEXP logarithm);
SEXP dsyMatrix_solve(SEXP a);
SEXP dsyMatrix_matrix_solve(SEXP a, SEXP b);

SEXP dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rt);

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0
SEXP dsyMatrix_as_dspMatrix(SEXP from);
#endif

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0
SEXP dsyMatrix_as_matrix(SEXP from, SEXP keep_dimnames);
#endif

#endif
