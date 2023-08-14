#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include "Mutils.h"

SEXP sparse_drop0(SEXP from, const char *class, double tol);
SEXP R_sparse_drop0(SEXP from, SEXP tol);

SEXP sparse_band(SEXP from, const char *class, int a, int b);
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2);

SEXP sparse_diag_get(SEXP obj, const char *class, int names);
SEXP R_sparse_diag_get(SEXP obj, SEXP names);

SEXP sparse_diag_set(SEXP from, const char *class, SEXP value);
SEXP R_sparse_diag_set(SEXP from, SEXP value);

SEXP sparse_diag_U2N(SEXP from, const char *class);
SEXP R_sparse_diag_U2N(SEXP from);

SEXP sparse_diag_N2U(SEXP from, const char *class);
SEXP R_sparse_diag_N2U(SEXP from);

SEXP sparse_transpose(SEXP from, const char *class, int lazy);
SEXP R_sparse_transpose(SEXP from, SEXP lazy);

SEXP sparse_force_symmetric(SEXP from, const char *class, char ul);
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo);

SEXP sparse_symmpart(SEXP from, const char *class);
SEXP R_sparse_symmpart(SEXP from);

SEXP sparse_skewpart(SEXP from, const char *class);
SEXP R_sparse_skewpart(SEXP from);

SEXP Tsparse_aggregate(SEXP from);

SEXP Csparse_is_diagonal(SEXP obj);
SEXP Rsparse_is_diagonal(SEXP obj);
SEXP Tsparse_is_diagonal(SEXP obj);
SEXP Csparse_is_triangular(SEXP obj, SEXP upper);
SEXP Rsparse_is_triangular(SEXP obj, SEXP upper);
SEXP Tsparse_is_triangular(SEXP obj, SEXP upper);
SEXP Csparse_is_symmetric(SEXP obj, SEXP checkDN);
SEXP Rsparse_is_symmetric(SEXP obj, SEXP checkDN);
#if 0 /* unimplemented ... currently going via CsparseMatrix */
SEXP Tsparse_is_symmetric(SEXP obj, SEXP checkDN);
#endif

SEXP CRsparse_colSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse);
SEXP CRsparse_rowSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse);

#endif
