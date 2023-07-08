#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include "Mutils.h"

SEXP matrix_as_dense(SEXP from, const char *code, char uplo, char diag,
                     int new, int transpose_if_vector);
SEXP R_matrix_as_dense(SEXP from, SEXP code, SEXP uplo, SEXP diag);

SEXP R_dense_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP diag);
SEXP R_dense_as_matrix(SEXP from);
SEXP R_geMatrix_as_matrix(SEXP from, SEXP pattern);
SEXP R_dense_as_vector(SEXP from);
SEXP R_geMatrix_as_vector(SEXP from, SEXP pattern);
SEXP R_dense_as_kind(SEXP from, SEXP kind);
SEXP dense_as_general(SEXP from, char kind, int new, int transpose_if_vector);
SEXP R_dense_as_general(SEXP from, SEXP kind);

SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2);
SEXP R_dense_colSums(SEXP obj, SEXP narm, SEXP mean);
SEXP R_dense_rowSums(SEXP obj, SEXP narm, SEXP mean);

#endif /* MATRIX_DENSE_H */
