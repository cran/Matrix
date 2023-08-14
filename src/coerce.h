#ifndef MATRIX_COERCE_H
#define MATRIX_COERCE_H

#include "Mutils.h"

SEXP matrix_as_dense(SEXP from, const char *zzz, char ul, char di,
                     int transpose_if_vector, int new);

SEXP R_matrix_as_dense(SEXP from, SEXP class, SEXP uplo, SEXP diag);

SEXP sparse_as_dense(SEXP from, const char *class, int packed);

SEXP R_sparse_as_dense(SEXP from, SEXP packed);

SEXP diagonal_as_dense(SEXP from, const char *class,
                       char shape, int packed, char ul);

SEXP R_diagonal_as_dense(SEXP from, SEXP shape, SEXP packed, SEXP uplo);

SEXP index_as_dense(SEXP from, const char *class, char kind);

SEXP R_index_as_dense(SEXP from, SEXP kind);

SEXP matrix_as_sparse(SEXP from, const char *zzz, char ul, char di,
                      int transpose_if_vector);

SEXP R_matrix_as_sparse(SEXP from, SEXP class, SEXP uplo, SEXP diag);

SEXP dense_as_sparse(SEXP from, const char *class, char repr);

SEXP R_dense_as_sparse(SEXP from, SEXP repr);

SEXP diagonal_as_sparse(SEXP from, const char *class,
                        char shape, char repr, char ul);

SEXP R_diagonal_as_sparse(SEXP from, SEXP shape, SEXP repr, SEXP uplo);

SEXP index_as_sparse(SEXP from, const char *class, char kind, char repr);

SEXP R_index_as_sparse(SEXP from, SEXP kind, SEXP repr);

SEXP dense_as_kind(SEXP from, const char *class, char kind);

SEXP R_dense_as_kind(SEXP from, SEXP kind);

SEXP sparse_as_kind(SEXP from, const char *class, char kind);

SEXP R_sparse_as_kind(SEXP from, SEXP kind);

SEXP diagonal_as_kind(SEXP from, const char *class, char kind);

SEXP R_diagonal_as_kind(SEXP from, SEXP kind);

SEXP index_as_kind(SEXP from, const char *class, char kind);

SEXP R_index_as_kind(SEXP from, SEXP kind);

SEXP dense_as_general(SEXP from, const char *class, int new);

SEXP R_dense_as_general(SEXP from);

SEXP sparse_as_general(SEXP from, const char *class);

SEXP R_sparse_as_general(SEXP from);

SEXP dense_as_unpacked(SEXP from, const char *class);

SEXP R_dense_as_unpacked(SEXP from);

SEXP dense_as_packed(SEXP from, const char *class, char ul, char di);

SEXP R_dense_as_packed(SEXP from, SEXP uplo, SEXP diag);

SEXP sparse_as_Csparse(SEXP from, const char *class);

SEXP R_sparse_as_Csparse(SEXP from);

SEXP sparse_as_Rsparse(SEXP from, const char *class);

SEXP R_sparse_as_Rsparse(SEXP from);

SEXP sparse_as_Tsparse(SEXP from, const char *class);

SEXP R_sparse_as_Tsparse(SEXP from);

SEXP R_Matrix_as_vector(SEXP from);

SEXP R_Matrix_as_matrix(SEXP from);

SEXP R_Matrix_as_unpacked(SEXP from);

SEXP R_Matrix_as_packed(SEXP from);

SEXP R_Matrix_as_Csparse(SEXP from);

SEXP R_Matrix_as_Rsparse(SEXP from);

SEXP R_Matrix_as_Tsparse(SEXP from);

SEXP R_Matrix_as_kind(SEXP from, SEXP kind, SEXP sparse);

SEXP R_Matrix_as_general(SEXP from, SEXP kind);

#endif
