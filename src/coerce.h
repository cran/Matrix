#ifndef MATRIX_COERCE_H
#define MATRIX_COERCE_H

#include <Rinternals.h>

SEXP vector_as_dense(SEXP, const char *, char, char, int, int, int, SEXP);
SEXP R_vector_as_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP matrix_as_dense(SEXP, const char *, char, char, int, int);
SEXP R_matrix_as_dense(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP sparse_as_dense(SEXP, const char *, int);
SEXP R_sparse_as_dense(SEXP, SEXP);

SEXP diagonal_as_dense(SEXP, const char *, char, char, int, char);
SEXP R_diagonal_as_dense(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP index_as_dense(SEXP, const char *, char);
SEXP R_index_as_dense(SEXP, SEXP);

SEXP vector_as_sparse(SEXP, const char *, char, char, int, int, int, SEXP);
SEXP R_vector_as_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP matrix_as_sparse(SEXP, const char *, char, char, int);
SEXP R_matrix_as_sparse(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP dense_as_sparse(SEXP, const char *, char);
SEXP R_dense_as_sparse(SEXP, SEXP);

SEXP diagonal_as_sparse(SEXP, const char *, char, char, char, char);
SEXP R_diagonal_as_sparse(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP index_as_sparse(SEXP, const char *, char, char);
SEXP R_index_as_sparse(SEXP, SEXP, SEXP);

SEXP dense_as_kind(SEXP, const char *, char, int);
SEXP R_dense_as_kind(SEXP, SEXP);

SEXP sparse_as_kind(SEXP, const char *, char);
SEXP R_sparse_as_kind(SEXP, SEXP);

SEXP diagonal_as_kind(SEXP, const char *, char);
SEXP R_diagonal_as_kind(SEXP, SEXP);

SEXP index_as_kind(SEXP, const char *, char);
SEXP R_index_as_kind(SEXP, SEXP);

SEXP dense_as_general(SEXP, const char *, int);
SEXP R_dense_as_general(SEXP);

SEXP sparse_as_general(SEXP, const char *);
SEXP R_sparse_as_general(SEXP);

SEXP dense_as_unpacked(SEXP, const char *);
SEXP R_dense_as_unpacked(SEXP);

SEXP dense_as_packed(SEXP, const char *, char, char);
SEXP R_dense_as_packed(SEXP, SEXP, SEXP);

SEXP sparse_as_Csparse(SEXP, const char *);
SEXP R_sparse_as_Csparse(SEXP);

SEXP sparse_as_Rsparse(SEXP, const char *);
SEXP R_sparse_as_Rsparse(SEXP);

SEXP sparse_as_Tsparse(SEXP, const char *);
SEXP R_sparse_as_Tsparse(SEXP);

SEXP R_Matrix_as_vector(SEXP);
SEXP R_Matrix_as_matrix(SEXP);
SEXP R_Matrix_as_unpacked(SEXP);
SEXP R_Matrix_as_packed(SEXP);
SEXP R_Matrix_as_Csparse(SEXP);
SEXP R_Matrix_as_Rsparse(SEXP);
SEXP R_Matrix_as_Tsparse(SEXP);
SEXP R_Matrix_as_kind(SEXP, SEXP, SEXP);
SEXP R_Matrix_as_general(SEXP, SEXP);

#endif /* MATRIX_COERCE_H */
