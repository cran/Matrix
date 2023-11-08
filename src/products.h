#ifndef MATRIX_PRODUCTS_H
#define MATRIX_PRODUCTS_H

#include <Rinternals.h>

SEXP R_dense_matmult(SEXP, SEXP, SEXP, SEXP);
SEXP R_sparse_matmult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_diagonal_matmult(SEXP, SEXP, SEXP, SEXP, SEXP);

#endif /* MATRIX_PRODUCTS_H */
