#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include <Rinternals.h>

SEXP dense_band(SEXP, const char *, int, int);
SEXP R_dense_band(SEXP, SEXP, SEXP);

SEXP dense_diag_get(SEXP, const char *, int);
SEXP R_dense_diag_get(SEXP, SEXP);

SEXP dense_diag_set(SEXP, const char *, SEXP, int);
SEXP R_dense_diag_set(SEXP, SEXP);

SEXP dense_transpose(SEXP, const char *);
SEXP R_dense_transpose(SEXP);

SEXP dense_force_symmetric(SEXP, const char *, char);
SEXP R_dense_force_symmetric(SEXP, SEXP);

SEXP dense_symmpart(SEXP, const char *);
SEXP R_dense_symmpart(SEXP);

SEXP dense_skewpart(SEXP, const char *);
SEXP R_dense_skewpart(SEXP);

int dense_is_symmetric(SEXP, const char *, int);
SEXP R_dense_is_symmetric(SEXP, SEXP);

int dense_is_triangular(SEXP, const char *, int);
SEXP R_dense_is_triangular(SEXP, SEXP);

int dense_is_diagonal(SEXP, const char *);
SEXP R_dense_is_diagonal(SEXP);

SEXP dense_marginsum(SEXP, const char *, int, int, int);
SEXP R_dense_marginsum(SEXP, SEXP, SEXP, SEXP);

SEXP dense_sum(SEXP, const char *, int);
SEXP R_dense_sum(SEXP, SEXP);

SEXP dense_prod(SEXP, const char *, int);
SEXP R_dense_prod(SEXP, SEXP);

#endif /* MATRIX_DENSE_H */
