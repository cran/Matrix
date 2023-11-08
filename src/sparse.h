#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include <Rinternals.h>

SEXP sparse_drop0(SEXP, const char *, double);
SEXP R_sparse_drop0(SEXP, SEXP);

SEXP sparse_diag_U2N(SEXP, const char *);
SEXP R_sparse_diag_U2N(SEXP);

SEXP sparse_diag_N2U(SEXP, const char *);
SEXP R_sparse_diag_N2U(SEXP);

SEXP sparse_band(SEXP, const char *, int, int);
SEXP R_sparse_band(SEXP, SEXP, SEXP);

SEXP sparse_diag_get(SEXP, const char *, int);
SEXP R_sparse_diag_get(SEXP, SEXP);

SEXP sparse_diag_set(SEXP, const char *, SEXP);
SEXP R_sparse_diag_set(SEXP, SEXP);

SEXP sparse_transpose(SEXP, const char *, int);
SEXP R_sparse_transpose(SEXP, SEXP);

SEXP sparse_force_symmetric(SEXP, const char *, char);
SEXP R_sparse_force_symmetric(SEXP, SEXP);

SEXP sparse_symmpart(SEXP, const char *);
SEXP R_sparse_symmpart(SEXP);

SEXP sparse_skewpart(SEXP, const char *);
SEXP R_sparse_skewpart(SEXP);

int sparse_is_symmetric(SEXP, const char *, int);
SEXP R_sparse_is_symmetric(SEXP, SEXP);

int sparse_is_triangular(SEXP, const char *, int);
SEXP R_sparse_is_triangular(SEXP, SEXP);

int sparse_is_diagonal(SEXP, const char *);
SEXP R_sparse_is_diagonal(SEXP);

SEXP sparse_marginsum(SEXP, const char *, int, int, int, int);
SEXP R_sparse_marginsum(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP sparse_sum(SEXP, const char *, int);
SEXP R_sparse_sum(SEXP, SEXP);

SEXP sparse_prod(SEXP, const char *, int);
SEXP R_sparse_prod(SEXP, SEXP);

SEXP Tsparse_aggregate(SEXP);

#endif /* MATRIX_SPARSE_H */
