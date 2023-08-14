#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include "Mutils.h"

SEXP dense_band(SEXP from, const char *class, int a, int b, int new);
SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2);

SEXP R_dense_colSums(SEXP obj, SEXP narm, SEXP mean);
SEXP R_dense_rowSums(SEXP obj, SEXP narm, SEXP mean);

#endif /* MATRIX_DENSE_H */
