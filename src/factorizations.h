#ifndef MATRIX_FACTORS_H
#define MATRIX_FACTORS_H

#include "Mutils.h"

SEXP denseLU_expand(SEXP obj);
SEXP BunchKaufman_expand(SEXP obj);

SEXP denseLU_determinant(SEXP obj, SEXP logarithm);
SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm);

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0
SEXP LU_validate(SEXP obj);
#endif /* MJ */

/* MJ: no longer needed ... prefer denseLU_expand() */
#if 0
SEXP LU_expand(SEXP x);
#endif /* MJ */

#endif
