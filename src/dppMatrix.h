#ifndef MATRIX_PPMATRIX_H
#define MATRIX_PPMATRIX_H

#include "dspMatrix.h"

SEXP dppMatrix_trf_(SEXP obj,  int warn);
SEXP dppMatrix_trf (SEXP obj, SEXP warn);

SEXP dppMatrix_rcond(SEXP obj);
SEXP dppMatrix_solve(SEXP a);
SEXP dppMatrix_matrix_solve(SEXP a, SEXP b);

#endif
