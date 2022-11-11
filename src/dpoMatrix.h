#ifndef MATRIX_POMATRIX_H
#define MATRIX_POMATRIX_H

#include "dsyMatrix.h"

SEXP dpoMatrix_trf_(SEXP obj,  int warn);
SEXP dpoMatrix_trf (SEXP obj, SEXP warn);

SEXP dpoMatrix_rcond(SEXP obj);
SEXP dpoMatrix_solve(SEXP a);
SEXP dpoMatrix_matrix_solve(SEXP a, SEXP b);

#endif
