#ifndef MATRIX_SSC_H
#define MATRIX_SSC_H

#include "taucs_utils.h"
#include "Metis_utils.h"

SEXP sscMatrix_validate(SEXP x);
SEXP sscMatrix_chol(SEXP x, SEXP pivot);
SEXP sscMatrix_inverse_factor(SEXP A);
SEXP sscMatrix_matrix_solve(SEXP a, SEXP b);
SEXP ssc_transpose(SEXP x);

#endif
