#ifndef MATRIX_TRS_H
#define MATRIX_TRS_H

#include "Mutils.h"

SEXP dsTMatrix_validate(SEXP x);
SEXP dsTMatrix_as_dsyMatrix(SEXP x);
SEXP dsTMatrix_as_dgTMatrix(SEXP x);

#endif
