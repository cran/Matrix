#ifndef MATRIX_LGCMATRIX_H
#define MATRIX_LGCMATRIX_H

#include "Mutils.h"

SEXP lgCMatrix_validate(SEXP x);

SEXP lcsc_to_matrix(SEXP x);
SEXP ncsc_to_matrix(SEXP x);

#endif
