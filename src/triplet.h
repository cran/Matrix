#ifndef MATRIX_TRIPLET_H
#define MATRIX_TRIPLET_H

#include "Mutils.h"

SEXP triplet_validate(SEXP x);
SEXP triplet_to_geMatrix(SEXP x);
SEXP triplet_to_matrix(SEXP x);

#endif
