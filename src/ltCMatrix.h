#ifndef MATRIX_LTCMATRIX_H
#define MATRIX_LTCMATRIX_H

#include <Rdefines.h>
#include "Mutils.h"
#include "triplet_to_col.h"

SEXP ltCMatrix_validate(SEXP x);
SEXP ltCMatrix_trans(SEXP x);

#endif
