#ifndef MATRIX_TSC_H
#define MATRIX_TSC_H

#include "Mutils.h"
#include "cscMatrix.h"

SEXP tsc_validate(SEXP x);
SEXP tsc_transpose(SEXP x);
SEXP tsc_to_triplet(SEXP x);

#endif
