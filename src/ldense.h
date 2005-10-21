#ifndef MATRIX_LDENSE_H
#define MATRIX_LDENSE_H

#include "Mutils.h"

SEXP lspMatrix_as_lsyMatrix(SEXP from);
SEXP lsyMatrix_as_lspMatrix(SEXP from);
SEXP ltpMatrix_as_ltrMatrix(SEXP from);
SEXP ltrMatrix_as_ltpMatrix(SEXP from);

#endif


