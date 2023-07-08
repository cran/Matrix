#ifndef MATRIX_SPARSEVECTOR_H
#define MATRIX_SPARSEVECTOR_H

#include "Mutils.h"

/* defined in ./sparse.c : */
SEXP R_sparse_as_general(SEXP);

SEXP v2spV(SEXP from);
SEXP CR2spV(SEXP from);

#endif /* MATRIX_SPARSEVECTOR_H */
