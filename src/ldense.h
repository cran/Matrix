#ifndef MATRIX_LDENSE_H
#define MATRIX_LDENSE_H

#include "Mutils.h"

/* MJ: no longer needed ... prefer more general (un)?pack() */
#if 0 
SEXP lspMatrix_as_lsyMatrix(SEXP from, SEXP kind);
SEXP lsyMatrix_as_lspMatrix(SEXP from, SEXP kind);
SEXP ltpMatrix_as_ltrMatrix(SEXP from, SEXP kind);
SEXP ltrMatrix_as_ltpMatrix(SEXP from, SEXP kind);
#endif /* MJ */

/* MJ: no longer needed ... prefer more general dense_as_general() */
#if 0
SEXP lsyMatrix_as_lgeMatrix(SEXP from, SEXP kind);
SEXP ltrMatrix_as_lgeMatrix(SEXP from, SEXP kind);
#endif /* MJ */

#endif


