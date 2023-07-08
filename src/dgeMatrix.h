#ifndef MATRIX_DGEMATRIX_H
#define MATRIX_DGEMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP dgeMatrix_Schur(SEXP x, SEXP vectors, SEXP isDGE);
SEXP dgeMatrix_exp(SEXP x);

/* MJ: unused */
#if 0
SEXP dgeMatrix_svd(SEXP x, SEXP nu, SEXP nv);
#endif /* MJ */

#endif
