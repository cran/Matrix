#ifndef MATRIX_DGCMATRIX_H
#define MATRIX_DGCMATRIX_H

#include "Mutils.h"

SEXP compressed_non_0_ij(SEXP x, SEXP colP);
SEXP dgCMatrix_lusol(SEXP x, SEXP y);
SEXP dgCMatrix_qrsol(SEXP x, SEXP y, SEXP ord);
SEXP dgCMatrix_cholsol(SEXP x, SEXP y);

#endif
