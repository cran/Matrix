#ifndef MATRIX_DGCMATRIX_H
#define MATRIX_DGCMATRIX_H

#include <R_ext/BLAS.h>
#include "Mutils.h"
#include "cs_utils.h"

SEXP xCMatrix_validate(SEXP x);
SEXP compressed_to_TMatrix(SEXP x, SEXP colP);
SEXP compressed_non_0_ij(SEXP x, SEXP colP);
SEXP R_to_CMatrix(SEXP x);
SEXP dgCMatrix_lusol(SEXP x, SEXP y);
SEXP dgCMatrix_qrsol(SEXP x, SEXP y);
SEXP dgCMatrix_QR(SEXP Ap, SEXP order);
SEXP dgCMatrix_LU(SEXP Ap, SEXP orderp, SEXP tolp);
SEXP dgCMatrix_matrix_solve(SEXP Ap, SEXP bp);

#endif
