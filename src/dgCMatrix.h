#ifndef MATRIX_DGCMATRIX_H
#define MATRIX_DGCMATRIX_H

#include <Rinternals.h>

SEXP dgCMatrix_lusol(SEXP, SEXP);
SEXP dgCMatrix_qrsol(SEXP, SEXP, SEXP);
SEXP dgCMatrix_cholsol(SEXP, SEXP);

#endif /* MATRIX_DGCMATRIX_H */
