#ifndef MATRIX_GEMATRIX_H
#define MATRIX_GEMATRIX_H

#include <R_ext/Boolean.h>
#include "Lapack-etc.h"
#include "Mutils.h"

SEXP dMatrix_validate(SEXP obj);

SEXP dgeMatrix_validate(SEXP obj);
SEXP dgeMatrix_norm(SEXP obj, SEXP norm);
SEXP dgeMatrix_rcond(SEXP obj, SEXP type);
/* for crossprod() and tcrossprod() -- dge*() and the generalized versions: */
SEXP dgeMatrix_crossprod(SEXP x, SEXP trans);
SEXP  geMatrix_crossprod(SEXP x, SEXP trans);
SEXP dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP  geMatrix_geMatrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP  geMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans);
//  %*% :
SEXP dgeMatrix_matrix_mm(SEXP a, SEXP b, SEXP right);
SEXP  geMatrix_matrix_mm(SEXP a, SEXP b, SEXP right);

SEXP dgeMatrix_addDiag(SEXP x, SEXP d);
SEXP dgeMatrix_getDiag(SEXP x);
SEXP lgeMatrix_getDiag(SEXP x);
SEXP dgeMatrix_setDiag(SEXP x, SEXP d);
SEXP lgeMatrix_setDiag(SEXP x, SEXP d);
SEXP dgeMatrix_LU (SEXP x, SEXP warn_singularity);
SEXP dgeMatrix_LU_(SEXP x, Rboolean warn_sing);
SEXP dgeMatrix_determinant(SEXP x, SEXP logarithm);
SEXP dgeMatrix_Schur(SEXP x, SEXP vectors, SEXP isDGE);
SEXP dgeMatrix_solve(SEXP a);
SEXP dgeMatrix_matrix_solve(SEXP a, SEXP b);
SEXP dgeMatrix_svd(SEXP x, SEXP nu, SEXP nv);
SEXP dgeMatrix_exp(SEXP x);
SEXP dgeMatrix_colsums(SEXP x, SEXP naRmP, SEXP cols, SEXP mean);

#endif
