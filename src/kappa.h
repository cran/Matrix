#ifndef MATRIX_KAPPA_H
#define MATRIX_KAPPA_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP dgeMatrix_norm(SEXP obj, SEXP type);
SEXP dgeMatrix_rcond(SEXP obj, SEXP trf, SEXP type);

SEXP dtrMatrix_norm(SEXP obj, SEXP type);
SEXP dtrMatrix_rcond(SEXP obj, SEXP type);

SEXP dtpMatrix_norm(SEXP obj, SEXP type);
SEXP dtpMatrix_rcond(SEXP obj, SEXP type);

SEXP dsyMatrix_norm(SEXP obj, SEXP type);
SEXP dsyMatrix_rcond(SEXP obj, SEXP trf, SEXP type);

SEXP dspMatrix_norm(SEXP obj, SEXP type);
SEXP dspMatrix_rcond(SEXP obj, SEXP trf, SEXP type);

SEXP dpoMatrix_rcond(SEXP obj, SEXP trf, SEXP type);

SEXP dppMatrix_rcond(SEXP obj, SEXP trf, SEXP type);

#endif /* MATRIX_KAPPA_H */
