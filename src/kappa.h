#ifndef MATRIX_KAPPA_H
#define MATRIX_KAPPA_H

#include <Rinternals.h>

SEXP dgeMatrix_norm(SEXP, SEXP);
SEXP dsyMatrix_norm(SEXP, SEXP);
SEXP dspMatrix_norm(SEXP, SEXP);
SEXP dtrMatrix_norm(SEXP, SEXP);
SEXP dtpMatrix_norm(SEXP, SEXP);

SEXP dgeMatrix_rcond(SEXP, SEXP, SEXP);
SEXP dsyMatrix_rcond(SEXP, SEXP, SEXP);
SEXP dspMatrix_rcond(SEXP, SEXP, SEXP);
SEXP dpoMatrix_rcond(SEXP, SEXP, SEXP);
SEXP dppMatrix_rcond(SEXP, SEXP, SEXP);
SEXP dtrMatrix_rcond(SEXP, SEXP);
SEXP dtpMatrix_rcond(SEXP, SEXP);

#endif /* MATRIX_KAPPA_H */
