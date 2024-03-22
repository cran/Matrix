#ifndef MATRIX_FACTOR_H
#define MATRIX_FACTOR_H

#include <Rinternals.h>

SEXP dgeMatrix_trf(SEXP, SEXP);
SEXP dsyMatrix_trf(SEXP, SEXP);
SEXP dspMatrix_trf(SEXP, SEXP);
SEXP dpoMatrix_trf(SEXP, SEXP, SEXP, SEXP);
SEXP dppMatrix_trf(SEXP, SEXP);
SEXP dgeMatrix_sch(SEXP, SEXP, SEXP);

SEXP dgCMatrix_trf(SEXP, SEXP, SEXP, SEXP);
SEXP dgCMatrix_orf(SEXP, SEXP, SEXP);
SEXP dpCMatrix_trf(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP BunchKaufman_expand(SEXP, SEXP);

SEXP CHMfactor_diag_get(SEXP, SEXP);
SEXP CHMfactor_update(SEXP, SEXP, SEXP);
SEXP CHMfactor_updown(SEXP, SEXP, SEXP);

#endif /* MATRIX_FACTOR_H */
