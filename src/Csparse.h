#ifndef MATRIX_CSPARSE_H
#define MATRIX_CSPARSE_H

#include <Rinternals.h>

SEXP CsparseMatrix_validate_maybe_sorting(SEXP);

SEXP dgCMatrix_lusol(SEXP, SEXP);
SEXP dgCMatrix_qrsol(SEXP, SEXP, SEXP);
SEXP dgCMatrix_cholsol(SEXP, SEXP);

SEXP dtCMatrix_diag(SEXP, SEXP);

SEXP Csparse_dmperm(SEXP, SEXP, SEXP);
SEXP Csparse_writeMM(SEXP, SEXP);

#endif /* MATRIX_CSPARSE_H */
