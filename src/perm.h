#ifndef MATRIX_PERM_H
#define MATRIX_PERM_H

#include <Rinternals.h>

SEXP R_isPerm(SEXP, SEXP);
SEXP R_signPerm(SEXP, SEXP);
SEXP R_invertPerm(SEXP, SEXP, SEXP);
SEXP R_asPerm(SEXP, SEXP, SEXP, SEXP);

#endif /* MATRIX_PERM_H */
