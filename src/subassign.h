#ifndef MATRIX_SUBASSIGN_H
#define MATRIX_SUBASSIGN_H

#include <Rinternals.h>

SEXP nCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP lCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP iCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP dCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP zCsparse_subassign(SEXP, SEXP, SEXP, SEXP);

#endif /* MATRIX_SUBASSIGN_H */
