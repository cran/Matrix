#ifndef MATRIX_CSPARSE_H
#define MATRIX_CSPARSE_H

#include <Rinternals.h>

SEXP CsparseMatrix_validate_maybe_sorting(SEXP);

SEXP tCsparse_diag(SEXP, SEXP);

SEXP nCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP lCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP iCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP dCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP zCsparse_subassign(SEXP, SEXP, SEXP, SEXP);

SEXP Csparse_dmperm(SEXP, SEXP, SEXP);

SEXP Csparse_MatrixMarket(SEXP, SEXP);

#endif /* MATRIX_CSPARSE_H */
