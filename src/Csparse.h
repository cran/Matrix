#ifndef MATRIX_CSPARSE_H
#define MATRIX_CSPARSE_H

#include "Mutils.h"

SEXP Csparse_validate_(SEXP x, Rboolean maybe_modify);
SEXP Csparse_validate2(SEXP x, SEXP maybe_modify);
SEXP Csparse_sort(SEXP x);
SEXP Csparse_horzcat(SEXP x, SEXP y);
SEXP Csparse_vertcat(SEXP x, SEXP y);

SEXP dCsparse_subassign(SEXP x, SEXP i_, SEXP j_, SEXP value);
SEXP lCsparse_subassign(SEXP x, SEXP i_, SEXP j_, SEXP value);
SEXP iCsparse_subassign(SEXP x, SEXP i_, SEXP j_, SEXP value);
SEXP nCsparse_subassign(SEXP x, SEXP i_, SEXP j_, SEXP value);
SEXP zCsparse_subassign(SEXP x, SEXP i_, SEXP j_, SEXP value);

SEXP Csparse_MatrixMarket(SEXP x, SEXP fname);
SEXP Csparse_dmperm(SEXP mat, SEXP seed, SEXP nAns);

SEXP diag_tC(SEXP obj, SEXP resultKind);

#endif
