#include "Mutils.h"

SEXP Csparse_Csparse_prod(SEXP a, SEXP b);
SEXP Csparse_crossprod(SEXP x, SEXP trans, SEXP triplet);
SEXP Csparse_dense_crossprod(SEXP a, SEXP b);
SEXP Csparse_dense_prod(SEXP a, SEXP b);
SEXP Csparse_horzcat(SEXP x, SEXP y);
SEXP Csparse_to_Tsparse(SEXP x);
SEXP Csparse_to_dense(SEXP x);
SEXP Csparse_transpose(SEXP x);
SEXP Csparse_validate(SEXP x);
SEXP Csparse_vertcat(SEXP x, SEXP y);
