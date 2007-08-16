#ifndef MATRIX_CSPARSE_H
#define MATRIX_CSPARSE_H

#include "Mutils.h"

SEXP Csparse_Csparse_prod(SEXP a, SEXP b);
SEXP Csparse_band(SEXP x, SEXP k1, SEXP k2);
SEXP Csparse_crossprod(SEXP x, SEXP trans, SEXP triplet);
SEXP Csparse_Csparse_crossprod(SEXP a, SEXP b, SEXP trans);
SEXP Csparse_dense_crossprod(SEXP a, SEXP b);
SEXP Csparse_dense_prod(SEXP a, SEXP b);
SEXP Csparse_diagU2N(SEXP x);
SEXP Csparse_drop(SEXP x, SEXP tol);
SEXP Csparse_horzcat(SEXP x, SEXP y);
SEXP Csparse_submatrix(SEXP x, SEXP i, SEXP j);
SEXP Csparse_symmetric_to_general(SEXP x);
SEXP Csparse_general_to_symmetric(SEXP x, SEXP uplo);
SEXP Csparse_MatrixMarket(SEXP x, SEXP fname);
SEXP Csparse_to_Tsparse(SEXP x, SEXP tri);
SEXP Csparse_to_dense(SEXP x);
SEXP Csparse_to_nz_pattern(SEXP x, SEXP tri);
SEXP Csparse_to_matrix(SEXP x);
SEXP Csparse_transpose(SEXP x, SEXP tri);
SEXP Csparse_validate(SEXP x);
SEXP Csparse_vertcat(SEXP x, SEXP y);

SEXP Rsparse_validate(SEXP x);

#endif
