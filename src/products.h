#ifndef MATRIX_PRODUCTS_H
#define MATRIX_PRODUCTS_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP dgeMatrix_crossprod(SEXP x, SEXP trans);
SEXP  geMatrix_crossprod(SEXP x, SEXP trans);
SEXP dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP   geMatrix_geMatrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans);
SEXP  geMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans);

SEXP dgeMatrix_matrix_mm(SEXP a, SEXP b, SEXP right);
SEXP  geMatrix_matrix_mm(SEXP a, SEXP b, SEXP right);

SEXP dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans);
SEXP dtrMatrix_matrix_mm   (SEXP a, SEXP b, SEXP right, SEXP trans);

SEXP dtpMatrix_matrix_mm(SEXP x, SEXP y, SEXP right, SEXP trans);
SEXP dgeMatrix_dtpMatrix_mm(SEXP x, SEXP y);

SEXP dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP right);

SEXP dspMatrix_matrix_mm(SEXP a, SEXP b);

SEXP Csp_dense_products(SEXP a, SEXP b,
                        Rboolean trans_a,
                        Rboolean trans_b,
                        Rboolean trans_ans);

SEXP Csparse_Csparse_prod(SEXP a, SEXP b, SEXP boolArith);
SEXP Csparse_Csparse_crossprod(SEXP a, SEXP b, SEXP trans, SEXP boolArith);
SEXP Csparse_crossprod(SEXP x, SEXP trans, SEXP boolArith);

SEXP Csparse_dense_prod     (SEXP a, SEXP b, SEXP trans);
SEXP Csparse_dense_crossprod(SEXP a, SEXP b, SEXP trans);

#endif
