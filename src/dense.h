#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include "Rdefines.h"
#include "R_ext/Lapack.h"
/* #include "flame.h" */

SEXP lsq_dense_Chol(SEXP X, SEXP y);
SEXP lsq_dense_QR(SEXP X, SEXP y);
SEXP lapack_qr(SEXP Xin, SEXP tl);

#endif
