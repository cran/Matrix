#ifndef MATRIX_TRMATRIX_H
#define MATRIX_TRMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

double get_norm_dtr(SEXP obj, const char *typstr);
SEXP dtrMatrix_norm(SEXP obj, SEXP type);
SEXP dtrMatrix_rcond(SEXP obj, SEXP type);
SEXP dtrMatrix_solve(SEXP a);
SEXP dtrMatrix_matrix_solve(SEXP a, SEXP b);

SEXP dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans);
SEXP dtrMatrix_matrix_mm   (SEXP a, SEXP b, SEXP right, SEXP trans);

SEXP dtrMatrix_chol2inv(SEXP a);

SEXP dtrMatrix_addDiag(SEXP x, SEXP d);

/* MJ: no longer needed ... prefer more general unpackedMatrix_diag_[gs]et() */
#if 0
SEXP dtrMatrix_getDiag(SEXP x);
SEXP ltrMatrix_getDiag(SEXP x);
SEXP dtrMatrix_setDiag(SEXP x, SEXP d);
SEXP ltrMatrix_setDiag(SEXP x, SEXP d);
#endif /* MJ */

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0
SEXP dtrMatrix_as_dtpMatrix(SEXP from);
#endif

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0
SEXP dtrMatrix_as_matrix(SEXP from, SEXP keep_dimnames);
#endif

#endif
