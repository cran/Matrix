#ifndef MATRIX_TPMATRIX_H
#define MATRIX_TPMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

double get_norm_dtp(SEXP obj, const char *typstr);
SEXP dtpMatrix_norm(SEXP obj, SEXP type);
SEXP dtpMatrix_rcond(SEXP obj, SEXP type);
SEXP dtpMatrix_solve(SEXP a);
SEXP dtpMatrix_matrix_solve(SEXP a, SEXP b);

SEXP dtpMatrix_matrix_mm(SEXP x, SEXP y, SEXP right, SEXP trans);
SEXP dgeMatrix_dtpMatrix_mm(SEXP x, SEXP y);

/* MJ: no longer needed ... prefer more general packedMatrix_diag_[gs]et() */
#if 0
SEXP dtpMatrix_getDiag(SEXP x);
SEXP ltpMatrix_getDiag(SEXP x);
SEXP dtpMatrix_setDiag(SEXP x, SEXP d);
SEXP ltpMatrix_setDiag(SEXP x, SEXP d);
/* was unused, not replaced: */
SEXP dtpMatrix_addDiag(SEXP x, SEXP d);
#endif

/* MJ: no longer needed ... prefer more general packedMatrix_unpack() */
#if 0
SEXP dtpMatrix_as_dtrMatrix(SEXP from);
#endif

#endif
