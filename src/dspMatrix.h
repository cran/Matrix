#ifndef MATRIX_SPMATRIX_H
#define MATRIX_SPMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

double get_norm_sp(SEXP obj, const char *typstr);
SEXP dspMatrix_norm(SEXP obj, SEXP type);
SEXP dspMatrix_rcond(SEXP obj, SEXP type);
SEXP dspMatrix_solve(SEXP a);
SEXP dspMatrix_matrix_solve(SEXP a, SEXP b);
SEXP dspMatrix_matrix_mm(SEXP a, SEXP b);
SEXP dspMatrix_trf(SEXP x);

/* MJ: no longer needed ... prefer more general packedMatrix_diag_[gs]et() */
#if 0
SEXP dspMatrix_getDiag(SEXP x);
SEXP lspMatrix_getDiag(SEXP x);
SEXP dspMatrix_setDiag(SEXP x, SEXP d);
SEXP lspMatrix_setDiag(SEXP x, SEXP d);
#endif

/* MJ: no longer needed ... prefer more general packedMatrix_unpack() */
#if 0
SEXP dspMatrix_as_dsyMatrix(SEXP from);
#endif
/* MJ */

#endif
