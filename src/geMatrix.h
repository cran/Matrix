#ifndef MATRIX_GEMATRIX_H
#define MATRIX_GEMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP geMatrix_validate(SEXP obj);
SEXP geMatrix_norm(SEXP obj, SEXP norm);
SEXP geMatrix_crossprod(SEXP x);
SEXP geMatrix_geMatrix_crossprod(SEXP x, SEXP y);
SEXP geMatrix_matrix_crossprod(SEXP x, SEXP y);
SEXP geMatrix_getDiag(SEXP x);
SEXP geMatrix_LU(SEXP x);
SEXP geMatrix_determinant(SEXP x, SEXP logarithm);
SEXP geMatrix_solve(SEXP a);
SEXP geMatrix_geMatrix_mm(SEXP a, SEXP b);
SEXP geMatrix_svd(SEXP x, SEXP nu, SEXP nv);

/* DGESDD - compute the singular value decomposition (SVD); of a   */
/* real M-by-N matrix A, optionally computing the left and/or      */
/* right singular vectors.  If singular vectors are desired, it uses a */
/* divide-and-conquer algorithm.                                   */
void F77_NAME(dgesdd)(const char *jobz,
		      const int *m, const int *n,
		      double *a, const int *lda, double *s,
		      double *u, const int *ldu,
		      double *vt, const int *ldvt,
		      double *work, const int *lwork, int *iwork, int *info);


#endif
