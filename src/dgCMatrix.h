#ifndef MATRIX_DGCMATRIX_H
#define MATRIX_DGCMATRIX_H

#include <R_ext/BLAS.h>
#include "Mutils.h"
#include "cs_utils.h"

SEXP xCMatrix_validate(SEXP x);
SEXP xRMatrix_validate(SEXP x);
SEXP compressed_non_0_ij(SEXP x, SEXP colP);

/* MJ: no longer needed ... prefer CRsparse_as_Tsparse() */
#if 0
SEXP compressed_to_TMatrix(SEXP x, SEXP colP);
#endif /* MJ */

/* MJ: no longer needed ... 
   now done via R_sparse_transpose(), tCRsparse_as_RCsparse() */
#if 0
SEXP R_to_CMatrix(SEXP x);
#endif /* MJ */

SEXP dgCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means);
SEXP igCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means);
SEXP lgCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means);
SEXP ngCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means);

/* SEXP dgCMatrix_lusol(SEXP x, SEXP y); */
SEXP dgCMatrix_qrsol(SEXP x, SEXP y, SEXP ord);
SEXP dgCMatrix_cholsol(SEXP x, SEXP y);
SEXP dgCMatrix_QR(SEXP Ap, SEXP order, SEXP keep_dimnames);

#ifdef Matrix_with_SPQR
SEXP dgCMatrix_SPQR(SEXP Ap, SEXP ordering, SEXP econ, SEXP tol);
#endif

SEXP dgCMatrix_LU(SEXP Ap, SEXP orderp, SEXP tolp, SEXP error_on_sing,
		  SEXP keep_dimnames);
SEXP dgCMatrix_matrix_solve(SEXP Ap, SEXP bp, SEXP give_sparse);

#endif
