#ifndef MATRIX_VALIDITY_H
#define MATRIX_VALIDITY_H

#include <Rinternals.h>

char* Dim_validate(SEXP);
SEXP R_Dim_validate(SEXP);

char* DimNames_validate(SEXP, int[]);
SEXP R_DimNames_validate(SEXP, SEXP);

SEXP R_DimNames_fixup(SEXP);

SEXP Matrix_validate(SEXP);
SEXP MatrixFactorization_validate(SEXP);

SEXP nMatrix_validate(SEXP);
SEXP lMatrix_validate(SEXP);
SEXP iMatrix_validate(SEXP);
SEXP dMatrix_validate(SEXP);
SEXP zMatrix_validate(SEXP);

SEXP compMatrix_validate(SEXP);
SEXP symmetricMatrix_validate(SEXP);
SEXP triangularMatrix_validate(SEXP);

SEXP diagonalMatrix_validate(SEXP);
SEXP indMatrix_validate(SEXP);
SEXP pMatrix_validate(SEXP);

SEXP CsparseMatrix_validate(SEXP);
SEXP RsparseMatrix_validate(SEXP);
SEXP TsparseMatrix_validate(SEXP);

SEXP sCMatrix_validate(SEXP);
SEXP tCMatrix_validate(SEXP);
SEXP sRMatrix_validate(SEXP);
SEXP tRMatrix_validate(SEXP);
SEXP sTMatrix_validate(SEXP);
SEXP tTMatrix_validate(SEXP);

SEXP xgCMatrix_validate(SEXP);
SEXP xsCMatrix_validate(SEXP);
SEXP xtCMatrix_validate(SEXP);
SEXP xgRMatrix_validate(SEXP);
SEXP xsRMatrix_validate(SEXP);
SEXP xtRMatrix_validate(SEXP);
SEXP xgTMatrix_validate(SEXP);
SEXP xsTMatrix_validate(SEXP);
SEXP xtTMatrix_validate(SEXP);

SEXP unpackedMatrix_validate(SEXP);
SEXP packedMatrix_validate(SEXP);

SEXP dpoMatrix_validate(SEXP);
SEXP dppMatrix_validate(SEXP);
SEXP corMatrix_validate(SEXP);
SEXP pcorMatrix_validate(SEXP);

SEXP sparseVector_validate(SEXP);
SEXP lsparseVector_validate(SEXP);
SEXP isparseVector_validate(SEXP);
SEXP dsparseVector_validate(SEXP);
SEXP zsparseVector_validate(SEXP);

SEXP denseLU_validate(SEXP);
SEXP sparseLU_validate(SEXP);
SEXP sparseQR_validate(SEXP);
SEXP BunchKaufman_validate(SEXP);
SEXP pBunchKaufman_validate(SEXP);
SEXP Cholesky_validate(SEXP);
SEXP pCholesky_validate(SEXP);
SEXP CHMfactor_validate(SEXP);
SEXP CHMsimpl_validate(SEXP);
SEXP CHMsuper_validate(SEXP);
SEXP dCHMsimpl_validate(SEXP);
SEXP dCHMsuper_validate(SEXP);
SEXP Schur_validate(SEXP);

void validObject(SEXP, const char *);

#endif /* MATRIX_VALIDITY_H */
