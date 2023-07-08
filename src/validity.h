#ifndef MATRIX_VALIDITY_H
#define MATRIX_VALIDITY_H

#include "Mutils.h"

char* Dim_validate(SEXP dim);
SEXP R_Dim_validate(SEXP dim);

char* DimNames_validate(SEXP dimnames, int pdim[]);
SEXP R_DimNames_validate(SEXP dimnames, SEXP dim);

SEXP R_DimNames_fixup(SEXP dn);

SEXP Matrix_validate(SEXP obj);
SEXP MatrixFactorization_validate(SEXP obj);

SEXP dMatrix_validate(SEXP obj);
SEXP lMatrix_validate(SEXP obj);
SEXP ndenseMatrix_validate(SEXP obj);
SEXP iMatrix_validate(SEXP obj);
SEXP zMatrix_validate(SEXP obj);

SEXP compMatrix_validate(SEXP obj);
SEXP symmetricMatrix_validate(SEXP obj);
SEXP triangularMatrix_validate(SEXP obj);

SEXP diagonalMatrix_validate(SEXP obj);
SEXP indMatrix_validate(SEXP obj);
SEXP pMatrix_validate(SEXP obj);

SEXP CsparseMatrix_validate(SEXP obj);
SEXP RsparseMatrix_validate(SEXP obj);
SEXP TsparseMatrix_validate(SEXP obj);

SEXP sCMatrix_validate(SEXP obj);
SEXP tCMatrix_validate(SEXP obj);
SEXP sRMatrix_validate(SEXP obj);
SEXP tRMatrix_validate(SEXP obj);
SEXP sTMatrix_validate(SEXP obj);
SEXP tTMatrix_validate(SEXP obj);

SEXP xgCMatrix_validate(SEXP obj);
SEXP xsCMatrix_validate(SEXP obj);
SEXP xtCMatrix_validate(SEXP obj);
SEXP xgRMatrix_validate(SEXP obj);
SEXP xsRMatrix_validate(SEXP obj);
SEXP xtRMatrix_validate(SEXP obj);
SEXP xgTMatrix_validate(SEXP obj);
SEXP xsTMatrix_validate(SEXP obj);
SEXP xtTMatrix_validate(SEXP obj);

SEXP unpackedMatrix_validate(SEXP obj);
SEXP packedMatrix_validate(SEXP obj);

SEXP dpoMatrix_validate(SEXP obj);
SEXP dppMatrix_validate(SEXP obj);
SEXP corMatrix_validate(SEXP obj);
SEXP pcorMatrix_validate(SEXP obj);

SEXP denseLU_validate(SEXP obj);
SEXP sparseLU_validate(SEXP obj);
SEXP sparseQR_validate(SEXP obj);
SEXP BunchKaufman_validate(SEXP obj);
SEXP pBunchKaufman_validate(SEXP obj);
SEXP Cholesky_validate(SEXP obj);
SEXP pCholesky_validate(SEXP obj);
SEXP CHMfactor_validate(SEXP obj);
SEXP CHMsimpl_validate(SEXP obj);
SEXP CHMsuper_validate(SEXP obj);
SEXP dCHMsimpl_validate(SEXP obj);
SEXP dCHMsuper_validate(SEXP obj);
SEXP Schur_validate(SEXP obj);

void validObject(SEXP obj, const char* cl);

#define SNPRINTF(_BUFFER_, _FORMAT_, ...) \
do { \
	_BUFFER_ = R_alloc(Matrix_ErrorBufferSize, sizeof(char)); \
	snprintf(_BUFFER_, Matrix_ErrorBufferSize, _FORMAT_, __VA_ARGS__); \
} while (0)

#define UPRET(_N_, _S_) \
do { \
	UNPROTECT(_N_); \
	return mkString(_(_S_)); \
} while (0)

#define FRUPRET(_PTR_, _M_, _N_, _S_) \
do { \
	Matrix_Free(_PTR_, _M_); \
	UNPROTECT(_N_); \
	return mkString(_(_S_)); \
} while (0)

#endif /* MATRIX_VALIDITY_H */
