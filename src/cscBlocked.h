#ifndef MATRIX_CSCBLOCKED_H
#define MATRIX_CSCBLOCKED_H

#include "Mutils.h"
#include "R_ldl.h"
#include "triplet_to_col.h"
#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP cscBlocked_validate(SEXP x);
void cscb_tri(enum CBLAS_UPLO upper, enum CBLAS_DIAG unit,
	      SEXP A, const int Parent[], SEXP AI);
void cscb_syrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans,
	       double alpha, SEXP A,
	       double beta, SEXP C);
int cscb_ldl(SEXP A, const int Parent[], SEXP L, SEXP D);
void cscb_trmm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
	       enum CBLAS_TRANSPOSE transa, enum CBLAS_DIAG diag,
	       double alpha, SEXP A, double B[], int m, int n, int ldb);
void cscb_trsm(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE transa, enum CBLAS_DIAG diag,
	       double alpha, SEXP A, double B[], int m, int n, int ldb);
void cscb_trcbm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
		enum CBLAS_TRANSPOSE transa, enum CBLAS_DIAG diag,
		double alpha, SEXP A, SEXP B);
void cscb_trcbsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
		 enum CBLAS_TRANSPOSE transa, enum CBLAS_DIAG diag,
		 double alpha, SEXP A, const int Parent[], SEXP B);
void cscb_cscbm(enum CBLAS_TRANSPOSE transa, enum CBLAS_TRANSPOSE transb,
		double alpha, SEXP A, SEXP B, double beta, SEXP C);
void cscb_mm(enum CBLAS_SIDE side, enum CBLAS_TRANSPOSE transa,
	     int m, int n, int k, double alpha, SEXP A,
	     const double B[], int ldb, double beta, double C[], int ldc);
#endif
