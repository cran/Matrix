#ifndef MATRIX_CSCBLOCKED_H
#define MATRIX_CSCBLOCKED_H

#include "Mutils.h"
#include <R_ext/Lapack.h>

SEXP cscBlocked_validate(SEXP x);
void cscBlocked_mm(char side, char transa, int m, int n, int k,
		   double alpha, int nr, int nc,
		   const int ap[], const int ai[],
		   const double ax[],
		   const double b[], int ldb,
		   double beta, double c[], int ldc);
void cscBlocked_tri(char upper, char unit, int n, int nr, int nc,
		    const int ap[], const int ai[], const double ax[],
		    int aip[], int aii[], double aix[]);

#endif
