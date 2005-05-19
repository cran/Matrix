#ifndef MATRIX_LCHOLCMATRIX_H
#define MATRIX_LCHOLCMATRIX_H

#include <Rdefines.h>
#include "Mutils.h"
#include "triplet_to_col.h"
#include "R_ldl.h"

SEXP lCholCMatrix_validate(SEXP x);
SEXP lCholCMatrix_solve(SEXP x);
SEXP lCholClgCsm(enum CBLAS_SIDE side, enum CBLAS_TRANSPOSE transa, int m,
		 int n, const int Parent[], SEXP BIP, int bp[]);
SEXP lCholCMatrix_lgCMatrix_solve(SEXP a, SEXP b);

#endif
