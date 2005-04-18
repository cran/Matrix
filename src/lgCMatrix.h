#ifndef MATRIX_LGCMATRIX_H
#define MATRIX_LCGMATRIX_H

#include <Rdefines.h>
#include "Mutils.h"
#include "triplet_to_col.h"

SEXP lgCMatrix_validate(SEXP x);
void Matrix_lgClgCmm(int tra, int trb, int m, int n, int k,
		     const int ai[], const int ap[],
		     const int bi[], const int bp[],
		     int beta, SEXP *CIP, int cp[]);
SEXP lgCMatrix_lgCMatrix_mm(SEXP a, SEXP b, SEXP transa, SEXP transb);

#endif
