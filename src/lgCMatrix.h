#ifndef MATRIX_LGCMATRIX_H
#define MATRIX_LGCMATRIX_H

#include <Rdefines.h>
#include "Mutils.h"
#include "triplet_to_col.h"

SEXP lgCMatrix_validate(SEXP x);

SEXP lcsc_to_matrix(SEXP x);

SEXP Matrix_lgClgCmm(int tra, int trb, int m, int n, int k,
		     const int ai[], const int ap[],
		     const int bi[], const int bp[],
		     int beta, SEXP CIP, int cp[]);
SEXP lgCMatrix_lgCMatrix_mm(SEXP a, SEXP b);
SEXP lgCMatrix_trans(SEXP x);
SEXP Matrix_lgCsyrk(int up, int tra, int n, int k, const int ai[],
		    const int ap[], int beta, SEXP CIP, int cp[]);
SEXP lgCMatrix_crossprod(SEXP x, SEXP trans, SEXP C);
SEXP lgCMatrix_picky_column(SEXP x);

#endif
