#ifndef MATRIX_TSC_H
#define MATRIX_TSC_H

#include "Mutils.h"
#include "dgCMatrix.h"

SEXP tsc_validate(SEXP x);
SEXP tsc_transpose(SEXP x);
SEXP tsc_to_dgTMatrix(SEXP x);
SEXP Parent_inverse(SEXP par, SEXP unitdiag);
int parent_inv_ap(int n, int countDiag, const int pr[], int ap[]);
void parent_inv_ai(int n, int countDiag, const int pr[], int ai[]);

#endif
