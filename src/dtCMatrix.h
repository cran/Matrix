#ifndef MATRIX_TSC_H
#define MATRIX_TSC_H

#include "Mutils.h"
#include "dgCMatrix.h"

SEXP dtCMatrix_validate(SEXP x);
SEXP Parent_inverse(SEXP par, SEXP unitdiag);
int parent_inv_ap(int n, int countDiag, const int pr[], int ap[]);
void parent_inv_ai(int n, int countDiag, const int pr[], int ai[]);
SEXP dtCMatrix_solve(SEXP a);
SEXP dtCMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed);
SEXP dtCMatrix_upper_solve(SEXP a);

#endif
