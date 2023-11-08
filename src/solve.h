#ifndef MATRIX_SOLVE_H
#define MATRIX_SOLVE_H

#include <Rinternals.h>

SEXP      denseLU_solve(SEXP, SEXP);
SEXP BunchKaufman_solve(SEXP, SEXP);
SEXP     Cholesky_solve(SEXP, SEXP);
SEXP    dtrMatrix_solve(SEXP, SEXP);
SEXP     sparseLU_solve(SEXP, SEXP, SEXP);
#if 0
/* MJ: we use 'sparseQR_matmult' instead */
SEXP     sparseQR_solve(SEXP, SEXP, SEXP);
#endif
SEXP    CHMfactor_solve(SEXP, SEXP, SEXP, SEXP);
SEXP    dtCMatrix_solve(SEXP, SEXP, SEXP);

SEXP sparseQR_matmult(SEXP, SEXP, SEXP, SEXP, SEXP);

#endif /* MATRIX_SOLVE_H */
