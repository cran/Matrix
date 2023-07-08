#ifndef MATRIX_FACTORIZATIONS_H
#define MATRIX_FACTORIZATIONS_H

#include "cs.h"
#include "chm_common.h"
#include "Lapack-etc.h"
#include "Mutils.h"

SEXP dgeMatrix_trf_(SEXP obj,  int warn);
SEXP dsyMatrix_trf_(SEXP obj,  int warn);
SEXP dspMatrix_trf_(SEXP obj,  int warn);
SEXP dpoMatrix_trf_(SEXP obj,  int warn,  int pivot, double tol);
SEXP dppMatrix_trf_(SEXP obj,  int warn);

SEXP dgeMatrix_trf (SEXP obj, SEXP warn);
SEXP dsyMatrix_trf (SEXP obj, SEXP warn);
SEXP dspMatrix_trf (SEXP obj, SEXP warn);
SEXP dpoMatrix_trf (SEXP obj, SEXP warn, SEXP pivot,   SEXP tol);
SEXP dppMatrix_trf (SEXP obj, SEXP warn);

int dgCMatrix_trf_(const cs *A, css **S, csn **N, int order, double tol);
SEXP dgCMatrix_trf(SEXP obj, SEXP order, SEXP tol, SEXP doError);

int dgCMatrix_orf_(const cs *A, css **S, csn **N, int order);
SEXP dgCMatrix_orf(SEXP obj, SEXP order, SEXP doError);

int dpCMatrix_trf_(cholmod_sparse *A, cholmod_factor **L,
                   int perm, int ldl, int super, double mult);
SEXP dpCMatrix_trf(SEXP obj,
                   SEXP perm, SEXP ldl, SEXP super, SEXP mult);

SEXP BunchKaufman_expand(SEXP obj, SEXP packed);

SEXP      denseLU_determinant(SEXP obj, SEXP logarithm);
SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm, SEXP packed);
SEXP     Cholesky_determinant(SEXP obj, SEXP logarithm, SEXP packed);
SEXP     sparseLU_determinant(SEXP obj, SEXP logarithm);
SEXP     sparseQR_determinant(SEXP obj, SEXP logarithm);
SEXP    CHMfactor_determinant(SEXP obj, SEXP logarithm, SEXP sqrt);

SEXP      denseLU_solve(SEXP a, SEXP b);
SEXP BunchKaufman_solve(SEXP a, SEXP b, SEXP packed);
SEXP     Cholesky_solve(SEXP a, SEXP b, SEXP packed);
SEXP     sparseLU_solve(SEXP a, SEXP b, SEXP sparse);
/* MJ: not needed since we have 'sparseQR_matmult' : */
#if 0
SEXP     sparseQR_solve(SEXP a, SEXP b, SEXP sparse);
#endif /* MJ */
SEXP    CHMfactor_solve(SEXP a, SEXP b, SEXP sparse, SEXP system);
SEXP    dtrMatrix_solve(SEXP a, SEXP b, SEXP packed);
SEXP    dtCMatrix_solve(SEXP a, SEXP b, SEXP sparse);

SEXP sparseQR_matmult(SEXP qr, SEXP y, SEXP op, SEXP complete, SEXP yxjj);

SEXP CHMfactor_diag_get(SEXP obj, SEXP square);
SEXP CHMfactor_update(SEXP obj, SEXP parent, SEXP mult);
SEXP CHMfactor_updown(SEXP obj, SEXP parent, SEXP update);

#endif /* MATRIX_FACTORIZATIONS_H */
