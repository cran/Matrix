#include "LU.h"
#include "Mutils.h"
#include "cscBlocked.h"
#include "cscMatrix.h"
#include "dense.h"
#include "factorizations.h"
#include "geMatrix.h"
#include "lmeRep.h"
#include "poMatrix.h"
#include "sscCrosstab.h"
#include "sscMatrix.h"
#include "ssclme.h"
#include "syMatrix.h"
#include "trMatrix.h"
#include "triplet.h"
#include "tscMatrix.h"
#include <R_ext/Rdynload.h>

#include "Syms.h"

static R_CallMethodDef CallEntries[] = {
    {"csc_check_column_sorting", (DL_FUNC) &csc_check_column_sorting, 1},
    {"nlme_replaceSlot", (DL_FUNC) &nlme_replaceSlot, 3},
    {"nlme_weight_matrix_list", (DL_FUNC) &nlme_weight_matrix_list, 4},
    {"LU_expand", (DL_FUNC) &LU_expand, 1},
    {"cscBlocked_validate", (DL_FUNC) &cscBlocked_validate, 1},
    {"csc_crossprod", (DL_FUNC) &csc_crossprod, 1},
    {"csc_matrix_crossprod", (DL_FUNC) &csc_matrix_crossprod, 2},
    {"csc_validate", (DL_FUNC) &csc_validate, 1},
    {"csc_to_triplet", (DL_FUNC) &csc_to_triplet, 1},
    {"csc_to_matrix", (DL_FUNC) &csc_to_matrix, 1},
    {"csc_to_geMatrix", (DL_FUNC) &csc_to_geMatrix, 1},
/*     {"csc_to_imagemat", (DL_FUNC) &csc_to_imagemat, 1}, */
    {"matrix_to_csc", (DL_FUNC) &matrix_to_csc, 1},
    {"triplet_to_csc", (DL_FUNC) &triplet_to_csc, 1},
    {"csc_getDiag", (DL_FUNC) &csc_getDiag, 1},
    {"csc_transpose", (DL_FUNC) &csc_transpose, 1},
    {"csc_matrix_mm", (DL_FUNC) &csc_matrix_mm, 2},
    {"lsq_dense_Chol", (DL_FUNC) &lsq_dense_Chol, 2},
    {"lsq_dense_QR", (DL_FUNC) &lsq_dense_QR, 2},
    {"lapack_qr", (DL_FUNC) &lapack_qr, 2},
    {"LU_validate", (DL_FUNC) &LU_validate, 1},
    {"Cholesky_validate", (DL_FUNC) &Cholesky_validate, 1},
    {"SVD_validate", (DL_FUNC) &SVD_validate, 1},
    {"geMatrix_validate", (DL_FUNC) &geMatrix_validate, 1},
    {"geMatrix_norm", (DL_FUNC) &geMatrix_norm, 2},
    {"geMatrix_crossprod", (DL_FUNC) &geMatrix_crossprod, 1},
    {"geMatrix_geMatrix_crossprod", (DL_FUNC) &geMatrix_geMatrix_crossprod, 2},
    {"geMatrix_matrix_crossprod", (DL_FUNC) &geMatrix_matrix_crossprod, 2},
    {"geMatrix_getDiag", (DL_FUNC) &geMatrix_getDiag, 1},
    {"geMatrix_LU", (DL_FUNC) &geMatrix_LU, 1},
    {"geMatrix_determinant", (DL_FUNC) &geMatrix_determinant, 2},
    {"geMatrix_solve", (DL_FUNC) &geMatrix_solve, 1},
    {"geMatrix_geMatrix_mm", (DL_FUNC) &geMatrix_geMatrix_mm, 2},
    {"lmeRep_validate", (DL_FUNC) &lmeRep_validate, 1},
    {"lmeRep_crosstab", (DL_FUNC) &lmeRep_crosstab, 1},
    {"lmeRep_create", (DL_FUNC) &lmeRep_create, 2},
    {"lmeRep_update_mm", (DL_FUNC) &lmeRep_update_mm, 3},
    {"lmeRep_initial", (DL_FUNC) &lmeRep_initial, 1},
    {"lmeRep_factor", (DL_FUNC) &lmeRep_factor, 1},
    {"lmeRep_invert", (DL_FUNC) &lmeRep_invert, 1},
    {"lmeRep_sigma", (DL_FUNC) &lmeRep_sigma, 2},
    {"lmeRep_coef", (DL_FUNC) &lmeRep_coef, 2},
    {"lmeRep_coefGets", (DL_FUNC) &lmeRep_coefGets, 3},
    {"lmeRep_fixef", (DL_FUNC) &lmeRep_fixef, 1},
    {"lmeRep_ranef", (DL_FUNC) &lmeRep_ranef, 1},
    {"lmeRep_ECMEsteps", (DL_FUNC) &lmeRep_ECMEsteps, 4},
    {"lmeRep_gradient", (DL_FUNC) &lmeRep_gradient, 3},
    {"lmeRep_variances", (DL_FUNC) &lmeRep_variances, 1},
    {"poMatrix_rcond", (DL_FUNC) &poMatrix_rcond, 2},
    {"poMatrix_solve", (DL_FUNC) &poMatrix_solve, 1},
    {"poMatrix_matrix_solve", (DL_FUNC) &poMatrix_matrix_solve, 2},
    {"poMatrix_geMatrix_solve", (DL_FUNC) &poMatrix_geMatrix_solve, 2},
    {"poMatrix_chol", (DL_FUNC) &poMatrix_chol, 1},
    {"sscCrosstab", (DL_FUNC) &sscCrosstab, 2},
    {"sscCrosstab_groupedPerm", (DL_FUNC) &sscCrosstab_groupedPerm, 1},
    {"sscCrosstab_project2", (DL_FUNC) &sscCrosstab_project2, 1},
    {"sscMatrix_validate", (DL_FUNC) &sscMatrix_validate, 1},
    {"sscMatrix_chol", (DL_FUNC) &sscMatrix_chol, 2},
    {"sscMatrix_inverse_factor", (DL_FUNC) &sscMatrix_inverse_factor, 1},
    {"sscMatrix_matrix_solve", (DL_FUNC) &sscMatrix_matrix_solve, 2},
    {"ssc_transpose", (DL_FUNC) &ssc_transpose, 1},
    {"sscMatrix_to_triplet", (DL_FUNC) &sscMatrix_to_triplet, 1},
    {"sscMatrix_ldl_symbolic", (DL_FUNC) &sscMatrix_ldl_symbolic, 2},
    {"ssclme_create", (DL_FUNC) &ssclme_create, 2},
    {"ssclme_transfer_dimnames", (DL_FUNC) &ssclme_transfer_dimnames, 3},
    {"ssclme_update_mm", (DL_FUNC) &ssclme_update_mm, 3},
    {"ssclme_inflate_and_factor", (DL_FUNC) &ssclme_inflate_and_factor, 1},
    {"ssclme_factor", (DL_FUNC) &ssclme_factor, 1},
    {"ssclme_invert", (DL_FUNC) &ssclme_invert, 1},
    {"ssclme_initial", (DL_FUNC) &ssclme_initial, 1},
    {"ssclme_fixef", (DL_FUNC) &ssclme_fixef, 1},
    {"ssclme_ranef", (DL_FUNC) &ssclme_ranef, 1},
    {"ssclme_sigma", (DL_FUNC) &ssclme_sigma, 2},
    {"ssclme_coef", (DL_FUNC) &ssclme_coef, 2},
    {"ssclme_coefUnc", (DL_FUNC) &ssclme_coefUnc, 1},
    {"ssclme_coefGetsUnc", (DL_FUNC) &ssclme_coefGetsUnc, 2},
    {"ssclme_coefGets", (DL_FUNC) &ssclme_coefGets, 3},
    {"ssclme_EMsteps", (DL_FUNC) &ssclme_EMsteps, 4},
    {"ssclme_fitted", (DL_FUNC) &ssclme_fitted, 4},
    {"ssclme_variances", (DL_FUNC) &ssclme_variances, 1},
    {"ssclme_grad", (DL_FUNC) &ssclme_grad, 4},
    {"ssclme_gradient", (DL_FUNC) &ssclme_gradient, 3},
    {"ssclme_Hessian", (DL_FUNC) &ssclme_Hessian, 3},
    {"ssclme_collapse", (DL_FUNC) &ssclme_collapse, 1},
    {"ssclme_to_lme", (DL_FUNC) &ssclme_to_lme, 8},
    {"syMatrix_validate", (DL_FUNC) &syMatrix_validate, 1},
    {"syMatrix_norm", (DL_FUNC) &syMatrix_norm, 2},
    {"syMatrix_rcond", (DL_FUNC) &syMatrix_rcond, 2},
    {"syMatrix_solve", (DL_FUNC) &syMatrix_solve, 1},
    {"syMatrix_matrix_solve", (DL_FUNC) &syMatrix_matrix_solve, 2},
    {"syMatrix_as_geMatrix", (DL_FUNC) &syMatrix_as_geMatrix, 1},
    {"syMatrix_as_matrix", (DL_FUNC) &syMatrix_as_matrix, 1},
    {"syMatrix_geMatrix_mm", (DL_FUNC) &syMatrix_geMatrix_mm, 2},
    {"syMatrix_geMatrix_mm_R", (DL_FUNC) &syMatrix_geMatrix_mm_R, 2},
    {"trMatrix_validate", (DL_FUNC) &trMatrix_validate, 1},
    {"trMatrix_norm", (DL_FUNC) &trMatrix_norm, 2},
    {"trMatrix_rcond", (DL_FUNC) &trMatrix_rcond, 2},
    {"trMatrix_solve", (DL_FUNC) &trMatrix_solve, 1},
    {"trMatrix_matrix_solve", (DL_FUNC) &trMatrix_matrix_solve, 2},
    {"trMatrix_as_geMatrix", (DL_FUNC) &trMatrix_as_geMatrix, 1},
    {"trMatrix_as_matrix", (DL_FUNC) &trMatrix_as_matrix, 1},
    {"trMatrix_getDiag", (DL_FUNC) &trMatrix_getDiag, 1},
    {"trMatrix_geMatrix_mm", (DL_FUNC) &trMatrix_geMatrix_mm, 2},
    {"trMatrix_geMatrix_mm_R", (DL_FUNC) &trMatrix_geMatrix_mm_R, 2},
    {"triplet_validate", (DL_FUNC) &triplet_validate, 1},
    {"triplet_to_geMatrix", (DL_FUNC) &triplet_to_geMatrix, 1},
    {"tsc_validate", (DL_FUNC) &tsc_validate, 1},
    {"tsc_transpose", (DL_FUNC) &tsc_transpose, 1},
    {"tsc_to_triplet", (DL_FUNC) &tsc_to_triplet, 1},
    {NULL, NULL, 0}
};

void R_init_Matrix(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    Matrix_DSym = install("D");
    Matrix_DIsqrtSym = install("DIsqrt");
    Matrix_DimSym = install("Dim");
    Matrix_GpSym = install("Gp");
    Matrix_LSym = install("L");
    Matrix_LiSym = install("Li");
    Matrix_LinvSym = install("Linv");
    Matrix_LpSym = install("Lp");
    Matrix_LxSym = install("Lx");
    Matrix_OmegaSym = install("Omega");
    Matrix_ParentSym = install("Parent");
    Matrix_RXXSym = install("RXX");
    Matrix_RZXSym = install("RZX");
    Matrix_XtXSym = install("XtX");
    Matrix_ZtXSym = install("ZtX");
    Matrix_ZZxSym = install("ZZx");
    Matrix_bVarSym = install("bVar");
    Matrix_cnamesSym = install("cnames");
    Matrix_devianceSym = install("deviance");
    Matrix_devCompSym = install("devComp");
    Matrix_diagSym = install("diag");
    Matrix_factorization = install("factorization");
    Matrix_iSym = install("i");
    Matrix_ipermSym = install("iperm");
    Matrix_jSym = install("j");
    Matrix_matSym = install("mat");
    Matrix_ncSym = install("nc");
    Matrix_pSym = install("p");
    Matrix_permSym = install("perm");
    Matrix_rcondSym = install("rcond");
    Matrix_statusSym = install("status");
    Matrix_uploSym = install("uplo");
    Matrix_xSym = install("x");
    Matrix_zSym = install("z");
}
