#include "Mutils.h"
#include "HBMM.h"
#include "dense.h"
#include "dgBCMatrix.h"
#include "dgCMatrix.h"
#include "dgTMatrix.h"
#include "dgeMatrix.h"
#include "dpoMatrix.h"
#include "dppMatrix.h"
#include "dsCMatrix.h"
#include "dsTMatrix.h"
#include "dspMatrix.h"
#include "dsyMatrix.h"
#include "dtCMatrix.h"
#include "dtTMatrix.h"
#include "dtrMatrix.h"
#include "dtpMatrix.h"
#include "factorizations.h"
#include "lCholCMatrix.h"
#include "lgCMatrix.h"
#include "lgTMatrix.h"
#include "lsCMatrix.h"
#include "ltCMatrix.h"
#include "lmer.h"
#include <R_ext/Rdynload.h>

#include "Syms.h"

static R_CallMethodDef CallEntries[] = {
    {"BunchKaufman_validate", (DL_FUNC) &BunchKaufman_validate, 1},
    {"pBunchKaufman_validate", (DL_FUNC) &pBunchKaufman_validate, 1},
    {"Cholesky_validate", (DL_FUNC) &Cholesky_validate, 1},
    {"pCholesky_validate", (DL_FUNC) &pCholesky_validate, 1},
    {"LU_expand", (DL_FUNC) &LU_expand, 1},
    {"LU_validate", (DL_FUNC) &LU_validate, 1},
    {"Matrix_expand_pointers", (DL_FUNC) &Matrix_expand_pointers, 1},
    {"Matrix_readHarwellBoeing", (DL_FUNC) &Matrix_readHarwellBoeing, 1},
    {"Matrix_readMatrixMarket", (DL_FUNC) &Matrix_readMatrixMarket, 1},
    {"Matrix_rWishart", (DL_FUNC) &Matrix_rWishart, 3},
    {"Matrix_writeHarwellBoeing", (DL_FUNC) &Matrix_writeHarwellBoeing, 3},
    {"Matrix_writeMatrixMarket", (DL_FUNC) &Matrix_writeMatrixMarket, 3},
    {"SVD_validate", (DL_FUNC) &SVD_validate, 1},
    {"csc_check_column_sorting", (DL_FUNC) &csc_check_column_sorting, 1},
    {"csc_crossprod", (DL_FUNC) &csc_crossprod, 1},
    {"csc_getDiag", (DL_FUNC) &csc_getDiag, 1},
    {"csc_matrix_crossprod", (DL_FUNC) &csc_matrix_crossprod, 3},
    {"csc_matrix_mm", (DL_FUNC) &csc_matrix_mm, 4},
    {"csc_tcrossprod", (DL_FUNC) &csc_tcrossprod, 1},
    {"compressed_to_dgTMatrix", (DL_FUNC) &compressed_to_dgTMatrix, 2},
    {"csc_to_dgeMatrix", (DL_FUNC) &csc_to_dgeMatrix, 1},
    {"csc_to_matrix", (DL_FUNC) &csc_to_matrix, 1},
    {"csc_transpose", (DL_FUNC) &csc_transpose, 1},
    {"dgBCMatrix_to_dgCMatrix", (DL_FUNC) &dgBCMatrix_to_dgCMatrix, 1},
    {"dgBCMatrix_to_dgTMatrix", (DL_FUNC) &dgBCMatrix_to_dgTMatrix, 1},
    {"dgBCMatrix_validate", (DL_FUNC) &dgBCMatrix_validate, 1},
    {"dgCMatrix_validate", (DL_FUNC) &dgCMatrix_validate, 1},
    {"dgTMatrix_to_csc", (DL_FUNC) &dgTMatrix_to_csc, 1},
    {"dgTMatrix_to_dgCMatrix", (DL_FUNC) &dgTMatrix_to_dgCMatrix, 1},
    {"dgTMatrix_to_dgeMatrix", (DL_FUNC) &dgTMatrix_to_dgeMatrix, 1},
    {"dgTMatrix_to_matrix", (DL_FUNC) &dgTMatrix_to_matrix, 1},
    {"dgTMatrix_validate", (DL_FUNC) &dgTMatrix_validate, 1},
    {"dgeMatrix_LU", (DL_FUNC) &dgeMatrix_LU, 1},
    {"dgeMatrix_Schur", (DL_FUNC) &dgeMatrix_Schur, 2},
    {"dgeMatrix_colsums", (DL_FUNC) &dgeMatrix_colsums, 4},
    {"dgeMatrix_crossprod", (DL_FUNC) &dgeMatrix_crossprod, 2},
    {"dgeMatrix_determinant", (DL_FUNC) &dgeMatrix_determinant, 2},
    {"dgeMatrix_dgeMatrix_crossprod", (DL_FUNC) &dgeMatrix_dgeMatrix_crossprod, 2},
    {"dgeMatrix_matrix_mm", (DL_FUNC) &dgeMatrix_matrix_mm, 4},
    {"dgeMatrix_matrix_solve", (DL_FUNC) &dgeMatrix_matrix_solve, 3},
    {"dgeMatrix_dtpMatrix_mm", (DL_FUNC) &dgeMatrix_dtpMatrix_mm, 2},
    {"dgeMatrix_exp", (DL_FUNC) &dgeMatrix_exp, 1},
    {"dgeMatrix_getDiag", (DL_FUNC) &dgeMatrix_getDiag, 1},
    {"dgeMatrix_matrix_crossprod", (DL_FUNC) &dgeMatrix_matrix_crossprod, 2},
    {"dgeMatrix_norm", (DL_FUNC) &dgeMatrix_norm, 2},
    {"dgeMatrix_rcond", (DL_FUNC) &dgeMatrix_rcond, 2},
    {"dgeMatrix_solve", (DL_FUNC) &dgeMatrix_solve, 1},
    {"dgeMatrix_validate", (DL_FUNC) &dgeMatrix_validate, 1},
    {"dpoMatrix_chol", (DL_FUNC) &dpoMatrix_chol, 1},
    {"dpoMatrix_dgeMatrix_solve", (DL_FUNC) &dpoMatrix_dgeMatrix_solve, 2},
    {"dpoMatrix_matrix_solve", (DL_FUNC) &dpoMatrix_matrix_solve, 2},
    {"dpoMatrix_rcond", (DL_FUNC) &dpoMatrix_rcond, 2},
    {"dpoMatrix_solve", (DL_FUNC) &dpoMatrix_solve, 1},
    {"dpoMatrix_validate", (DL_FUNC) &dpoMatrix_validate, 1},
    {"dppMatrix_chol", (DL_FUNC) &dppMatrix_chol, 1},
    {"dppMatrix_matrix_solve", (DL_FUNC) &dppMatrix_matrix_solve, 3},
    {"dppMatrix_rcond", (DL_FUNC) &dppMatrix_rcond, 2},
    {"dppMatrix_solve", (DL_FUNC) &dppMatrix_solve, 1},
    {"dppMatrix_validate", (DL_FUNC) &dppMatrix_validate, 1},
    {"dsCMatrix_chol", (DL_FUNC) &dsCMatrix_chol, 2},
    {"dsCMatrix_inverse_factor", (DL_FUNC) &dsCMatrix_inverse_factor, 1},
    {"dsCMatrix_ldl_symbolic", (DL_FUNC) &dsCMatrix_ldl_symbolic, 2},
    {"dsCMatrix_matrix_solve", (DL_FUNC) &dsCMatrix_matrix_solve, 3},
    {"dsCMatrix_to_dgTMatrix", (DL_FUNC) &dsCMatrix_to_dgTMatrix, 1},
    {"dsCMatrix_validate", (DL_FUNC) &dsCMatrix_validate, 1},
    {"dsTMatrix_as_dsCMatrix", (DL_FUNC) &dsTMatrix_as_dsCMatrix, 1},
    {"dsTMatrix_as_dsyMatrix", (DL_FUNC) &dsTMatrix_as_dsyMatrix, 1},
    {"dsTMatrix_validate", (DL_FUNC) &dsTMatrix_validate, 1},
    {"dsyMatrix_as_dgeMatrix", (DL_FUNC) &dsyMatrix_as_dgeMatrix, 1},
    {"dsyMatrix_as_dspMatrix", (DL_FUNC) &dsyMatrix_as_dspMatrix, 1},
    {"dsyMatrix_as_matrix", (DL_FUNC) &dsyMatrix_as_matrix, 1},
    {"dsyMatrix_dgeMatrix_mm", (DL_FUNC) &dsyMatrix_dgeMatrix_mm, 2},
    {"dsyMatrix_dgeMatrix_mm_R", (DL_FUNC) &dsyMatrix_dgeMatrix_mm_R, 2},
    {"dsyMatrix_matrix_solve", (DL_FUNC) &dsyMatrix_matrix_solve, 2},
    {"dsyMatrix_dgeMatrix_solve", (DL_FUNC) &dsyMatrix_dgeMatrix_solve, 2},
    {"dsyMatrix_norm", (DL_FUNC) &dsyMatrix_norm, 2},
    {"dsyMatrix_rcond", (DL_FUNC) &dsyMatrix_rcond, 2},
    {"dsyMatrix_solve", (DL_FUNC) &dsyMatrix_solve, 1},
    {"dsyMatrix_validate", (DL_FUNC) &dsyMatrix_validate, 1},
    {"dspMatrix_trf", (DL_FUNC) &dspMatrix_trf, 1},
    {"dspMatrix_as_dsyMatrix", (DL_FUNC) &dspMatrix_as_dsyMatrix, 1},
    {"dspMatrix_matrix_mm", (DL_FUNC) &dspMatrix_matrix_mm, 3},
    {"dspMatrix_matrix_solve", (DL_FUNC) &dspMatrix_matrix_solve, 3},
    {"dspMatrix_norm", (DL_FUNC) &dspMatrix_norm, 2},
    {"dspMatrix_rcond", (DL_FUNC) &dspMatrix_rcond, 2},
    {"dspMatrix_solve", (DL_FUNC) &dspMatrix_solve, 1},
    {"dspMatrix_trf", (DL_FUNC) &dspMatrix_trf, 1},
    {"dspMatrix_validate", (DL_FUNC) &dspMatrix_validate, 1},
    {"dtTMatrix_as_dtrMatrix", (DL_FUNC) &dtTMatrix_as_dtrMatrix, 1},
    {"dtTMatrix_validate", (DL_FUNC) &dtTMatrix_validate, 1},
    {"dtpMatrix_as_dtrMatrix", (DL_FUNC) &dtpMatrix_as_dtrMatrix, 1},
    {"dtpMatrix_dgeMatrix_mm", (DL_FUNC) &dtpMatrix_dgeMatrix_mm, 2},
    {"dtpMatrix_getDiag", (DL_FUNC) &dtpMatrix_getDiag, 1},
    {"dtpMatrix_matrix_mm", (DL_FUNC) &dtpMatrix_matrix_mm, 2},
    {"dtpMatrix_matrix_solve", (DL_FUNC) &dtpMatrix_matrix_solve, 2},
    {"dtpMatrix_norm", (DL_FUNC) &dtpMatrix_norm, 2},
    {"dtpMatrix_rcond", (DL_FUNC) &dtpMatrix_rcond, 2},
    {"dtpMatrix_solve", (DL_FUNC) &dtpMatrix_solve, 1},
    {"dtpMatrix_validate", (DL_FUNC) &dtpMatrix_validate, 1},
    {"dtrMatrix_as_dgeMatrix", (DL_FUNC) &dtrMatrix_as_dgeMatrix, 1},
    {"dtrMatrix_as_dtpMatrix", (DL_FUNC) &dtrMatrix_as_dtpMatrix, 1},
    {"dtrMatrix_as_matrix", (DL_FUNC) &dtrMatrix_as_matrix, 1},
    {"dtrMatrix_matrix_mm", (DL_FUNC) &dtrMatrix_matrix_mm, 4},
    {"dtrMatrix_getDiag", (DL_FUNC) &dtrMatrix_getDiag, 1},
    {"dtrMatrix_matrix_solve", (DL_FUNC) &dtrMatrix_matrix_solve, 3},
    {"dtrMatrix_norm", (DL_FUNC) &dtrMatrix_norm, 2},
    {"dtrMatrix_rcond", (DL_FUNC) &dtrMatrix_rcond, 2},
    {"dtrMatrix_solve", (DL_FUNC) &dtrMatrix_solve, 1},
    {"dtrMatrix_validate", (DL_FUNC) &dtrMatrix_validate, 1},
    {"glmer_MCMCsamp", (DL_FUNC) &glmer_MCMCsamp, 6},
    {"glmer_PQL", (DL_FUNC) &glmer_PQL, 1},
    {"glmer_devAGQ", (DL_FUNC) &glmer_devAGQ, 3},
    {"glmer_finalize", (DL_FUNC) &glmer_finalize, 1},
    {"glmer_fixed_update", (DL_FUNC) &glmer_fixed_update, 3},
    {"glmer_init", (DL_FUNC) &glmer_init, 1},
    {"glmer_ranef_update", (DL_FUNC) &glmer_ranef_update, 4},
    {"lapack_qr", (DL_FUNC) &lapack_qr, 2},
    {"lCholCMatrix_solve", (DL_FUNC) &lCholCMatrix_solve, 1},
    {"lCholCMatrix_lgCMatrix_solve", (DL_FUNC) &lCholCMatrix_lgCMatrix_solve, 2},
    {"lCholCMatrix_validate", (DL_FUNC) &lCholCMatrix_validate, 1},
    {"lgCMatrix_crossprod", (DL_FUNC) &lgCMatrix_crossprod, 3},
    {"lgCMatrix_lgCMatrix_mm", (DL_FUNC) &lgCMatrix_lgCMatrix_mm, 2},
    {"lgCMatrix_picky_column", (DL_FUNC) &lgCMatrix_picky_column, 1},
    {"lgCMatrix_trans", (DL_FUNC) &lgCMatrix_trans, 1},
    {"lgCMatrix_validate", (DL_FUNC) &lgCMatrix_validate, 1},
    {"lgTMatrix_as_lgCMatrix", (DL_FUNC) &lgTMatrix_as_lgCMatrix, 1},
    {"lgTMatrix_validate", (DL_FUNC) &lgTMatrix_validate, 1},
    {"lmer_Crosstab", (DL_FUNC) &lmer_Crosstab, 1},
    {"lmer_MCMCsamp", (DL_FUNC) &lmer_MCMCsamp, 4},
    {"lmer_ECMEsteps", (DL_FUNC) &lmer_ECMEsteps, 3},
    {"lmer_coef", (DL_FUNC) &lmer_coef, 2},
    {"lmer_coefGets", (DL_FUNC) &lmer_coefGets, 3},
    {"lmer_create", (DL_FUNC) &lmer_create, 3},
    {"lmer_factor", (DL_FUNC) &lmer_factor, 1},
    {"lmer_firstDer", (DL_FUNC) &lmer_firstDer, 2},
    {"lmer_fitted", (DL_FUNC) &lmer_fitted, 3},
    {"lmer_fixef", (DL_FUNC) &lmer_fixef, 1},
    {"lmer_gradient", (DL_FUNC) &lmer_gradient, 2},
    {"lmer_inflate", (DL_FUNC) &lmer_inflate, 1},
    {"lmer_initial", (DL_FUNC) &lmer_initial, 1},
    {"lmer_invert", (DL_FUNC) &lmer_invert, 1},
    {"lmer_ranef", (DL_FUNC) &lmer_ranef, 1},
    {"lmer_secondDer", (DL_FUNC) &lmer_secondDer, 1},
    {"lmer_sigma", (DL_FUNC) &lmer_sigma, 2},
    {"lmer_update_mm", (DL_FUNC) &lmer_update_mm, 2},
    {"lmer_validate", (DL_FUNC) &lmer_validate, 1},
    {"lmer_variances", (DL_FUNC) &lmer_variances, 1},
    {"lsCMatrix_chol", (DL_FUNC) &lsCMatrix_chol, 2},
    {"lsCMatrix_trans", (DL_FUNC) &lsCMatrix_trans, 1},
    {"lsCMatrix_validate", (DL_FUNC) &lsCMatrix_validate, 1},
    {"ltCMatrix_trans", (DL_FUNC) &ltCMatrix_trans, 1},
    {"ltCMatrix_validate", (DL_FUNC) &ltCMatrix_validate, 1},
    {"lsq_dense_Chol", (DL_FUNC) &lsq_dense_Chol, 2},
    {"lsq_dense_QR", (DL_FUNC) &lsq_dense_QR, 2},
    {"matrix_to_csc", (DL_FUNC) &matrix_to_csc, 1},
    {"ssc_transpose", (DL_FUNC) &ssc_transpose, 1},
    {"tsc_to_dgTMatrix", (DL_FUNC) &tsc_to_dgTMatrix, 1},
    {"tsc_transpose", (DL_FUNC) &tsc_transpose, 1},
    {"tsc_validate", (DL_FUNC) &tsc_validate, 1},
    {NULL, NULL, 0}
};

void R_init_Matrix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE); 
    Matrix_DIsqrtSym = install("DIsqrt");
    Matrix_DSym = install("D");
    Matrix_DimSym = install("Dim");
    Matrix_DimNamesSym = install("Dimnames");
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
    Matrix_ZZxSym = install("ZZx");
    Matrix_ZZpOSym = install("ZZpO");
    Matrix_ZtXSym = install("ZtX");
    Matrix_ZtZSym = install("ZtZ");
    Matrix_bVarSym = install("bVar");
    Matrix_cnamesSym = install("cnames");
    Matrix_devCompSym = install("devComp");
    Matrix_devianceSym = install("deviance");
    Matrix_diagSym = install("diag");
    Matrix_factorSym = install("factors");
    Matrix_flistSym = install("flist");
    Matrix_iSym = install("i");
    Matrix_ipermSym = install("iperm");
    Matrix_jSym = install("j");
    Matrix_matSym = install("mat");
    Matrix_methodSym = install("method");
    Matrix_ncSym = install("nc");
    Matrix_pSym = install("p");
    Matrix_permSym = install("perm");
    Matrix_rcondSym = install("rcond");
    Matrix_statusSym = install("status");
    Matrix_uploSym = install("uplo");
    Matrix_xSym = install("x");
    Matrix_zSym = install("z");
}
