#include <Rinternals.h> /* SEXP, Rcomplex */
#include "Csparse.h"
#include "CHMfactor.h"
#include "abIndex.h"
#include "chm_common.h"
#include "dense.h"
#include "dgCMatrix.h"
#include "dgeMatrix.h"
#include "factorizations.h"
#include "kappa.h"
#include "packedMatrix.h"
#include "products.h"
#include "sparse.h"
#include "sparseVector.h"
#include "subscript.h"
#include "unpackedMatrix.h"
#include "validity.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "Syms.h"
Rcomplex Matrix_zzero, Matrix_zone, Matrix_zna;

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define EXTDEF(name, n)   {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
	CALLDEF(CHM_set_common_env, 1),
	CALLDEF(get_SuiteSparse_version, 0),

	CALLDEF(m_encodeInd,  4),
	CALLDEF(m_encodeInd2, 5),

	CALLDEF(Matrix_rle_i, 2),
	CALLDEF(Matrix_rle_d, 2),

	CALLDEF(Csparse_Csparse_prod, 3),
	CALLDEF(Csparse_Csparse_crossprod, 4),
	CALLDEF(Csparse_MatrixMarket, 2),
	CALLDEF(Csparse_crossprod, 4),
	CALLDEF(Csparse_dense_crossprod, 3),
	CALLDEF(Csparse_dense_prod, 3),
	CALLDEF(Csparse_drop, 2),
	CALLDEF(Csparse_horzcat, 2),
	CALLDEF(Csparse_sort, 1),

	CALLDEF(dCsparse_subassign, 4),
	CALLDEF(lCsparse_subassign, 4),
	CALLDEF(iCsparse_subassign, 4),
	CALLDEF(nCsparse_subassign, 4),
	CALLDEF(zCsparse_subassign, 4),

	CALLDEF(Csparse_validate2, 2),
	CALLDEF(Csparse_vertcat, 2),
	CALLDEF(Csparse_dmperm, 3),
	CALLDEF(diag_tC, 2),

	CALLDEF(Matrix_expand_pointers, 1),
	CALLDEF(R_rbind2_vector, 2),
	CALLDEF(R_all0, 1),
	CALLDEF(R_any0, 1),

	CALLDEF(compressed_non_0_ij, 2),
	CALLDEF(dgCMatrix_lusol, 2),
	CALLDEF(dgCMatrix_qrsol, 3),
	CALLDEF(dgCMatrix_cholsol, 2),

	CALLDEF(dgeMatrix_Schur, 3),
	CALLDEF(dgeMatrix_exp, 1),

	CALLDEF(dgeMatrix_crossprod, 2),
	CALLDEF (geMatrix_crossprod, 2),
	CALLDEF(dgeMatrix_dgeMatrix_crossprod, 3),
	CALLDEF (geMatrix_geMatrix_crossprod, 3),
	CALLDEF(dgeMatrix_matrix_crossprod, 3),
	CALLDEF (geMatrix_matrix_crossprod, 3),
	CALLDEF(dgeMatrix_matrix_mm, 3),
	CALLDEF (geMatrix_matrix_mm, 3),
	CALLDEF(dtrMatrix_dtrMatrix_mm, 4),
	CALLDEF(dtrMatrix_matrix_mm, 4),
	CALLDEF(dtpMatrix_matrix_mm, 4),
	CALLDEF(dgeMatrix_dtpMatrix_mm, 2),
	CALLDEF(dsyMatrix_matrix_mm, 3),
	CALLDEF(dspMatrix_matrix_mm, 2),

	CALLDEF(Matrix_validate, 1),
	CALLDEF(MatrixFactorization_validate, 1),

	CALLDEF(dMatrix_validate, 1),
	CALLDEF(lMatrix_validate, 1),
	CALLDEF(ndenseMatrix_validate, 1),
	CALLDEF(iMatrix_validate, 1),
	CALLDEF(zMatrix_validate, 1),

	CALLDEF(compMatrix_validate, 1),
	CALLDEF(symmetricMatrix_validate, 1),
	CALLDEF(triangularMatrix_validate, 1),

	CALLDEF(diagonalMatrix_validate, 1),
	CALLDEF(indMatrix_validate, 1),
	CALLDEF(pMatrix_validate, 1),

	CALLDEF(CsparseMatrix_validate, 1),
	CALLDEF(RsparseMatrix_validate, 1),
	CALLDEF(TsparseMatrix_validate, 1),

	CALLDEF(sCMatrix_validate, 1),
	CALLDEF(tCMatrix_validate, 1),
	CALLDEF(sRMatrix_validate, 1),
	CALLDEF(tRMatrix_validate, 1),
	CALLDEF(sTMatrix_validate, 1),
	CALLDEF(tTMatrix_validate, 1),

	CALLDEF(xgCMatrix_validate, 1),
	CALLDEF(xsCMatrix_validate, 1),
	CALLDEF(xtCMatrix_validate, 1),
	CALLDEF(xgRMatrix_validate, 1),
	CALLDEF(xsRMatrix_validate, 1),
	CALLDEF(xtRMatrix_validate, 1),
	CALLDEF(xgTMatrix_validate, 1),
	CALLDEF(xsTMatrix_validate, 1),
	CALLDEF(xtTMatrix_validate, 1),

	CALLDEF(unpackedMatrix_validate, 1),
	CALLDEF(packedMatrix_validate, 1),

	CALLDEF(dpoMatrix_validate, 1),
	CALLDEF(dppMatrix_validate, 1),
	CALLDEF(corMatrix_validate, 1),
	CALLDEF(pcorMatrix_validate, 1),

	CALLDEF(denseLU_validate, 1),
	CALLDEF(sparseLU_validate, 1),
	CALLDEF(sparseQR_validate, 1),
	CALLDEF(BunchKaufman_validate, 1),
	CALLDEF(pBunchKaufman_validate, 1),
	CALLDEF(Cholesky_validate, 1),
	CALLDEF(pCholesky_validate, 1),
	CALLDEF(CHMfactor_validate, 1),
	CALLDEF(CHMsimpl_validate, 1),
	CALLDEF(CHMsuper_validate, 1),
	CALLDEF(dCHMsimpl_validate, 1),
	CALLDEF(dCHMsuper_validate, 1),
	CALLDEF(Schur_validate, 1),

	CALLDEF(R_Dim_validate, 1),
	CALLDEF(R_DimNames_validate, 2),

	CALLDEF(R_DimNames_fixup, 1),
	CALLDEF(R_DimNames_is_symmetric, 1),
	CALLDEF(R_symmDN, 1),
	CALLDEF(R_revDN, 1),
	CALLDEF(R_Matrix_nonvirtual, 2),
	CALLDEF(R_Matrix_kind, 2),
	CALLDEF(R_Matrix_shape, 1),
	CALLDEF(R_Matrix_repr, 1),
	CALLDEF(R_index_triangle, 4),
	CALLDEF(R_index_diagonal, 3),
	CALLDEF(R_nnz, 3),
	CALLDEF(R_isPerm, 2),
	CALLDEF(R_signPerm, 2),
	CALLDEF(R_invertPerm, 3),
	CALLDEF(R_asPerm, 4),
	CALLDEF(R_set_factor, 4),
	CALLDEF(R_empty_factors, 2),

	CALLDEF(R_subscript_1ary, 2),
	CALLDEF(R_subscript_1ary_mat, 2),
	CALLDEF(R_subscript_2ary, 3),

	CALLDEF(R_sparse_as_dense, 2),
	CALLDEF(R_sparse_as_matrix, 1),
	CALLDEF(R_sparse_as_vector, 1),
	CALLDEF(R_sparse_as_kind, 3),
	CALLDEF(R_sparse_as_general, 1),

	CALLDEF(R_diagonal_as_sparse, 4),
	CALLDEF(R_diagonal_as_dense, 3),
	CALLDEF(R_diagonal_as_kind, 2),

	CALLDEF(R_sparse_drop0, 1),
	CALLDEF(R_sparse_band, 3),
	CALLDEF(R_sparse_diag_get, 2),
	CALLDEF(R_sparse_diag_set, 2),
	CALLDEF(R_sparse_diag_U2N, 1),
	CALLDEF(R_sparse_diag_N2U, 1),
	CALLDEF(R_sparse_transpose, 1),
	CALLDEF(R_sparse_force_symmetric, 2),
	CALLDEF(R_sparse_symmpart, 1),
	CALLDEF(R_sparse_skewpart, 1),

	CALLDEF(CRsparse_as_Tsparse, 1),
	CALLDEF(Tsparse_as_CRsparse, 2),
	CALLDEF(Tsparse_aggregate, 1),
	CALLDEF(tCRsparse_as_RCsparse, 1),
	CALLDEF(Csparse_is_diagonal, 1),
	CALLDEF(Rsparse_is_diagonal, 1),
	CALLDEF(Tsparse_is_diagonal, 1),
	CALLDEF(Csparse_is_triangular, 2),
	CALLDEF(Rsparse_is_triangular, 2),
	CALLDEF(Tsparse_is_triangular, 2),
	CALLDEF(Csparse_is_symmetric, 2),
	CALLDEF(Rsparse_is_symmetric, 2),
	CALLDEF(CRsparse_colSums, 4),
	CALLDEF(CRsparse_rowSums, 4),

	CALLDEF(R_matrix_as_dense, 4),
	CALLDEF(R_dense_as_general, 2),
	CALLDEF(R_dense_as_sparse, 4),
	CALLDEF(R_dense_as_kind, 2),
	CALLDEF(R_dense_as_matrix, 1),
	CALLDEF(R_geMatrix_as_matrix, 2),
	CALLDEF(R_dense_as_vector, 1),
	CALLDEF(R_geMatrix_as_vector, 2),
	CALLDEF(R_dense_band, 3),
	CALLDEF(R_dense_colSums, 3),
	CALLDEF(R_dense_rowSums, 3),

	CALLDEF(matrix_is_symmetric, 2),
	CALLDEF(matrix_is_triangular, 2),
	CALLDEF(matrix_is_diagonal, 1),
	CALLDEF(matrix_symmpart, 1),
	CALLDEF(matrix_skewpart, 1),

	CALLDEF(unpackedMatrix_pack, 4),
	CALLDEF(unpackedMatrix_force_symmetric, 2),
	CALLDEF(unpackedMatrix_is_symmetric, 2),
	CALLDEF(unpackedMatrix_is_triangular, 2),
	CALLDEF(unpackedMatrix_is_diagonal, 1),
	CALLDEF(unpackedMatrix_transpose, 1),
	CALLDEF(unpackedMatrix_diag_get, 2),
	CALLDEF(unpackedMatrix_diag_set, 2),
	CALLDEF(unpackedMatrix_symmpart, 1),
	CALLDEF(unpackedMatrix_skewpart, 1),

	CALLDEF(packedMatrix_unpack, 2),
	CALLDEF(packedMatrix_force_symmetric, 2),
	CALLDEF(packedMatrix_is_symmetric, 2),
	CALLDEF(packedMatrix_is_triangular, 2),
	CALLDEF(packedMatrix_is_diagonal, 1),
	CALLDEF(packedMatrix_transpose, 1),
	CALLDEF(packedMatrix_diag_get, 2),
	CALLDEF(packedMatrix_diag_set, 2),
	CALLDEF(packedMatrix_symmpart, 1),
	CALLDEF(packedMatrix_skewpart, 1),

	CALLDEF(dgeMatrix_trf, 2),
	CALLDEF(dsyMatrix_trf, 2),
	CALLDEF(dspMatrix_trf, 2),
	CALLDEF(dpoMatrix_trf, 4),
	CALLDEF(dppMatrix_trf, 2),
	CALLDEF(dgCMatrix_trf, 4),
	CALLDEF(dgCMatrix_orf, 3),
	CALLDEF(dpCMatrix_trf, 5),

	CALLDEF(BunchKaufman_expand, 2),

	CALLDEF(denseLU_determinant, 2),
	CALLDEF(BunchKaufman_determinant, 3),
	CALLDEF(Cholesky_determinant, 3),
	CALLDEF(sparseLU_determinant, 2),
	CALLDEF(sparseQR_determinant, 2),
	CALLDEF(CHMfactor_determinant, 3),

	CALLDEF(denseLU_solve, 2),
	CALLDEF(BunchKaufman_solve, 3),
	CALLDEF(Cholesky_solve, 3),
	CALLDEF(sparseLU_solve, 3),
/* MJ: not needed since we have 'sparseQR_matmult' : */
#if 0
	CALLDEF(sparseQR_solve, 3),
#endif
	CALLDEF(CHMfactor_solve, 4),
	CALLDEF(dtrMatrix_solve, 3),
	CALLDEF(dtCMatrix_solve, 3),

	CALLDEF(sparseQR_matmult, 5),

	CALLDEF(CHMfactor_diag_get, 2),
	CALLDEF(CHMfactor_update, 3),
	CALLDEF(CHMfactor_updown, 3),

	CALLDEF(dgeMatrix_norm, 2),
	CALLDEF(dtrMatrix_norm, 2),
	CALLDEF(dtpMatrix_norm, 2),
	CALLDEF(dsyMatrix_norm, 2),
	CALLDEF(dspMatrix_norm, 2),

	CALLDEF(dgeMatrix_rcond, 3),
	CALLDEF(dtrMatrix_rcond, 2),
	CALLDEF(dtpMatrix_rcond, 2),
	CALLDEF(dsyMatrix_rcond, 3),
	CALLDEF(dspMatrix_rcond, 3),
	CALLDEF(dpoMatrix_rcond, 3),
	CALLDEF(dppMatrix_rcond, 3),

	CALLDEF(v2spV, 1),
	CALLDEF(CR2spV, 1),

	{NULL, NULL, 0}
};

static const R_ExternalMethodDef ExtEntries[] = {
	EXTDEF(Mmatrix, 7),
	{NULL, NULL, 0}
};

void attribute_visible R_init_Matrix(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, CallEntries, NULL, ExtEntries);
	R_useDynamicSymbols(dll, FALSE);

/* These are callable from other packages' C code: */

#define RREGDEF(name)  R_RegisterCCallable("Matrix", #name, (DL_FUNC) name)

	/* Matrix: SEXP -> CHOLMOD */
	RREGDEF(as_cholmod_dense);
	RREGDEF(as_cholmod_factor);
	RREGDEF(as_cholmod_sparse);
	RREGDEF(as_cholmod_triplet);

	/* Matrix: CHOLMOD -> SEXP */
	RREGDEF(chm_factor_to_SEXP);
	RREGDEF(chm_sparse_to_SEXP);
	RREGDEF(chm_triplet_to_SEXP);

	/* Matrix: miscellaneous */
	RREGDEF(chm_factor_ldetL2);
	RREGDEF(chm_factor_update);
	RREGDEF(numeric_as_chm_dense);
#if 0
	RREGDEF(Csparse_diagU2N);
	RREGDEF(dpoMatrix_chol);
#else
	R_RegisterCCallable(
	"Matrix", "Csparse_diagU2N", (DL_FUNC) R_sparse_diag_U2N);
	R_RegisterCCallable(
	"Matrix", "dpoMatrix_chol", (DL_FUNC) dpoMatrix_trf);
#endif

	/* CHOLMOD: */
	RREGDEF(cholmod_aat);
	RREGDEF(cholmod_add);
	RREGDEF(cholmod_allocate_dense);
	RREGDEF(cholmod_allocate_sparse);
	RREGDEF(cholmod_allocate_triplet);
	RREGDEF(cholmod_analyze);
	RREGDEF(cholmod_analyze_p);
	RREGDEF(cholmod_band_inplace);
	RREGDEF(cholmod_change_factor);
	RREGDEF(cholmod_copy);
	RREGDEF(cholmod_copy_dense);
	RREGDEF(cholmod_copy_factor);
	RREGDEF(cholmod_copy_sparse);
	RREGDEF(cholmod_dense_to_sparse);
	RREGDEF(cholmod_factor_to_sparse);
	RREGDEF(cholmod_factorize);
	RREGDEF(cholmod_factorize_p);
	RREGDEF(cholmod_finish);
	RREGDEF(cholmod_free_dense);
	RREGDEF(cholmod_free_factor);
	RREGDEF(cholmod_free_sparse);
	RREGDEF(cholmod_free_triplet);
	RREGDEF(cholmod_nnz);
	RREGDEF(cholmod_scale);
	RREGDEF(cholmod_sdmult);
	RREGDEF(cholmod_solve);
	RREGDEF(cholmod_solve2);
	RREGDEF(cholmod_sort);
	RREGDEF(cholmod_sparse_to_dense);
	RREGDEF(cholmod_sparse_to_triplet);
	RREGDEF(cholmod_speye);
	RREGDEF(cholmod_spsolve);
	RREGDEF(cholmod_ssmult);
	RREGDEF(cholmod_start);
	RREGDEF(cholmod_submatrix);
	RREGDEF(cholmod_transpose);
	RREGDEF(cholmod_triplet_to_sparse);
	RREGDEF(cholmod_vertcat);
	RREGDEF(cholmod_updown);

	R_cholmod_start(&c);
#if 0
	/* TODO: needs more work in ./chm_common.c, etc. */
	R_cholmod_start(&cl);
#endif

	Matrix_DimNamesSym = install("Dimnames");
	Matrix_DimSym      = install("Dim");
	Matrix_LSym        = install("L");
	Matrix_QSym        = install("Q");
	Matrix_RSym        = install("R");
	Matrix_TSym        = install("T");
	Matrix_USym        = install("U");
	Matrix_VSym        = install("V");
	Matrix_betaSym     = install("beta");
	Matrix_diagSym     = install("diag");
	Matrix_factorSym   = install("factors");
	Matrix_iSym        = install("i");
	Matrix_jSym        = install("j");
	Matrix_lengthSym   = install("length");
	Matrix_marginSym   = install("margin");
	Matrix_pSym        = install("p");
	Matrix_permSym     = install("perm");
	Matrix_qSym        = install("q");
	Matrix_sdSym       = install("sd");
	Matrix_uploSym     = install("uplo");
	Matrix_xSym        = install("x");

	MatrixNamespace = R_FindNamespace(mkString("Matrix"));
	if (MatrixNamespace == R_UnboundValue)
		error(_("missing 'Matrix' namespace; should never happen"));
#ifdef Matrix_Debug
	if(isEnvironment(MatrixNamespace))
		Rprintf("MatrixNamespace: %s\n",
		        CHAR(asChar(eval(lang2(install("format"), MatrixNamespace),
		                         R_GlobalEnv))));
	else
#else
	if (!isEnvironment(MatrixNamespace))
#endif
		error(_("'Matrix' namespace not determined correctly"));

	Matrix_zzero.r = 0.0; Matrix_zone.r = 1.0; Matrix_zna.r = NA_REAL;
	Matrix_zzero.i = 0.0; Matrix_zone.i = 0.0; Matrix_zna.i = NA_REAL;
}

void R_unload_Matrix(DllInfo *dll)
{
	cholmod_finish(&c);
}
