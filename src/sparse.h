#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include "Mutils.h"

/* defined in ./sparseVector.c : */
SEXP v2spV(SEXP);

SEXP sparse_as_dense(SEXP from, int packed);
SEXP R_sparse_as_dense(SEXP from, SEXP packed);
SEXP R_sparse_as_matrix(SEXP from);
SEXP R_sparse_as_vector(SEXP from);
SEXP sparse_as_kind(SEXP from, char kind, int drop0);
SEXP R_sparse_as_kind(SEXP from, SEXP kind, SEXP drop0);
SEXP R_sparse_as_general(SEXP from);

SEXP R_diagonal_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP drop0);
SEXP R_diagonal_as_dense(SEXP from, SEXP code, SEXP uplo);
SEXP R_diagonal_as_kind(SEXP from, SEXP kind);

SEXP R_sparse_drop0(SEXP from);
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2);
SEXP R_sparse_diag_get(SEXP obj, SEXP nms);
SEXP R_sparse_diag_set(SEXP obj, SEXP val);
SEXP R_sparse_diag_U2N(SEXP obj);
SEXP R_sparse_diag_N2U(SEXP obj);
SEXP R_sparse_transpose(SEXP from);
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo_to);
SEXP R_sparse_symmpart(SEXP from);
SEXP R_sparse_skewpart(SEXP from);

SEXP CRsparse_as_Tsparse(SEXP from);
SEXP Tsparse_as_CRsparse(SEXP from, SEXP Csparse);
SEXP Tsparse_aggregate(SEXP from);
SEXP tCRsparse_as_RCsparse(SEXP from);

SEXP Csparse_is_diagonal(SEXP obj);
SEXP Rsparse_is_diagonal(SEXP obj);
SEXP Tsparse_is_diagonal(SEXP obj);
SEXP Csparse_is_triangular(SEXP obj, SEXP upper);
SEXP Rsparse_is_triangular(SEXP obj, SEXP upper);
SEXP Tsparse_is_triangular(SEXP obj, SEXP upper);
SEXP Csparse_is_symmetric(SEXP obj, SEXP checkDN);
SEXP Rsparse_is_symmetric(SEXP obj, SEXP checkDN);
#if 0 /* unimplemented ... currently going via CsparseMatrix */
SEXP Tsparse_is_symmetric(SEXP obj, SEXP checkDN);
#endif

SEXP CRsparse_colSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse);
SEXP CRsparse_rowSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse);

#define SPARSE_CASES(_SEXPTYPE_, _DO_) \
do { \
	switch (_SEXPTYPE_) { \
	case LGLSXP: \
		_DO_(int, LOGICAL, 0, 1, ISNZ_LOGICAL); \
		break; \
	case INTSXP: \
		_DO_(int, INTEGER, 0, 1, ISNZ_INTEGER); \
		break; \
	case REALSXP: \
		_DO_(double, REAL, 0.0, 1.0, ISNZ_REAL); \
		break; \
	case CPLXSXP: \
		_DO_(Rcomplex, COMPLEX, Matrix_zzero, Matrix_zone, ISNZ_COMPLEX); \
		break; \
	default: \
		break; \
	} \
} while (0)

#endif
