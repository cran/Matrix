#include "Lapack-etc.h"
#include "cs-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "factorizations.h"

/* defined in ./attrib.c : */
SEXP get_factor(SEXP, const char *);
void set_factor(SEXP, const char *, SEXP);

static
SEXP dgeMatrix_trf_(SEXP obj, int warn)
{
	SEXP val = PROTECT(newObject("denseLU")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	SET_SLOT(val, Matrix_DimSym, dim);
	SET_SLOT(val, Matrix_DimNamesSym, dimnames);
	if (r > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, r)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
		int *pperm = INTEGER(perm), info;
#ifdef MATRIX_ENABLE_ZMATRIX
		if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
		Matrix_memcpy(py, px, XLENGTH(y), sizeof(Rcomplex));
		F77_CALL(zgetrf)(&m, &n, py, &m, pperm, &info);
		ERROR_LAPACK_2(zgetrf, info, warn, U);
		} else {
#endif
		double *px = REAL(x), *py = REAL(y);
		Matrix_memcpy(py, px, XLENGTH(y), sizeof(double));
		F77_CALL(dgetrf)(&m, &n, py, &m, pperm, &info);
		ERROR_LAPACK_2(dgetrf, info, warn, U);
#ifdef MATRIX_ENABLE_ZMATRIX
		}
#endif
		SET_SLOT(val, Matrix_permSym, perm);
		SET_SLOT(val, Matrix_xSym, y);
		UNPROTECT(3); /* y, x, perm */
	}
	UNPROTECT(3); /* dimnames, dim, val */
	return val;
}

static
SEXP dsyMatrix_trf_(SEXP obj, int warn)
{
	SEXP val = PROTECT(newObject("BunchKaufman")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int n = INTEGER(dim)[1];
	char ul = *CHAR(STRING_ELT(uplo, 0));
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, n)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
		int *pperm = INTEGER(perm), info, lwork = -1;
#ifdef MATRIX_ENABLE_ZMATRIX
		if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x), *py = COMPLEX(y), tmp, *work;
		Matrix_memset(py, 0, XLENGTH(y), sizeof(Rcomplex));
		F77_CALL(zlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
		F77_CALL(zsytrf)(&ul, &n, py, &n, pperm, &tmp, &lwork, &info FCONE);
		lwork = (int) tmp.r;
		Matrix_Calloc(work, lwork, Rcomplex);
		F77_CALL(zsytrf)(&ul, &n, py, &n, pperm, work, &lwork, &info FCONE);
		Matrix_Free(work, lwork);
		ERROR_LAPACK_2(zsytrf, info, warn, D);
		} else {
#endif
		double *px = REAL(x), *py = REAL(y), tmp, *work;
		Matrix_memset(py, 0, XLENGTH(y), sizeof(double));
		F77_CALL(dlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
		F77_CALL(dsytrf)(&ul, &n, py, &n, pperm, &tmp, &lwork, &info FCONE);
		lwork = (int) tmp;
		Matrix_Calloc(work, lwork, double);
		F77_CALL(dsytrf)(&ul, &n, py, &n, pperm, work, &lwork, &info FCONE);
		Matrix_Free(work, lwork);
		ERROR_LAPACK_2(dsytrf, info, warn, D);
#ifdef MATRIX_ENABLE_ZMATRIX
		}
#endif
		SET_SLOT(val, Matrix_permSym, perm);
		SET_SLOT(val, Matrix_xSym, y);
		UNPROTECT(3); /* y, x, perm */
	}
	UNPROTECT(4); /* uplo, dimnames, dim, val */
	return val;
}

static
SEXP dspMatrix_trf_(SEXP obj, int warn)
{
	SEXP val = PROTECT(newObject("pBunchKaufman")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int n = INTEGER(dim)[1];
	char ul = *CHAR(STRING_ELT(uplo, 0));
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, n)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
		int *pperm = INTEGER(perm), info;
#ifdef MATRIX_ENABLE_ZMATRIX
		if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
		Matrix_memcpy(py, px, XLENGTH(y), sizeof(Rcomplex));
		F77_CALL(zsptrf)(&ul, &n, py, pperm, &info FCONE);
		ERROR_LAPACK_2(zsptrf, info, warn, D);
		} else {
#endif
		double *px = REAL(x), *py = REAL(y);
		Matrix_memcpy(py, px, XLENGTH(y), sizeof(double));
		F77_CALL(dsptrf)(&ul, &n, py, pperm, &info FCONE);
		ERROR_LAPACK_2(dsptrf, info, warn, D);
#ifdef MATRIX_ENABLE_ZMATRIX
		}
#endif
		SET_SLOT(val, Matrix_permSym, perm);
		SET_SLOT(val, Matrix_xSym, y);
		UNPROTECT(3); /* y, x, perm */
	}
	UNPROTECT(4); /* uplo, dimnames, dim, val */
	return val;
}

static
SEXP dpoMatrix_trf_(SEXP obj, int warn, int pivot, double tol)
{
	SEXP val = PROTECT(newObject("Cholesky")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int n = INTEGER(dim)[1];
	char ul = *CHAR(STRING_ELT(uplo, 0));
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
		int info;
#ifdef MATRIX_ENABLE_ZMATRIX
		if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
		Matrix_memset(py, 0, XLENGTH(y), sizeof(Rcomplex));
		F77_CALL(zlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
		if (!pivot) {
		F77_CALL(zpotrf)(&ul, &n, py, &n, &info FCONE);
		ERROR_LAPACK_3(zpotrf, info, warn, 6);
		} else {
		SEXP perm = PROTECT(allocVector(INTSXP, n));
		int *pperm = INTEGER(perm), rank;
		Rcomplex *work = (Rcomplex *) R_alloc((size_t) 2 * n, sizeof(Rcomplex));
		F77_CALL(zpstrf)(&ul, &n, py, &n, pperm, &rank, &tol, work, &info FCONE);
		ERROR_LAPACK_4(zpstrf, info, rank, warn);
		if (info > 0) {
			int j, d = n - rank;
			py += (R_xlen_t) rank * n + rank;
			for (j = rank; j < n; ++j) {
				Matrix_memset(py, 0, d, sizeof(Rcomplex));
				py += n;
			}
		}
		SET_SLOT(val, Matrix_permSym, perm);
		UNPROTECT(1); /* perm */
		}
		} else {
#endif
		double *px = REAL(x), *py = REAL(y);
		Matrix_memset(py, 0, XLENGTH(y), sizeof(double));
		F77_CALL(dlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
		if (!pivot) {
		F77_CALL(dpotrf)(&ul, &n, py, &n, &info FCONE);
		ERROR_LAPACK_3(dpotrf, info, warn, 6);
		} else {
		SEXP perm = PROTECT(allocVector(INTSXP, n));
		int *pperm = INTEGER(perm), rank;
		double *work = (double *) R_alloc((size_t) 2 * n, sizeof(double));
		F77_CALL(dpstrf)(&ul, &n, py, &n, pperm, &rank, &tol, work, &info FCONE);
		ERROR_LAPACK_4(dpstrf, info, rank, warn);
		if (info > 0) {
			int j, d = n - rank;
			py += (R_xlen_t) rank * n + rank;
			for (j = rank; j < n; ++j) {
				Matrix_memset(py, 0, d, sizeof(double));
				py += n;
			}
		}
		SET_SLOT(val, Matrix_permSym, perm);
		UNPROTECT(1); /* perm */
		}
#ifdef MATRIX_ENABLE_ZMATRIX
		}
#endif
		SET_SLOT(val, Matrix_xSym, y);
		UNPROTECT(2); /* y, x */
	}
	UNPROTECT(4); /* uplo, dimnames, dim, val */
	return val;
}

static
SEXP dppMatrix_trf_(SEXP obj, int warn)
{
	SEXP val = PROTECT(newObject("pCholesky")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int n = INTEGER(dim)[1];
	char ul = *CHAR(STRING_ELT(uplo, 0));
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
		int info;
#ifdef MATRIX_ENABLE_ZMATRIX
		if (TYPEOF(x) == CPLXSXP) {
			Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
			Matrix_memcpy(py, px, XLENGTH(y), sizeof(Rcomplex));
			F77_CALL(zpptrf)(&ul, &n, py, &info FCONE);
			ERROR_LAPACK_3(zpptrf, info, warn, 6);
		} else {
#endif
			double *px = REAL(x), *py = REAL(y);
			Matrix_memcpy(py, px, XLENGTH(y), sizeof(double));
			F77_CALL(dpptrf)(&ul, &n, py, &info FCONE);
			ERROR_LAPACK_3(dpptrf, info, warn, 6);
#ifdef MATRIX_ENABLE_ZMATRIX
		}
#endif
		SET_SLOT(val, Matrix_xSym, y);
		UNPROTECT(2); /* y, x */
	}
	UNPROTECT(4); /* uplo, dimnames, dim, val */
	return val;
}

SEXP dgeMatrix_trf(SEXP obj, SEXP warn)
{
	SEXP val = get_factor(obj, "denseLU");
	if (isNull(val)) {
		PROTECT(val = dgeMatrix_trf_(obj, asInteger(warn)));
		set_factor(obj, "denseLU", val);
		UNPROTECT(1);
	}
	return val;
}

SEXP dsyMatrix_trf(SEXP obj, SEXP warn)
{
	SEXP val = get_factor(obj, "BunchKaufman");
	if (isNull(val)) {
		PROTECT(val = dsyMatrix_trf_(obj, asInteger(warn)));
		set_factor(obj, "BunchKaufman", val);
		UNPROTECT(1);
	}
	return val;
}

SEXP dspMatrix_trf(SEXP obj, SEXP warn)
{
	SEXP val = get_factor(obj, "pBunchKaufman");
	if (isNull(val)) {
		PROTECT(val = dspMatrix_trf_(obj, asInteger(warn)));
		set_factor(obj, "pBunchKaufman", val);
		UNPROTECT(1);
	}
	return val;
}

SEXP dpoMatrix_trf(SEXP obj, SEXP warn, SEXP pivot, SEXP tol)
{
	int pivot_ = asLogical(pivot);
	SEXP val = get_factor(obj, (pivot_) ? "Cholesky~" : "Cholesky");
	if (isNull(val)) {
		double tol_ = asReal(tol);
		PROTECT(val = dpoMatrix_trf_(obj, asInteger(warn), pivot_, tol_));
		set_factor(obj, (pivot_) ? "Cholesky~" : "Cholesky", val);
		UNPROTECT(1);
	}
	return val;
}

SEXP dppMatrix_trf(SEXP obj, SEXP warn)
{
	SEXP val = get_factor(obj, "pCholesky");
	if (isNull(val)) {
		PROTECT(val = dppMatrix_trf_(obj, asInteger(warn)));
		set_factor(obj, "pCholesky", val);
		UNPROTECT(1);
	}
	return val;
}

#define DO_FREE(_A_, _S_, _N_) \
do { \
	if (!(_A_)) \
		_A_ = Matrix_cs_spfree(_A_); \
	if (!(_S_)) \
		_S_ = Matrix_cs_sfree (_S_); \
	if (!(_N_)) \
		_N_ = Matrix_cs_nfree (_N_); \
} while (0)

#define DO_SORT(_A_) \
do { \
	Matrix_cs_dropzeros(_A_); \
	T = Matrix_cs_transpose(_A_, 1); \
	if (!T) { \
		DO_FREE(T, *S, *N); \
		return 0; \
	} \
	_A_ = Matrix_cs_spfree(_A_); \
	_A_ = Matrix_cs_transpose(T, 1); \
	if (!(_A_)) { \
		DO_FREE(T, *S, *N); \
		return 0; \
	} \
	T = Matrix_cs_spfree(T); \
} while (0)

static
int dgCMatrix_trf_(const Matrix_cs *A, Matrix_css **S, Matrix_csn **N,
                   int order, double tol)
{
	Matrix_cs *T = NULL;
	if (!(*S = Matrix_cs_sqr(order, A, 0)) ||
	    !(*N = Matrix_cs_lu(A, *S, tol))) {
		DO_FREE(T, *S, *N);
		return 0;
	}
	DO_SORT((*N)->L);
	DO_SORT((*N)->U);
	return 1;
}

SEXP dgCMatrix_trf(SEXP obj, SEXP order, SEXP tol, SEXP doError)
{
	double tol_ = asReal(tol);
	if (ISNAN(tol_))
		error(_("'%s' is not a number"), "tol");

	int order_ = asInteger(order);
	if (order_ == NA_INTEGER)
		order_ = (tol_ == 1.0) ? 2 : 1;
	else if (order_ < 0 || order_ > 3)
		order_ = 0;

	SEXP val = get_factor(obj, (order_) ? "sparseLU~" : "sparseLU");
	if (!isNull(val))
		return val;
	PROTECT(val = newObject("sparseLU"));

	Matrix_cs *A = M2CXS(obj, 1);
	MCS_XTYPE_SET(A->xtype);

	Matrix_css *S = NULL;
	Matrix_csn *N = NULL;
	int *P = NULL;

	if (A->m != A->n)
		error(_("LU factorization of m-by-n %s requires m == n"),
		      ".gCMatrix");
	if (!dgCMatrix_trf_(A, &S, &N, order_, tol_) ||
	    !(P = Matrix_cs_pinv(N->pinv, A->m))) {
		if (!P) {
			S = Matrix_cs_sfree(S);
			N = Matrix_cs_nfree(N);
		}
		if (asLogical(doError))
			error(_("LU factorization of %s failed: out of memory or near-singular"),
			      ".gCMatrix");
		/* Defensive code will check with is(., "sparseLU") : */
		UNPROTECT(1); /* val */
		return ScalarLogical(NA_LOGICAL);
	}

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	SET_SLOT(val, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(val, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP L = PROTECT(CXS2M(N->L, 1, 't')),
		U = PROTECT(CXS2M(N->U, 1, 't')),
		uplo = PROTECT(mkString("L"));
	SET_SLOT(L, Matrix_uploSym, uplo);
	SET_SLOT(val, Matrix_LSym, L);
	SET_SLOT(val, Matrix_USym, U);
	UNPROTECT(3); /* uplo, U, L */

	SEXP p = PROTECT(allocVector(INTSXP, A->m));
	Matrix_memcpy(INTEGER(p), P, A->m, sizeof(int));
	SET_SLOT(val, Matrix_pSym, p);
	UNPROTECT(1); /* p */
	if (order_ > 0) {
		SEXP q = PROTECT(allocVector(INTSXP, A->n));
		Matrix_memcpy(INTEGER(q), S->q, A->n, sizeof(int));
		SET_SLOT(val, Matrix_qSym, q);
		UNPROTECT(1); /* q */
	}

	S = Matrix_cs_sfree(S);
	N = Matrix_cs_nfree(N);
	P = Matrix_cs_free(P);

	set_factor(obj, (order_) ? "sparseLU~" : "sparseLU", val);
	UNPROTECT(1); /* val */
	return val;
}

static
int dgCMatrix_orf_(const Matrix_cs *A, Matrix_css **S, Matrix_csn **N,
                   int order)
{
	Matrix_cs *T = NULL;
	if (!(*S = Matrix_cs_sqr(order, A, 1)) ||
	    !(*N = Matrix_cs_qr(A, *S))) {
		DO_FREE(T, *S, *N);
		return 0;
	}
	DO_SORT((*N)->L);
	DO_SORT((*N)->U);
	return 1;
}

SEXP dgCMatrix_orf(SEXP obj, SEXP order, SEXP doError)
{
	int order_ = asInteger(order);
	if (order_ < 0 || order_ > 3)
		order_ = 0;

	SEXP val = get_factor(obj, (order_) ? "sparseQR~" : "sparseQR");
	if (!isNull(val))
		return val;
	PROTECT(val = newObject("sparseQR"));

	Matrix_cs *A = M2CXS(obj, 1);
	MCS_XTYPE_SET(A->xtype);

	Matrix_css *S = NULL;
	Matrix_csn *N = NULL;
	int *P = NULL;

	if (A->m < A->n)
		error(_("QR factorization of m-by-n %s requires m >= n"),
		      ".gCMatrix");
	if (!dgCMatrix_orf_(A, &S, &N, order_) ||
	    !(P = Matrix_cs_pinv(S->pinv, S->m2))) {
		if (!P) {
			S = Matrix_cs_sfree(S);
			N = Matrix_cs_nfree(N);
		}
		if (asLogical(doError))
			error(_("QR factorization of %s failed: out of memory"),
			      ".gCMatrix");
		/* Defensive code will check with is(., "sparseQR") : */
		UNPROTECT(1); /* val */
		return ScalarLogical(NA_LOGICAL);
	}

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	SET_SLOT(val, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(val, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP V = PROTECT(CXS2M(N->L, 1, 'g')),
		R = PROTECT(CXS2M(N->U, 1, 'g'));
	SET_SLOT(val, Matrix_VSym, V);
	SET_SLOT(val, Matrix_RSym, R);
	UNPROTECT(2); /* R, V */

	SEXP beta = PROTECT(allocVector(REALSXP, A->n));
	Matrix_memcpy(REAL(beta), N->B, A->n, sizeof(double));
	SET_SLOT(val, Matrix_betaSym, beta);
	UNPROTECT(1); /* beta */

	SEXP p = PROTECT(allocVector(INTSXP, S->m2));
	Matrix_memcpy(INTEGER(p), P, S->m2, sizeof(int));
	SET_SLOT(val, Matrix_pSym, p);
	UNPROTECT(1); /* p */
	if (order_ > 0) {
		SEXP q = PROTECT(allocVector(INTSXP, A->n));
		Matrix_memcpy(INTEGER(q), S->q, A->n, sizeof(int));
		SET_SLOT(val, Matrix_qSym, q);
		UNPROTECT(1); /* q */
	}

	S = Matrix_cs_sfree(S);
	N = Matrix_cs_nfree(N);
	P = Matrix_cs_free(P);

	set_factor(obj, (order_) ? "sparseQR~" : "sparseQR", val);
	UNPROTECT(1); /* val */
	return val;
}

#undef DO_FREE
#undef DO_SORT

static
int dpCMatrix_trf_(cholmod_sparse *A, cholmod_factor **L,
                   int perm, int ldl, int super, double mult)
{
	/* defined in ./chm_common.c : */
	void R_cholmod_common_envget(void);
	void R_cholmod_common_envset(void);

	R_cholmod_common_envset();

	if (*L == NULL) {
		if (perm == 0) {
			c.nmethods = 1;
			c.method[0].ordering = CHOLMOD_NATURAL;
			c.postorder = 0;
		}

		c.supernodal = (super == NA_LOGICAL) ? CHOLMOD_AUTO :
			((super != 0) ? CHOLMOD_SUPERNODAL : CHOLMOD_SIMPLICIAL);

		*L = cholmod_analyze(A, &c);
	}

	if (super == NA_LOGICAL)
		super = (*L)->is_super;
	if (super != 0)
		ldl = 0;

	c.final_asis = 0;
	c.final_super = super != 0;
	c.final_ll = ldl == 0;
	c.final_pack = 1;
	c.final_monotonic = 1;

	double beta[2];
	beta[0] = mult;
	beta[1] = 0.0;
	int res = cholmod_factorize_p(A, beta, NULL, 0, *L, &c);

	R_cholmod_common_envget();

	return res;
}

SEXP dpCMatrix_trf(SEXP obj,
                   SEXP perm, SEXP ldl, SEXP super, SEXP mult)
{
	int perm_ = asLogical(perm), ldl_ = asLogical(ldl),
		super_ = asLogical(super);
	double mult_ = asReal(mult);
	if (!R_FINITE(mult_))
		error(_("'%s' is not a number or not finite"), "mult");

	SEXP trf = R_NilValue;
	char nm[] = "spdCholesky";
	if (perm_)
		nm[1] = 'P';
	if (super_ != NA_LOGICAL && super_ != 0)
		ldl_ = 0;
	if (super_ == NA_LOGICAL || super_ == 0) {
		if (ldl_)
			nm[2] = 'D';
		trf = get_factor(obj, nm);
	}
	if (isNull(trf) && (super_ == NA_LOGICAL || super_ != 0)) {
		nm[0] = 'S';
		nm[2] = 'd';
		trf = get_factor(obj, nm);
	}

	int cached = !isNull(trf);
	if (cached && mult_ == 0.0)
		return trf;

	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(trf, &pid);
	cholmod_sparse *A = M2CHS(obj, 1);
	cholmod_factor *L = NULL;

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = *CHAR(STRING_ELT(uplo, 0));
	A->stype = (ul == 'U') ? 1 : -1;

	if (cached) {
		L = M2CHF(trf, 1);
		L = cholmod_copy_factor(L, &c);
		dpCMatrix_trf_(A, &L, perm_, ldl_, super_, mult_);
	} else {
		dpCMatrix_trf_(A, &L, perm_, ldl_, super_, mult_);
		if (super_ == NA_LOGICAL) {
			nm[0] = (L->is_super) ? 'S' : 's';
			nm[2] = (L->is_ll   ) ? 'd' : 'D';
		}
	}
	REPROTECT(trf = CHF2M(L, 1), pid);
	cholmod_free_factor(&L, &c);

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	set_symmetrized_DimNames(trf, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (!cached && mult_ == 0.0)
		set_factor(obj, nm, trf);
	UNPROTECT(1); /* trf */
	return trf;
}

SEXP BunchKaufman_expand(SEXP obj, SEXP packed)
{
	SEXP P_ = PROTECT(newObject("pMatrix")),
		T_ = PROTECT(newObject("dtCMatrix")),
		D_ = PROTECT(newObject("dsCMatrix")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int i, j, s, n = INTEGER(dim)[0];
	R_xlen_t n1a = (R_xlen_t) n + 1;
	if (n > 0) {
		SET_SLOT(P_, Matrix_DimSym, dim);
		SET_SLOT(T_, Matrix_DimSym, dim);
		SET_SLOT(D_, Matrix_DimSym, dim);
	}
	UNPROTECT(1); /* dim */

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	if (!upper) {
		SET_SLOT(T_, Matrix_uploSym, uplo);
		SET_SLOT(D_, Matrix_uploSym, uplo);
	}
	UNPROTECT(1); /* uplo */

	SEXP diag = PROTECT(mkString("U"));
	SET_SLOT(T_, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
		D_p = PROTECT(allocVector(INTSXP, n1a));
	int *ppivot = INTEGER(pivot), *D_pp = INTEGER(D_p),
		b = n, dp = (upper) ? 1 : 2;
	D_pp[0] = 0;
	j = 0;
	while (j < n) {
		if (ppivot[j] > 0) {
			D_pp[j+1] = D_pp[j] + 1;
			j += 1;
		} else {
			D_pp[j+1] = D_pp[j] + dp;
			D_pp[j+2] = D_pp[j] + 3;
			j += 2;
			--b;
		}
	}
	SET_SLOT(D_, Matrix_pSym, D_p);
	UNPROTECT(1); /* D_p */

	SEXP P, P_perm, T, T_p, T_i, T_x,
		D_i = PROTECT(allocVector(INTSXP, D_pp[n])),
		D_x = PROTECT(allocVector(REALSXP, D_pp[n])),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int *P_pperm, *T_pp, *T_pi, *D_pi = INTEGER(D_i);
	double *T_px, *D_px = REAL(D_x), *px = REAL(x);

	int unpacked = !asLogical(packed);

	R_xlen_t len = (R_xlen_t) 2 * b + 1, k = (upper) ? len - 2 : 0;
	SEXP res = PROTECT(allocVector(VECSXP, len));

	j = 0;
	while (b--) {
		s = (ppivot[j] > 0) ? 1 : 2;
		dp = (upper) ? j : n - j - s;

		PROTECT(P = duplicate(P_));
		PROTECT(P_perm = allocVector(INTSXP, n));
		PROTECT(T = duplicate(T_));
		PROTECT(T_p = allocVector(INTSXP, n1a));
		PROTECT(T_i = allocVector(INTSXP, (R_xlen_t) s * dp));
		PROTECT(T_x = allocVector(REALSXP, (R_xlen_t) s * dp));

		P_pperm = INTEGER(P_perm);
		T_pp = INTEGER(T_p);
		T_pi = INTEGER(T_i);
		T_px = REAL(T_x);
		T_pp[0] = 0;

		for (i = 0; i < j; ++i) {
			T_pp[i+1] = 0;
			P_pperm[i] = i + 1;
		}
		for (i = j; i < j+s; ++i) {
			T_pp[i+1] = T_pp[i] + dp;
			P_pperm[i] = i + 1;
		}
		for (i = j+s; i < n; ++i) {
			T_pp[i+1] = T_pp[i];
			P_pperm[i] = i + 1;
		}

		if (s == 1) {
			P_pperm[j] = ppivot[j];
			P_pperm[ppivot[j]-1] = j + 1;
		} else if (upper) {
			P_pperm[j] = -ppivot[j];
			P_pperm[-ppivot[j]-1] = j + 1;
		} else {
			P_pperm[j+1] = -ppivot[j];
			P_pperm[-ppivot[j]-1] = j + 2;
		}

		if (upper) {
			for (i = 0; i < j; ++i) {
				*(T_pi++) = i;
				*(T_px++) = *(px++);
			}
			*(D_pi++) = j;
			*(D_px++) = *(px++);
			++j;
			if (unpacked)
				px += n - j;
			if (s == 2) {
				for (i = 0; i < j-1; ++i) {
					*(T_pi++) = i;
					*(T_px++) = *(px++);
				}
				*(D_pi++) = j - 1;
				*(D_pi++) = j;
				*(D_px++) = *(px++);
				*(D_px++) = *(px++);
				++j;
				if (unpacked)
					px += n - j;
			}
		} else {
			if (s == 2) {
				*(D_pi++) = j;
				*(D_pi++) = j + 1;
				*(D_px++) = *(px++);
				*(D_px++) = *(px++);
				for (i = j+2; i < n; ++i) {
					*(T_pi++) = i;
					*(T_px++) = *(px++);
				}
				++j;
				if (unpacked)
					px += j;
			}
			*(D_pi++) = j;
			*(D_px++) = *(px++);
			for (i = j+1; i < n; ++i) {
				*(T_pi++) = i;
				*(T_px++) = *(px++);
			}
			++j;
			if (unpacked)
				px += j;
		}

		SET_SLOT(P, Matrix_permSym, P_perm);
		SET_SLOT(T, Matrix_pSym, T_p);
		SET_SLOT(T, Matrix_iSym, T_i);
		SET_SLOT(T, Matrix_xSym, T_x);

		if (upper) {
			SET_VECTOR_ELT(res, k-1, P);
			SET_VECTOR_ELT(res, k  , T);
			k -= 2;
		} else {
			SET_VECTOR_ELT(res, k  , P);
			SET_VECTOR_ELT(res, k+1, T);
			k += 2;
		}
		UNPROTECT(6); /* T_x, T_i, T_p, T, P_perm, P */
	}

	SET_SLOT(D_, Matrix_iSym, D_i);
	SET_SLOT(D_, Matrix_xSym, D_x);
	SET_VECTOR_ELT(res, len-1, D_);

	UNPROTECT(8); /* res, x, D_x, D_i, pivot, D_, T_, P_ */
	return res;
}

SEXP CHMfactor_diag_get(SEXP obj, SEXP square)
{
	cholmod_factor *L = M2CHF(obj, 1);
	int n = (int) L->n, square_ = asLogical(square);
	SEXP y = PROTECT(allocVector(REALSXP, n));
	double *py = REAL(y);
	if (L->is_super) {
		int k, j, nc,
			nsuper = (int) L->nsuper,
			*psuper = (int *) L->super,
			*ppi = (int *) L->pi,
			*ppx = (int *) L->px;
		double *px = (double *) L->x, *px_;
		R_xlen_t nr1a;
		for (k = 0; k < nsuper; ++k) {
			nc = psuper[k+1] - psuper[k];
			nr1a = (R_xlen_t) (ppi[k+1] - ppi[k]) + 1;
			px_ = px + ppx[k];
			for (j = 0; j < nc; ++j) {
				*py = *px_;
				if (square_)
					*py *= *py;
				++py;
				px_ += nr1a;
			}
		}
	} else {
		square_ = square_ && L->is_ll;
		int j, *pp = (int *) L->p;
		double *px = (double *) L->x;
		for (j = 0; j < n; ++j) {
			*py = px[pp[j]];
			if (square_)
				*py *= *py;
			++py;
		}
	}
	UNPROTECT(1);
	return y;
}

SEXP CHMfactor_update(SEXP obj, SEXP parent, SEXP mult)
{
	/* defined in ./objects.c : */
	char Matrix_shape(SEXP);

	double mult_ = asReal(mult);
	if (!R_FINITE(mult_))
		error(_("'%s' is not a number or not finite"), "mult");

	cholmod_factor *L = cholmod_copy_factor(M2CHF(obj, 1), &c);
	cholmod_sparse *A = M2CHS(parent, 1);
	if (Matrix_shape(parent) == 's') {
		SEXP uplo = GET_SLOT(parent, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		A->stype = (ul == 'U') ? 1 : -1;
	}

	dpCMatrix_trf_(A, &L, 0, !L->is_ll, L->is_super, mult_);

	SEXP res = PROTECT(CHF2M(L, 1));
	cholmod_free_factor(&L, &c);

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1);

	UNPROTECT(1);
	return res;
}

SEXP CHMfactor_updown(SEXP obj, SEXP parent, SEXP update)
{
	/* defined in ./objects.c : */
	char Matrix_shape(SEXP);

	cholmod_factor *L = cholmod_copy_factor(M2CHF(obj, 1), &c);
	cholmod_sparse *A = M2CHS(parent, 1);
	if (Matrix_shape(parent) == 's') {
		SEXP uplo = GET_SLOT(parent, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		A->stype = (ul == 'U') ? 1 : -1;
	}

	cholmod_updown(asLogical(update) != 0, A, L, &c);

	SEXP res = PROTECT(CHF2M(L, 1));
	cholmod_free_factor(&L, &c);

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1);

	UNPROTECT(1);
	return res;
}
