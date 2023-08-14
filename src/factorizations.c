#include <Rmath.h> /* math.h, logspace_add, logspace_sub */
#include "factorizations.h"

static cs *dgC2cs(SEXP obj)
{
	cs *A = (cs *) R_alloc(1, sizeof(cs));
	memset(A, 0, sizeof(cs));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	A->nzmax = LENGTH(i);
	A->m = INTEGER(dim)[0];
	A->n = INTEGER(dim)[1];
	A->p = INTEGER(p);
	A->i = INTEGER(i);
	A->x = REAL(x);
	A->nz = -1;
	UNPROTECT(4);
	return A;
}

static SEXP cs2dgC(const cs *A, const char *cl)
{
	int nnz = ((int *) A->p)[A->n];
	R_xlen_t np1 = (R_xlen_t) A->n + 1;
	SEXP obj = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(allocVector(INTSXP, np1)),
		i = PROTECT(allocVector(INTSXP, nnz)),
		x = PROTECT(allocVector(REALSXP, nnz));
	INTEGER(dim)[0] = A->m;
	INTEGER(dim)[1] = A->n;
	Matrix_memcpy(INTEGER(p), A->p, np1, sizeof(int));
	Matrix_memcpy(INTEGER(i), A->i, nnz, sizeof(int));
	Matrix_memcpy(REAL(x), A->x, nnz, sizeof(double));
	SET_SLOT(obj, Matrix_pSym, p);
	SET_SLOT(obj, Matrix_iSym, i);
	SET_SLOT(obj, Matrix_xSym, x);
	UNPROTECT(5);
	return obj;
}

static cholmod_sparse *dgC2cholmod(SEXP obj)
{
	cholmod_sparse *A = (cholmod_sparse *) R_alloc(1, sizeof(cholmod_sparse));
	memset(A, 0, sizeof(cholmod_sparse));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	A->nrow = (size_t) INTEGER(dim)[0];
	A->ncol = (size_t) INTEGER(dim)[1];
	A->p = INTEGER(p);
	A->i = INTEGER(i);
	A->x = REAL(x);
	A->nzmax = (size_t) ((int *) A->p)[A->ncol];
	A->stype = 0;
	A->itype = CHOLMOD_INT;
	A->xtype = CHOLMOD_REAL;
	A->dtype = CHOLMOD_DOUBLE;
	A->sorted = 1;
	A->packed = 1;
	UNPROTECT(4);
	return A;
}

static SEXP cholmod2dgC(cholmod_sparse *A, const char *cl)
{
	if (A->itype != CHOLMOD_INT ||
	    A->xtype != CHOLMOD_REAL ||
	    A->dtype != CHOLMOD_DOUBLE)
		error(_("wrong '%s' or '%s' or '%s'"), "itype", "xtype", "dtype");
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		error(_("dimensions cannot exceed %s"), "2^31-1");
	if (!A->sorted || !A->packed || A->stype != 0)
		cholmod_sort(A, &c);
	int m = (int) A->nrow, n = (int) A->ncol,
		nnz = ((int *) A->p)[A->ncol];
	R_xlen_t n1a = (R_xlen_t) n + 1;
	SEXP obj = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(allocVector(INTSXP, n1a)),
		i = PROTECT(allocVector(INTSXP, nnz)),
		x = PROTECT(allocVector(REALSXP, nnz));
	INTEGER(dim)[0] = m;
	INTEGER(dim)[1] = n;
	Matrix_memcpy(INTEGER(p), A->p, n1a, sizeof(int));
	Matrix_memcpy(INTEGER(i), A->i, nnz, sizeof(int));
	Matrix_memcpy(REAL(x), A->x, nnz, sizeof(double));
	SET_SLOT(obj, Matrix_pSym, p);
	SET_SLOT(obj, Matrix_iSym, i);
	SET_SLOT(obj, Matrix_xSym, x);
	UNPROTECT(5);
	return obj;
}

static cholmod_dense *dge2cholmod(SEXP obj)
{
	cholmod_dense *A = (cholmod_dense *) R_alloc(1, sizeof(cholmod_dense));
	memset(A, 0, sizeof(cholmod_dense));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	A->nzmax = (size_t) XLENGTH(x);
	A->nrow = (size_t) INTEGER(dim)[0];
	A->ncol = (size_t) INTEGER(dim)[1];
	A->d = A->nrow;
	A->x = REAL(x);
	A->xtype = CHOLMOD_REAL;
	A->dtype = CHOLMOD_DOUBLE;
	UNPROTECT(2);
	return A;
}

static SEXP cholmod2dge(const cholmod_dense *A, const char *cl)
{
	if (A->xtype != CHOLMOD_REAL || A->dtype != CHOLMOD_DOUBLE)
		error(_("wrong '%s' or '%s'"), "xtype", "dtype");
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		error(_("dimensions cannot exceed %s"), "2^31-1");
	int m = (int) A->nrow, n = (int) A->ncol;
	if ((double) m * n > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	size_t d = A->d;
	R_xlen_t mn = (R_xlen_t) m * n;
	SEXP obj = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		x = PROTECT(allocVector(REALSXP, mn));
	double *px = REAL(x), *py = (double *) A->x;
	INTEGER(dim)[0] = m;
	INTEGER(dim)[1] = n;
	if (d == m)
		Matrix_memcpy(px, py, mn, sizeof(double));
	else {
		int j;
		for (j = 0; j < n; ++j) {
			Matrix_memcpy(px, py, m, sizeof(double));
			px += m;
			py += d;
		}
	}
	SET_SLOT(obj, Matrix_xSym, x);
	UNPROTECT(3);
	return obj;
}

static cholmod_factor *mf2cholmod(SEXP obj)
{
	static const char *valid[] = { "dCHMsimpl", "dCHMsuper", "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		error(_("expected %s or %s"), "dCHMsimpl", "dCHMsuper");
	SEXP type = PROTECT(GET_SLOT(obj, install("type"))),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		colcount = PROTECT(GET_SLOT(obj, install("colcount"))),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int *ptype = INTEGER(type);
	cholmod_factor *L = (cholmod_factor *) R_alloc(1, sizeof(cholmod_factor));
	memset(L, 0, sizeof(cholmod_factor));

	L->ordering = ptype[0];
	L->is_super = ptype[2];

	L->n = (size_t) INTEGER(dim)[0];
	L->minor = L->n; /* FIXME: could be wrong for obj <- new(...) */
	L->ColCount = INTEGER(colcount);

	if (L->ordering != CHOLMOD_NATURAL) {
		SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
		L->Perm = INTEGER(perm);
		UNPROTECT(1);
	} else {
		/* cholmod_check_factor allows L->Perm == NULL,
		   but cholmod_copy_factor does not test, so it segfaults ...
		*/
		int j, n = (int) L->n, *Perm = (int *) R_alloc(L->n, sizeof(int));
		for (j = 0; j < n; ++j)
			Perm[j] = j;
		L->Perm = Perm;
	}

	L->itype = CHOLMOD_INT;
	L->xtype = CHOLMOD_REAL;
	L->dtype = CHOLMOD_DOUBLE;
	L->x = REAL(x);

	if (L->is_super) {
		SEXP super = PROTECT(GET_SLOT(obj, install("super"))),
			pi = PROTECT(GET_SLOT(obj, install("pi"))),
			px = PROTECT(GET_SLOT(obj, install("px"))),
			s = PROTECT(GET_SLOT(obj, install("s")));
		L->super = INTEGER(super);
		L->pi = INTEGER(pi);
		L->px = INTEGER(px);
		L->s = INTEGER(s);
		L->nsuper = (size_t) LENGTH(super) - 1;
		L->ssize = (size_t) ((int *) L->pi)[L->nsuper];
		L->xsize = (size_t) ((int *) L->px)[L->nsuper];
		L->maxcsize = (size_t) ptype[4];
		L->maxesize = (size_t) ptype[5];
		L->is_ll = 1;
		L->is_monotonic = 1;
		UNPROTECT(4);
	} else {
		SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			nz = PROTECT(GET_SLOT(obj, install("nz"))),
			nxt = PROTECT(GET_SLOT(obj, install("nxt"))),
			prv = PROTECT(GET_SLOT(obj, install("prv")));
		L->p = INTEGER(p);
		L->i = INTEGER(i);
		L->nz = INTEGER(nz);
		L->next = INTEGER(nxt);
		L->prev = INTEGER(prv);
		L->nzmax = (size_t) ((int *) L->p)[L->n];
		L->is_ll = ptype[1];
		L->is_monotonic = ptype[3];
		UNPROTECT(5);
	}

	UNPROTECT(4);
	return L;
}

static SEXP cholmod2mf(const cholmod_factor *L)
{
	if (L->itype != CHOLMOD_INT ||
	    L->xtype != CHOLMOD_REAL ||
	    L->dtype != CHOLMOD_DOUBLE)
		error(_("wrong '%s' or '%s' or '%s"), "itype", "xtype", "dtype");
	if (L->n > INT_MAX)
		error(_("dimensions cannot exceed %s"), "2^31-1");
	if (L->super) {
		if (L->maxcsize > INT_MAX)
			error(_("'%s' would overflow \"%s\""), "maxcsize", "integer");
	} else {
		if (L->n == INT_MAX)
			error(_("n+1 would overflow \"%s\""), "integer");
	}
	if (L->minor < L->n) {
		if (L->is_ll)
			error(_("leading principal minor of order %d is not positive"),
			      (int) L->minor + 1);
		else
			error(_("leading principal minor of order %d is zero"),
			      (int) L->minor + 1);
	}

	SEXP obj = PROTECT(NEW_OBJECT_OF_CLASS((L->is_super) ? "dCHMsuper" : "dCHMsimpl")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		type = PROTECT(allocVector(INTSXP, 6)),
		colcount = PROTECT(allocVector(INTSXP, L->n)),
		x = PROTECT(allocVector(REALSXP, (L->is_super) ? L->xsize : L->nzmax));
	int *pdim = INTEGER(dim), *ptype = INTEGER(type);
	pdim[0] = pdim[1] = (int) L->n;
	ptype[0] = L->ordering;
	ptype[1] = L->is_ll;
	ptype[2] = L->is_super;
	ptype[3] = L->is_monotonic;
	ptype[4] = (int) L->maxcsize;
	ptype[5] = (int) L->maxesize;
	Matrix_memcpy(INTEGER(colcount), L->ColCount, L->n, sizeof(int));
	Matrix_memcpy(REAL(x), L->x, XLENGTH(x), sizeof(double));
	SET_SLOT(obj, Matrix_DimSym, dim);
	SET_SLOT(obj, install("type"), type);
	SET_SLOT(obj, install("colcount"), colcount);
	SET_SLOT(obj, Matrix_xSym, x);

	if (L->ordering != CHOLMOD_NATURAL) {
		SEXP perm = PROTECT(allocVector(INTSXP, L->n));
		Matrix_memcpy(INTEGER(perm), L->Perm, L->n, sizeof(int));
		SET_SLOT(obj, Matrix_permSym, perm);
		UNPROTECT(1);
	}

	if (L->is_super) {
		SEXP super = PROTECT(allocVector(INTSXP, L->nsuper + 1)),
			pi = PROTECT(allocVector(INTSXP, L->nsuper + 1)),
			px = PROTECT(allocVector(INTSXP, L->nsuper + 1)),
			s = PROTECT(allocVector(INTSXP, L->ssize));
		Matrix_memcpy(INTEGER(super), L->super, L->nsuper + 1, sizeof(int));
		Matrix_memcpy(INTEGER(pi), L->pi, L->nsuper + 1, sizeof(int));
		Matrix_memcpy(INTEGER(px), L->px, L->nsuper + 1, sizeof(int));
		Matrix_memcpy(INTEGER(s), L->s, L->ssize, sizeof(int));
		SET_SLOT(obj, install("super"), super);
		SET_SLOT(obj, install("pi"), pi);
		SET_SLOT(obj, install("px"), px);
		SET_SLOT(obj, install("s"), s);
		UNPROTECT(4);
	} else {
		SEXP p = PROTECT(allocVector(INTSXP, L->n + 1)),
			i = PROTECT(allocVector(INTSXP, L->nzmax)),
			nz = PROTECT(allocVector(INTSXP, L->n)),
			nxt = PROTECT(allocVector(INTSXP, L->n + 2)),
			prv = PROTECT(allocVector(INTSXP, L->n + 2));
		Matrix_memcpy(INTEGER(p), L->p, L->n + 1, sizeof(int));
		Matrix_memcpy(INTEGER(i), L->i, L->nzmax, sizeof(int));
		Matrix_memcpy(INTEGER(nz), L->nz, L->n, sizeof(int));
		Matrix_memcpy(INTEGER(nxt), L->next, L->n + 2, sizeof(int));
		Matrix_memcpy(INTEGER(prv), L->prev, L->n + 2, sizeof(int));
		SET_SLOT(obj, Matrix_pSym, p);
		SET_SLOT(obj, Matrix_iSym, i);
		SET_SLOT(obj, install("nz"), nz);
		SET_SLOT(obj, install("nxt"), nxt);
		SET_SLOT(obj, install("prv"), prv);
		UNPROTECT(5);
	}

	UNPROTECT(5);
	return obj;
}

#define ERROR_LAPACK_1(_ROUTINE_, _INFO_) \
do { \
	if ((_INFO_) < 0) \
		error(_("LAPACK routine '%s': argument %d had illegal value"), \
		      #_ROUTINE_, -(_INFO_)); \
} while (0)

#define ERROR_LAPACK_2(_ROUTINE_, _INFO_, _WARN_, _LETTER_) \
do { \
	ERROR_LAPACK_1(_ROUTINE_, _INFO_); \
	if ((_INFO_) > 0 && (_WARN_) > 0) { \
		if (_WARN_ > 1) \
			error  (_("LAPACK routine '%s': matrix is exactly singular, %s[i,i]=0, i=%d"), \
			        #_ROUTINE_, #_LETTER_, (_INFO_)); \
		else \
			warning(_("LAPACK routine '%s': matrix is exactly singular, %s[i,i]=0, i=%d"), \
			        #_ROUTINE_, #_LETTER_, (_INFO_)); \
	} \
} while (0)

#define ERROR_LAPACK_3(_ROUTINE_, _INFO_, _WARN_, _NPROTECT_) \
do { \
	ERROR_LAPACK_1(_ROUTINE_, _INFO_); \
	if ((_INFO_) > 0 && (_WARN_) > 0) { \
		if (_WARN_ > 1) \
			error  (_("LAPACK routine '%s': leading principal minor of order %d is not positive"), \
			        #_ROUTINE_, (_INFO_)); \
		else { \
			warning(_("LAPACK routine '%s': leading principal minor of order %d is not positive"), \
			        #_ROUTINE_, (_INFO_)); \
			UNPROTECT(_NPROTECT_); \
			return ScalarInteger(_INFO_); \
		} \
	} \
} while (0)

#define ERROR_LAPACK_4(_ROUTINE_, _INFO_, _RANK_, _WARN_) \
	do { \
		ERROR_LAPACK_1(_ROUTINE_, _INFO_); \
		if ((_INFO_) > 0 && (_WARN_) > 0) { \
			if (_WARN_ > 1) \
				error  (_("LAPACK routine '%s': matrix is rank deficient or not positive definite, the _computed_ rank is %d"), \
				        #_ROUTINE_, (_RANK_)); \
			else \
				warning(_("LAPACK routine '%s': matrix is rank deficient or not positive definite, the _computed_ rank is %d"), \
				        #_ROUTINE_, (_RANK_)); \
		} \
	} while (0)

SEXP dgeMatrix_trf_(SEXP obj, int warn)
{
	SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("denseLU")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	int *pdim = INTEGER(dim), r = (pdim[0] < pdim[1]) ? pdim[0] : pdim[1];
	SET_SLOT(val, Matrix_DimSym, dim);
	SET_SLOT(val, Matrix_DimNamesSym, dimnames);
	if (r > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, r)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		x = duplicate(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
		int *pperm = INTEGER(perm), info;
		double *px = REAL(x);

		F77_CALL(dgetrf)(pdim, pdim + 1, px, pdim, pperm, &info);
		ERROR_LAPACK_2(dgetrf, info, warn, U);

		SET_SLOT(val, Matrix_permSym, perm);
		SET_SLOT(val, Matrix_xSym, x);
		UNPROTECT(2); /* x, perm */
	}
	UNPROTECT(3); /* dimnames, dim, val */
	return val;
}

SEXP dsyMatrix_trf_(SEXP obj, int warn)
{
	SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("BunchKaufman")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, n)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			y = PROTECT(allocVector(REALSXP, XLENGTH(x)));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		int *pperm = INTEGER(perm), lwork = -1, info;
		double *px = REAL(x), *py = REAL(y), tmp, *work;

		Matrix_memset(py, 0, XLENGTH(y), sizeof(double));
		F77_CALL(dlacpy)(&ul, pdim, pdim, px, pdim, py, pdim FCONE);
		F77_CALL(dsytrf)(&ul, pdim, py, pdim, pperm, &tmp, &lwork,
		                 &info FCONE);
		lwork = (int) tmp;
		Matrix_Calloc(work, lwork, double);
		F77_CALL(dsytrf)(&ul, pdim, py, pdim, pperm, work, &lwork,
		                 &info FCONE);
		Matrix_Free(work, lwork);
		ERROR_LAPACK_2(dsytrf, info, warn, D);

		SET_SLOT(val, Matrix_permSym, perm);
		SET_SLOT(val, Matrix_xSym, y);
		UNPROTECT(3); /* y, x, perm */
	}
	UNPROTECT(4); /* uplo, dimnames, dim, val */
	return val;
}

SEXP dspMatrix_trf_(SEXP obj, int warn)
{
	SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("pBunchKaufman")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, n)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		x = duplicate(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		int *pperm = INTEGER(perm), info;
		double *px = REAL(x);

		F77_CALL(dsptrf)(&ul, pdim, px, pperm, &info FCONE);
		ERROR_LAPACK_2(dsptrf, info, warn, D);

		SET_SLOT(val, Matrix_permSym, perm);
		SET_SLOT(val, Matrix_xSym, x);
		UNPROTECT(2); /* x, perm */
	}
	UNPROTECT(4); /* uplo, dimnames, dim, val */
	return val;
}

SEXP dpoMatrix_trf_(SEXP obj, int warn, int pivot, double tol)
{
	SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("Cholesky")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			y = PROTECT(allocVector(REALSXP, XLENGTH(x)));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		int info;
		double *px = REAL(x), *py = REAL(y);

		Matrix_memset(py, 0, XLENGTH(y), sizeof(double));
		F77_CALL(dlacpy)(&ul, pdim, pdim, px, pdim, py, pdim FCONE);

		if (pivot) {
			SEXP perm = PROTECT(allocVector(INTSXP, n));
			int *pperm = INTEGER(perm), rank;
			double *work = (double *) R_alloc((size_t) 2 * n, sizeof(double));

			F77_CALL(dpstrf)(&ul, pdim, py, pdim,
			                 pperm, &rank, &tol, work, &info FCONE);
			ERROR_LAPACK_4(dpstrf, info, rank, warn);

			if (info > 0) {
				/* Zero the trailing (n-rank)-by-(n-rank) principal
				   submatrix to guarantee that the result is valid
				   under the assumption that the factorized matrix
				   is positive semidefinite */
				int j;
				R_xlen_t len = (R_xlen_t) (n - rank);
				py += (R_xlen_t) rank * n + rank;
				for (j = rank; j < n; ++j) {
					Matrix_memset(py, 0, len, sizeof(double));
					py += n;
				}
			}

			SET_SLOT(val, Matrix_permSym, perm);
			UNPROTECT(1); /* perm */
		} else {
			F77_CALL(dpotrf)(&ul, pdim, py, pdim, &info FCONE);
			ERROR_LAPACK_3(dpotrf, info, warn, 6);
		}

		SET_SLOT(val, Matrix_xSym, y);
		UNPROTECT(2); /* y, x */
	}
	UNPROTECT(4); /* uplo, dimnames, dim, val */
	return val;
}

SEXP dppMatrix_trf_(SEXP obj, int warn)
{
	SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("pCholesky")),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		x = duplicate(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		int info;
		double *px = REAL(x);

		F77_CALL(dpptrf)(&ul, pdim, px, &info FCONE);
		ERROR_LAPACK_3(dpptrf, info, warn, 5);

		SET_SLOT(val, Matrix_xSym, x);
		UNPROTECT(1); /* x */
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

int dgCMatrix_trf_(const cs *A, css **S, csn **N, int order, double tol)
{
	cs *T = NULL;

#define FREE_AND_RETURN_ZERO \
	do { \
		if (*S) \
			*S = cs_sfree(*S); \
		if (*N) \
			*N = cs_nfree(*N); \
		if (T) \
			T = cs_spfree(T); \
		return 0; \
	} while (0)

	/* Symbolic analysis and numeric factorization : */
	if (!(*S = cs_sqr(order, A, 0)) ||
	    !(*N = cs_lu(A, *S, tol)))
		FREE_AND_RETURN_ZERO;
	/* Drop zeros from L and sort it : */
	cs_dropzeros((*N)->L);
	T = cs_transpose((*N)->L, 1);
	if (!T)
		FREE_AND_RETURN_ZERO;
	(*N)->L = cs_spfree((*N)->L);
	(*N)->L = cs_transpose(T, 1);
	if (!(*N)->L)
		FREE_AND_RETURN_ZERO;
	T = cs_spfree(T);
	/* Drop zeros from U and sort it : */
	T = cs_transpose((*N)->U, 1);
	if (!T)
		FREE_AND_RETURN_ZERO;
	(*N)->U = cs_spfree((*N)->U);
	(*N)->U = cs_transpose(T, 1);
	if (!(*N)->U)
		FREE_AND_RETURN_ZERO;
	T = cs_spfree(T);
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
	PROTECT(val = NEW_OBJECT_OF_CLASS("sparseLU"));

	const cs *A = dgC2cs(obj);
	css *S = NULL;
	csn *N = NULL;
	int *pp = NULL;
	if (A->m != A->n)
		error(_("LU factorization of m-by-n %s requires m == n"),
		      "dgCMatrix");
	if (!dgCMatrix_trf_(A, &S, &N, order_, tol_) ||
	    !(pp = cs_pinv(N->pinv, A->m))) {
		if (!pp) {
			S = cs_sfree(S);
			N = cs_nfree(N);
		}
		if (asLogical(doError))
			error(_("LU factorization of %s failed: out of memory or near-singular"),
			      "dgCMatrix");
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

	SEXP L = PROTECT(cs2dgC(N->L, "dtCMatrix")),
	U = PROTECT(cs2dgC(N->U, "dtCMatrix")),
	uplo = PROTECT(mkString("L"));
	SET_SLOT(L, Matrix_uploSym, uplo);
	SET_SLOT(val, Matrix_LSym, L);
	SET_SLOT(val, Matrix_USym, U);
	UNPROTECT(3); /* uplo, U, L */

	SEXP p = PROTECT(allocVector(INTSXP, A->m));
	Matrix_memcpy(INTEGER(p), pp, A->m, sizeof(int));
	SET_SLOT(val, Matrix_pSym, p);
	UNPROTECT(1); /* p */
	if (order_ > 0) {
		SEXP q = PROTECT(allocVector(INTSXP, A->n));
		int *pq = S->q;
		Matrix_memcpy(INTEGER(q), pq, A->n, sizeof(int));
		SET_SLOT(val, Matrix_qSym, q);
		UNPROTECT(1); /* q */
	}

	S = cs_sfree(S);
	N = cs_nfree(N);
	pp = cs_free(pp);

	set_factor(obj, (order_) ? "sparseLU~" : "sparseLU", val);
	UNPROTECT(1); /* val */
	return val;
}

int dgCMatrix_orf_(const cs *A, css **S, csn **N, int order)
{
	cs *T = NULL;

	/* Symbolic analysis and numeric factorization : */
	if (!(*S = cs_sqr(order, A, 1)) ||
	    !(*N = cs_qr(A, *S)))
		FREE_AND_RETURN_ZERO;
	/* Drop zeros from V and sort it : */
	cs_dropzeros((*N)->L);
	T = cs_transpose((*N)->L, 1);
	if (!T)
		FREE_AND_RETURN_ZERO;
	(*N)->L = cs_spfree((*N)->L);
	(*N)->L = cs_transpose(T, 1);
	if (!(*N)->L)
		FREE_AND_RETURN_ZERO;
	T = cs_spfree(T);
	/* Drop zeros from R and sort it : */
	T = cs_transpose((*N)->U, 1);
	if (!T)
		FREE_AND_RETURN_ZERO;
	(*N)->U = cs_spfree((*N)->U);
	(*N)->U = cs_transpose(T, 1);
	if (!(*N)->U)
		FREE_AND_RETURN_ZERO;
	T = cs_spfree(T);
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
	PROTECT(val = NEW_OBJECT_OF_CLASS("sparseQR"));

	const cs *A = dgC2cs(obj);
	css *S = NULL;
	csn *N = NULL;
	int *pp = NULL;
	if (A->m < A->n)
		error(_("QR factorization of m-by-n %s requires m >= n"),
		      "dgCMatrix");
	if (!dgCMatrix_orf_(A, &S, &N, order_) ||
	    !(pp = cs_pinv(S->pinv, S->m2))) {
		if (!pp) {
			S = cs_sfree(S);
			N = cs_nfree(N);
		}
		if (asLogical(doError))
			error(_("QR factorization of %s failed: out of memory"),
			      "dgCMatrix");
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

	SEXP V = PROTECT(cs2dgC(N->L, "dgCMatrix")),
	R = PROTECT(cs2dgC(N->U, "dgCMatrix"));
	SET_SLOT(val, Matrix_VSym, V);
	SET_SLOT(val, Matrix_RSym, R);
	UNPROTECT(2); /* R, V */

	SEXP beta = PROTECT(allocVector(REALSXP, A->n));
	double *pbeta = N->B;
	Matrix_memcpy(REAL(beta), pbeta, A->n, sizeof(double));
	SET_SLOT(val, Matrix_betaSym, beta);
	UNPROTECT(1); /* beta */

	SEXP p = PROTECT(allocVector(INTSXP, S->m2));
	Matrix_memcpy(INTEGER(p), pp, S->m2, sizeof(int));
	SET_SLOT(val, Matrix_pSym, p);
	UNPROTECT(1); /* p */
	if (order_ > 0) {
		SEXP q = PROTECT(allocVector(INTSXP, A->n));
		int *pq = S->q;
		Matrix_memcpy(INTEGER(q), pq, A->n, sizeof(int));
		SET_SLOT(val, Matrix_qSym, q);
		UNPROTECT(1); /* q */
	}

	S = cs_sfree(S);
	N = cs_nfree(N);
	pp = cs_free(pp);

	set_factor(obj, (order_) ? "sparseQR~" : "sparseQR", val);
	UNPROTECT(1); /* val */
	return val;
}

int dpCMatrix_trf_(cholmod_sparse *A, cholmod_factor **L,
                   int perm, int ldl, int super, double mult)
{
	CHM_store_common();

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
	int res = cholmod_factorize_p(A, beta, (int *) NULL, 0, *L, &c);

	CHM_restore_common();

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
	cholmod_sparse *A = dgC2cholmod(obj);
	cholmod_factor *L = NULL;

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */
	A->stype = (ul == 'U') ? 1 : -1;

	if (cached) {
		L = mf2cholmod(trf);
		L = cholmod_copy_factor(L, &c);
		dpCMatrix_trf_(A, &L, perm_, ldl_, super_, mult_);
	} else {
		dpCMatrix_trf_(A, &L, perm_, ldl_, super_, mult_);
		if (super_ == NA_LOGICAL) {
			nm[0] = (L->is_super) ? 'S' : 's';
			nm[2] = (L->is_ll   ) ? 'd' : 'D';
		}
	}
	REPROTECT(trf = cholmod2mf(L), pid);
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
	SEXP P_ = PROTECT(NEW_OBJECT_OF_CLASS("pMatrix")),
		T_ = PROTECT(NEW_OBJECT_OF_CLASS("dtCMatrix")),
		D_ = PROTECT(NEW_OBJECT_OF_CLASS("dsCMatrix")),
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

static SEXP mkDet(double modulus, int logarithm, int sign)
{
	SEXP nms = PROTECT(allocVector(STRSXP, 2)),
		cl = PROTECT(mkString("det")),
		det = PROTECT(allocVector(VECSXP, 2)),
		det0 = PROTECT(ScalarReal((logarithm) ? modulus : exp(modulus))),
		det1 = PROTECT(ScalarInteger(sign)),
		det0a = PROTECT(ScalarLogical(logarithm));
	SET_STRING_ELT(nms, 0, mkChar("modulus"));
	SET_STRING_ELT(nms, 1, mkChar("sign"));
	setAttrib(det, R_NamesSymbol, nms);
	setAttrib(det, R_ClassSymbol, cl);
	setAttrib(det0, install("logarithm"), det0a);
	SET_VECTOR_ELT(det, 0, det0);
	SET_VECTOR_ELT(det, 1, det1);
	UNPROTECT(6);
	return det;
}

SEXP denseLU_determinant(SEXP obj, SEXP logarithm)
{

#define DETERMINANT_START(_MAYBE_NOT_SQUARE_) \
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), n = pdim[0]; \
	if ((_MAYBE_NOT_SQUARE_) && pdim[1] != n) \
		error(_("determinant of non-square matrix is undefined")); \
	UNPROTECT(1); /* dim */ \
	int givelog = asLogical(logarithm) != 0, sign = 1; \
	double modulus = 0.0; /* result for n == 0 */

	DETERMINANT_START(1);
	if (n > 0) {
		SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		int j, *ppivot = INTEGER(pivot);
		R_xlen_t n1a = (R_xlen_t) n + 1;
		double *px = REAL(x);

		for (j = 0; j < n; ++j, px += n1a, ++ppivot) {
			if (*px < 0.0) {
				modulus += log(-(*px));
				if (*ppivot == j + 1)
					sign = -sign;
			} else {
				/* incl. 0, NaN cases */
				modulus += log(*px);
				if (*ppivot != j + 1)
					sign = -sign;
			}
		}
		UNPROTECT(2); /* x, pivot */
	}
	return mkDet(modulus, givelog, sign);
}

SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm, SEXP packed)
{
	DETERMINANT_START(0);
	if (n > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
		UNPROTECT(1); /* uplo */

		SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		int j = 0, unpacked = !asLogical(packed),
			*ppivot = INTEGER(pivot);
		R_xlen_t n1a = (R_xlen_t) n + 1;
		double *px = REAL(x), a, b, c, logab, logcc;
		while (j < n) {
			if (ppivot[j] > 0) {
				if (*px < 0.0) {
					modulus += log(-(*px));
					sign = -sign;
				} else {
					/* incl. 0, NaN cases */
					modulus += log(*px);
				}
				px += (unpacked) ? n1a : ((upper) ? j + 2 : n - j);
				j += 1;
			} else {
				a = *px;
				if (upper) {
					px += (unpacked) ? n1a : j + 2;
					b = *px;
					c = *(px - 1);
					px += (unpacked) ? n1a : j + 3;
				} else {
					c = *(px + 1);
					px += (unpacked) ? n1a : n - j;
					b = *px;
					px += (unpacked) ? n1a : n - j - 1;
				}
				logab = log(fabs(a)) + log(fabs(b));
				logcc = 2.0 * log(fabs(c));
				if ((a < 0.0) != (b < 0.0)) {
					/* det = ab - cc = -(abs(ab) + cc) < 0 */
					modulus += logspace_add(logab, logcc);
					sign = -sign;
				} else if (logab < logcc) {
					/* det = ab - cc = -(cc - ab) < 0 */
					modulus += logspace_sub(logcc, logab);
					sign = -sign;
				} else {
					/* det = ab - cc > 0 */
					modulus += logspace_sub(logab, logcc);
				}
				j += 2;
			}
		}
		UNPROTECT(2); /* x, pivot */
	}
	return mkDet(modulus, givelog, sign);
}

SEXP Cholesky_determinant(SEXP obj, SEXP logarithm, SEXP packed)
{
	DETERMINANT_START(0);
	if (n > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
		UNPROTECT(1); /* uplo */

		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		int j, unpacked = !asLogical(packed);
		R_xlen_t n1a = (R_xlen_t) n + 1;
		double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			if (*px < 0.0) {
				modulus += log(-(*px));
				sign = -sign;
			} else {
				/* incl. 0, NaN cases */
				modulus += log(*px);
			}
			px += (unpacked) ? n1a : ((upper) ? j + 2 : n - j);
		}
		modulus *= 2.0;
		UNPROTECT(1); /* x */
	}
	return mkDet(modulus, givelog, sign);
}

SEXP sparseLU_determinant(SEXP obj, SEXP logarithm)
{
	DETERMINANT_START(0);
	if (n > 0) {
		SEXP U = PROTECT(GET_SLOT(obj, Matrix_USym)),
			p = PROTECT(GET_SLOT(U, Matrix_pSym)),
			i = PROTECT(GET_SLOT(U, Matrix_iSym)),
			x = PROTECT(GET_SLOT(U, Matrix_xSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k = 0, kend;
		double *px = REAL(x);

		for (j = 0; j < n; ++j) {
			kend = *(++pp);
			if (kend > k && pi[kend - 1] == j) {
				if (px[kend - 1] < 0.0) {
					modulus += log(-px[kend - 1]);
					sign = -sign;
				} else {
					/* incl. 0, NaN cases */
					modulus += log(px[kend - 1]);
				}
			} else {
				UNPROTECT(4); /* x, i, p, U */
				return mkDet(R_NegInf, givelog, 1);
			}
			k = kend;
		}
		UNPROTECT(4); /* x, i, p, U */

		PROTECT(p = GET_SLOT(obj, Matrix_pSym));
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		UNPROTECT(1); /* p */
		PROTECT(p = GET_SLOT(obj, Matrix_qSym));
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		UNPROTECT(1); /* p */
	}
	return mkDet(modulus, givelog, sign);
}

SEXP sparseQR_determinant(SEXP obj, SEXP logarithm)
{
	DETERMINANT_START(1);
	if (n > 0) {
		SEXP R = PROTECT(GET_SLOT(obj, Matrix_RSym));
		PROTECT(dim = GET_SLOT(R, Matrix_DimSym));
		if (INTEGER(dim)[0] > n)
			error(_("%s(<%s>) does not support structurally rank deficient case"),
			      "determinant", "sparseQR");
		UNPROTECT(1); /* dim */

		SEXP p = PROTECT(GET_SLOT(R, Matrix_pSym)),
			i = PROTECT(GET_SLOT(R, Matrix_iSym)),
			x = PROTECT(GET_SLOT(R, Matrix_xSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k = 0, kend;
		double *px = REAL(x);

		for (j = 0; j < n; ++j) {
			kend = *(++pp);
			if (kend > k && pi[kend - 1] == j) {
				if (px[kend - 1] < 0.0) {
					modulus += log(-px[kend - 1]);
					sign = -sign;
				} else {
					/* incl. 0, NaN cases */
					modulus += log(px[kend - 1]);
				}
			} else {
				UNPROTECT(4); /* x, i, p, R */
				return mkDet(R_NegInf, givelog, 1);
			}
			k = kend;
		}
		UNPROTECT(4); /* x, i, p, U */

		PROTECT(p = GET_SLOT(obj, Matrix_pSym));
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		UNPROTECT(1); /* p */
		PROTECT(p = GET_SLOT(obj, Matrix_qSym));
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		UNPROTECT(1); /* p */
		if (n % 2)
			sign = -sign;
	}
	return mkDet(modulus, givelog, sign);
}

SEXP CHMfactor_determinant(SEXP obj, SEXP logarithm, SEXP sqrt)
{
	DETERMINANT_START(0);
	if (n > 0) {
		int sqrt_ = asLogical(sqrt);
		cholmod_factor *L = mf2cholmod(obj);
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
					modulus += log(*px_);
					px_ += nr1a;
				}
			}
			modulus *= 2.0;
		} else {
			int j, *pp = (int *) L->p;
			double *px = (double *) L->x;
			if (L->is_ll) {
				for (j = 0; j < n; ++j)
					modulus += log(px[pp[j]]);
				modulus *= 2.0;
			} else {
				for (j = 0; j < n; ++j) {
					if (px[pp[j]] < 0.0) {
						if (sqrt_)
							return mkDet(R_NaN, givelog, 1);
						modulus += log(-px[pp[j]]);
						sign = -sign;
					} else {
						/* incl. 0, NaN cases */
						modulus += log(px[pp[j]]);
					}
				}
			}
		}
		if (sqrt_)
			modulus *= 0.5;
	}
	return mkDet(modulus, givelog, sign);
}

static void solveDN(SEXP rdn, SEXP adn, SEXP bdn)
{
	SEXP s;
	if (!isNull(s = VECTOR_ELT(adn, 1)))
		SET_VECTOR_ELT(rdn, 0, s);
	if (!isNull(s = VECTOR_ELT(bdn, 1)))
		SET_VECTOR_ELT(rdn, 1, s);
	PROTECT(adn = getAttrib(adn, R_NamesSymbol));
	PROTECT(bdn = getAttrib(bdn, R_NamesSymbol));
	if(!isNull(adn) || !isNull(bdn)) {
		PROTECT(s = allocVector(STRSXP, 2));
		if (!isNull(adn))
			SET_STRING_ELT(s, 0, STRING_ELT(adn, 1));
		if (!isNull(bdn))
			SET_STRING_ELT(s, 1, STRING_ELT(bdn, 1));
		setAttrib(rdn, R_NamesSymbol, s);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return;
}

SEXP denseLU_solve(SEXP a, SEXP b)
{

#define SOLVE_START(_MAYBE_NOT_SQUARE_) \
	SEXP adim = PROTECT(GET_SLOT(a, Matrix_DimSym)); \
	int *padim = INTEGER(adim), m = padim[0], n = m; \
	if ((_MAYBE_NOT_SQUARE_) && padim[1] != m) \
		error(_("'%s' is not square"), "a"); \
	UNPROTECT(1); /* adim */ \
	if (!isNull(b)) { \
	SEXP bdim = PROTECT(GET_SLOT(b, Matrix_DimSym)); \
	int *pbdim = INTEGER(bdim); \
	if (pbdim[0] != m) \
		error(_("dimensions of '%s' and '%s' are inconsistent"), \
		      "a", "b"); \
	n = pbdim[1]; \
	UNPROTECT(1); /* bdim */ \
	}

#define SOLVE_FINISH \
	SEXP rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym)), \
	adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)); \
	if (isNull(b)) \
		revDN(rdimnames, adimnames); \
	else { \
		SEXP bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)); \
		solveDN(rdimnames, adimnames, bdimnames); \
		UNPROTECT(1); /* bdimnames */ \
	} \
	UNPROTECT(2); /* adimnames, rdimnames */

	SOLVE_START(1);
	SEXP r = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),
		rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
	int *prdim = INTEGER(rdim);
	prdim[0] = m;
	prdim[1] = n;
	if (m > 0) {
		SEXP rx, ax = PROTECT(GET_SLOT(a, Matrix_xSym)),
			apivot = PROTECT(GET_SLOT(a, Matrix_permSym));
		int info;
		if (isNull(b)) {
			PROTECT(rx = duplicate(ax));
			int lwork = -1;
			double work0, *work = &work0;
			F77_CALL(dgetri)(&m, REAL(rx), &m, INTEGER(apivot),
			                 work, &lwork, &info);
			ERROR_LAPACK_1(dgetri, info);
			lwork = (int) work0;
			work = (double *) R_alloc((size_t) lwork, sizeof(double));
			F77_CALL(dgetri)(&m, REAL(rx), &m, INTEGER(apivot),
			                 work, &lwork, &info);
			ERROR_LAPACK_2(dgetri, info, 2, U);
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			rx = duplicate(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
			F77_CALL(dgetrs)("N", &m, &n, REAL(ax), &m, INTEGER(apivot),
			                 REAL(rx), &m, &info FCONE);
			ERROR_LAPACK_1(dgetrs, info);
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(3); /* rx, apivot, ax */
	}
	SOLVE_FINISH;
	UNPROTECT(2); /* rdim, r */
	return r;
}

SEXP BunchKaufman_solve(SEXP a, SEXP b, SEXP packed)
{
	SOLVE_START(0);
	int unpacked = !asLogical(packed);
	const char *cl = (!isNull(b)) ? "dgeMatrix" :
		((unpacked) ? "dsyMatrix" : "dspMatrix");
	SEXP r = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
		rdim = PROTECT(GET_SLOT(r, Matrix_DimSym)),
		auplo = PROTECT(GET_SLOT(a, Matrix_uploSym));
	int *prdim = INTEGER(rdim);
	prdim[0] = m;
	prdim[1] = n;
	if (isNull(b))
		SET_SLOT(r, Matrix_uploSym, auplo);
	if (m > 0) {
		SEXP rx, ax = PROTECT(GET_SLOT(a, Matrix_xSym)),
			apivot = PROTECT(GET_SLOT(a, Matrix_permSym));
		char ul = *CHAR(STRING_ELT(auplo, 0));
		int info;
		if (isNull(b)) {
			PROTECT(rx = duplicate(ax));
			double *work = (double *) R_alloc((size_t) m, sizeof(double));
			if (unpacked) {
				F77_CALL(dsytri)(&ul, &m, REAL(rx), &m, INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(dsytri, info, 2, D);
			} else {
				F77_CALL(dsptri)(&ul, &m, REAL(rx),     INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(dsptri, info, 2, D);
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			rx = duplicate(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
			if (unpacked) {
				F77_CALL(dsytrs)(&ul, &m, &n, REAL(ax), &m,
				                 INTEGER(apivot),
				                 REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dsytrs, info);
			} else {
				F77_CALL(dsptrs)(&ul, &m, &n, REAL(ax),
				                 INTEGER(apivot),
				                 REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dsptrs, info);
			}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(3); /* rx, apivot, ax */
	}
	SOLVE_FINISH;
	UNPROTECT(3); /* auplo, rdim, r */
	return r;
}

SEXP Cholesky_solve(SEXP a, SEXP b, SEXP packed)
{
	SOLVE_START(0);
	int unpacked = !asLogical(packed);
	const char *cl = (!isNull(b)) ? "dgeMatrix" :
		((unpacked) ? "dpoMatrix" : "dppMatrix");
	SEXP r = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
		rdim = PROTECT(GET_SLOT(r, Matrix_DimSym)),
		auplo = PROTECT(GET_SLOT(a, Matrix_uploSym));
	int *prdim = INTEGER(rdim);
	prdim[0] = m;
	prdim[1] = n;
	if (isNull(b))
		SET_SLOT(r, Matrix_uploSym, auplo);
	if (m > 0) {
		SEXP rx, ax = PROTECT(GET_SLOT(a, Matrix_xSym));
		char ul = *CHAR(STRING_ELT(auplo, 0));
		int info, nprotect = 2;

		SEXP aperm = NULL;
		if (HAS_SLOT(a, Matrix_permSym)) {
			SEXP tmp = GET_SLOT(a, Matrix_permSym);
			if (LENGTH(tmp) > 0) {
				PROTECT(aperm = tmp);
				++nprotect;
			}
		} /* else 'a' is a dtrMatrix or dtpMatrix, as in chol2inv */

		if (isNull(b)) {
			PROTECT(rx = duplicate(ax));
			SET_SLOT(r, Matrix_uploSym, auplo);
			if (unpacked) {
				F77_CALL(dpotri)(&ul, &m, REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_2(dpotri, info, 2, L);
				if (aperm)
					symPerm(REAL(rx), n, ul, INTEGER(aperm), 1, 1);
			} else {
				F77_CALL(dpptri)(&ul, &m, REAL(rx),     &info FCONE);
				ERROR_LAPACK_2(dpptri, info, 2, L);
				if (aperm) {
					/* FIXME: symPerm() supporting packed matrices */
					double *work;
					size_t lwork = (size_t) n * n;
					Matrix_Calloc(work, lwork, double);
					ddense_unpack(work, REAL(rx), n, ul, 'N');
					symPerm(work, n, ul, INTEGER(aperm), 1, 1);
					ddense_pack  (REAL(rx), work, n, ul, 'N');
					Matrix_Free(work, lwork);
				}
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			rx = duplicate(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
			if (aperm)
				rowPerm(REAL(rx), m, n, INTEGER(aperm), 1, 0);
			if (unpacked) {
				F77_CALL(dpotrs)(&ul, &m, &n, REAL(ax), &m,
				                 REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dpotrs, info);
			} else {
				F77_CALL(dpptrs)(&ul, &m, &n, REAL(ax),
				                 REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dpptrs, info);
			}
			if (aperm)
				rowPerm(REAL(rx), m, n, INTEGER(aperm), 1, 1);
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(nprotect); /* rx, aperm, ax */
	}
	SOLVE_FINISH;
	UNPROTECT(3); /* auplo, rdim, r */
	return r;
}

SEXP sparseLU_solve(SEXP a, SEXP b, SEXP sparse)
{

#define ERROR_SOLVE_OOM(_A_, _B_) \
	error(_("%s(<%s>, <%s>) failed: out of memory"), "solve", #_A_, #_B_)

	SOLVE_START(0);
	SEXP r,
	aL = PROTECT(GET_SLOT(a, Matrix_LSym)),
	aU = PROTECT(GET_SLOT(a, Matrix_USym)),
	ap = PROTECT(GET_SLOT(a, Matrix_pSym)),
	aq = PROTECT(GET_SLOT(a, Matrix_qSym));
	int j,
		*pap = INTEGER(ap),
		*paq = (LENGTH(aq)) ? INTEGER(aq) : (int *) NULL;
	double *work = (double *) R_alloc((size_t) m, sizeof(double));
	cs *L = dgC2cs(aL), *U = dgC2cs(aU);
	if (!asLogical(sparse)) {
		PROTECT(r = NEW_OBJECT_OF_CLASS("dgeMatrix"));
		SEXP rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
		int *prdim = INTEGER(rdim);
		prdim[0] = m;
		prdim[1] = n;
		UNPROTECT(1); /* rdim */
		R_xlen_t mn = (R_xlen_t) m * n;
		SEXP rx = PROTECT(allocVector(REALSXP, mn));
		double *prx = REAL(rx);
		if (isNull(b)) {
			Matrix_memset(prx, 0, mn, sizeof(double));
			for (j = 0; j < n; ++j) {
				prx[j] = 1.0;
				cs_pvec(pap, prx, work, m);
				cs_lsolve(L, work);
				cs_usolve(U, work);
				cs_ipvec(paq, work, prx, m);
				prx += m;
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			double *pbx = REAL(bx);
			for (j = 0; j < n; ++j) {
				cs_pvec(pap, pbx, work, m);
				cs_lsolve(L, work);
				cs_usolve(U, work);
				cs_ipvec(paq, work, prx, m);
				prx += m;
				pbx += m;
			}
			UNPROTECT(1); /* bx */
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
	} else {
		cs *B, *X;
		int *papinv = cs_pinv(pap, m);
		if (!papinv)
			ERROR_SOLVE_OOM(sparseLU, dgCMatrix);
		if (isNull(b)) {
			B = cs_spalloc(m, n, n, 1, 0);
			if (!B)
				ERROR_SOLVE_OOM(sparseLU, dgCMatrix);
			for (j = 0; j < n; ++j) {
				B->p[j] = j;
				B->i[j] = j;
				B->x[j] = 1.0;
			}
			B->p[n] = n;
			X = cs_permute(        B, papinv, (int *) NULL, 1);
			B = cs_spfree(B);
		} else
			X = cs_permute(dgC2cs(b), papinv, (int *) NULL, 1);
		papinv = cs_free(papinv);
		if (!X)
			ERROR_SOLVE_OOM(sparseLU, dgCMatrix);
		B = X;

		int i, k, top, nz, nzmax,
			*iwork = (int *) R_alloc((size_t) 2 * m, sizeof(int));

#define DO_TRIANGULAR_SOLVE(_A_, _LOA_, _FRB_, _CLA_, _CLB_) \
		do { \
			X = cs_spalloc(m, n, B->nzmax, 1, 0); \
			if (!X) { \
				if (_FRB_) \
					B = cs_spfree(B); \
				ERROR_SOLVE_OOM(_CLA_, _CLB_); \
			} \
			X->p[0] = nz = 0; \
			nzmax = X->nzmax; \
			for (j = 0, k = 0; j < n; ++j) { \
				top = cs_spsolve(_A_, B, j, iwork, work, (int *) NULL, _LOA_); \
				if (m - top > INT_MAX - nz) { \
					if (_FRB_) \
						B = cs_spfree(B); \
					X = cs_spfree(X); \
					error(_("attempt to construct sparse matrix with more than %s nonzero elements"), \
					      "2^31-1"); \
				} \
				nz += m - top; \
				if (nz > nzmax) { \
					nzmax = (nz <= INT_MAX / 2) ? 2 * nz : INT_MAX; \
					if (!cs_sprealloc(X, nzmax)) { \
						if (_FRB_) \
							B = cs_spfree(B); \
						X = cs_spfree(X); \
						ERROR_SOLVE_OOM(_CLA_, _CLB_); \
					} \
				} \
				X->p[j + 1] = nz; \
				if (_LOA_) { \
					for (i = top; i <    m; ++i) { \
						X->i[k] =      iwork[i]; \
						X->x[k] = work[iwork[i]]; \
						++k; \
					} \
				} else { \
					for (i = m - 1; i >= top; --i) { \
						X->i[k] =      iwork[i]; \
						X->x[k] = work[iwork[i]]; \
						++k; \
					} \
				} \
			} \
			if (_FRB_) \
				B = cs_spfree(B); \
			B = X; \
		} while (0)

		DO_TRIANGULAR_SOLVE(L, 1, 1, sparseLU, dgCMatrix);
		DO_TRIANGULAR_SOLVE(U, 0, 1, sparseLU, dgCMatrix);

		if (paq) {
			X = cs_permute(B, paq, (int *) NULL, 1);
			B = cs_spfree(B);
			if (!X)
				ERROR_SOLVE_OOM(sparseLU, dgCMatrix);
			B = X;
		}

		/* Drop zeros from B and sort it : */
		cs_dropzeros(B);
		X = cs_transpose(B, 1);
		B = cs_spfree(B);
		if (!X)
			ERROR_SOLVE_OOM(sparseLU, dgCMatrix);
		B = cs_transpose(X, 1);
		X = cs_spfree(X);
		if (!B)
			ERROR_SOLVE_OOM(sparseLU, dgCMatrix);

		PROTECT(r = cs2dgC(B, "dgCMatrix"));
		B = cs_spfree(B);
	}
	SOLVE_FINISH;
	UNPROTECT(5); /* r, aq, ap, aU, aL */
	return r;
}

/* MJ: not needed since we have 'sparseQR_matmult' : */
#if 0

SEXP sparseQR_solve(SEXP a, SEXP b, SEXP sparse)
{
	return R_NilValue;
}

#endif /* MJ */

SEXP CHMfactor_solve(SEXP a, SEXP b, SEXP sparse, SEXP system)
{
	/* see top of :
	   ./CHOLMOD/Cholesky/cholmod_solve.c
	   ./CHOLMOD/Cholesky/cholmod_spsolve.c
	*/
	static const char *valid[] = {
		"A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt", "" };
	int ivalid = -1;
	if (TYPEOF(system) != STRSXP || LENGTH(system) < 1 ||
	    (system = STRING_ELT(system, 0)) == NA_STRING ||
	    (ivalid = strmatch(CHAR(system), valid)) < 0)
		error(_("invalid '%s' to %s()"), "system", __func__);
	SOLVE_START(0);
	SEXP r;
	int j;
	cholmod_factor *L = mf2cholmod(a);
	if (!asLogical(sparse)) {
		cholmod_dense *B, *X;
		if (isNull(b)) {
			B = cholmod_allocate_dense(m, n, m, CHOLMOD_REAL, &c);
			if (!B)
				ERROR_SOLVE_OOM(CHMfactor, dgeMatrix);
			R_xlen_t m1a = (R_xlen_t) m + 1;
			double *px = (double *) B->x;
			Matrix_memset(px, 0, (R_xlen_t) m * n, sizeof(double));
			for (j = 0; j < n; ++j) {
				*px = 1.0;
				px += m1a;
			}
			X = cholmod_solve(ivalid, L, B, &c);
			if (!X)
				ERROR_SOLVE_OOM(CHMfactor, dgeMatrix);
			cholmod_free_dense(&B, &c);
			const char *cl = (ivalid < 2) ? "dpoMatrix" :
				((ivalid < 7) ? "dtrMatrix" : "dgeMatrix");
			PROTECT(r = cholmod2dge(X, cl));
		} else {
			B = dge2cholmod(b);
			X = cholmod_solve(ivalid, L, B, &c);
			if (!X)
				ERROR_SOLVE_OOM(CHMfactor, dgeMatrix);
			PROTECT(r = cholmod2dge(X, "dgeMatrix"));
		}
		cholmod_free_dense(&X, &c);
	} else {
		cholmod_sparse *B, *X;
		if (isNull(b)) {
			B = cholmod_allocate_sparse(m, n, n, 1, 1, 0, CHOLMOD_REAL, &c);
			if (!B)
				ERROR_SOLVE_OOM(CHMfactor, dgCMatrix);
			int *pp = (int *) B->p, *pi = (int *) B->i;
			double *px = (double *) B->x;
			for (j = 0; j < n; ++j) {
				pp[j] = j;
				pi[j] = j;
				px[j] = 1.0;
			}
			pp[n] = n;
			X = cholmod_spsolve(ivalid, L, B, &c);
			if (!X)
				ERROR_SOLVE_OOM(CHMfactor, dgCMatrix);
			cholmod_free_sparse(&B, &c);
			if (ivalid < 7) {
				X->stype = (ivalid == 2 || ivalid == 4) ? -1 : 1;
				cholmod_sort(X, &c);
				if (!X)
					ERROR_SOLVE_OOM(CHMfactor, dgCMatrix);
			}
			const char *cl = (ivalid < 2) ? "dsCMatrix" :
				((ivalid < 7) ? "dtCMatrix" : "dgCMatrix");
			PROTECT(r = cholmod2dgC(X, cl));
		} else {
			B = dgC2cholmod(b);
			X = cholmod_spsolve(ivalid, L, B, &c);
			if (!X)
				ERROR_SOLVE_OOM(CHMfactor, dgCMatrix);
			PROTECT(r = cholmod2dgC(X, "dgCMatrix"));
		}
		cholmod_free_sparse(&X, &c);
	}
	if (isNull(b) && (ivalid == 2 || ivalid == 4)) {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(r, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	SOLVE_FINISH;
	UNPROTECT(1); /* r */
	return r;
}

SEXP dtrMatrix_solve(SEXP a, SEXP b, SEXP packed)
{
	SOLVE_START(0);
	int unpacked = !asLogical(packed);
	const char *cl = (!isNull(b)) ? "dgeMatrix" :
		((unpacked) ? "dtrMatrix" : "dtpMatrix");
	SEXP r = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
		rdim = PROTECT(GET_SLOT(r, Matrix_DimSym)),
		auplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
		adiag = PROTECT(GET_SLOT(a, Matrix_diagSym));
	int *prdim = INTEGER(rdim);
	prdim[0] = m;
	prdim[1] = n;
	if (isNull(b)) {
		SET_SLOT(r, Matrix_uploSym, auplo);
		SET_SLOT(r, Matrix_diagSym, adiag);
	}
	if (m > 0) {
		SEXP rx, ax = PROTECT(GET_SLOT(a, Matrix_xSym));
		char ul = *CHAR(STRING_ELT(auplo, 0)),
			di = *CHAR(STRING_ELT(adiag, 0));
		int info;
		if (isNull(b)) {
			PROTECT(rx = duplicate(ax));
			SET_SLOT(r, Matrix_uploSym, auplo);
			SET_SLOT(r, Matrix_diagSym, adiag);
			if (unpacked) {
				F77_CALL(dtrtri)(&ul, &di, &m, REAL(rx), &m,
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(dtrtri, info, 2, A);
			} else {
				F77_CALL(dtptri)(&ul, &di, &m, REAL(rx),
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(dtptri, info, 2, A);
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			rx = duplicate(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
			if (unpacked) {
				F77_CALL(dtrtrs)(&ul, "N", &di, &m, &n, REAL(ax), &m,
				                 REAL(rx), &m, &info FCONE FCONE FCONE);
				ERROR_LAPACK_1(dtrtrs, info);
			} else {
				// https://bugs.r-project.org/show_bug.cgi?id=18534
				F77_CALL(dtptrs)(&ul, "N", &di, &m, &n, REAL(ax),
				                 REAL(rx), &m, &info
#ifdef usePR18534fix
				                 FCONE FCONE FCONE);
#else
				                 FCONE FCONE);
#endif
				ERROR_LAPACK_1(dtptrs, info);
			}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(2); /* rx, ax */
	}
	SOLVE_FINISH;
	UNPROTECT(4); /* adiag, auplo, rdim, r */
	return r;
}

SEXP dtCMatrix_solve(SEXP a, SEXP b, SEXP sparse)
{
	SOLVE_START(0);
	SEXP r, auplo = PROTECT(GET_SLOT(a, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(auplo, 0));
	int j;
	cs *A = dgC2cs(a);
	if (!asLogical(sparse)) {
		const char *cl = (isNull(b)) ? "dtrMatrix" : "dgeMatrix";
		PROTECT(r = NEW_OBJECT_OF_CLASS(cl));

		SEXP rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
		int *prdim = INTEGER(rdim);
		prdim[0] = m;
		prdim[1] = n;
		UNPROTECT(1); /* rdim */

		R_xlen_t mn = (R_xlen_t) m * n;
		SEXP rx = PROTECT(allocVector(REALSXP, mn));
		double *prx = REAL(rx);
		if (isNull(b)) {
			Matrix_memset(prx, 0, mn, sizeof(double));
			for (j = 0; j < n; ++j) {
				prx[j] = 1.0;
				if (ul == 'U')
					cs_usolve(A, prx);
				else
					cs_lsolve(A, prx);
				prx += m;
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			double *pbx = REAL(bx);
			Matrix_memcpy(prx, pbx, mn, sizeof(double));
			UNPROTECT(1); /* bx */
			for (j = 0; j < n; ++j) {
				if (ul == 'U')
					cs_usolve(A, prx);
				else
					cs_lsolve(A, prx);
				prx += m;
			}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
	} else {
		const char *cl = (isNull(b)) ? "dtCMatrix" : "dgCMatrix";
		cs *B, *X;

		if (isNull(b)) {
			B = cs_spalloc(m, n, n, 1, 0);
			if (!B)
				ERROR_SOLVE_OOM(dtCMatrix, dgCMatrix);
			for (j = 0; j < n; ++j) {
				B->p[j] = j;
				B->i[j] = j;
				B->x[j] = 1.0;
			}
			B->p[n] = n;
		} else
			B = dgC2cs(b);

		int i, k, top, nz, nzmax,
			*iwork = (int *) R_alloc((size_t) 2 * m, sizeof(int));
		double *work = (double *) R_alloc((size_t) m, sizeof(double));

		DO_TRIANGULAR_SOLVE(A, ul != 'U', isNull(b), dtCMatrix, dgCMatrix);

		/* Drop zeros from B and sort it : */
		cs_dropzeros(B);
		X = cs_transpose(B, 1);
		B = cs_spfree(B);
		if (!X)
			ERROR_SOLVE_OOM(dtCMatrix, dgCMatrix);
		B = cs_transpose(X, 1);
		X = cs_spfree(X);
		if (!B)
			ERROR_SOLVE_OOM(dtCMatrix, dgCMatrix);

		PROTECT(r = cs2dgC(B, cl));
		B = cs_spfree(B);
	}
	if (isNull(b))
		SET_SLOT(r, Matrix_uploSym, auplo);
	SOLVE_FINISH;
	UNPROTECT(2); /* r, auplo */
	return r;
}

SEXP sparseQR_matmult(SEXP qr, SEXP y, SEXP op, SEXP complete, SEXP yxjj)
{
	SEXP V = PROTECT(GET_SLOT(qr, Matrix_VSym)),
		beta = PROTECT(GET_SLOT(qr, Matrix_betaSym)),
		p = PROTECT(GET_SLOT(qr, Matrix_pSym));
	const cs *V_ = dgC2cs(V);
	double *pbeta = REAL(beta);
	int m = V_->m, r = V_->n, n, i, j, op_ = asInteger(op),
		*pp = INTEGER(p), nprotect = 6;

	SEXP yx;
	double *pyx;
	if (isNull(y)) {
		n = (asLogical(complete)) ? m : r;

		R_xlen_t mn = (R_xlen_t) m * n, m1a = (R_xlen_t) m + 1;
		PROTECT(yx = allocVector(REALSXP, mn));
		pyx = REAL(yx);
		Matrix_memset(pyx, 0, mn, sizeof(double));

	if (isNull(yxjj)) {
		for (j = 0; j < n; ++j) {
			*pyx = 1.0;
			pyx += m1a;
		}
	} else if (TYPEOF(yxjj) == REALSXP && XLENGTH(yxjj) >= n) {
		double *pyxjj = REAL(yxjj);
		for (j = 0; j < n; ++j) {
			*pyx = *pyxjj;
			pyx += m1a;
			pyxjj += 1;
		}
	} else
		error(_("invalid '%s' to %s()"), "yxjj", __func__);
	} else {
		SEXP ydim = PROTECT(GET_SLOT(y, Matrix_DimSym));
		int *pydim = INTEGER(ydim);
		if (pydim[0] != m)
			error(_("dimensions of '%s' and '%s' are inconsistent"),
			     "qr", "y");
		n = pydim[1];
		UNPROTECT(1); /* ydim */

		PROTECT(yx = GET_SLOT(y, Matrix_xSym));
	}
	pyx = REAL(yx);

	SEXP a = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),
		adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
		ax = yx;
	int *padim = INTEGER(adim);
	padim[0] = (op_ != 0) ? m : r;
	padim[1] = n;
	if (!isNull(y) || padim[0] != m) {
		PROTECT(ax = allocVector(REALSXP, (R_xlen_t) padim[0] * padim[1]));
		++nprotect;
	}
	double *pax = REAL(ax), *work = NULL;
	if (op_ < 5)
		work = (double *) R_alloc((size_t) m, sizeof(double));

	switch (op_) {
	case 0: /* qr.coef : A = P2 R1^{-1} Q1' P1 y */
	{
		SEXP R = PROTECT(GET_SLOT(qr, Matrix_RSym)),
			q = PROTECT(GET_SLOT(qr, Matrix_qSym));
		const cs *R_ = dgC2cs(R);
		int *pq = (LENGTH(q)) ? INTEGER(q) : (int *) NULL;

		for (j = 0; j < n; ++j) {
			cs_pvec(pp, pyx, work, m);
			for (i = 0; i < r; ++i)
				cs_happly(V_, i, pbeta[i], work);
			cs_usolve(R_, work);
			cs_ipvec(pq, work, pax, r);
			pyx += m;
			pax += r;
		}

		UNPROTECT(2); /* q, R */
		break;
	}
	case 1: /* qr.fitted : A = P1' Q1 Q1' P1 y */
		for (j = 0; j < n; ++j) {
			cs_pvec(pp, pyx, work, m);
			for (i = 0; i < r; ++i)
				cs_happly(V_, i, pbeta[i], work);
			if (r < m)
				Matrix_memset(work + r, 0, m - r, sizeof(double));
			for (i = r - 1; i >= 0; --i)
				cs_happly(V_, i, pbeta[i], work);
			cs_ipvec(pp, work, pax, m);
			pyx += m;
			pax += m;
		}
		break;
	case 2: /* qr.resid : A = P1' Q2 Q2' P1 y */
		for (j = 0; j < n; ++j) {
			cs_pvec(pp, pyx, work, m);
			for (i = 0; i < r; ++i)
				cs_happly(V_, i, pbeta[i], work);
			if (r > 0)
				Matrix_memset(work, 0, r, sizeof(double));
			for (i = r - 1; i >= 0; --i)
				cs_happly(V_, i, pbeta[i], work);
			cs_ipvec(pp, work, pax, m);
			pyx += m;
			pax += m;
		}
		break;
	case 3: /* qr.qty {w/ perm.} : A = Q' P1 y */
		for (j = 0; j < n; ++j) {
			cs_pvec(pp, pyx, work, m);
			Matrix_memcpy(pax, work, m, sizeof(double));
			for (i = 0; i < r; ++i)
				cs_happly(V_, i, pbeta[i], pax);
			pyx += m;
			pax += m;
		}
		break;
	case 4: /* qr.qy {w/ perm.} : A = P1' Q y */
		for (j = 0; j < n; ++j) {
			Matrix_memcpy(work, pyx, m, sizeof(double));
			for (i = r - 1; i >= 0; --i)
				cs_happly(V_, i, pbeta[i], work);
			cs_ipvec(pp, work, pax, m);
			pyx += m;
			pax += m;
		}
	break;
	case 5: /* qr.qty {w/o perm.} : A = Q' y */
		if (ax != yx)
			Matrix_memcpy(pax, pyx, (R_xlen_t) m * n, sizeof(double));
		for (j = 0; j < n; ++j) {
			for (i = 0; i < r; ++i)
				cs_happly(V_, i, pbeta[i], pax);
			pax += m;
		}
	break;
	case 6: /* qr.qy {w/o perm.} : A = Q y */
		if (ax != yx)
			Matrix_memcpy(pax, pyx, (R_xlen_t) m * n, sizeof(double));
		for (j = 0; j < n; ++j) {
			for (i = r - 1; i >= 0; --i)
				cs_happly(V_, i, pbeta[i], pax);
			pax += m;
		}
	break;
	default:
		error(_("invalid '%s' to %s()"), "op", __func__);
		break;
	}

	SET_SLOT(a, Matrix_xSym, ax);
	UNPROTECT(nprotect); /* ax, adim, a, yx, p, beta, V */
	return a;
}

SEXP CHMfactor_diag_get(SEXP obj, SEXP square)
{
	cholmod_factor *L = mf2cholmod(obj);
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
	double mult_ = asReal(mult);
	if (!R_FINITE(mult_))
		error(_("'%s' is not a number or not finite"), "mult");

	cholmod_factor *L = cholmod_copy_factor(mf2cholmod(obj), &c);
	cholmod_sparse *A = dgC2cholmod(parent);
	if (Matrix_shape(parent) == 's') {
		SEXP uplo = PROTECT(GET_SLOT(parent, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		A->stype = (ul == 'U') ? 1 : -1;
		UNPROTECT(1);
	}

	dpCMatrix_trf_(A, &L, 0, !L->is_ll, L->is_super, mult_);

	SEXP res = PROTECT(cholmod2mf(L));
	cholmod_free_factor(&L, &c);

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_DimNamesSym, dimnames);

	UNPROTECT(2);
	return res;
}

SEXP CHMfactor_updown(SEXP obj, SEXP parent, SEXP update)
{
	cholmod_factor *L = cholmod_copy_factor(mf2cholmod(obj), &c);
	cholmod_sparse *A = dgC2cholmod(parent);
	if (Matrix_shape(parent) == 's') {
		SEXP uplo = PROTECT(GET_SLOT(parent, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		A->stype = (ul == 'U') ? 1 : -1;
		UNPROTECT(1);
	}

	cholmod_updown(asLogical(update) != 0, A, L, &c);

	SEXP res = PROTECT(cholmod2mf(L));
	cholmod_free_factor(&L, &c);

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_DimNamesSym, dimnames);

	UNPROTECT(2);
	return res;
}
