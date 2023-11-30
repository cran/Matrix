#include "Mdefines.h"
#include "chm_common.h"

/* NB: mostly parallel to CsparseMatrix_validate in ./validity.c */
SEXP checkpi(SEXP p, SEXP i, int m, int n)
{

#define MKMS(_FORMAT_, ...) mkString(Matrix_sprintf(_FORMAT_, __VA_ARGS__))

	if (TYPEOF(p) != INTSXP)
		return MKMS(_("'%s' slot is not of type \"%s\""),
		            "p", "integer");
	if (XLENGTH(p) - 1 != n)
		return MKMS(_("'%s' slot does not have length %s"),
		            "p", "Dim[2]+1");
	int *pp = INTEGER(p);
	if (pp[0] != 0)
		return MKMS(_("first element of '%s' slot is not 0"),
		            "p");
	int j;
	for (j = 1; j <= n; ++j) {
		if (pp[j] == NA_INTEGER)
			return MKMS(_("'%s' slot contains NA"),
			            "p");
		if (pp[j] < pp[j - 1])
			return MKMS(_("'%s' slot is not nondecreasing"),
			            "p");
		if (pp[j] - pp[j - 1] > m)
			return MKMS(_("first differences of '%s' slot exceed %s"),
			            "p", "Dim[1]");
	}

	if (TYPEOF(i) != INTSXP)
		return MKMS(_("'%s' slot is not of type \"%s\""),
		            "i", "integer");
	if (XLENGTH(i) < pp[n])
		return MKMS(_("'%s' slot has length less than %s"),
		            "i", "p[length(p)]");
	int *pi = INTEGER(i), k, kend, ik, i0, sorted = 1;
	for (j = 1, k = 0; j <= n; ++j) {
		kend = pp[j];
		i0 = -1;
		while (k < kend) {
			ik = pi[k];
			if (ik == NA_INTEGER)
				return MKMS(_("'%s' slot contains NA"),
				            "i");
			if (ik < 0 || ik >= m)
				return MKMS(_("'%s' slot has elements not in {%s}"),
				            "i", "0,...,Dim[1]-1");
			if (ik < i0)
				sorted = 0;
			else if (ik == i0)
				return MKMS(_("'%s' slot is not increasing within columns after sorting"),
				            "i");
			i0 = ik;
			++k;
		}
	}

	SEXP ans = allocVector(LGLSXP, 1);
	LOGICAL(ans)[0] = sorted;
	return ans;
}

/**
 * Coerce from CHMfactor to (cholmod_factor *)
 *
 * Sets the members of a pointed-to cholmod_factor struct, using "data"
 * obtained from slots of a CHMfactor.  The result should _not_ be
 * freed using cholmod_free_factor, as the resulting members point to
 * memory controlled by R, not by CHOLMOD.
 *
 * @param L a pointer to a cholmod_factor struct, to be modified in-place.
 * @param from an S4 object inheriting from virtual class CHMfactor.
 *
 * @return L.
 */
/* NB: mostly parallel to M2CHF in ./cholmod-etc.c */
cholmod_factor *sexp_as_cholmod_factor(cholmod_factor *L, SEXP from)
{
	static const char *valid[] = {
		"dCHMsuper", "dCHMsimpl", "nCHMsuper", "nCHMsimpl", "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *class = valid[ivalid];
	memset(L, 0, sizeof(cholmod_factor));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
		type = PROTECT(GET_SLOT(from, install("type"))),
		perm = PROTECT(GET_SLOT(from, Matrix_permSym)),
		colcount = PROTECT(GET_SLOT(from, install("colcount")));
	L->n = INTEGER(dim)[0];
	L->minor = L->n; /* FIXME: could be wrong for from <- new(...) */
	L->ordering = INTEGER(type)[0];
	if (L->ordering != CHOLMOD_NATURAL)
		L->Perm = INTEGER(perm);
	else {
		/* cholmod_check_factor allows L->Perm == NULL,
		   but cholmod_copy_factor does not test, so it segfaults ...
		*/
		int j, n = (int) L->n, *Perm = (int *) R_alloc(L->n, sizeof(int));
		for (j = 0; j < n; ++j)
			Perm[j] = j;
		L->Perm = Perm;
	}
	L->ColCount = INTEGER(colcount);
	L->is_super = INTEGER(type)[2];
	if (L->is_super) {
		L->is_ll = 1;
		L->is_monotonic = 1;
		SEXP super = PROTECT(GET_SLOT(from, install("super"))),
			pi = PROTECT(GET_SLOT(from, install("pi"))),
			px = PROTECT(GET_SLOT(from, install("px"))),
			s = PROTECT(GET_SLOT(from, install("s")));
		L->super = INTEGER(super);
		L->pi = INTEGER(pi);
		L->px = INTEGER(px);
		L->s = INTEGER(s);
		L->nsuper = LENGTH(super) - 1;
		L->ssize = ((int *) L->pi)[L->nsuper];
		L->xsize = ((int *) L->px)[L->nsuper];
		L->maxcsize = INTEGER(type)[4];
		L->maxesize = INTEGER(type)[5];
		UNPROTECT(4);
	} else {
		L->is_ll = INTEGER(type)[1];
		L->is_monotonic = INTEGER(type)[3];
		if (class[0] != 'n') {
			SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym)),
				i = PROTECT(GET_SLOT(from, Matrix_iSym)),
				nz = PROTECT(GET_SLOT(from, install("nz"))),
				nxt = PROTECT(GET_SLOT(from, install("nxt"))),
				prv = PROTECT(GET_SLOT(from, install("prv")));
			L->p = INTEGER(p);
			L->i = INTEGER(i);
			L->nz = INTEGER(nz);
			L->next = INTEGER(nxt);
			L->prev = INTEGER(prv);
			L->nzmax = ((int *) L->p)[L->n];
			UNPROTECT(5);
		}
	}
	L->itype = CHOLMOD_INT;
	L->dtype = CHOLMOD_DOUBLE;
	if (class[0] != 'n') {
		SEXP x = GET_SLOT(from, Matrix_xSym);
		switch (TYPEOF(x)) {
		case CPLXSXP:
			L->x = COMPLEX(x);
			L->xtype = CHOLMOD_COMPLEX;
			break;
		case REALSXP:
			L->x = REAL(x);
			L->xtype = CHOLMOD_REAL;
			break;
		default:
			ERROR_INVALID_TYPE(x, __func__);
			break;
		}
	}

	if (!cholmod_check_factor(L, &c))
		error(_("'%s' failed in '%s'"), "cholmod_check_factor", __func__);
	UNPROTECT(4);
	return L;
}

/**
 * Coerce from CsparseMatrix to (cholmod_sparse *)
 *
 * Sets the members of a pointed-to cholmod_sparse struct, using "data"
 * obtained from slots of a CsparseMatrix.  The result should _not_ be
 * freed using cholmod_free_sparse, as the resulting members point to
 * memory controlled by R, not by CHOLMOD.
 *
 * @param A a pointer to a cholmod_sparse struct, to be modified in-place.
 * @param from an S4 object inheriting from virtual class CsparseMatrix.
 * @param checkUnit a boolean indicating if the unit diagonal of formally
 *     unit triangular CsparseMatrix should be allocated.
 * @param sortInPlace a boolean indicating if unsorted CsparseMatrix
 *     should be sorted in place to avoid an allocation.
 *
 * @return A.
 */
/* NB: mostly parallel to M2CHS in ./cholmod-etc.c */
cholmod_sparse *sexp_as_cholmod_sparse(cholmod_sparse *A, SEXP from,
                                       Rboolean checkUnit, Rboolean sortInPlace)
{
	/* MJ: Do users really ever pass invalid 'from' ... ?
	       If not, then the code here could be simplified tremendously ...
	*/

	static const char *valid[] = {
		"dgCMatrix", "dsCMatrix", "dtCMatrix",
		"lgCMatrix", "lsCMatrix", "ltCMatrix",
		"ngCMatrix", "nsCMatrix", "ntCMatrix", "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *class = valid[ivalid];
	memset(A, 0, sizeof(cholmod_sparse));

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym)),
	     i = PROTECT(GET_SLOT(from, Matrix_iSym)),
	   cpi = PROTECT(checkpi(p, i, m, n));
	if (TYPEOF(cpi) != LGLSXP)
		error(_("'%s' failed in '%s': %s"),
		      "checkpi", __func__, CHAR(STRING_ELT(cpi, 0)));
	int *pp = INTEGER(p), *pi = INTEGER(i), sorted = LOGICAL(cpi)[0];
	size_t np = (size_t) XLENGTH(p), ni = (size_t) XLENGTH(i);
	if (!sorted && !sortInPlace) {
		int *tmp;
		tmp = (int *) R_alloc(np, sizeof(int));
		memcpy(tmp, pp, np * sizeof(int));
		pp = tmp;
		tmp = (int *) R_alloc(ni, sizeof(int));
		memcpy(tmp, pi, ni * sizeof(int));
		pi = tmp;
	}

	A->nrow = m;
	A->ncol = n;
	A->p = pp;
	A->i = pi;
	A->nzmax = ni;
	A->stype = 0;
	A->itype = CHOLMOD_INT;
	A->xtype = CHOLMOD_PATTERN;
	A->dtype = CHOLMOD_DOUBLE;
	A->sorted = LOGICAL(cpi)[0];
	A->packed = 1;

	if (ni > pp[n]) { /* overallocated */
		A->packed = 0;
		int j, *tmp = (int *) R_alloc(n, sizeof(int));
		for (j = 0; j < n; ++j)
			tmp[j] = pp[j + 1] - pp[j];
		A->nz = tmp;
	}
	if (class[1] == 's') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		A->stype = (ul == 'U') ? 1 : -1;
	}
	if (class[0] != 'n') {
		SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
		size_t nx = (size_t) XLENGTH(x);
		switch (class[0]) {
		case 'l':
		case 'i':
		{
			int *px = (TYPEOF(x) == LGLSXP) ? LOGICAL(x) : INTEGER(x);
			double *rtmp = (double *) R_alloc(nx, sizeof(double));
			for (size_t ix = 0; ix < nx; ++ix)
				rtmp[ix] = (px[ix] == NA_INTEGER)
					? NA_REAL : (double) px[ix];
			A->x = rtmp;
			A->xtype = CHOLMOD_REAL;
			break;
		}
		case 'd':
		{
			double *px = REAL(x);
			if (!sorted && !sortInPlace) {
				double *rtmp = (double *) R_alloc(nx, sizeof(double));
				memcpy(rtmp, px, nx * sizeof(double));
				px = rtmp;
			}
			A->x = px;
			A->xtype = CHOLMOD_REAL;
			break;
		}
		case 'z':
		{
			Rcomplex *px = COMPLEX(x);
			if (!sorted && !sortInPlace) {
				Rcomplex *rtmp = (Rcomplex *) R_alloc(nx, sizeof(Rcomplex));
				memcpy(rtmp, px, nx * sizeof(Rcomplex));
				px = rtmp;
			}
			A->x = px;
			A->xtype = CHOLMOD_COMPLEX;
			break;
		}
		default:
			break;
		}
		UNPROTECT(1); /* x */
	}
	if (!sorted && !cholmod_sort(A, &c))
		error(_("'%s' failed in '%s'"), "cholmod_sort", __func__);
	if (checkUnit && class[1] == 't' && n > 0) {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		char di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N') {
			double one[] = { 1.0, 0.0 };
			cholmod_sparse
				*eye = cholmod_speye(n, n, A->xtype, &c),
				*A1a = cholmod_add(A, eye, one, one, 1, 1, &c);
			memcpy(A, A1a, sizeof(cholmod_sparse));
			A->p = (int *) R_alloc(A1a->ncol + 1, sizeof(int));
			memcpy(A->p, A1a->p, (A1a->ncol + 1) * sizeof(int));
			A->i = (int *) R_alloc(A1a->nzmax, sizeof(int));
			memcpy(A->i, A1a->i, A1a->nzmax * sizeof(int));
			if (A1a->xtype != CHOLMOD_PATTERN) {
				size_t size = (A1a->xtype == CHOLMOD_REAL)
					? sizeof(double) : sizeof(Rcomplex);
				A->x = R_alloc(A1a->nzmax, size);
				memcpy(A->x, A1a->x, A1a->nzmax * size);
			}
			cholmod_free_sparse(&eye, &c);
			cholmod_free_sparse(&A1a, &c);
		}
	}

	UNPROTECT(3); /* cpi, i, p */
	return A;
}

/**
 * Coerce from TsparseMatrix to (cholmod_triplet *)
 *
 * Sets the members of a pointed-to cholmod_triplet struct, using "data"
 * obtained from slots of a TsparseMatrix.  The result should _not_ be
 * freed using cholmod_free_sparse, as the resulting members point to
 * memory controlled by R, not by CHOLMOD.
 *
 * @param A a pointer to a cholmod_triplet struct, to be modified in-place.
 * @param from an S4 object inheriting from virtual class TsparseMatrix.
 * @param checkUnit a boolean indicating if the unit diagonal of formally
 *     unit triangular TsparseMatrix should be allocated.
 *
 * @return A.
 */
cholmod_triplet *sexp_as_cholmod_triplet(cholmod_triplet *A, SEXP from,
                                         Rboolean checkUnit)
{
	static const char *valid[] = {
		"dgTMatrix", "dsTMatrix", "dtTMatrix",
		"lgTMatrix", "lsTMatrix", "ltTMatrix",
		"ngTMatrix", "nsTMatrix", "ntTMatrix", "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *class = valid[ivalid];
	memset(A, 0, sizeof(cholmod_triplet));

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	SEXP i = PROTECT(GET_SLOT(from, Matrix_pSym)),
		j = PROTECT(GET_SLOT(from, Matrix_iSym));
	int *pi = INTEGER(i), *pj = INTEGER(j);
	size_t nnz0 = (size_t) XLENGTH(i), nnz1 = nnz0;

	if (checkUnit && class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		char di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N')
			nnz1 += n;
	}

	if (nnz0 < nnz1) {
		int *tmp;
		tmp = (int *) R_alloc(nnz1, sizeof(int));
		memcpy(tmp, pi, nnz1 * sizeof(int));
		pi = tmp;
		tmp = (int *) R_alloc(nnz1, sizeof(int));
		memcpy(tmp, pj, nnz1 * sizeof(int));
		pj = tmp;

		pi += nnz0; pj += nnz0;
		for (int d = 0; d < n; ++d)
			*(pi++) = *(pj++) = d;
		pi -= nnz1; pj -= nnz1;
	}

	A->nrow = m;
	A->ncol = n;
	A->i = pi;
	A->j = pj;
	A->nzmax = nnz1;
	A->nnz = nnz1;
	A->stype = 0;
	A->itype = CHOLMOD_INT;
	A->xtype = CHOLMOD_PATTERN;
	A->dtype = CHOLMOD_DOUBLE;

	if (class[1] == 's') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		A->stype = (ul == 'U') ? 1 : -1;
	}
	if (class[0] != 'n') {
		SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
		switch (class[0]) {
		case 'l':
		case 'i':
		{
			int *px = (TYPEOF(x) == LGLSXP) ? LOGICAL(x) : INTEGER(x);
			double *rtmp = (double *) R_alloc(nnz1, sizeof(double));
			for (size_t k =    0; k < nnz0; ++k)
				rtmp[k] = (px[k] == NA_INTEGER)
					? NA_REAL : (double) px[k];
			for (size_t k = nnz0; k < nnz1; ++k)
				rtmp[k] = 1.0;
			A->x = rtmp;
			A->xtype = CHOLMOD_REAL;
			break;
		}
		case 'd':
		{
			double *px = REAL(x);
			if (nnz0 < nnz1) {
				double *rtmp = (double *) R_alloc(nnz1, sizeof(double));
				memcpy(rtmp, px, nnz0 * sizeof(double));
				for (size_t k = nnz0; k < nnz1; ++k)
					rtmp[k] = 1.0;
				px = rtmp;
			}
			A->x = px;
			A->xtype = CHOLMOD_REAL;
			break;
		}
		case 'z':
		{
			Rcomplex *px = COMPLEX(x);
			if (nnz0 < nnz1) {
				Rcomplex *rtmp = (Rcomplex *) R_alloc(nnz1, sizeof(Rcomplex));
				memcpy(rtmp, px, nnz0 * sizeof(Rcomplex));
				for (size_t k = nnz0; k < nnz1; ++k)
					rtmp[k] = Matrix_zone;
				px = rtmp;
			}
			A->x = px;
			A->xtype = CHOLMOD_COMPLEX;
			break;
		}
		default:
			break;
		}
		UNPROTECT(1); /* x */
	}

	UNPROTECT(2); /* j, i */
	return A;
}

/**
 * Coerce from [nlidz]geMatrix or vector to (cholmod_dense *)
 *
 * Sets the members of a pointed-to cholmod_dense struct, using "data"
 * obtained from slots of a [nlidz]geMatrix.  The result should _not_ be
 * freed using cholmod_free_dense, as the resulting members point to
 * memory controlled by R, not by CHOLMOD.
 *
 * @param A a pointer to a cholmod_dense struct, to be modified in-place.
 * @param from an S4 object inheriting from class [nlidz]geMatrix _or_
 *     a traditional vector of type "logical", "integer", "double", or
 *     "complex" (to be handled as a 1-column matrix if not a matrix).
 *
 * @return A.
 */
/* NB: mostly parallel to M2CHD in ./cholmod-etc.c */
cholmod_dense *sexp_as_cholmod_dense(cholmod_dense *A, SEXP from)
{
	static const char *valid[] = {
		"dgeMatrix", "lgeMatrix", "ngeMatrix", "" };
	int ivalid = R_check_class_etc(from, valid);
	memset(A, 0, sizeof(cholmod_dense));

	int m, n;
	if (ivalid < 0) {
		switch (TYPEOF(from)) {
		case LGLSXP:
		case INTSXP:
		case REALSXP:
		case CPLXSXP:
			break;
		default:
			ERROR_INVALID_TYPE(from, __func__);
			break;
		}
		SEXP dim = getAttrib(from, R_DimSymbol);
		if (TYPEOF(dim) == INTSXP && LENGTH(dim) == 2) {
			m = INTEGER(dim)[0];
			n = INTEGER(dim)[1];
		} else {
			m = LENGTH(from);
			n = 1;
		}
	} else {
		SEXP dim = GET_SLOT(from, Matrix_DimSym);
		m = INTEGER(dim)[0];
		n = INTEGER(dim)[1];
		from = GET_SLOT(from, Matrix_xSym);
	}

	A->nrow = m;
	A->ncol = n;
	A->nzmax = A->nrow * A->ncol;
	A->d = A->nrow;
	A->dtype = CHOLMOD_DOUBLE;

	size_t nx = (size_t) XLENGTH(from);
	switch (TYPEOF(from)) {
	case LGLSXP:
	case INTSXP:
	{
		int *px = (TYPEOF(from) == LGLSXP) ? LOGICAL(from) : INTEGER(from),
			pattern = ivalid == 2;
		double *rtmp = (double *) R_alloc(nx + 1, sizeof(double));
		for (size_t ix = 0; ix < nx; ++ix)
			rtmp[ix] = (px[ix] == NA_INTEGER)
				? ((pattern) ? 1.0 : NA_REAL) : (double) px[ix];
		A->x = rtmp;
		A->xtype = CHOLMOD_REAL;
		break;
	}
	case REALSXP:
		A->x = REAL(from);
		A->xtype = CHOLMOD_REAL;
		break;
	case CPLXSXP:
		A->x = COMPLEX(from);
		A->xtype = CHOLMOD_COMPLEX;
		break;
	default:
		ERROR_INVALID_TYPE(from, __func__);
		break;
	}

	return A;
}

/**
 * Coerce from (double *) to (cholmod_dense *) with given dimensions
 *
 * An analogue of base::matrix(data, nrow, ncol),
 * where typeof(data)=="double" and length(data)==nrow*ncol.
 *
 * @param A a pointer to a cholmod_dense struct, to be modified in-place.
 * @param data a pointer to an nrow*ncol*sizeof(double) block of memory.
 * @param nrow the desired number of rows.
 * @param ncol the desired number of columns.
 *
 * @return A.
 */
cholmod_dense *numeric_as_cholmod_dense(cholmod_dense *A,
                                        double *data, int nrow, int ncol)
{
	memset(A, 0, sizeof(cholmod_dense));
	A->nrow = nrow;
	A->ncol = ncol;
	A->nzmax = A->nrow * A->ncol;
	A->d = A->nrow;
	A->x = data;
	A->xtype = CHOLMOD_REAL;
	A->dtype = CHOLMOD_DOUBLE;
	return A;
}

/**
 * Coerce from (cholmod_factor *) to CHMfactor
 *
 * Allocates an S4 object inheriting from virtual class CHMfactor
 * and copies into the slots from members of a pointed-to cholmod_factor
 * struct.  The specific class of the result is determined by struct
 * members xtype and is_super.
 *
 * @param L a pointer to a cholmod_factor struct.
 * @param doFree a flag indicating if and how to free L before returning.
 *     (0) don't free, (>0) free with cholmod_free_factor, (<0) free with
 *     R_Free.
 *
 * @return A CHMfactor.
 */
/* NB: mostly parallel to CHF2M in ./cholmod-etc.c */
SEXP cholmod_factor_as_sexp(cholmod_factor *L, int doFree)
{

#define FREE_THEN(_EXPR_) \
	do { \
		if (doFree != 0) { \
			if (doFree < 0) \
				R_Free(L); \
			else if (L->itype == CHOLMOD_INT) \
				cholmod_free_factor(&L, &c); \
			else \
				cholmod_l_free_factor(&L, &cl); \
			_EXPR_; \
		} \
	} while (0)

	if (L->itype != CHOLMOD_INT)
		FREE_THEN(error(_("wrong '%s'"), "itype"));
	if (L->xtype != CHOLMOD_PATTERN &&
	    L->xtype != CHOLMOD_REAL && L->xtype != CHOLMOD_COMPLEX)
		FREE_THEN(error(_("wrong '%s'"), "xtype"));
	if (L->dtype != CHOLMOD_DOUBLE)
		FREE_THEN(error(_("wrong '%s'"), "dtype"));
	if (L->n > INT_MAX)
		FREE_THEN(error(_("dimensions cannot exceed %s"), "2^31-1"));
	if (L->super) {
		if (L->maxcsize > INT_MAX)
			FREE_THEN(error(_("'%s' would overflow type \"%s\""),
			                "maxcsize", "integer"));
	} else {
		if (L->n == INT_MAX)
			FREE_THEN(error(_("n+1 would overflow type \"%s\""),
			                "integer"));
	}
	if (L->minor < L->n) {
		if (L->is_ll)
			FREE_THEN(error(_("leading principal minor of order %d is not positive"),
			                (int) L->minor + 1));
		else
			FREE_THEN(error(_("leading principal minor of order %d is zero"),
			                (int) L->minor + 1));
	}
	char class[] = ".CHM.....";
	class[0] = (L->xtype == CHOLMOD_PATTERN)
		? 'n' : ((L->xtype == CHOLMOD_COMPLEX) ? 'z' : 'd');
	memcpy(class + 4, (L->is_super) ? "super" : "simpl", 5);
	SEXP to = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(to, Matrix_DimSym));
	INTEGER(dim)[0] = INTEGER(dim)[1] = (int) L->n;
	if (L->ordering != CHOLMOD_NATURAL) {
		SEXP perm = PROTECT(allocVector(INTSXP, L->n));
		memcpy(INTEGER(perm), L->Perm, L->n * sizeof(int));
		SET_SLOT(to, Matrix_permSym, perm);
		UNPROTECT(1);
	}
	SEXP type = PROTECT(allocVector(INTSXP, 6)),
		colcount = PROTECT(allocVector(INTSXP, L->n));
	INTEGER(type)[0] = L->ordering;
	INTEGER(type)[1] = (L->is_super) ? 1 : L->is_ll;
	INTEGER(type)[2] = (L->is_super) ? 1 : 0;
	INTEGER(type)[3] = (L->is_super) ? 1 : L->is_monotonic;
	INTEGER(type)[4] = (L->is_super) ? (int) L->maxcsize : 0;
	INTEGER(type)[5] = (L->is_super) ? (int) L->maxesize : 0;
	memcpy(INTEGER(colcount), L->ColCount, L->n * sizeof(int));
	SET_SLOT(to, install("type"), type);
	SET_SLOT(to, install("colcount"), colcount);
	if (L->is_super) {
		SEXP super = PROTECT(allocVector(INTSXP, L->nsuper + 1)),
			pi = PROTECT(allocVector(INTSXP, L->nsuper + 1)),
			px = PROTECT(allocVector(INTSXP, L->nsuper + 1)),
			s = PROTECT(allocVector(INTSXP, L->ssize));
		memcpy(INTEGER(super), L->super, (L->nsuper + 1) * sizeof(int));
		memcpy(INTEGER(pi), L->pi, (L->nsuper + 1) * sizeof(int));
		memcpy(INTEGER(px), L->px, (L->nsuper + 1) * sizeof(int));
		memcpy(INTEGER(s), L->s, L->ssize * sizeof(int));
		SET_SLOT(to, install("super"), super);
		SET_SLOT(to, install("pi"), pi);
		SET_SLOT(to, install("px"), px);
		SET_SLOT(to, install("s"), s);
		UNPROTECT(4);
	} else if (L->xtype != CHOLMOD_PATTERN) {
		SEXP p = PROTECT(allocVector(INTSXP, L->n + 1)),
			i = PROTECT(allocVector(INTSXP, L->nzmax)),
			nz = PROTECT(allocVector(INTSXP, L->n)),
			nxt = PROTECT(allocVector(INTSXP, L->n + 2)),
			prv = PROTECT(allocVector(INTSXP, L->n + 2));
		memcpy(INTEGER(p), L->p, (L->n + 1) * sizeof(int));
		memcpy(INTEGER(i), L->i, L->nzmax * sizeof(int));
		memcpy(INTEGER(nz), L->nz, L->n * sizeof(int));
		memcpy(INTEGER(nxt), L->next, (L->n + 2) * sizeof(int));
		memcpy(INTEGER(prv), L->prev, (L->n + 2) * sizeof(int));
		SET_SLOT(to, Matrix_pSym, p);
		SET_SLOT(to, Matrix_iSym, i);
		SET_SLOT(to, install("nz"), nz);
		SET_SLOT(to, install("nxt"), nxt);
		SET_SLOT(to, install("prv"), prv);
		UNPROTECT(5);
	}
	if (L->xtype != CHOLMOD_PATTERN) {
		SEXP x;
		R_xlen_t nx = (R_xlen_t) ((L->is_super) ? L->xsize : L->nzmax);
		if (L->xtype == CHOLMOD_COMPLEX) {
			PROTECT(x = allocVector(CPLXSXP, nx));
			memcpy(COMPLEX(x), L->x, (size_t) nx * sizeof(Rcomplex));
		} else {
			PROTECT(x = allocVector(REALSXP, nx));
			memcpy(REAL(x), L->x, (size_t) nx * sizeof(double));
		}
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1);
	}
	FREE_THEN();

#undef FREE_THEN

	UNPROTECT(4);
	return to;
}

/**
 * Coerce from (cholmod_sparse *) to CsparseMatrix
 *
 * Allocates an S4 object inheriting from virtual class CsparseMatrix
 * and copies into the slots from members of a pointed-to cholmod_sparse
 * struct.  The specific class of the result is determined by struct
 * members xtype and stype and by arguments ttype and doLogic.
 *
 * @param A a pointer to a cholmod_sparse struct.
 * @param doFree a flag indicating if and how to free A before returning.
 *     (0) don't free, (>0) free with cholmod_free_sparse, (<0) free with
 *     R_Free.
 * @param ttype a flag indicating if the result should be a .tCMatrix.
 *     (0) not .tCMatrix, (>0) .tCMatrix with uplo="U", (<0) .tCMatrix
 *     with uplo="L".  If ttype=0, then the result is a .gCMatrix or
 *     .sCMatrix depending on stype.  (0) .gCMatrix, (>0) .sCMatrix with
 *     uplo="U", (<0) .sCMatrix with uplo="L".
 * @param doLogic a flag indicating if the result should be a l.CMatrix
 *     if xtype=CHOLMOD_REAL.
 * @param diagString a null-terminated string or NULL.  The diag slot
 *     of a .tCMatrix result is "N" if and only if diagString is NULL
 *     or diagString[0] is 'N'.
 * @param dimnames an R object specifying the Dimnames slot of the result,
 *     unused if not a list of length 2.
 *
 * @return A CsparseMatrix.
 */
/* NB: mostly parallel to CHS2M in ./cholmod-etc.c */
SEXP cholmod_sparse_as_sexp(cholmod_sparse *A, int doFree,
                            int ttype, int doLogic, const char *diagString,
                            SEXP dimnames)
{

#define FREE_THEN(_EXPR_) \
	do { \
		if (doFree != 0) { \
			if (doFree < 0) \
				R_Free(A); \
			else if (A->itype == CHOLMOD_INT) \
				cholmod_free_sparse(&A, &c); \
			else \
				cholmod_l_free_sparse(&A, &cl); \
			_EXPR_; \
		} \
	} while (0)

	if (A->itype != CHOLMOD_INT)
		FREE_THEN(error(_("wrong '%s'"), "itype"));
	if (A->xtype != CHOLMOD_PATTERN &&
	    A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		FREE_THEN(error(_("wrong '%s'"), "xtype"));
	if (A->dtype != CHOLMOD_DOUBLE)
		FREE_THEN(error(_("wrong '%s'"), "dtype"));
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		FREE_THEN(error(_("dimensions cannot exceed %s"), "2^31-1"));
	if (A->stype != 0 || !A->sorted || !A->packed)
		cholmod_sort(A, &c);
	int m = (int) A->nrow, n = (int) A->ncol, nnz = ((int *) A->p)[A->ncol];
	R_xlen_t n1a = (R_xlen_t) n + 1;
	char class[] = "..CMatrix";
	class[0] = (A->xtype == CHOLMOD_PATTERN)
		? 'n' : ((A->xtype == CHOLMOD_COMPLEX) ? 'z' : ((doLogic) ? 'l' : 'd'));
	class[1] = (ttype != 0) ? 't' : ((A->stype != 0) ? 's' : 'g');
	SEXP to = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(to, Matrix_DimSym)),
		p = PROTECT(allocVector(INTSXP, n1a)),
		i = PROTECT(allocVector(INTSXP, nnz));
	INTEGER(dim)[0] = m;
	INTEGER(dim)[1] = n;
	memcpy(INTEGER(p), A->p, (size_t) n1a * sizeof(int));
	memcpy(INTEGER(i), A->i, (size_t) nnz * sizeof(int));
	SET_SLOT(to, Matrix_pSym, p);
	SET_SLOT(to, Matrix_iSym, i);
	if (A->xtype != CHOLMOD_PATTERN) {
		SEXP x;
		if (A->xtype == CHOLMOD_COMPLEX) {
			PROTECT(x = allocVector(CPLXSXP, nnz));
			memcpy(COMPLEX(x), A->x, (size_t) nnz * sizeof(Rcomplex));
		} else if (!doLogic) {
			PROTECT(x = allocVector(REALSXP, nnz));
			memcpy(REAL(x), A->x, (size_t) nnz * sizeof(double));
		} else {
			PROTECT(x = allocVector(LGLSXP, nnz));
			int k, *px = LOGICAL(x);
			double *py = (double *) A->x;
			for (k = 0; k < nnz; ++k)
				px[k] = (ISNAN(py[k])) ? NA_LOGICAL : (py[k] != 0.0);
		}
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1);
	}
	if (ttype < 0 || A->stype < 0) {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1);
	}
	if (ttype != 0 && diagString && diagString[0] != 'N') {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1);
	}
	if (TYPEOF(dimnames) == VECSXP && LENGTH(dimnames) == 2)
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	FREE_THEN();

#undef FREE_THEN

	UNPROTECT(4);
	return to;
}

/**
 * Coerce from (cholmod_triplet *) to TsparseMatrix
 *
 * Allocates an S4 object inheriting from virtual class TsparseMatrix
 * and copies into the slots from members of a pointed-to cholmod_triplet
 * struct.  The specific class of the result is determined by struct
 * members xtype and stype and by arguments ttype and doLogic.
 *
 * @param A a pointer to a cholmod_triplet struct.
 * @param doFree a flag indicating if and how to free A before returning.
 *     (0) don't free, (>0) free with cholmod_free_triplet, (<0) free with
 *     R_Free.
 * @param ttype a flag indicating if the result should be a .tTMatrix.
 *     (0) not .tTMatrix, (>0) .tTMatrix with uplo="U", (<0) .tTMatrix
 *     with uplo="L".  If ttype=0, then the result is a .gTMatrix or
 *     .sTMatrix depending on stype.  (0) .gTMatrix, (>0) .sTMatrix with
 *     uplo="U", (<0) .sTMatrix with uplo="L".
 * @param doLogic a flag indicating if the result should be an l.TMatrix
 *     if xtype=CHOLMOD_REAL.
 * @param diagString a null-terminated string or NULL.  The diag slot
 *     of a .tTMatrix result is "N" if and only if diagString is NULL
 *     or diagString[0] is 'N'.
 * @param dimnames an R object specifying the Dimnames slot of the result,
 *     unused if not a list of length 2.
 *
 * @return A TsparseMatrix.
 */
SEXP cholmod_triplet_as_sexp(cholmod_triplet *A, int doFree,
                             int ttype, int doLogic, const char *diagString,
                             SEXP dimnames)
{

#define FREE_THEN(_EXPR_) \
	do { \
		if (doFree != 0) { \
			if (doFree < 0) \
				R_Free(A); \
			else if (A->itype == CHOLMOD_INT) \
				cholmod_free_triplet(&A, &c); \
			else \
				cholmod_l_free_triplet(&A, &cl); \
			_EXPR_; \
		} \
	} while (0)

	if (A->itype != CHOLMOD_INT)
		FREE_THEN(error(_("wrong '%s'"), "itype"));
	if (A->xtype != CHOLMOD_PATTERN &&
	    A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		FREE_THEN(error(_("wrong '%s'"), "xtype"));
	if (A->dtype != CHOLMOD_DOUBLE)
		FREE_THEN(error(_("wrong '%s'"), "dtype"));
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		FREE_THEN(error(_("dimensions cannot exceed %s"), "2^31-1"));
	int m = (int) A->nrow, n = (int) A->ncol;
	R_xlen_t nnz = (R_xlen_t) A->nnz;
	char class[] = "..TMatrix";
	class[0] = (A->xtype == CHOLMOD_PATTERN)
		? 'n' : ((A->xtype == CHOLMOD_COMPLEX) ? 'z' : ((doLogic) ? 'l' : 'd'));
	class[1] = (ttype != 0) ? 't' : ((A->stype != 0) ? 's' : 'g');
	SEXP to = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(to, Matrix_DimSym)),
		i = PROTECT(allocVector(INTSXP, nnz)),
		j = PROTECT(allocVector(INTSXP, nnz));
	INTEGER(dim)[0] = m;
	INTEGER(dim)[1] = n;
	memcpy(INTEGER(i), A->i, (size_t) nnz * sizeof(int));
	memcpy(INTEGER(j), A->j, (size_t) nnz * sizeof(int));
	SET_SLOT(to, Matrix_iSym, i);
	SET_SLOT(to, Matrix_jSym, j);
	if (A->xtype != CHOLMOD_PATTERN) {
		SEXP x;
		if (A->xtype == CHOLMOD_COMPLEX) {
			PROTECT(x = allocVector(CPLXSXP, nnz));
			memcpy(COMPLEX(x), A->x, (size_t) nnz * sizeof(Rcomplex));
		} else if (!doLogic) {
			PROTECT(x = allocVector(REALSXP, nnz));
			memcpy(REAL(x), A->x, (size_t) nnz * sizeof(double));
		} else {
			PROTECT(x = allocVector(LGLSXP, nnz));
			int k, *px = LOGICAL(x);
			double *py = (double *) A->x;
			for (k = 0; k < nnz; ++k)
				px[k] = (ISNAN(py[k])) ? NA_LOGICAL : (py[k] != 0.0);
		}
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1);
	}
	if (ttype < 0 || A->stype < 0) {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1);
	}
	if (ttype != 0 && diagString && diagString[0] != 'N') {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1);
	}
	if (TYPEOF(dimnames) == VECSXP && LENGTH(dimnames) == 2)
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	FREE_THEN();

#undef FREE_THEN

	UNPROTECT(4);
	return to;
}

/**
 * Coerce from (cholmod_dense *) to [dz]geMatrix
 *
 * Allocates an S4 object of class [dz]geMatrix
 * and copies into the slots from members of a pointed-to cholmod_dense
 * struct.  The specific class of the result is determined by struct
 * member xtype.
 *
 * @param A a pointer to a cholmod_dense struct.
 * @param doFree a flag indicating if and how to free A before returning.
 *     (0) don't free, (>0) free with cholmod_free_dense, (<0) free with
 *     R_Free.
 *
 * @return A [dz]geMatrix.
 */
/* NB: mostly parallel to CHD2M in ./cholmod-etc.c */
SEXP cholmod_dense_as_sexp(cholmod_dense *A, int doFree)
{

#define FREE_THEN(_EXPR_) \
	do { \
		if (doFree != 0) { \
			if (doFree < 0) \
				R_Free(A); \
			else \
				cholmod_free_dense(&A, &c); \
			_EXPR_; \
		} \
	} while (0)

	if (A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		FREE_THEN(error(_("wrong '%s'"), "xtype"));
	if (A->dtype != CHOLMOD_DOUBLE)
		FREE_THEN(error(_("wrong '%s'"), "dtype"));
	if (A->d != A->nrow) /* MJ: currently no need to support this case */
		FREE_THEN(error(_("leading dimension not equal to number of rows")));
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		FREE_THEN(error(_("dimensions cannot exceed %s"), "2^31-1"));
	int m = (int) A->nrow, n = (int) A->ncol;
	if ((Matrix_int_fast64_t) m * n > R_XLEN_T_MAX)
		FREE_THEN(error(_("attempt to allocate vector of length exceeding %s"),
		                "R_XLEN_T_MAX"));
	char class[] = ".geMatrix";
	class[0] = (A->xtype == CHOLMOD_COMPLEX) ? 'z' : 'd';
	SEXP to = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(to, Matrix_DimSym));
	INTEGER(dim)[0] = m;
	INTEGER(dim)[1] = n;
	SEXP x;
	if (A->xtype == CHOLMOD_COMPLEX) {
		PROTECT(x = allocVector(CPLXSXP, (R_xlen_t) m * n));
		memcpy(COMPLEX(x), A->x, (size_t) m * n * sizeof(Rcomplex));
	} else {
		PROTECT(x = allocVector(REALSXP, (R_xlen_t) m * n));
		memcpy(REAL(x), A->x, (size_t) m * n * sizeof(double));
	}
	SET_SLOT(to, Matrix_xSym, x);
	FREE_THEN();

#undef FREE_THEN

	UNPROTECT(3);
	return to;
}

/**
 * Log determinant from Cholesky factorization
 *
 * Computes log(det(A)) given the Cholesky factorization of A as
 * P1 * A * P1' = L1 * D * L1' = L * L', L = L1 * sqrt(D).  The
 * result is computed as sum(log(diag(D))) or 2*sum(log(diag(L))),
 * depending on members is_super and is_ll of the supplied struct.
 * Note that CHOLMOD does not require diag(D) to be positive and
 * that this routine does not check (FIXME).
 *
 * @param L a pointer to a cholmod_factor struct.  It is assumed that
 *     L->xtype=CHOLMOD_REAL.
 */
double cholmod_factor_ldetA(cholmod_factor *L)
{
	int i, j, p;
	double ans = 0;
	if (L->is_super) {
		int *lpi = (int *) L->pi, *lsup = (int *) L->super;
		for (i = 0; i < L->nsuper; i++) {
			int nrp1 = 1 + lpi[i + 1] - lpi[i],
				nc = lsup[i + 1] - lsup[i];
			double *x = (double *) L->x + ((int *) L->px)[i];
			for (R_xlen_t jn = 0, j = 0; j < nc; j++, jn += nrp1)
				ans += 2.0 * log(fabs(x[jn]));
		}
	} else {
		int *li = (int *) L->i, *lp = (int *) L->p;
		double *lx = (double *) L->x;
		for (j = 0; j < L->n; j++) {
			for (p = lp[j]; li[p] != j && p < lp[j + 1]; p++)
				;
			if (li[p] != j) {
				error(_("invalid simplicial Cholesky factorization: structural zero on main diagonal in column %d"),
				      j);
				break;
			}
			ans += log(lx[p] * ((L->is_ll) ? lx[p] : 1.0));
		}
	}
	return ans;
}

/**
 * Update a Cholesky factorization
 *
 * Updates in-place the Cholesky factorization of a symmetric matrix
 * X+alpha*I with the Cholesky factorization of
 * (1) A+beta*I, where A is a symmetric matrix sharing the nonzero pattern
 * of X, or
 * (2) A*A'+beta*I, where A is a general matrix sharing the nonzero pattern
 * of Y, assuming that X = Y*Y'.
 *
 * @param L a pointer to a cholmod_factor struct, to be modified in-place.
 * @param A a pointer to a cholmod_sparse struct.
 * @param beta a multiplier, typically positive, to guarantee strict
 *     diagonal dominance.
 *
 * @return L.
 */
cholmod_factor *cholmod_factor_update(cholmod_factor *L, cholmod_sparse *A,
                                      double beta)
{
	int ll = L->is_ll;
	double z[2];
	z[0] = beta;
	z[1] = 0.0;
	if (!cholmod_factorize_p(A, z, NULL, 0, L, &c))
		error(_("'%s' failed in '%s'"), "cholmod_factorize_p", __func__);
	if (L->is_ll != ll &&
	    !cholmod_change_factor(L->xtype, ll, L->is_super, 1, 1, L, &c))
		error(_("'%s' failed in '%s'"), "cholmod_change_factor", __func__);
	return L;
}

#if 0
static
int R_cholmod_print_function(const char *fmt, ...)
{
	va_list(ap);
	va_start(ap, fmt);
	Rprintf((char *) fmt, ap);
	va_end(ap);
	return 0;
}
#endif

static
void R_cholmod_error_handler(int status, const char *file, int line,
                             const char *message)
{
	R_cholmod_common_envget();
	if (status < 0)
		error(_("CHOLMOD error '%s' at file '%s', line %d"),
		      message, file, line);
	else
		warning(_("CHOLMOD warning '%s' at file '%s', line %d"),
		        message, file, line);
}

int R_cholmod_start(cholmod_common *Common)
{
	int ans = cholmod_start(Common);
	if (!ans)
		error(_("'%s' failed in '%s'"), "cholmod_start", __func__);
#if 0
	/* No longer, with SuiteSparse 5.7.1 : */
	Common->print_function =
# if 0
		R_cholmod_print_function;
# else
		NULL;
# endif
#endif
	Common->error_handler = R_cholmod_error_handler;
	return ans;
}

int R_cholmod_finish(cholmod_common *Common)
{
	int ans = cholmod_finish(Common);
	if (!ans)
		error(_("'%s' failed in '%s'"), "cholmod_finish", __func__);
	return ans;
}

SEXP cholmod_common_env;

static
SEXP
	dboundSym,
	grow0Sym,
	grow1Sym,
	grow2Sym,
	maxrankSym,
	supernodal_switchSym,
	supernodalSym,
	final_asisSym,
	final_superSym,
	final_llSym,
	final_packSym,
	final_monotonicSym,
	final_resymbolSym,
	prefer_zomplexSym,
	prefer_upperSym,
	quick_return_if_not_posdefSym,
	nmethodsSym,
	postorderSym,
	m0_ordSym;

SEXP R_cholmod_common_envini(SEXP rho) {
	if (!isEnvironment(rho))
		ERROR_INVALID_TYPE(rho, __func__);
	cholmod_common_env = rho;
	dboundSym = install("dbound");
	grow0Sym = install("grow0");
	grow1Sym = install("grow1");
	grow2Sym = install("grow2");
	maxrankSym = install("maxrank");
	supernodal_switchSym = install("supernodal_switch");
	supernodalSym = install("supernodal");
	final_asisSym = install("final_asis");
	final_superSym = install("final_super");
	final_llSym = install("final_ll");
	final_packSym = install("final_pack");
	final_monotonicSym = install("final_monotonic");
	final_resymbolSym = install("final_resymbol");
	prefer_zomplexSym = install("final_zomplex");
	prefer_upperSym = install("final_upper");
	quick_return_if_not_posdefSym = install("quick_return_if_not_posdef");
	nmethodsSym = install("nmethods");
	postorderSym = install("postorder");
	m0_ordSym = install("m0.ord");
	R_cholmod_common_envset();
	return R_NilValue;
}

void R_cholmod_common_envset(void) {
	SEXP rho = cholmod_common_env, tmp;

#define SET_FRAME_FROM_MEMBER(_MEMBER_, _KIND_) \
	do { \
		PROTECT(tmp = Scalar ## _KIND_(c. _MEMBER_)); \
		defineVar(_MEMBER_ ## Sym, tmp, rho); \
		UNPROTECT(1); \
	} while (0)

	SET_FRAME_FROM_MEMBER(dbound,                     Real);
	SET_FRAME_FROM_MEMBER(grow0,                      Real);
	SET_FRAME_FROM_MEMBER(grow1,                      Real);
	SET_FRAME_FROM_MEMBER(grow2,                      Integer);
	SET_FRAME_FROM_MEMBER(maxrank,                    Integer);
	SET_FRAME_FROM_MEMBER(supernodal_switch,          Real);
	SET_FRAME_FROM_MEMBER(supernodal,                 Logical);
	SET_FRAME_FROM_MEMBER(final_asis,                 Logical);
	SET_FRAME_FROM_MEMBER(final_super,                Logical);
	SET_FRAME_FROM_MEMBER(final_ll,                   Logical);
	SET_FRAME_FROM_MEMBER(final_pack,                 Logical);
	SET_FRAME_FROM_MEMBER(final_monotonic,            Logical);
	SET_FRAME_FROM_MEMBER(final_resymbol,             Logical);
	SET_FRAME_FROM_MEMBER(prefer_zomplex,             Logical);
	SET_FRAME_FROM_MEMBER(prefer_upper,               Logical);
	SET_FRAME_FROM_MEMBER(quick_return_if_not_posdef, Logical);
	SET_FRAME_FROM_MEMBER(nmethods,                   Integer);
	SET_FRAME_FROM_MEMBER(postorder,                  Logical);

	PROTECT(tmp = ScalarInteger(c.method[0].ordering));
	defineVar(m0_ordSym, tmp, rho);
	UNPROTECT(1);

	return;
}

void R_cholmod_common_envget(void) {
	SEXP rho = cholmod_common_env, tmp;

#define GET_MEMBER_FROM_FRAME(_MEMBER_, _KIND_) \
	do { \
		PROTECT(tmp = findVarInFrame(rho, _MEMBER_ ## Sym)); \
		c. _MEMBER_ = as ## _KIND_(tmp); \
		UNPROTECT(1); \
	} while (0)

	GET_MEMBER_FROM_FRAME(dbound,                     Real);
	GET_MEMBER_FROM_FRAME(grow0,                      Real);
	GET_MEMBER_FROM_FRAME(grow1,                      Real);
	GET_MEMBER_FROM_FRAME(grow2,                      Integer);
	GET_MEMBER_FROM_FRAME(maxrank,                    Integer);
	GET_MEMBER_FROM_FRAME(supernodal_switch,          Real);
	GET_MEMBER_FROM_FRAME(supernodal,                 Logical);
	GET_MEMBER_FROM_FRAME(final_asis,                 Logical);
	GET_MEMBER_FROM_FRAME(final_super,                Logical);
	GET_MEMBER_FROM_FRAME(final_ll,                   Logical);
	GET_MEMBER_FROM_FRAME(final_pack,                 Logical);
	GET_MEMBER_FROM_FRAME(final_monotonic,            Logical);
	GET_MEMBER_FROM_FRAME(final_resymbol,             Logical);
	GET_MEMBER_FROM_FRAME(prefer_zomplex,             Logical);
	GET_MEMBER_FROM_FRAME(prefer_upper,               Logical);
	GET_MEMBER_FROM_FRAME(quick_return_if_not_posdef, Logical);
	GET_MEMBER_FROM_FRAME(nmethods,                   Integer);
	GET_MEMBER_FROM_FRAME(postorder,                  Logical);

	PROTECT(tmp = findVarInFrame(rho, m0_ordSym));
	c.method[0].ordering = asInteger(tmp);
	UNPROTECT(1);

	return;
}
