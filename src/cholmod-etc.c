#include "Mdefines.h"
#include "idz.h"
#include "cholmod-etc.h"

cholmod_common c ;
cholmod_common cl;

cholmod_factor *M2CHF(SEXP obj, int values)
{
	cholmod_factor *L = (cholmod_factor *) R_alloc(1, sizeof(cholmod_factor));
	memset(L, 0, sizeof(cholmod_factor));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		type = PROTECT(GET_SLOT(obj, install("type"))),
		perm = PROTECT(GET_SLOT(obj, Matrix_permSym)),
		colcount = PROTECT(GET_SLOT(obj, install("colcount"))),
		x = PROTECT(getAttrib(obj, Matrix_xSym));
	L->n = INTEGER(dim)[0];
	L->minor = L->n; /* FIXME: could be wrong for obj <- new(...) */
	L->ordering = INTEGER(type)[0];
	if (L->ordering != CHOLMOD_NATURAL)
		L->Perm = INTEGER(perm);
	else {
		/* cholmod_check_factor allows L->Perm == NULL,
		   but cholmod_copy_factor does not test, so it segfaults ...
		*/
		int n = (int) L->n, *Perm = (int *) R_alloc(L->n, sizeof(int));
		for (int j = 0; j < n; ++j)
			Perm[j] = j;
		L->Perm = Perm;
	}
	L->ColCount = INTEGER(colcount);
	L->is_super = INTEGER(type)[2];
	if (L->is_super) {
		L->is_ll = 1;
		L->is_monotonic = 1;
		SEXP super = PROTECT(GET_SLOT(obj, install("super"))),
			pi = PROTECT(GET_SLOT(obj, install("pi"))),
			px = PROTECT(GET_SLOT(obj, install("px"))),
			s = PROTECT(GET_SLOT(obj, install("s")));
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
		if (values && x != R_NilValue) {
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
			L->nzmax = ((int *) L->p)[L->n];
			UNPROTECT(5);
		}
	}
	L->itype = CHOLMOD_INT;
	L->dtype = CHOLMOD_DOUBLE;
	if (values && x != R_NilValue) {
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
	UNPROTECT(5);
	return L;
}

cholmod_sparse *M2CHS(SEXP obj, int values)
{
	cholmod_sparse *A = (cholmod_sparse *) R_alloc(1, sizeof(cholmod_sparse));
	memset(A, 0, sizeof(cholmod_sparse));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		x = PROTECT(getAttrib(obj, Matrix_xSym));
	A->nrow = INTEGER(dim)[0];
	A->ncol = INTEGER(dim)[1];
	A->p = INTEGER(p);
	A->i = INTEGER(i);
	A->nzmax = ((int *) A->p)[A->ncol];
	A->stype = 0;
	A->itype = CHOLMOD_INT;
	A->xtype = CHOLMOD_PATTERN;
	A->dtype = CHOLMOD_DOUBLE;
	A->sorted = 1;
	A->packed = 1;
	if (values && x != R_NilValue) {
		switch (TYPEOF(x)) {
		case CPLXSXP:
			A->x = COMPLEX(x);
			A->xtype = CHOLMOD_COMPLEX;
			break;
		case REALSXP:
			A->x = REAL(x);
			A->xtype = CHOLMOD_REAL;
			break;
		default:
			ERROR_INVALID_TYPE(x, __func__);
			break;
		}
	}
	UNPROTECT(4);
	return A;
}

cholmod_dense *M2CHD(SEXP obj, int trans)
{
	cholmod_dense *A = (cholmod_dense *) R_alloc(1, sizeof(cholmod_dense));
	memset(A, 0, sizeof(cholmod_dense));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int m = INTEGER(dim)[0], n = INTEGER(dim)[1];
	A->nrow = ((trans) ? n : m);
	A->ncol = ((trans) ? m : n);
	A->nzmax = A->nrow * A->ncol;
	A->d = A->nrow;
	A->dtype = CHOLMOD_DOUBLE;
	switch (TYPEOF(x)) {
	case CPLXSXP:
	{
		Rcomplex *px = COMPLEX(x);
		if (!trans)
			A->x = px;
		else {
			Rcomplex *py = R_Calloc(A->nzmax, Rcomplex);
			ztranspose2(py, px, m, n);
			A->x = py; /* NB: caller must do R_Free(A->x) */
		}
		A->xtype = CHOLMOD_COMPLEX;
		break;
	}
	case REALSXP:
	{
		double *px = REAL(x);
		if (!trans)
			A->x = px;
		else {
			double *py = R_Calloc(A->nzmax, double);
			dtranspose2(py, px, m, n);
			A->x = py; /* NB: caller must do R_Free(A->x) */
		}
		A->xtype = CHOLMOD_REAL;
		break;
	}
	default:
		ERROR_INVALID_TYPE(x, __func__);
		break;
	}
	UNPROTECT(2);
	return A;
}

SEXP CHF2M(cholmod_factor *L, int values)
{
	if (L->itype != CHOLMOD_INT)
		error(_("wrong '%s'"), "itype");
	if (values && L->xtype != CHOLMOD_REAL && L->xtype != CHOLMOD_COMPLEX)
		error(_("wrong '%s'"), "xtype");
	if (values && L->dtype != CHOLMOD_DOUBLE)
		error(_("wrong '%s'"), "dtype");
	if (L->n > INT_MAX)
		error(_("dimensions cannot exceed %s"), "2^31-1");
	if (L->super) {
		if (L->maxcsize > INT_MAX)
			error(_("'%s' would overflow type \"%s\""),
			      "maxcsize", "integer");
	} else {
		if (L->n == INT_MAX)
			error(_("n+1 would overflow type \"%s\""),
			      "integer");
	}
	if (L->minor < L->n) {
		if (L->is_ll)
			error(_("leading principal minor of order %d is not positive"),
			      (int) L->minor + 1);
		else
			error(_("leading principal minor of order %d is zero"),
			      (int) L->minor + 1);
	}
	char cl[] = ".CHM.....";
	cl[0] = (!values) ? 'n' : ((L->xtype == CHOLMOD_COMPLEX) ? 'z' : 'd');
	memcpy(cl + 4, (L->is_super) ? "super" : "simpl", 5);
	SEXP obj = PROTECT(newObject(cl)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	INTEGER(dim)[0] = INTEGER(dim)[1] = (int) L->n;
	if (L->ordering != CHOLMOD_NATURAL) {
		SEXP perm = PROTECT(allocVector(INTSXP, L->n));
		Matrix_memcpy(INTEGER(perm), L->Perm, L->n, sizeof(int));
		SET_SLOT(obj, Matrix_permSym, perm);
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
	Matrix_memcpy(INTEGER(colcount), L->ColCount, L->n, sizeof(int));
	SET_SLOT(obj, install("type"), type);
	SET_SLOT(obj, install("colcount"), colcount);
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
	} else if (values) {
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
	if (values) {
		SEXP x;
		R_xlen_t nx = (R_xlen_t) ((L->is_super) ? L->xsize : L->nzmax);
		if (L->xtype == CHOLMOD_COMPLEX) {
			PROTECT(x = allocVector(CPLXSXP, nx));
			Matrix_memcpy(COMPLEX(x), L->x, nx, sizeof(Rcomplex));
		} else {
			PROTECT(x = allocVector(REALSXP, nx));
			Matrix_memcpy(REAL(x), L->x, nx, sizeof(double));
		}
		SET_SLOT(obj, Matrix_xSym, x);
		UNPROTECT(1);
	}
	UNPROTECT(4);
	return obj;
}

SEXP CHS2M(cholmod_sparse *A, int values, char shape)
{
	cholmod_sparse *A_ = A;
	if (A->itype != CHOLMOD_INT)
		error(_("wrong '%s'"), "itype");
	if (values && A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		error(_("wrong '%s'"), "xtype");
	if (values && A->dtype != CHOLMOD_DOUBLE)
		error(_("wrong '%s'"), "dtype");
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		error(_("dimensions cannot exceed %s"), "2^31-1");
	if (!A->sorted)
		cholmod_sort(A, &c);
	if (!A->packed || A->stype != 0)
		A = cholmod_copy(A, A->stype, 1, &c);
	char cl[] = "..CMatrix";
	cl[0] = (!values) ? 'n' : ((A->xtype == CHOLMOD_COMPLEX) ? 'z' : 'd');
	cl[1] = shape;
	int m = (int) A->nrow, n = (int) A->ncol, nnz = ((int *) A->p)[A->ncol];
	R_xlen_t n1a = (R_xlen_t) n + 1;
	SEXP obj = PROTECT(newObject(cl)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(allocVector(INTSXP, n1a)),
		i = PROTECT(allocVector(INTSXP, nnz));
	INTEGER(dim)[0] = m;
	INTEGER(dim)[1] = n;
	Matrix_memcpy(INTEGER(p), A->p, n1a, sizeof(int));
	Matrix_memcpy(INTEGER(i), A->i, nnz, sizeof(int));
	SET_SLOT(obj, Matrix_pSym, p);
	SET_SLOT(obj, Matrix_iSym, i);
	if (values) {
		SEXP x;
		if (A->xtype == CHOLMOD_COMPLEX) {
			PROTECT(x = allocVector(CPLXSXP, nnz));
			Matrix_memcpy(COMPLEX(x), A->x, nnz, sizeof(Rcomplex));
		} else {
			PROTECT(x = allocVector(REALSXP, nnz));
			Matrix_memcpy(REAL(x), A->x, nnz, sizeof(double));
		}
		SET_SLOT(obj, Matrix_xSym, x);
		UNPROTECT(1);
	}
	if (A != A_)
		cholmod_free_sparse(&A, &c);
	UNPROTECT(4);
	return obj;
}

SEXP CHD2M(cholmod_dense *A, int trans, char shape)
{
	if (A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		error(_("wrong '%s'"), "xtype");
	if (A->dtype != CHOLMOD_DOUBLE)
		error(_("wrong '%s'"), "dtype");
	if (A->d != A->nrow) /* MJ: currently no need to support this case */
		error(_("leading dimension not equal to number of rows"));
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		error(_("dimensions cannot exceed %s"), "2^31-1");
	int m = (int) A->nrow, n = (int) A->ncol;
	if ((Matrix_int_fast64_t) m * n > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	char cl[] = "...Matrix";
	cl[0] = (A->xtype == CHOLMOD_COMPLEX) ? 'z' : 'd';
	cl[1] = shape;
	cl[2] = (shape == 'g')
		? 'e' : ((shape == 's') ? 'y' : ((shape == 'p') ? 'o' : 'r'));
	SEXP obj = PROTECT(newObject(cl)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	INTEGER(dim)[0] = (trans) ? n : m;
	INTEGER(dim)[1] = (trans) ? m : n;
	SEXP x;
	if (A->xtype == CHOLMOD_COMPLEX) {
		PROTECT(x = allocVector(CPLXSXP, (R_xlen_t) m * n));
		Rcomplex *px = COMPLEX(x), *py = (Rcomplex *) A->x;
		if (!trans)
			Matrix_memcpy(px, py, (R_xlen_t) m * n, sizeof(Rcomplex));
		else
			ztranspose2(px, py, m, n);
	} else {
		PROTECT(x = allocVector(REALSXP, (R_xlen_t) m * n));
		double *px = REAL(x), *py = (double *) A->x;
		if (!trans)
			Matrix_memcpy(px, py, (R_xlen_t) m * n, sizeof(double));
		else
			dtranspose2(px, py, m, n);
	}
	SET_SLOT(obj, Matrix_xSym, x);
	UNPROTECT(3);
	return obj;
}
