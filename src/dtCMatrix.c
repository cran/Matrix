				/* Sparse triangular numeric matrices */
#include "dtCMatrix.h"
#include "cs_utils.h"

/* This should be use for *BOTH* triangular and symmetric Csparse: */
SEXP tCMatrix_validate(SEXP x)
{
    SEXP val = xCMatrix_validate(x);/* checks x slot */
    if(isString(val))
	return(val);
    else {
	SEXP
	    islot = GET_SLOT(x, Matrix_iSym),
	    pslot = GET_SLOT(x, Matrix_pSym);
	int uploT = (*uplo_P(x) == 'U'),
	    k, nnz = length(islot),
	    *xi = INTEGER(islot),
	    *xj = INTEGER(PROTECT(allocVector(INTSXP, nnz)));

#define RETURN(_CH_)   UNPROTECT(1); return (_CH_);

	expand_cmprPt(length(pslot) - 1, INTEGER(pslot), xj);

	/* Maybe FIXME: ">" should be ">="	for diag = 'U' (uplo = 'U') */
	if(uploT) {
	    for (k = 0; k < nnz; k++)
		if(xi[k] > xj[k]) {
		    RETURN(mkString(_("uplo='U' must not have sparse entries in lower diagonal")));
		}
	}
	else {
	    for (k = 0; k < nnz; k++)
		if(xi[k] < xj[k]) {
		    RETURN(mkString(_("uplo='L' must not have sparse entries in upper diagonal")));
		}
	}

	RETURN(ScalarLogical(1));
    }
}
#undef RETURN

/**
 * Derive the column pointer vector for the inverse of L from the parent array
 *
 * @param n length of parent array
 * @param countDiag 0 for a unit triangular matrix with implicit diagonal, otherwise 1
 * @param pr parent vector describing the elimination tree
 * @param ap array of length n+1 to be filled with the column pointers
 *
 * @return the number of non-zero entries (ap[n])
 */
int parent_inv_ap(int n, int countDiag, const int pr[], int ap[])
{
    int *sz = Alloca(n, int), j;
    R_CheckStack();

    for (j = n - 1; j >= 0; j--) {
	int parent = pr[j];
	sz[j] = (parent < 0) ?  countDiag : (1 + sz[parent]);
    }
    ap[0] = 0;
    for (j = 0; j < n; j++)
	ap[j+1] = ap[j] + sz[j];
    return ap[n];
}

/**
 * Derive the row index array for the inverse of L from the parent array
 *
 * @param n length of parent array
 * @param countDiag 0 for a unit triangular matrix with implicit diagonal, otherwise 1
 * @param pr parent vector describing the elimination tree
 * @param ai row index vector of length ap[n]
 */
void parent_inv_ai(int n, int countDiag, const int pr[], int ai[])
{
    int j, k, pos = 0;
    for (j = 0; j < n; j++) {
	if (countDiag) ai[pos++] = j;
	for (k = pr[j]; k >= 0; k = pr[k]) ai[pos++] = k;
    }
}

SEXP Parent_inverse(SEXP par, SEXP unitdiag)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *ap, *ai, *dims, *pr = INTEGER(par),
	countDiag = 1 - asLogical(unitdiag),
	j, n = length(par), nnz;
    double *ax;

    if (!isInteger(par)) error(_("par argument must be an integer vector"));
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, n + 1));
    ap = INTEGER(GET_SLOT(ans, Matrix_pSym));
    nnz = parent_inv_ap(n, countDiag, pr, ap);
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    ai = INTEGER(GET_SLOT(ans, Matrix_iSym));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    ax = REAL(GET_SLOT(ans, Matrix_xSym));
    for (j = 0; j < nnz; j++) ax[j] = 1.;
    SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
    dims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    dims[0] = dims[1] = n;
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_diagSym, (countDiag ? mkString("N") : mkString("U")));
    parent_inv_ai(n, countDiag, pr, ai);
    UNPROTECT(1);
    return ans;
}

SEXP dtCMatrix_solve(SEXP a)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    CSP A = AS_CSP(a);
    int *bp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, (A->n) + 1)),
	bnz = 10 * A->n,	/* initial estimate of nnz in b */
	lo = uplo_P(a)[0] == 'L', top;
    /* These arrays must use Calloc because of possible Realloc */
    int *ti = Calloc(bnz, int), p, j, nz, pos = 0;
    double *tx = Calloc(bnz, double);
    cs *u = cs_spalloc(A->n, 1,1,1,0);	/* Sparse unit vector */
    double  *wrk = Alloca(A->n, double);
    int *xi = Alloca(2*A->n, int);	/* for cs_reach */
    R_CheckStack();

    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(a, Matrix_DimSym)));
    SET_DimNames(ans, a);
    SET_SLOT(ans, Matrix_uploSym, duplicate(GET_SLOT(a, Matrix_uploSym)));
    SET_SLOT(ans, Matrix_diagSym, duplicate(GET_SLOT(a, Matrix_diagSym)));
    /* initialize the "constant part" of the sparse unit vector */
    u->x[0] = 1.;
    u->p[0] = 0; u->p[1] = 1;
    bp[0] = 0;
    for (j = 0; j < A->n; j++) {
	u->i[0] = j;			/* u := j'th unit vector */
	/* (wrk[top:n],xi[top:n]) :=  A^{-1} u  : */
	top = cs_spsolve (A, u, 0, xi, wrk, 0, lo);
	nz = A->n - top;
	bp[j + 1] = nz + bp[j];
	if (bp[j + 1] > bnz) {
	    while (bp[j + 1] > bnz) bnz *= 2;
	    ti = Realloc(ti, bnz, int);
	    tx = Realloc(tx, bnz, double);
	}
	if (lo)
	    for(p = top; p < A->n; p++, pos++) {
		ti[pos] = xi[p];
		tx[pos] = wrk[xi[p]];
	    }
	else /* upper triagonal */
	    for(p = A->n - 1; p >= top; p--, pos++) {
		ti[pos] = xi[p];
		tx[pos] = wrk[xi[p]];
	    }
    }
    nz = bp[A->n];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP,  nz)), ti, nz);
    Memcpy(   REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)), tx, nz);

    Free(ti); Free(tx); cs_spfree(u);
    UNPROTECT(1);
    return ans;
}

SEXP dtCMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed)
{
    int cl = asLogical(classed);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    CSP A = AS_CSP(a);
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol));
    int j, n = bdims[0], nrhs = bdims[1], lo = (*uplo_P(a) == 'L');
    double *bx;
    R_CheckStack();

    if (*adims != n || nrhs < 1 || *adims < 1 || *adims != adims[1])
	error(_("Dimensions of system to be solved are inconsistent"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)), bdims, 2);
    /* FIXME: copy dimnames or Dimnames as well */
    bx = Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * nrhs)),
		REAL(cl ? GET_SLOT(b, Matrix_xSym):b), n * nrhs);
    for (j = 0; j < nrhs; j++)
	lo ? cs_lsolve(A, bx + n * j) : cs_usolve(A, bx + n * j);
    UNPROTECT(1);
    return ans;
}

SEXP dtCMatrix_upper_solve(SEXP a)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int lo = uplo_P(a)[0] == 'L', unit = diag_P(a)[0] == 'U',
	n = INTEGER(GET_SLOT(a, Matrix_DimSym))[0],
	*ai = INTEGER(GET_SLOT(a,Matrix_iSym)),
	*ap = INTEGER(GET_SLOT(a, Matrix_pSym)),
	*bp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1));
    int bnz = 10 * ap[n];	  /* initial estimate of nnz in b */
    int *ti = Alloca(bnz, int), j, nz;
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)), *tx = Alloca(bnz, double),
	*tmp = Alloca(n, double);
    R_CheckStack();

    if (lo || (!unit))
	error(_("Code written for unit upper triangular unit matrices"));
    bp[0] = 0;
    for (j = 0; j < n; j++) {
	int i, i1 = ap[j + 1];
	AZERO(tmp, n);
	for (i = ap[j]; i < i1; i++) {
	    int ii = ai[i], k;
	    if (unit) tmp[ii] -= ax[i];
	    for (k = bp[ii]; k < bp[ii + 1]; k++) tmp[ti[k]] -= ax[i] * tx[k];
	}
	for (i = 0, nz = 0; i < n; i++) if (tmp[i]) nz++;
	bp[j + 1] = bp[j] + nz;
	if (bp[j + 1] > bnz) {
	    while (bp[j + 1] > bnz) bnz *= 2;
	    ti = Realloc(ti, bnz, int);
	    tx = Realloc(tx, bnz, double);
	}
	i1 = bp[j];
	for (i = 0; i < n; i++) if (tmp[i]) {ti[i1] = i; tx[i1] = tmp[i]; i1++;}
    }
    nz = bp[n];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)), ti, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)), tx, nz);
    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(a, Matrix_DimSym)));
    SET_DimNames(ans, a);
    SET_SLOT(ans, Matrix_uploSym, duplicate(GET_SLOT(a, Matrix_uploSym)));
    SET_SLOT(ans, Matrix_diagSym, duplicate(GET_SLOT(a, Matrix_diagSym)));
    UNPROTECT(1);
    return ans;
}
