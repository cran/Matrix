				/* Sparse triangular numeric matrices */
#include "dtCMatrix.h"
#include "cs_utils.h"

SEXP tsc_validate(SEXP x)
{
    return triangularMatrix_validate(x);
    /* see ./dsCMatrix.c or ./dtpMatrix.c  on how to do more testing here */
}

#if 0    
SEXP tsc_transpose(SEXP x)
{
    cholmod_sparse *cx = as_cholmod_sparse(x);
    
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix"))),
	islot = GET_SLOT(x, Matrix_iSym);
    int nnz = length(islot),
	*adims, *xdims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int up = uplo_P(x)[0] == 'U';

    adims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    adims[0] = xdims[1]; adims[1] = xdims[0];

    if(*diag_P(x) == 'U')
	SET_SLOT(ans, Matrix_diagSym, duplicate(GET_SLOT(x, Matrix_diagSym)));
    SET_SLOT(ans, Matrix_uploSym, mkString(up ? "L" : "U"));

    csc_compTr(xdims[0], xdims[1], nnz,
	       INTEGER(GET_SLOT(x, Matrix_pSym)), INTEGER(islot),
	       REAL(GET_SLOT(x, Matrix_xSym)),
	       INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, xdims[0] + 1)),
	       INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nnz)),
	       REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nnz)));
    UNPROTECT(1);
    return ans;
}
#endif

SEXP tsc_to_dgTMatrix(SEXP x)
{
    SEXP ans;
    if (*diag_P(x) != 'U')
	ans = compressed_to_dgTMatrix(x, ScalarLogical(1));
    else {			/* unit triangular matrix */
	SEXP islot = GET_SLOT(x, Matrix_iSym),
	    pslot = GET_SLOT(x, Matrix_pSym);
	int *ai, *aj, j,
	    n = length(pslot) - 1,
	    nod = length(islot),
	    nout = n + nod,
	    *p = INTEGER(pslot);
	double *ax;

	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgTMatrix")));
	SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
	SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nout));
	ai = INTEGER(GET_SLOT(ans, Matrix_iSym));
	Memcpy(ai, INTEGER(islot), nod);
	SET_SLOT(ans, Matrix_jSym, allocVector(INTSXP, nout));
	aj = INTEGER(GET_SLOT(ans, Matrix_jSym));
	SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nout));
	ax = REAL(GET_SLOT(ans, Matrix_xSym));
	Memcpy(ax, REAL(GET_SLOT(x, Matrix_xSym)), nod);
	for (j = 0; j < n; j++) {
	    int jj, npj = nod + j,  p2 = p[j+1];
	    aj[npj] = ai[npj] = j;
	    ax[npj] = 1.;
	    for (jj = p[j]; jj < p2; jj++) aj[jj] = j;
	}
	UNPROTECT(1);
    }
    return ans;
}

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
    int *sz = Calloc(n, int), j;

    for (j = n - 1; j >= 0; j--) {
	int parent = pr[j];
	sz[j] = (parent < 0) ?  countDiag : (1 + sz[parent]);
    }
    ap[0] = 0;
    for (j = 0; j < n; j++)
	ap[j+1] = ap[j] + sz[j];
    Free(sz);
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

#if 0
SEXP dtCMatrix_solve(SEXP a)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int lo = uplo_P(a)[0] == 'L', unit = diag_P(a)[0] == 'U',
	n = INTEGER(GET_SLOT(a, Matrix_DimSym))[0], nnz,
	*ai = INTEGER(GET_SLOT(a,Matrix_iSym)),
	*ap = INTEGER(GET_SLOT(a, Matrix_pSym)), *bi,
	*bp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1));
    int bnz = 10 * ap[n];	  /* initial estimate of nnz in b */
    int *ri = Calloc(bnz, int), *ind = Calloc(n, int), j;
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)), *bx;

    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(a, Matrix_DimSym)));
    SET_SLOT(ans, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(a, Matrix_DimNamesSym)));
    SET_SLOT(ans, Matrix_uploSym, duplicate(GET_SLOT(a, Matrix_uploSym)));
    SET_SLOT(ans, Matrix_diagSym, duplicate(GET_SLOT(a, Matrix_diagSym)));
    
    if (!(lo && unit))
	error("code for non-unit or upper triangular not yet written");
    /* Initially bp will contain increasing negative values ending at zero. */
    /* Later we add the negative of bp[0] to all values. */
    bp[n] = 0;
    for (j = n - 1; j >= 0; j--) { /* columns in reverse order */
	int i, i1 = ap[j], i2 = ap[j + 1], k, nr;
	if (i1 < i2) AZERO(ind, n);
	for (i = i1; i < i2; i++) {
	    ind[ai[i]] = 1;
	    for (k = -bp[ai[i] + 1]; k < -bp[ai[i]]; k++) ind[ri[k]] = 1;
	}
	for (k = 0, nr = 0; k < n; k++) if (ind[k]) nr++;
	if ((nr - bp[j + 1]) > bnz) {
	    while (nr > (bnz + bp[j + 1])) bnz *= 2;
	    ri = Realloc(ri, bnz, int);
	}
	bp[j] = bp[j + 1] - nr;
	for (k = 0, i = -bp[j + 1]; k < n; k++) if (ind[k]) ri[i++] = k;
    }
    bnz = -bp[0];
    bi = INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, bnz));
    bx = REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, bnz));
    for (j = 0; j < n; j++) {
	int bpnew = bp[j] + bnz;
	Memcpy(bi + bpnew, ri - bp[j], bp[j + 1] - bp[j]);
	bp[j] = bpnew;
    }
    /* insert code for calculating the actual values here */
    for (j = 0; j < bnz; j++) bx[j] = 1;

    Free(ind); Free(ri);
    UNPROTECT(1);
    return ans;
}
#else
SEXP dtCMatrix_solve(SEXP a)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    cs *A = Matrix_as_cs(a);
    int *bp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, (A->n) + 1)),
	lo = uplo_P(a)[0] == 'L',
	bnz = 10 * A->n;	/* initial estimate of nnz in b */
    int *ti = Calloc(bnz, int), i, j, nz, pos = 0;
    double *tx = Calloc(bnz, double), *wrk = Calloc(A->n, double);

    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(a, Matrix_DimSym)));
    SET_SLOT(ans, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(a, Matrix_DimNamesSym)));
    SET_SLOT(ans, Matrix_uploSym, duplicate(GET_SLOT(a, Matrix_uploSym)));
    SET_SLOT(ans, Matrix_diagSym, duplicate(GET_SLOT(a, Matrix_diagSym)));
    bp[0] = 0;
    for (j = 0; j < A->n; j++) {
	AZERO(wrk, A->n);
	wrk[j] = 1;
	lo ? cs_lsolve(A, wrk) : cs_usolve(A, wrk);
	for (i = 0, nz = 0; i < A->n; i++) if (wrk[i]) nz++;
	bp[j + 1] = nz + bp[j];
	if (bp[j + 1] > bnz) {
	    while (bp[j + 1] > bnz) bnz *= 2;
	    ti = Realloc(ti, bnz, int);
	    tx = Realloc(tx, bnz, double);
	}
	for (i = 0; i < A->n; i++)
	    if (wrk[i]) {ti[pos] = i; tx[pos] = wrk[i]; pos++;}
    }
    nz = bp[A->n];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)), ti, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)), tx, nz);

    Free(A); Free(ti); Free(tx);
    UNPROTECT(1);
    return ans;
}
#endif

SEXP dtCMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed)
{
    int cl = asLogical(classed);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    cs *A = Matrix_as_cs(a);    
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol));
    int j, n = bdims[0], nrhs = bdims[1], lo = (*uplo_P(a) == 'L');
    double *bx;

    if (*adims != n || nrhs < 1 || *adims < 1 || *adims != adims[1])
	error(_("Dimensions of system to be solved are inconsistent"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)), bdims, 2);
    /* copy dimnames or Dimnames as well */
    bx = Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * nrhs)),
		REAL(cl ? GET_SLOT(b, Matrix_xSym):b), n * nrhs);
    for (j = 0; j < nrhs; j++)
	lo ? cs_lsolve(A, bx + n * j) : cs_usolve(A, bx + n * j);
    Free(A);
    UNPROTECT(1);
    return ans;
}

SEXP dtCMatrix_upper_solve(SEXP a)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int lo = uplo_P(a)[0] == 'L', unit = diag_P(a)[0] == 'U',
	n = INTEGER(GET_SLOT(a, Matrix_DimSym))[0], nnz,
	*ai = INTEGER(GET_SLOT(a,Matrix_iSym)),
	*ap = INTEGER(GET_SLOT(a, Matrix_pSym)), *bi,
	*bp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1));
    int bnz = 10 * ap[n];	  /* initial estimate of nnz in b */
    int *ti = Calloc(bnz, int), j, nz;
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)), *tx = Calloc(bnz, double),
	*tmp = Calloc(n, double);
    
    if (lo || (!unit)) error(_("Code written for unit upper triangular unit matrices"));
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
    Free(tmp); Free(tx); Free(ti);
    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(a, Matrix_DimSym)));
    SET_SLOT(ans, Matrix_DimNamesSym, duplicate(GET_SLOT(a, Matrix_DimNamesSym)));
    SET_SLOT(ans, Matrix_uploSym, duplicate(GET_SLOT(a, Matrix_uploSym)));
    SET_SLOT(ans, Matrix_diagSym, duplicate(GET_SLOT(a, Matrix_diagSym)));
    UNPROTECT(1);
    return ans;
}
