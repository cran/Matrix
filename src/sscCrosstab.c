#include "sscCrosstab.h"

SEXP sscCrosstab(SEXP flist, SEXP upper)
{
    int
	**fpt,
	*Ap,
	*Gp,
	*TTi,
	*Ti,
	*Tj,
	*dims,
	i,
	ncol = 0,
	nfac = length(flist),
	nfc2 = (nfac * (nfac - 1))/2, /* nfac choose 2 */
	nobs = length(VECTOR_ELT(flist, 0)),
	ntrpl,
	nz,
	pos,
	up = asLogical(upper);
    double
	*TTx,
	*Tx;
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("sscCrosstab")));

    if (!isNewList(flist) || nfac < 1)
	error("flist must be a non-empty list");
    SET_SLOT(val, Matrix_GpSym, allocVector(INTSXP, nfac + 1));
    Gp = INTEGER(GET_SLOT(val, Matrix_GpSym));
    fpt = (int **) R_alloc(nfac, sizeof(int *));
    for (i = 0; i < nfac; i++) {
	SEXP el = VECTOR_ELT(flist, i);
	if (!inherits(el, "factor"))
	    error("flist must be a non-empty list of factors");
	if (length(el) != nobs)
	    error("All elements of flist must have the same length");
	Gp[i] = ncol;
	ncol += length(getAttrib(el, R_LevelsSymbol));
	fpt[i] = INTEGER(el);
    }
    Gp[nfac] = ncol;
    SET_SLOT(val, Matrix_uploSym, ScalarString(mkChar(up ? "U" : "L")));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    dims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    dims[0] = dims[1] = ncol;
    ntrpl = nfc2 * nobs + ncol;
    Ti = Calloc(ntrpl, int); Tj = Calloc(ntrpl, int); TTi = Calloc(ntrpl, int);
    Tx = Calloc(ntrpl, double); TTx = Calloc(ntrpl, double);
				/* Generate the triplet form of the result */
    for (i = 0; i < ncol; i++) {
	Ti[i] = Tj[i] = i;	/* The diagonals - these will store counts */
	Tx[i] = 0.0;
    }
    pos = ncol;
    for (i = 0; i < nobs ; i++) {
	int j, jcol, k;
	for (j = 0; j < nfac; j++) {
	    jcol = Gp[j] + fpt[j][i] - 1;
	    Tx[jcol] += 1.;	/* increment diagonal count */
	    for (k = j + 1; k < nfac; k++) { /* off-diagonals */
		int irow = Gp[k] + fpt[k][i] - 1;
		if (up) {
		    Ti[pos] = jcol; Tj[pos] = irow;
		} else {
		    Tj[pos] = jcol; Ti[pos] = irow;
		}
		Tx[pos] = 1.;
		pos++;
	    }
	}
    }
    SET_SLOT(val, Matrix_pSym, allocVector(INTSXP, ncol + 1));
    Ap = INTEGER(GET_SLOT(val, Matrix_pSym));
    triplet_to_col(ncol, ncol, ntrpl, Ti, Tj, Tx, Ap, TTi, TTx);
    nz = Ap[ncol];		/* non-zeros in Z'Z crosstab */
    SET_SLOT(val, Matrix_iSym, allocVector(INTSXP, nz));
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, nz));
    Memcpy(INTEGER(GET_SLOT(val, Matrix_iSym)), TTi, nz);
    Memcpy(REAL(GET_SLOT(val, Matrix_xSym)), TTx, nz);
    Free(Ti); Free(Tj); Free(Tx);
    Free(TTi); Free(TTx);
    
    UNPROTECT(1);
    return val;
}

extern void ssclme_fill_LIp(int n, const int Parent[], int LIp[]);

SEXP sscCrosstab_L_LI_sizes(SEXP ctab, SEXP permexp)
{
    SEXP ans = PROTECT(allocVector(INTSXP, 4));
    int *Ai = INTEGER(GET_SLOT(ctab, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(ctab, Matrix_pSym)),
	*aa = INTEGER(ans),
	*perm = INTEGER(permexp),
	n = INTEGER(GET_SLOT(ctab, Matrix_DimSym))[1],
	*Lp = Calloc(n + 1, int),
	*Parent = Calloc(n, int),
	*Lnz = Calloc(n, int),
	*Flag = Calloc(n, int);

    ldl_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag,
		 (int *) NULL, (int *) NULL); /* P & Pinv */
    aa[0] = Lp[n];
    ssclme_fill_LIp(n, Parent, Lp);
    aa[1] = Lp[n];
    ssc_symbolic_permute(n, 1, perm, Ap, Ai);
    ldl_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag,
		 (int *) NULL, (int *) NULL); /* P & Pinv */
    aa[2] = Lp[n];
    ssclme_fill_LIp(n, Parent, Lp);
    aa[3] = Lp[n];
    Free(Flag); Free(Lnz); Free(Parent); Free(Lp); 
    UNPROTECT(1);
    return ans;
}
    
SEXP sscCrosstab_groupedPerm(SEXP ctab)
{
    SEXP
	GpSlot = GET_SLOT(ctab, Matrix_GpSym),
	iSlot = GET_SLOT(ctab, Matrix_iSym),
	pSlot = GET_SLOT(ctab, Matrix_pSym);
    int *Ai = INTEGER(iSlot),
	*Ap = INTEGER(pSlot),
	*Gp = INTEGER(GpSlot),
	nl1 = Gp[1],		/* number of levels of first factor */
	*icounts,		/* counts for rows */
	*jcounts = Calloc(nl1, int), /* counts for columns */
	nl2,			/* delay initializing until nf is known */
	ii, j, jj,
	n = length(pSlot) - 1,	/* number of columns */
	nf = length(GpSlot) - 1, /* number of factors */
	nz = length(iSlot),	/* number of non-zeros */
	nod = nz - n,		/* number of off-diagonals */
	*np = Calloc(nl1 + 1, int), /* new array pointers */
	*ni = Calloc(nod, int),	/* new array row indices */
	p1, p3;
    SEXP ans = PROTECT(allocVector(INTSXP, n));
    int *perm = INTEGER(ans);

    if (toupper(*CHAR(STRING_ELT(GET_SLOT(ctab, Matrix_uploSym), 0))) != 'L')
	error("Lower triangle required in sscCrosstab object");

    for (j = 0; j < n; j++) perm[j] = j; /* initialize permutation to identity */

    if (nf > 1) {
	nl2 = Gp[2];
	icounts = Calloc(nl2 - nl1, int);
	np[0] = 0;
	for (j = 0; j < nl1; j++) jcounts[j] = 0;
	p1 = 0;			/* copy off-diagonals from first nl1 cols */
	for (j = 0; j < nl1; j++) { /* and nl1 <= row < nl2 */
	    int p3 = Ap[j + 1];
	    for (jj = Ap[j]; jj < p3; jj++) {
		int i = Ai[jj];
		if (nl1 <= i && i < nl2) {
		    ni[p1++] = i - nl1;	
		    jcounts[j]++;
		}
	    }
	    np[j+1] = p1;
	}

	for (ii = 1; ii <= nl2; ii++) {
	    int maxrc, rr;
	    int i, minjc = nl2+1;	/* find minimum positive jcount */
	    for (j = 0; j < nl1; j++) {
		if (jcounts[i] > 0 && jcounts[j] < minjc) minjc = jcounts[j];
	    }
				/* accumulate the row counts on cols where jcount=minjc */
	    for (i = 0; i < nl2; i++) icounts[i] = 0;
	    for (j = 0; j < nl1; j++) {
		if (jcounts[j] == minjc) {
		    int p2 = np[j+1];
		    for (i = np[j]; i < p2; i++) icounts[ni[i]]++;
		}
	    }
				/* find the last row whose count in max(icount) */
	    maxrc = -1;
	    for (i = 0; i < nl2; i++) {
		int ic = icounts[i];
		if (ic >= maxrc) {
		    maxrc = ic;
		    rr = i;
		}
	    }
	    perm[nl1 + nl2 - ii] = rr;
				/* update icounts, np and ni */
	    p1 = p3 = 0; 
	    for (j = 0; j < nl1; j++) {
		int p2 = np[j+1];
		for (i = np[j]; i < p2; i++) {
		    if (ni[i] != rr) {
			if (i != p1) ni[p1] = ni[i];
			p1++;
		    }
		}
		np[j] = p3;
		p3 = p1;	/* save the count */
	    }
	}
	Free(icounts);
    }
    Free(np); Free(ni); Free(jcounts);
    UNPROTECT(1);
    return ans;
}
