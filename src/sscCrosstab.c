#include "sscCrosstab.h"

SEXP sscCrosstab(SEXP flist, SEXP upper)
{
    int **fpt, *Ap, *Gp, *TTi, *Ti, *Tj, *dims, i,
	ncol = 0,
	nfac = length(flist),
	nfc2 = (nfac * (nfac - 1))/2, /* nfac choose 2 */
	nobs = length(VECTOR_ELT(flist, 0)),
	ntrpl, nz, pos, up = asLogical(upper);
    double *TTx, *Tx;
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

static
void col_metis_order(int j0, int j1, int i2,
		     const int Tp[], const int Ti[], int ans[])
{
    int j, nz = 0;		/* count off-diagonal pairs */
    for (j = j0; j < j1; j++) {	/* columns of interest */
	int ii, nr = 0, p2 = Tp[j + 1];
	for (ii = Tp[j]; ii < p2; ii++) {
	    int i = Ti[ii];
	    if (j1 <= i && i < i2) nr++; /* verify row index */
	}
	nz += (nr * (nr - 1))/2; /* add number of pairs of rows */
    }
    if (nz > 0) {		/* Form an ssc Matrix */
	int j, n = i2 - j1,	/* number of rows */
	    nnz = n + nz, pos;
	int *Ap = Calloc(n + 1, int),
	    *Ai = Calloc(nnz, int),
	    *Tj = Calloc(nnz, int),
	    *TTi = Calloc(nnz, int),
	    *perm = Calloc(n, int),
	    *iperm = Calloc(n, int);

	for (j = 0; j < n; j++) { /* diagonals */
	    TTi[j] = Tj[j] = j;
	}
	pos = n;
	for (j = j0; j < j1; j++) { /* create the pairs */
	    int ii, p2 = Tp[j + 1];
	    for (ii = Tp[j]; ii < p2; ii++) {
		int r1 = Ti[ii], i1;
		if (j1 <= r1 && r1 < i2) {
		    for (i1 = ii + 1; i1 < p2; i1++) {
			int r2 = Ti[i1];
			if (r2 < i2) {
			    TTi[pos] = r2 - j1;
			    Tj[pos] = r1 - j1;
			    pos++;
			}
		    }
		}
	    }
	}
	triplet_to_col(n, n, nnz, TTi, Tj, (double *) NULL,
		       Ap, Ai, (double *) NULL);
	ssc_metis_order(n, Ap, Ai, perm, iperm);
	for (j = j1; j < i2; j++) ans[j] = j1 + iperm[j - j1];
	Free(TTi); Free(Tj); Free(Ai); Free(Ap);
	Free(perm); Free(iperm);
    }
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
	i,
	n = length(pSlot) - 1,	/* number of columns */
	nf = length(GpSlot) - 1, /* number of factors */
	up;
    SEXP ans = PROTECT(allocVector(INTSXP, n));

    up = *CHAR(STRING_ELT(GET_SLOT(ctab, Matrix_uploSym), 0)) != 'L';
    if (nf > 1 && up) {			/* transpose */
	int nz = length(iSlot);
	int *ai = Calloc(nz, int),
	    *ap = Calloc(n + 1, int);
	double *ax = Calloc(nz, double);

	csc_components_transpose(n, n, nz, Ap, Ai,
				 REAL(GET_SLOT(ctab, Matrix_xSym)),
				 ap, ai, ax);
	Ap = ap;
	Ai = ai;
	Free(ax);		/* don't need values, only positions */
    }
    for (i = 0; i < n; i++) {
	INTEGER(ans)[i] = i;    /* initialize permutation to identity */
    }
    for (i = 1; i < nf; i++) {
	col_metis_order(Gp[i - 1], Gp[i], Gp[i+1], Ap, Ai, INTEGER(ans));
    }
    if (nf > 1 && up) {Free(Ap); Free(Ai);}
    UNPROTECT(1);
    return ans;
}

/**
 * Project the (2,1) component of an sscCrosstab object into the (2,2)
 * component (for illustration only)
 *
 * @param ctab pointer to a sscCrosstab object
 *
 * @return a pointer to an dsCMatrix giving the projection of the 2,1 component
 */
SEXP sscCrosstab_project(SEXP ctab)
{
    SEXP
	GpSlot = GET_SLOT(ctab, Matrix_GpSym),
	iSlot = GET_SLOT(ctab, Matrix_iSym),
	pSlot = GET_SLOT(ctab, Matrix_pSym);
    int *Ai = INTEGER(iSlot),
	*Ap = INTEGER(pSlot),
	*Gp = INTEGER(GpSlot),
	j, j0, j1, i2,
	n = length(pSlot) - 1,	/* number of columns */
	nf = length(GpSlot) - 1, /* number of factors */
	nz, up;
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix")));

    up = *CHAR(STRING_ELT(GET_SLOT(ctab, Matrix_uploSym), 0)) != 'L';
    if (nf > 1 && up) {			/* transpose */
	int nz = length(iSlot);
	int *ai = Calloc(nz, int),
	    *ap = Calloc(n + 1, int);
	double *ax = Calloc(nz, double);

	csc_components_transpose(n, n, nz, Ap, Ai,
				 REAL(GET_SLOT(ctab, Matrix_xSym)),
				 ap, ai, ax);
	Ap = ap;
	Ai = ai;
	Free(ax);		/* don't need values, only positions */
    }

    nz = 0;			/* count of off-diagonal pairs */
    j0 = 0; j1 = Gp[1]; i2 = Gp[2];
    for (j = j0; j < j1; j++) {	/* columns of interest */
	int ii, nr = 0, p2 = Ap[j + 1];
	for (ii = Ap[j]; ii < p2; ii++) {
	    int i = Ai[ii];
	    if (j1 <= i && i < i2) nr++; /* verify row index */
	}
	nz += (nr * (nr - 1))/2; /* add number of pairs of rows */
    }
    if (nz > 0) {		/* Form an ssc Matrix */
	int j, n = i2 - j1,	/* number of rows */
	    nnz = n + nz, pos;
	int *AAp,
	    *AAi = Calloc(nnz, int),
	    *Tj = Calloc(nnz, int),
	    *TTi = Calloc(nnz, int);
	double *Ax;

	SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, n + 1));
	AAp = INTEGER(GET_SLOT(ans, Matrix_pSym));
	for (j = 0; j < n; j++) { /* diagonals */
	    TTi[j] = Tj[j] = j;
	}
	pos = n;
	for (j = j0; j < j1; j++) { /* create the pairs */
	    int ii, p2 = Ap[j + 1];
	    for (ii = Ap[j]; ii < p2; ii++) {
		int r1 = Ai[ii], i1;
		if (j1 <= r1 && r1 < i2) {
		    for (i1 = ii + 1; i1 < p2; i1++) {
			int r2 = Ai[i1];
			if (r2 < i2) {
			    TTi[pos] = r2 - j1;
			    Tj[pos] = r1 - j1;
			    pos++;
			}
		    }
		}
	    }
	}
	triplet_to_col(n, n, nnz, TTi, Tj, (double *) NULL,
		       AAp, AAi, (double *) NULL);
	nz = AAp[n];
	SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nz));
	Memcpy(INTEGER(GET_SLOT(ans, Matrix_iSym)), AAi, nz);
	SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nz));
	Ax = REAL(GET_SLOT(ans, Matrix_xSym));
	for (j = 0; j < nz; j++) Ax[j] = 1.;
	SET_SLOT(ans, Matrix_uploSym, mkString("L"));
	SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
	AAp = INTEGER(GET_SLOT(ans, Matrix_DimSym));
	AAp[0] = AAp[1] = n;
	Free(TTi); Free(Tj); Free(AAi);
    }
    if (nf > 1 && up) {Free(Ap); Free(Ai);}
    UNPROTECT(1);
    return ans;
}

/**
 * Project the first group of columns in an sscCrosstab object onto the
 * remaining columns.
 *
 * @param ctab pointer to a sscCrosstab object
 *
 * @return a pointer to an dsCMatrix with the projection
 */
SEXP sscCrosstab_project2(SEXP ctab)
{
    SEXP
	GpSlot = GET_SLOT(ctab, Matrix_GpSym),
	iSlot = GET_SLOT(ctab, Matrix_iSym),
	pSlot = GET_SLOT(ctab, Matrix_pSym);
    int *Ai = INTEGER(iSlot),
	*Ap = INTEGER(pSlot),
	*Gp = INTEGER(GpSlot),
	i, i2, j, j1, k, k2,
	nf = length(GpSlot) - 1, /* number of factors */
	nz, pos, up, *AAp, *Ti, *Cp, *ind;
    double *Ax = REAL(GET_SLOT(ctab, Matrix_xSym));
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix")));


    if (nf < 2) error("sscCrosstab_project2 requires more than one group");
    up = *CHAR(STRING_ELT(GET_SLOT(ctab, Matrix_uploSym), 0)) != 'L';
    if (up) {			/* tranpose */
	int n = length(pSlot) - 1;	/* number of columns */
	int nnz = length(iSlot);
	int *ai = Calloc(nnz, int);
	int *ap = Calloc(n + 1, int);
	double *ax = Calloc(nnz, double);

	csc_components_transpose(n, n, nnz, Ap, Ai, Ax, ap, ai, ax);
	Ap = ap; Ai = ai; Ax = ax;
    }

    j1 = Gp[1]; i2 = Gp[nf];	/* group boundaries */
    Cp = Calloc(j1, int);	/* first pos in col with row > i */
    for (j = 0; j < j1; j++) Cp[j] = Ap[j] + 1;	/* Ap[j] is the diagonal */

    nz = Ap[i2] - Ap[j1];	/* upper bound on nonzeros in result */
    for (j = 0; j < j1; j++) {
	int nr = Ap[j + 1] - Ap[j] - 1;
	nz += (nr * (nr - 1))/2; /* add number of pairs of rows below diag */
    }

    Ti = Calloc(nz, int);	/* temporary row indices */
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, i2 - j1 + 1));
    AAp = INTEGER(GET_SLOT(ans, Matrix_pSym)); /* column pointers */

    AAp[0] = 0; pos = 0;
    ind = Calloc(i2 - j1, int);	/* indicator of rows in same column */
    for (i = j1; i < i2; i++) {
	for (k = j1; k < i2; k++) ind[k - j1] = 0;
	for (j = 0; j < j1; j++) {
	    if (Ai[Cp[j]] == i) { /* go down the column */
		k2 = Ap[j+1];
		for(k = Cp[j] + 1; k < k2; k++) ind[Ai[k] - j1] = 1;
		Cp[j]++;
	    }
	}

	Ti[pos++] = i - j1;	/* diagonal element */
	for (k = i+1; k < i2; k++) { /* projected pairs */
	    int ii = k - j1;
	    if (ind[ii]) Ti[pos++] = ii;
	}
	k2 = Ap[i+1];
	for (k = Ap[i] + 1; k < k2; k++) { /* previous off-diagonals */
	    Ti[pos++] = Ai[k] - j1;
	}
	AAp[i - j1 + 1] = pos;
    }
    nz = AAp[i2 - j1];
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nz));
    Memcpy(INTEGER(GET_SLOT(ans, Matrix_iSym)), Ti, nz);
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nz));
    Ax = REAL(GET_SLOT(ans, Matrix_xSym));
    for (j = 0; j < nz; j++) Ax[j] = 1.;
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
    AAp = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    AAp[0] = AAp[1] = i2 - j1;
    Free(Ti); Free(Cp); Free(ind);
    if (up) {Free(Ap); Free(Ai); free(Ax);}
    UNPROTECT(1);
    return ans;
}
