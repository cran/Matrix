#include "bCrosstab.h"
/* TODO:
 * - Only do a fill-reducing permutation on the first non-nested factor
 * - Alternatively: change the algorithm for the fill-reducing
 *   permutation to a greedy or a picky algorithm.
 * - Rewrite this whole section using the lgTMatrix and lgCMatrix classes.
 *   A symbolic manipulation of a sparse matrix is equivalent to
 *   manipulating the sparse matrix with logical entries.
 * - Create coercion methods for dgCMatrix -> lgCMatrix and vice
 *   versa.  Same for dgTMatrix -> lgTMatrix.
 */

/** 
 * Replace the structure of C by the structure of CL^{-1} where L is the
 * unit lower triangular sparse matrix from an LDL' Cholesky decomposition
 * 
 * @param anc number of columns in A
 * @param Parent parent array for A
 * @param C a dgBCMatrix object to be updated
 */
static void
symbolic_right_unit_sm(int anc, const int Parent[], SEXP C)
{
    SEXP cip = GET_SLOT(C, Matrix_iSym),
	cpp = GET_SLOT(C, Matrix_pSym);
    int *Flag,
	*ci = INTEGER(cip),
	*cp = INTEGER(cpp),
	*ncp,
	cnr, cnz = length(cip),
	i, j;

    if ((length(cpp) - 1) != anc) /* A is square so can compare no of cols */
	error(_("No. of rows in A (%d) does not match no. of cols in C (%d)"),
	      anc, length(cpp) - 1); 
    i = 1;			/* check for A being the identity */
    for (j = 0; j < anc; j++) {
	if (Parent[j] >= 0) {
	    i = 0;
	    break;
	}
    }
    if (i) return;		/* A is the identity */
    cnr = 0;			/* number of rows in C (= max(ci + 1)) */
    for (i = 0; i < cnz; i++) {
	int ri = ci[i] + 1;
	if (cnr < ri) cnr = ri;
    }
    Flag = Calloc(cnr, int);
    ncp = Calloc(anc + 1, int);	/* new column pointers */

    ncp[0] = 0;
    for (j = 0; j < anc; j++) {
	int cj2 = cp[j + 1], kk, kc;
	for (i = 0; i < cnr; i++) Flag[i] = 0;
	ncp[j+1] = ncp[j] + cj2 - cp[j];
				/* positions of current column j of C */
	for (kc = cp[j]; kc < cj2; kc++) Flag[ci[kc]] = 1;
				/* other positions in column j of product */
	for (kk = Parent[j]; kk >= 0; kk = Parent[kk]) {
	    int kk2 = cp[kk + 1];
	    for (kc = cp[kk]; kc < kk2; kc++) {
		if (!Flag[ci[kc]]) {
		    ncp[j+1]++;
		    Flag[ci[kc]] = 1;
		}
	    }
	}
    }
    if (ncp[anc] > cp[anc]) {
	int *dims, *nci, nnz = ncp[anc], pos = 0;
	double *ncx;

	SET_SLOT(C, Matrix_iSym, allocVector(INTSXP,  nnz));
	nci = INTEGER(GET_SLOT(C, Matrix_iSym));
	dims = INTEGER(getAttrib(GET_SLOT(C, Matrix_xSym), R_DimSymbol));
	SET_SLOT(C, Matrix_xSym, alloc3Darray(REALSXP, dims[0], dims[1], nnz));
	ncx = REAL(GET_SLOT(C, Matrix_xSym));
	for (i = 0; i < nnz; i++) ncx[i] = 1.;
				/* As Diana Krall said, "Just do it again." */
	for (j = 0; j < anc; j++) {
	    int cj2 = cp[j + 1], kc, kk;
	    for (i = 0; i < cnr; i++) Flag[i] = 0;
	    for (kc = cp[j]; kc < cj2; kc++) Flag[ci[kc]] = 1;
	    for (kk = Parent[j]; kk >= 0; kk = Parent[kk]) {
		int kk2 = cp[kk + 1];
		for (kc = cp[kk]; kc < kk2; kc++) Flag[ci[kc]] = 1;
	    }
	    for (i = 0; i < cnr; i++) if (Flag[i]) nci[pos++] = i;
	}
	Memcpy(cp, ncp, anc + 1);
    }	
    Free(Flag); Free(ncp);
}

/** 
 * Replace the structure of C by the structure of CA^{-T}
 * 
 * @param anc number of column blocks in A
 * @param Parent parent array for column blocks of A
 * @param C a dgBCMatrix object to be updated
 */
static void
symbolic_right_unit_mm_trans(int anc, const int Parent[], SEXP C)
{
    SEXP cip = GET_SLOT(C, Matrix_iSym),
	cpp = GET_SLOT(C, Matrix_pSym);
    int *ci = INTEGER(cip),
	*cp = INTEGER(cpp),
	cnz = length(cip),
	i, j, nextra = 0;

    if ((length(cpp) - 1) != anc)
	error(_("No. of cols in A (%d) does not match no. of cols in C (%d)"),
	      anc, length(cpp) - 1);
    
    i = 1;			/* check for A being the identity */
    for (j = 0; j < anc; j++) {
	if (Parent[j] >= 0) {
	    i = 0;
	    break;
	}
    }
    if (i) return;		/* A is the identity */

    for (j = 0; j < anc; j++) { /* bound the number of extra triplets */
	int cj2 = cp[j + 1], ka, kc;
	for (ka = Parent[j]; ka >= 0; ka = Parent[ka]) {
	    for (kc = cp[j]; kc < cj2; kc++) {
		if (check_csc_index(cp, ci, ci[kc], ka, -1) < 0) nextra++;
	    }
	}
    }
    if (nextra) {
	int cnr, ntot = cnz + nextra, pos;
	int *dims = INTEGER(getAttrib(GET_SLOT(C, Matrix_xSym), R_DimSymbol)),
	    *Ti = Memcpy((int *) Calloc(ntot, int), ci, cnz),
	    *Tj = expand_column_pointers(anc, cp, Calloc(ntot, int)),
	    *Ci = Calloc(ntot, int);

	for (j = 0, pos = cnz; j < anc; j++) {
	    int cj2 = cp[j + 1], ka, kc;
	    for (ka = Parent[j]; ka >= 0; ka = Parent[ka]) {
		for (kc = cp[j]; kc < cj2; kc++) {
		    if (check_csc_index(cp, ci, ci[kc], ka, -1) < 0) {
			Tj[pos] = ka;
			Ti[pos] = ci[kc];
			pos++;
		    }
		}
	    }
	}
	for (j = 0, cnr = 0; j < cnz; j++) { /* determine number of rows in C */
	    int rr = ci[j] + 1;
	    if (rr > cnr) cnr = rr;
	}
	triplet_to_col(cnr, anc, ntot, Ti, Tj, (double *) NULL,
		       INTEGER(cpp), Ci, (double *) NULL);
	cnz = cp[anc];
	SET_SLOT(C, Matrix_iSym, allocVector(INTSXP, cnz));
	SET_SLOT(C, Matrix_xSym, alloc3Darray(REALSXP, dims[0], dims[1], cnz));
	Free(Ti); Free(Tj); Free(Ci);
    }
}
    
/** 
 * Update a block of L in the blocked crosstabulation
 * 
 * @param L pointer to a unit lower triangular list of logical
 *          compressed sparse column-oriented matrices
 * @param ZZpO pointer to a list of upper triangular diagonal blocks
 *          stored as compressed sparse column-oriented matrices.
 * @param j index of updating column block
 * @param k column index of block to be updated 
 * @param i row index of block to be updated (j < k <= i)
 */
static void
block_update(SEXP L, SEXP ZZpO, int j, int k, int i)
{
    SEXP tb = (i == k) ? VECTOR_ELT(ZZpO, i) : VECTOR_ELT(L, Lind(i, k)),
	ib = VECTOR_ELT(L, Lind(i, j)),
	kb = VECTOR_ELT(L, Lind(k, j));
    SEXP tpp = GET_SLOT(tb, Matrix_pSym),
	kpp = GET_SLOT(kb, Matrix_pSym);
    int *ti = INTEGER(GET_SLOT(tb, Matrix_iSym)),
	*tp = INTEGER(tpp),
	*ii = INTEGER(GET_SLOT(ib, Matrix_iSym)),
	*ip = INTEGER(GET_SLOT(ib, Matrix_pSym)),
	*ki = INTEGER(GET_SLOT(kb, Matrix_iSym)),
	*kp = INTEGER(kpp),
	tnc = length(tpp) - 1,
	knc = length(kpp) - 1;
    int jj, extra;

    if (k > i || j >= k)
	error(_("i,j,k values of %d,%d,%d do not satisfy j < k <= i"),
	      i, j, k);
				/* bound the number of extra elements */
    extra = 0;
    for (jj = 0; jj < knc; jj++) {
	int i1, kk, i2 = ip[jj + 1], k2 = kp[jj + 1];
	for (kk = kp[jj]; kk < k2; kk++) {
	    for (i1 = ip[jj]; i1 < i2; i1++) {
		    if ((check_csc_index(tp, ti, ii[i1], ki[kk], -1) < 0) &&
				/* only update upper triangle of
				 * diagonal blocks */
			((k != i) || (ii[i1] <= ki[kk]))) extra++;
	    }
	}
    }
    if (!extra) return;
    {
	int pos, nnz = tp[tnc];
	int ntot = nnz + extra, tnr;
	int *Ai = Calloc(ntot, int),
	    *Ti = Calloc(ntot, int),
	    *Tj = Calloc(ntot, int),
	    *dims;
	double *Ax;

	Memcpy(Ti, ti, nnz);	/* make a copy of the row indices */
	for (pos = 0, jj = 0; jj < tnc; jj++) {	/* fill in the column indices */
	    int j2 = tp[jj + 1];
	    for (; pos < j2; pos++) Tj[pos] = jj;
	}
				/* add the extra elements */
	for (jj = 0; jj < knc; jj++) {
	    int i1, kk, i2 = ip[jj + 1], k2 = kp[jj + 1];
	    for (kk = kp[jj]; kk < k2; kk++) {
		for (i1 = ip[jj]; i1 < i2; i1++) {
		    if ((check_csc_index(tp, ti, ii[i1], ki[kk], -1) < 0) &&
			((k != i) || (ii[i1] <= ki[kk]))) { 
			Ti[pos] = ii[i1];
			Tj[pos] = ki[kk];
			pos++;
		    }
		}
	    }
	}
	/* FIXME: Pass nlev instead -  dimensions are nlev[i], nlev[k] */
	/* Determine maximum row index in T */
	tnr = -1; for (jj = 0; jj < ntot; jj++) if (Ti[jj] > tnr) tnr = Ti[jj];
	tnr++;			/* increment by 1 to get number of rows */
	triplet_to_col(tnr, tnc, ntot, Ti, Tj, (double *) NULL,
		       tp, Ai, (double *) NULL);
	nnz = tp[tnc];
	SET_SLOT(tb, Matrix_iSym, allocVector(INTSXP, nnz));
	Memcpy(INTEGER(GET_SLOT(tb, Matrix_iSym)), Ai, nnz);
	dims = INTEGER(getAttrib(GET_SLOT(tb, Matrix_xSym), R_DimSymbol));
	SET_SLOT(tb, Matrix_xSym,
		 alloc3Darray(REALSXP, dims[0], dims[1], nnz));
	Ax = REAL(GET_SLOT(tb, Matrix_xSym));
	for (j = 0; j < nnz; j++) Ax[j] = 1.;
	Free(Ai); Free(Ti); Free(Tj);
	return;
    }
}

/** 
 * Permute the levels of one of the grouping factors in a bCrosstab object
 * 
 * @param ctab Pointer to a bCrosstab object
 * @param nf number of factors in ctab
 * @param jj index (0-based) of the factor levels to permute
 * @param nlev number of levels of the grouping factors
 * @param iperm inverse of the permutation
 */
static void
bCrosstab_permute(SEXP ctab, int nf, int jj,
		  const int nlev[], const int iperm[])
{
    int j;
    for (j = 0; j < nf; j++) {
	int ind = (j < jj ? Lind(jj, j) : Lind(j, jj)),
	    ncol = (j < jj ? nlev[j] : nlev[jj]),
	    nrow = (j < jj ? nlev[jj] : nlev[j]);
	SEXP cscb = VECTOR_ELT(ctab, ind),
	    cscbi = GET_SLOT(cscb, Matrix_iSym);
	int *cp = INTEGER(GET_SLOT(cscb, Matrix_pSym)),
	    nnz = length(cscbi);
	double *cx = REAL(GET_SLOT(cscb, Matrix_xSym));
	int *mj = expand_column_pointers(ncol, cp, Calloc(nnz, int));
	int *mi = Memcpy(Calloc(nnz, int), INTEGER(cscbi), nnz);
	double *mx = Memcpy(Calloc(nnz, double), cx, nnz);

	if (j <= jj) int_permute(mi, nnz, iperm);
	if (j >= jj) int_permute(mj, nnz, iperm);
	if (j == jj) make_upper_triangular(mi, mj, nnz);
	triplet_to_col(nrow, ncol, nnz, mi, mj, mx, cp, INTEGER(cscbi), cx);
	Free(mi); Free(mj); Free(mx);
    }
}

static void
symmetric_permute(SEXP A, int nlev, const int iperm[])
{
    SEXP AiP = GET_SLOT(A, Matrix_iSym);
    int *Ap = INTEGER(GET_SLOT(A, Matrix_pSym)),
	nnz = length(AiP);
    double *Ax = REAL(GET_SLOT(A, Matrix_xSym));
    int *mj = expand_column_pointers(nlev, Ap, Calloc(nnz, int));
    int *mi = Memcpy(Calloc(nnz, int), INTEGER(AiP), nnz);
    double *mx = Memcpy(Calloc(nnz, double), Ax, nnz);

    int_permute(mi, nnz, iperm);
    int_permute(mj, nnz, iperm);
    make_upper_triangular(mi, mj, nnz);
    triplet_to_col(nlev, nlev, nnz, mi, mj, mx, Ap, INTEGER(AiP), Ax);
    Free(mi); Free(mj); Free(mx);
}

/** 
 * Apply a permutation vector to the levels of a factor.
 *
 * The dest pointer is assumed to point to a copy of the src pointer's
 * contents.
 * 
 * @param dest pointer to the destination factor
 * @param src pointer to the source factor
 * @param perm permutation vector (0-based)
 * @param iperm inverse permutation vector (0-based)
 */
static void
factor_levels_permute(SEXP dest, SEXP src, const int perm[],
		      const int iperm[])
{
    SEXP dlev = getAttrib(dest, R_LevelsSymbol),
	slev = getAttrib(src, R_LevelsSymbol);
    int nlev = length(dlev), flen = length(dest);
    int *d = INTEGER(dest), *s = INTEGER(src), i;

    if (length(slev) != nlev)
	error(_("number of levels in src and dest must match"));
    if (length(src) != flen)
	error(_("length of src and dest must match"));
    for (i = 0; i < nlev; i++)
	SET_STRING_ELT(dlev, i, STRING_ELT(slev, perm[i]));
    for (i = 0; i < flen; i++)
	d[i] = 1 + iperm[s[i]-1];
}

/** 
 * Create and populate slots in an lmer object from the blocked crosstabulation.
 * 
 * @param val Pointer to an lmer object
 */
void
lmer_populate(SEXP val)
{
    SEXP D, L, Parent, ZZpO, flist = GET_SLOT(val, Matrix_flistSym),
	perm, Omega, ZtZ = GET_SLOT(val, Matrix_ZtZSym);
    SEXP fnms = getAttrib(flist, R_NamesSymbol);
    int j, k, nf = length(flist);
    int *nc = INTEGER(GET_SLOT(val, Matrix_ncSym)), *Gp,
	*nlev = Calloc(nf, int), npairs = (nf * (nf + 1))/2;
    char *statnms[] = {"factored", "inverted", ""},
	*devnms[] = {"ML", "REML", ""},
	*pnms[] = {"index", "block", ""};
	
    /* Allocate fixed-sized slots */
    SET_SLOT(val, Matrix_statusSym, Matrix_make_named(LGLSXP, statnms));
    SET_SLOT(val, Matrix_devianceSym, Matrix_make_named(REALSXP, devnms));
    SET_SLOT(val, Matrix_devCompSym, allocVector(REALSXP, 4));
    /* Allocate slots that are lists of length nf */
    ZZpO = ALLOC_SLOT(val, Matrix_ZZpOSym, VECSXP, nf);
    setAttrib(ZZpO, R_NamesSymbol, duplicate(fnms));
    D = ALLOC_SLOT(val, Matrix_DSym, VECSXP, nf);
    setAttrib(D, R_NamesSymbol, duplicate(fnms));
    perm = ALLOC_SLOT(val, Matrix_permSym, VECSXP, nf);
    setAttrib(perm, R_NamesSymbol, duplicate(fnms));    
    Parent = ALLOC_SLOT(val, Matrix_ParentSym, VECSXP, nf);
    setAttrib(Parent, R_NamesSymbol, duplicate(fnms));
    Omega = ALLOC_SLOT(val, Matrix_OmegaSym, VECSXP, nf);
    setAttrib(Omega, R_NamesSymbol, duplicate(fnms));
    
    /* Allocate peculiar length slots */
    SET_SLOT(val, Matrix_LSym, allocVector(VECSXP, npairs));
    L = GET_SLOT(val, Matrix_LSym);
    SET_SLOT(val, Matrix_GpSym, allocVector(INTSXP, nf + 1));
    Gp = INTEGER(GET_SLOT(val, Matrix_GpSym));
    Gp[0] = 0;
    for (j = 0; j < nf; j++) {
	nlev[j] = length(getAttrib(VECTOR_ELT(flist, j), R_LevelsSymbol));
	Gp[j + 1] = Gp[j] + nc[j] * nlev[j];
	SET_VECTOR_ELT(D, j, alloc3Darray(REALSXP, nc[j], nc[j], nlev[j]));
	SET_VECTOR_ELT(Omega, j, allocMatrix(REALSXP, nc[j], nc[j]));
	SET_VECTOR_ELT(ZZpO, j, duplicate(VECTOR_ELT(ZtZ, Lind(j, j))));
	for (k = j; k < nf; k++)
	    SET_VECTOR_ELT(L, Lind(k, j),
			   duplicate(VECTOR_ELT(ZtZ, Lind(k, j))));
    }
    SET_SLOT(val, Matrix_XtXSym, allocMatrix(REALSXP, nc[nf], nc[nf]));
    AZERO(REAL(GET_SLOT(val, Matrix_XtXSym)), nc[nf] * nc[nf]);
    SET_SLOT(val, Matrix_RXXSym, allocMatrix(REALSXP, nc[nf], nc[nf]));
    AZERO(REAL(GET_SLOT(val, Matrix_RXXSym)), nc[nf] * nc[nf]);
    SET_SLOT(val, Matrix_ZtXSym, allocMatrix(REALSXP, Gp[nf], nc[nf]));
    SET_SLOT(val, Matrix_RZXSym, allocMatrix(REALSXP, Gp[nf], nc[nf]));
    for (j = 0; j < nf; j++) {
	int dind = Lind(j, j), i;
	SEXP ctd = VECTOR_ELT(ZZpO, j); /* diagonal in crosstab */
	SEXP Ljj = VECTOR_ELT(L, dind),
	    cpp = GET_SLOT(ctd, Matrix_pSym),
	    cip = GET_SLOT(ctd, Matrix_iSym), parent;
	int *Lp = INTEGER(GET_SLOT(Ljj, Matrix_pSym)), *Perm,
	    *cp = INTEGER(cpp),
	    *ci = INTEGER(cip), *dims,
	    ncj = length(cpp) - 1,
	    nnz = length(cip);
				
	SET_VECTOR_ELT(Parent, j, Matrix_make_named(VECSXP, pnms));
	parent = VECTOR_ELT(Parent, j);
	SET_VECTOR_ELT(parent, 0, allocVector(INTSXP, ncj));
	SET_VECTOR_ELT(parent, 1, allocVector(INTSXP, ncj));
	SET_VECTOR_ELT(perm, j, allocVector(INTSXP, ncj));
	Perm = INTEGER(VECTOR_ELT(perm, j));
	dims = INTEGER(getAttrib(GET_SLOT(ctd, Matrix_xSym), R_DimSymbol));
	if (nnz > ncj) {	/* calculate fill-reducing permutation */
	    SEXP fac = VECTOR_ELT(flist, j);
	    SEXP fcp = PROTECT(duplicate(fac));
	    int *iPerm = Calloc(ncj, int);

	    ssc_metis_order(ncj, cp, ci, Perm, iPerm);
				/* apply to the crosstabulation, L, and ZZpO */
	    bCrosstab_permute(ZtZ, nf, j, nlev, iPerm);
	    bCrosstab_permute(L, nf, j, nlev, iPerm);
	    symmetric_permute(VECTOR_ELT(ZZpO, j), nlev[j], iPerm);
				/* apply to the factor */
	    factor_levels_permute(fac, fcp, Perm, iPerm);
				/* symbolic analysis to get Parent */
	    R_ldl_symbolic(ncj, cp, ci, Lp, INTEGER(VECTOR_ELT(parent, 0)), 
			 (int *) NULL, (int *) NULL);
	    for (i = 0; i < ncj; i++)
		INTEGER(VECTOR_ELT(parent, 1))[i] =
		    (INTEGER(VECTOR_ELT(parent, 0))[i] < 0) ? -1 : j;
	    nnz = Lp[ncj];
	    SET_SLOT(Ljj, Matrix_iSym, allocVector(INTSXP, nnz));
	    SET_SLOT(Ljj, Matrix_xSym,
		     alloc3Darray(REALSXP, dims[0], dims[1], nnz));
	    Free(iPerm); UNPROTECT(1);
	} else {
	    for (i = 0; i < ncj; i++) {
		Lp[i] = 0;
		INTEGER(VECTOR_ELT(parent,0))[i] = -1;
		INTEGER(VECTOR_ELT(parent,1))[i] = -1;
		Perm[i] = i;
	    }
	    Lp[ncj] = 0;
	    SET_SLOT(Ljj, Matrix_iSym, allocVector(INTSXP, 0));
	    SET_SLOT(Ljj, Matrix_xSym,
		     alloc3Darray(REALSXP, dims[0], dims[1], 0));
	}
	for (k = j+1; k < nf; k++) { /* Update other blocks in this column */
	    symbolic_right_unit_mm_trans(ncj, INTEGER(VECTOR_ELT(parent, 0)),
					 VECTOR_ELT(L, Lind(k,j)));
	}
	for (k = j+1; k < nf; k++) { /* Update remaining columns */
	    for (i = k; i < nf; i++) block_update(L, ZZpO, j, k, i);
	}
    }
    /* Convert blockwise Parent arrays to extended Parent arrays */
    for (j = 0; j < (nf - 1); j++) { /* Parent[nf] does not need conversion */
	SEXP Ljp1j = VECTOR_ELT(L, Lind(j + 1, j)),
	    LpP = GET_SLOT(Ljp1j, Matrix_pSym);
	int *Li = INTEGER(GET_SLOT(Ljp1j, Matrix_iSym)),
	    *Lp = INTEGER(LpP),
	    *block = INTEGER(VECTOR_ELT(VECTOR_ELT(Parent, j), 1)),
	    *parent = INTEGER(VECTOR_ELT(VECTOR_ELT(Parent, j), 0)),
	    i, nlev = length(LpP) - 1;
	for (i = 0; i < nlev; i++) {
	    if (block[i] < 0) {
		block[i] = j + 1;
		parent[i] = Li[Lp[i]];
	    }
	}	    
    }
    Free(nlev);
}
