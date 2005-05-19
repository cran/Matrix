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
	int tnr, nnz = tp[tnc];
	int ntot = nnz + extra, pos = nnz;
	int *Ai = Calloc(ntot, int),
	    *Ti = Calloc(ntot, int),
	    *Tj = Calloc(ntot, int),
	    *dims;
	double *Ax;

	Memcpy(Ti, ti, nnz);	/* make a copy of the row indices */
	expand_cmprPt(tnc, tp, Tj); /* fill in existing column indices */
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
	Memcpy(INTEGER(ALLOC_SLOT(tb, Matrix_iSym, INTSXP, nnz)), Ai, nnz);
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
 * Convert the parent array to the row indices (But it doesn't work.)
 * 
 * @param n order of the triangular matrix
 * @param p array of length n + 1 of column pointers (for checking)
 * @param parent parent array (length n)
 * @param i array of row indices to be written
 *
 * @return i
 */
int *
Parent2rows(int n, const int p[], const int parent[], int i[])
{
    int ii, j, par;

    for (j = 0; j < n; j++) {
	for(ii = p[j], par = parent[j]; par >= 0; ii++, par = parent[par]) {}
	if (ii != p[j + 1])
	    error(_("p[%d] = %d but parent array gives %d"),
		  j + 1, p[j + 1], ii);
    }
    return i;
}

/** 
 * Replace the structure of C by the structure of CA^{-T}
 * 
 * @param anc number of column blocks in A
 * @param Parent parent array for column blocks of A
 * @param C a dgBCMatrix object to be updated
 */
static void
symbolic_right_unit_sm_trans(int anc, const int Parent[], SEXP C)
{
    SEXP cip = GET_SLOT(C, Matrix_iSym),
	cpp = GET_SLOT(C, Matrix_pSym), tmp;
    int *ci = INTEGER(cip), *cj, *Ti, *Tp,
	*cp = INTEGER(cpp), cnr,
	cnz = length(cip),
	i, j;

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

				/* determine number of rows in C */
    cnr = -1;
    for (i = 0; i < cnz; i++) {
	int cii = ci[i];
	if (cii > cnr) cnr = cii;
    }
    cnr++;			/* max 0-based index is one less the
				 * no. of rows */
				
    cj = expand_cmprPt(anc, cp, Calloc(cnz, int));
    Ti = Calloc(cnz, int);
    Tp = Calloc(cnr, int);	/* transpose C */
    triplet_to_col(anc, cnr, cnz, cj, ci, (double *) NULL,
		   Tp, Ti, (double *) NULL);
    cj = Realloc(cj, cnr, int);	/* Create column pointers for empty matrix */
    for (i = 0; i < cnr; i++) cj[i] = 0;
    PROTECT(tmp = lCholClgCsm(LFT, NTR, anc, cnr, Parent, Ti, Tp, cip, cj));
				/* transpose the result */
    cnz = cj[cnr];
    Free(Tp);
    Tp = expand_cmprPt(cnr, cj, Calloc(cnz, int));
    triplet_to_col(cnr, anc, cnz, Tp, INTEGER(tmp), (double *) NULL,
		   cp, INTEGER(ALLOC_SLOT(C, Matrix_iSym, INTSXP, cnz)),
		   (double *) NULL);
    UNPROTECT(1);
}

    
/** 
 * Update a diagonal block of ZZpO in the blocked crosstabulation
 * 
 * @param db pointer to the diagonal block
 * @param odb pointer to the off-diagonal block
 */
static R_INLINE void
diag_update(SEXP db, SEXP odb, int n, int k)
{
    SET_SLOT(db, Matrix_iSym,
	     Matrix_lgCsyrk(1, 0, n, k,
			    INTEGER(GET_SLOT(odb, Matrix_iSym)),
			    INTEGER(GET_SLOT(odb, Matrix_pSym)),
			    1,
			    GET_SLOT(db, Matrix_iSym),
			    INTEGER(GET_SLOT(db, Matrix_pSym))));
}

/** 
 * Update an off-diagonal block of L from the blocked crosstabulation
 * 
 * @param A lower block
 * @param B upper block
 * @param C product block
 * @param nrA number of rows in A
 */
static R_INLINE void
offdiag_update(SEXP A, SEXP B, SEXP C, int m, int n, int k)
{
    SET_SLOT(C, Matrix_iSym,
	     Matrix_lgClgCmm(0, 1, m, n, k,
			     INTEGER(GET_SLOT(A, Matrix_iSym)),
			     INTEGER(GET_SLOT(A, Matrix_pSym)),
			     INTEGER(GET_SLOT(B, Matrix_iSym)),
			     INTEGER(GET_SLOT(B, Matrix_pSym)),
			     1,
			     GET_SLOT(C, Matrix_iSym),
			     INTEGER(GET_SLOT(C, Matrix_pSym))));
}
