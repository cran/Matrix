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
