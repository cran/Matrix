#include "ssclme.h"

#define slot_dup(dest, src, sym)  SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))

/** 
 * Using the sscCrosstab object from the grouping factors, generate
 * the slots in an ssclme object related to the symmetric sparse
 * matrix representation of Z'Z.  If the model matrices for the
 * grouping factors have only one column each then the structure can
 * be copied, otherwise it must be generated from the sscCrosstab and
 * the number of columns per grouping factor.
 * 
 * @param nf number of factors
 * @param nc vector of length nf+2 with number of columns in model matrices
 * @param ctab pointer to the sscCrosstab object
 * @param ssc pointer to an ssclme object to be filled out
 */
static
void ssclme_copy_ctab(int nf, const int nc[], SEXP ctab, SEXP ssc)
{
    int *snc, i, copyonly = 1;

    SET_SLOT(ssc, Matrix_ncSym, allocVector(INTSXP, nf + 2));
    snc = INTEGER(GET_SLOT(ssc, Matrix_ncSym));
    for (i = 0; i <= nf; i++) {
	snc[i] = nc[i];
	if (nc[i] > 1 && i < nf) copyonly = 0;
    }
    if (copyonly) {
	slot_dup(ssc, ctab, Matrix_pSym);
	slot_dup(ssc, ctab, Matrix_iSym);
	slot_dup(ssc, ctab, Matrix_xSym);
	slot_dup(ssc, ctab, Matrix_DimSym);
	slot_dup(ssc, ctab, Matrix_GpSym);
	return;
    }
    {
	int
	    *GpIn = INTEGER(GET_SLOT(ctab, Matrix_GpSym)),
	    *GpOut,
	    *AiIn = INTEGER(GET_SLOT(ctab, Matrix_iSym)),
	    *AiOut,
	    *ApIn = INTEGER(GET_SLOT(ctab, Matrix_pSym)),
	    *ApOut,
	    nIn = GpIn[nf], nOut, nzOut,
	    *dims,
	    *map = Calloc(nIn + 1, int), /* column number map */
	    *ncc = Calloc(nIn, int); /* number of columns out for this
				      * col in */

	SET_SLOT(ssc, Matrix_GpSym, allocVector(INTSXP, nf + 1));
	GpOut = INTEGER(GET_SLOT(ssc, Matrix_GpSym));
	map[0] = GpOut[0] = 0;
	for (i = 0; i < nf; i++) {
	    int j, GpIni = GpIn[i], GpInip1 = GpIn[i+1], nci = nc[i];
	    GpOut[i+1] = GpOut[i] + (GpInip1 - GpIni) * nci;
	    for (j = GpIni; j < GpInip1; j++) {
		ncc[j] = nci;
		map[j+1] = map[j] + nci;
	    }
	}
	nOut = GpOut[nf];	/* size of output matrix */
	SET_SLOT(ssc, Matrix_DimSym, allocVector(INTSXP, 2));
	dims = INTEGER(GET_SLOT(ssc, Matrix_DimSym));
	dims[0] = dims[1] = nOut;
	SET_SLOT(ssc, Matrix_pSym, allocVector(INTSXP, nOut + 1));
	ApOut = INTEGER(GET_SLOT(ssc, Matrix_pSym));
	ApOut[0] = 0;
	for (i = 0; i < nf; i++) { /* determine the column pointers */
	    int j, jout = GpOut[i], nci = nc[i], p2 = GpIn[i+1];
	    for (j = GpIn[i]; j < p2; j++) {
		int k, nj = 0, p3 = ApIn[j+1];
		for (k = ApIn[j]; k < p3; k++) {
		    nj += ncc[AiIn[k]];
		}
		nj -= nci - 1;
		ApOut[jout+1] = ApOut[jout] + nj;
		jout++;
		for (k = 1; k < nci; k++) {
		    ApOut[jout+1] = ApOut[jout] + nj + k;
		    jout++;
		}
	    }
	}
	nzOut = ApOut[nOut];	/* number of non-zeros in output */
	SET_SLOT(ssc, Matrix_xSym, allocVector(REALSXP, nzOut));
	memset(REAL(GET_SLOT(ssc, Matrix_xSym)), 0,
	       sizeof(double) * nzOut);
	SET_SLOT(ssc, Matrix_iSym, allocVector(INTSXP, nzOut));
	AiOut = INTEGER(GET_SLOT(ssc, Matrix_iSym));
	for (i = 0; i < nf; i++) { /* fill in the rows */
	    int j, jj, nci = nc[i], p2 = GpIn[i+1];
	    for (j = GpIn[i]; j < p2; j++) { /* first col for input col */
		int k, mj = map[j], p3 = ApIn[j+1], pos = ApOut[mj];
		for (k = ApIn[j]; k < p3; k++) {
		    int ii = AiIn[k], ncci = ncc[ii];
		    AiOut[pos++] = map[ii];
		    if (ii < j) {	/* above the diagonal */
			for (jj = 1; jj < ncci; jj++) {
			    AiOut[pos+1] = AiOut[pos] + 1;
			    pos++;
			}
		    }
		    for (jj = 1; jj < nci; jj++) { /* repeat the column adding 
						    * another diagonal element */
			int mjj = mj + jj, pj = ApOut[mjj], pjm1 = ApOut[mjj-1];
			Memcpy(AiOut + pj, AiOut + pjm1, pj - pjm1);
			AiOut[ApOut[mjj + 1] - 1] = mjj;
		    }
		}
	    }
	}
	Free(map); Free(ncc);
    }
}

/** 
 * Calculate and store the maximum number of off-diagonal elements in
 * the inverse of L, based on the elimination tree.  The maximum is
 * itself stored in the Parent array.  (FIXME: come up with a better design.)
 * 
 * @param n number of columns in the matrix
 * @param Parent elimination tree for the matrix
 */
static void ssclme_calc_maxod(int n, int Parent[])
{
    int *sz = Calloc(n, int), i, mm = -1;
    for (i = n - 1; i >= 0; i--) {
	sz[i] = (Parent[i] < 0) ? 0 : (1 + sz[Parent[i]]);
	if (sz[i] > mm) mm = sz[i];
    }
    Parent[n] = mm;
    Free(sz);
}

/** 
 * Create an ssclme object from a list of grouping factors, sorted in
 * order of non-increasing numbers of levels, and an integer vector of
 * the number of columns in the model matrices.  There is one more
 * element in ncv than in facs.  The last element is the number of
 * columns in the model matrix for the fixed effects plus the
 * response.  (i.e. p+1)
 * 
 * @param facs pointer to a list of grouping factors
 * @param ncv pointer to an integer vector of number of columns per model matrix
 * 
 * @return pointer to an ssclme object
 */
SEXP
ssclme_create(SEXP facs, SEXP ncv)
{
    SEXP ctab, nms, ssc, tmp,
	val = PROTECT(allocVector(VECSXP, 2)),
	dd = PROTECT(allocVector(INTSXP, 3));	/* dimensions of 3-D arrays */
    int *Ai, *Ap, *Gp, *Lp, *Parent,
	*nc, Lnz, i, nf = length(facs), nzcol, pp1,
	*dims = INTEGER(dd);

    if (length(ncv) != (nf + 1))
	error("length of nc (%d) should be length of facs (%d) + 1",
	      length(ncv), nf);
    SET_VECTOR_ELT(val, 0, NEW_OBJECT(MAKE_CLASS("ssclme")));
    ssc = VECTOR_ELT(val, 0);
				/* Pairwise cross-tabulation */
    ctab = PROTECT(sscCrosstab(facs, ScalarLogical(1)));
    SET_VECTOR_ELT(val, 1, sscCrosstab_groupedPerm(ctab));
    if (length(VECTOR_ELT(val, 1)) > 0) {/* Fill-reducing permutation */
	ssc_symbolic_permute(INTEGER(GET_SLOT(ctab, Matrix_DimSym))[1],
			     1, INTEGER(VECTOR_ELT(val, 1)),
			     INTEGER(GET_SLOT(ctab, Matrix_pSym)),
			     INTEGER(GET_SLOT(ctab, Matrix_iSym)));
    }
    ssclme_copy_ctab(nf, INTEGER(ncv), ctab, ssc);
    UNPROTECT(1);		/* ctab */

    nzcol = INTEGER(GET_SLOT(ssc, Matrix_DimSym))[1];
    Gp = INTEGER(GET_SLOT(ssc, Matrix_GpSym));
    Ap = INTEGER(GET_SLOT(ssc, Matrix_pSym));
    Ai = INTEGER(GET_SLOT(ssc, Matrix_iSym));
    nc = INTEGER(GET_SLOT(ssc, Matrix_ncSym));
    nc[nf + 1] = length(VECTOR_ELT(facs, 0)); /* number of observations */
				/* Create slots */
    pp1 = nc[nf];
    SET_SLOT(ssc, Matrix_XtXSym, allocMatrix(REALSXP, pp1, pp1));
    SET_SLOT(ssc, Matrix_RXXSym, allocMatrix(REALSXP, pp1, pp1));
    SET_SLOT(ssc, Matrix_ZtXSym, allocMatrix(REALSXP, nzcol, pp1));
    SET_SLOT(ssc, Matrix_RZXSym, allocMatrix(REALSXP, nzcol, pp1));
				/* Zero symmetric matrices (cosmetic) */
    memset(REAL(GET_SLOT(ssc, Matrix_XtXSym)), 0,
	   sizeof(double) * pp1 * pp1); 
    memset(REAL(GET_SLOT(ssc, Matrix_RXXSym)), 0,
	   sizeof(double) * pp1 * pp1);
    SET_SLOT(ssc, Matrix_LpSym, allocVector(INTSXP, nzcol + 1));
    Lp = INTEGER(GET_SLOT(ssc, Matrix_LpSym));
    SET_SLOT(ssc, Matrix_ParentSym, allocVector(INTSXP, nzcol + 1));
    Parent = INTEGER(GET_SLOT(ssc, Matrix_ParentSym));
    SET_SLOT(ssc, Matrix_DSym, allocVector(REALSXP, nzcol));
    SET_SLOT(ssc, Matrix_DIsqrtSym, allocVector(REALSXP, nzcol));
    ldl_symbolic(nzcol, Ap, Ai, Lp, Parent,
		 (int *) R_alloc(nzcol, sizeof(int)), /* Lnz */
		 (int *) R_alloc(nzcol, sizeof(int)), /* Flag */
		 (int *) NULL, (int *) NULL); /* P & Pinv */
    ssclme_calc_maxod(nzcol, Parent);
    Lnz = Lp[nzcol];
    SET_SLOT(ssc, Matrix_LiSym, allocVector(INTSXP, Lnz));
    SET_SLOT(ssc, Matrix_LxSym, allocVector(REALSXP, Lnz));
    SET_SLOT(ssc, Matrix_OmegaSym, allocVector(VECSXP, nf));
    tmp = GET_SLOT(ssc, Matrix_OmegaSym);
    setAttrib(tmp, R_NamesSymbol, getAttrib(facs, R_NamesSymbol));
    for (i = 0; i < nf; i++) {
	SET_VECTOR_ELT(tmp, i, allocMatrix(REALSXP, nc[i], nc[i]));
	memset(REAL(VECTOR_ELT(tmp, i)), 0,
	       sizeof(double) * nc[i] * nc[i]);
    }
    SET_SLOT(ssc, Matrix_devianceSym, allocVector(REALSXP, 2));
    tmp = GET_SLOT(ssc, Matrix_devianceSym);
    setAttrib(tmp, R_NamesSymbol, allocVector(STRSXP, 2));
    nms = getAttrib(tmp, R_NamesSymbol);
    SET_STRING_ELT(nms, 0, mkChar("ML"));
    SET_STRING_ELT(nms, 1, mkChar("REML"));
    SET_SLOT(ssc, Matrix_devCompSym, allocVector(REALSXP, 4));
    SET_SLOT(ssc, Matrix_statusSym, allocVector(LGLSXP, 2));
    tmp = GET_SLOT(ssc, Matrix_statusSym);
    LOGICAL(tmp)[0] = LOGICAL(tmp)[1] = 0;
    setAttrib(tmp, R_NamesSymbol, allocVector(STRSXP, 2));
    nms = getAttrib(tmp, R_NamesSymbol);
    SET_STRING_ELT(nms, 0, mkChar("factored"));
    SET_STRING_ELT(nms, 1, mkChar("inverted"));
    SET_SLOT(ssc, Matrix_bVarSym, allocVector(VECSXP, nf));
    tmp = GET_SLOT(ssc, Matrix_bVarSym);
    setAttrib(tmp, R_NamesSymbol, getAttrib(facs, R_NamesSymbol));
    for (i = 0; i < nf; i++) {
	int nci = nc[i], mi = (Gp[i+1] - Gp[i])/nc[i];

	dims[0] = dims[1] = nci;
	dims[2] = mi;
	SET_VECTOR_ELT(tmp, i, allocArray(REALSXP, dd));
	memset(REAL(VECTOR_ELT(tmp, i)), 0,
	       sizeof(double) * nci * nci * mi);
    }
    UNPROTECT(2);
    return val;
}

/** 
 * Copy information on Z'Z accumulated in the bVar array to Z'Z
 * 
 * @param ncj number of columns in this block
 * @param Gpj initial column for this group
 * @param Gpjp initial column for the next group
 * @param bVj pointer to the ncj x ncj x mj array to be filled
 * @param Ap column pointer array for Z'Z
 * @param Ai row indices for Z'Z
 * @param Ax elements of Z'Z
 */
static
void bVj_to_A(int ncj, int Gpj, int Gpjp, const double bVj[],
	      const int Ap[], const int Ai[], double Ax[])
{
    int i, diag, k;
    for (i = Gpj; i < Gpjp; i += ncj) {
	for (k = 0; k < ncj; k++) {
	    diag = Ap[i + k + 1] - 1;
	    if (Ai[diag] != i+k)
		error("Expected Ai[%d] to be %d (i.e on diagonal) not %d",
		      diag, i+k, Ai[diag]);
	    Memcpy(Ax + diag - k, bVj + (i+k-Gpj)*ncj, k + 1);
	}
    }
}

/** 
 * Copy the dimnames from the list of grouping factors and the model
 * matrices for the grouping factors into the appropriate parts of the
 * ssclme object.
 * 
 * @param x pointer to an ssclme object
 * @param facs pointer to a list of factors
 * @param mmats pointer to a list of model matrices
 * 
 * @return NULL
 */
SEXP
ssclme_transfer_dimnames(SEXP x, SEXP facs, SEXP mmats)
{
    SEXP bVar = GET_SLOT(x, Matrix_bVarSym),
	nms2 = PROTECT(allocVector(VECSXP, 2)),
	nms3 = PROTECT(allocVector(VECSXP, 3));
    int i, nf = length(mmats) - 1;
    SEXP xcols = VECTOR_ELT(GetArrayDimnames(VECTOR_ELT(mmats, nf)), 1);

    for (i = 0; i < nf; i++) {
	SEXP cnms = VECTOR_ELT(GetArrayDimnames(VECTOR_ELT(mmats, i)), 1);
	SET_VECTOR_ELT(nms3, 0, cnms);
	SET_VECTOR_ELT(nms3, 1, cnms);
	SET_VECTOR_ELT(nms3, 2,
		       getAttrib(VECTOR_ELT(facs, i), R_LevelsSymbol));
	dimnamesgets(VECTOR_ELT(bVar, i), duplicate(nms3));
    }
    SET_VECTOR_ELT(nms2, 0, xcols);
    SET_VECTOR_ELT(nms2, 1, xcols);
    dimnamesgets(GET_SLOT(x, Matrix_XtXSym), nms2);
    dimnamesgets(GET_SLOT(x, Matrix_RXXSym), nms2);
    UNPROTECT(2);
    return R_NilValue;
}

/** 
 * Update the numerical entries x, ZtX, and XtX in an ssclme object
 * according to a set of model matrices.
 * 
 * @param x pointer to an ssclme object
 * @param facs pointer to a list of grouping factors
 * @param mmats pointer to a list of model matrices
 * 
 * @return NULL
 */
SEXP
ssclme_update_mm(SEXP x, SEXP facs, SEXP mmats)
{
    SEXP bVar = GET_SLOT(x, Matrix_bVarSym);
    int
	*Ai = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(x, Matrix_pSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	i, j, k,
	ione = 1,
	nf = length(mmats) - 1,
	nobs = nc[nf + 1],
	nzcol = Gp[nf],
	nz = Ap[nzcol],
	pp1 = nc[nf];
    double
	**Z = Calloc(nf + 1, double *),
	*Ax = REAL(GET_SLOT(x, Matrix_xSym)),
	*XtX = REAL(GET_SLOT(x, Matrix_XtXSym)),
	*ZtX = REAL(GET_SLOT(x, Matrix_ZtXSym)),
	one = 1.0,
	zero = 0.0;

    for (i = 0; i <= nf; i++) {
	int *dims = INTEGER(getAttrib(VECTOR_ELT(mmats, i), R_DimSymbol)),
	    nci = nc[i];
	if (nobs != dims[0])
	    error("Expected %d rows in the %d'th model matrix. Got %d",
		  nobs, i+1, dims[0]);
	if (nci != dims[1])
	    error("Expected %d columns in the %d'th model matrix. Got %d",
		  nci, i+1, dims[1]);
	Z[i] = REAL(VECTOR_ELT(mmats, i));
    }
				/* Create XtX - X is Z[nf] */
    F77_CALL(dsyrk)("U", "T", nc+nf, &nobs, &one,
		    Z[nf], &nobs, &zero, XtX, nc + nf);
				/* Zero the accumulators */
    memset((void *) ZtX, 0, sizeof(double) * nzcol * pp1);
    memset((void *) Ax, 0, sizeof(double) * nz);
    for (j = 0; j < nf; j++) { /* Create ZtX */
	int *fpj = INTEGER(VECTOR_ELT(facs, j)), ncj = nc[j],
	    Ncj = ncj > 1;
	double
	    *bVj = REAL(VECTOR_ELT(bVar, j)),
	    *Zj = Z[j],
	    *zxj = ZtX + Gp[j];

	if (Ncj) {		/* bVj will accumulate Z'Z blocks */
	    memset(bVj, 0, sizeof(double) * ncj * (Gp[j+1]-Gp[j]));
	}
	for (i = 0; i < nobs; i++) { /* accumulate diagonal of ZtZ */
	    int fpji = fpj[i] - 1, /* factor indices are 1-based */
		dind = Ap[Gp[j] + fpji * ncj + 1] - 1;
	    if (Ai[dind] != (Gp[j] + fpji * ncj))
		error("logic error in ssclme_update_mm");
	    if (Ncj) {		/* use bVar to accumulate */
		F77_CALL(dsyrk)("U", "T", &ncj, &ione, &one, Zj+i,
				&nobs, &one, bVj + fpji*ncj*ncj, &ncj);
	    } else {		/* update scalars directly */
		Ax[dind] += Zj[i] * Zj[i];
	    }
				/* update rows of Z'X */
	    F77_CALL(dgemm)("T", "N", &ncj, &pp1, &ione, &one,
			    Zj + i, &nobs, Z[nf] + i, &nobs,
			    &one, zxj + fpji * ncj, &nzcol);
	}
	if (Ncj) bVj_to_A(ncj, Gp[j], Gp[j+1], bVj, Ap, Ai, Ax);
	for (k = j+1; k < nf; k++) { /* off-diagonals */
	    int *fpk = INTEGER(VECTOR_ELT(facs, k)),
		*Apk = Ap + Gp[k],
		nck = nc[k],
		scalar = ncj == 1 && nck == 1;
	    double
		*Zk = Z[k], *work;
	    if (!scalar) work = Calloc(ncj * nck, double);
	    for (i = 0; i < nobs; i++) {
		int ii, ind = -1, fpji = fpj[i] - 1,
		    row = Gp[j] + fpji * ncj,
		    fpki = fpk[i] - 1,
		    lastind = Apk[fpki*nck + 1];
		for (ii = Apk[fpki*nck]; ii < lastind; ii++) {
		    if (Ai[ii] == row) {
			ind = ii;
			break;
		    }
		}
		if (ind < 0) error("logic error in ssclme_update_mm");
		if (scalar) {	/* update scalars directly */
		    Ax[ind] += Zj[i] * Zk[i];
		} else {
		    int jj, offset = ind - Apk[fpki * nck];
		    F77_CALL(dgemm)("T", "N", &ncj, &nck, &ione, &one,
				    Zj + i, &nobs, Zk + i, &nobs,
				    &zero, work, &ncj);
		    for (jj = 0; jj < nck; jj++) {
			ind = Apk[fpki * nck + jj] + offset;
			if (Ai[ind] != row)
			    error("logic error in ssclme_update_mm");
			for (ii = 0; ii < ncj; ii++) {
			    Ax[ind++] += work[jj * ncj + ii];
			}
		    }
		}
	    }
	    if (!scalar) Free(work);
	}
    }
    Free(Z);
    ssclme_transfer_dimnames(x, facs, mmats);
    status[0] = status[1] = 0;
    return R_NilValue;
}

/** 
 * Inflate Z'Z according to Omega and create the factorization LDL'
 * 
 * @param x pointer to an ssclme object
 * 
 * @return NULL
 */
SEXP ssclme_inflate_and_factor(SEXP x)
{
    SEXP
	GpSlot = GET_SLOT(x, Matrix_GpSym),
	Omega = GET_SLOT(x, Matrix_OmegaSym);
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[1];
    int
	*Ai = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(x, Matrix_pSym)),
	*Flag = Calloc(n, int),
	*Gp = INTEGER(GpSlot),
	*Lnz = Calloc(n, int),
	*Pattern = Calloc(n, int),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	j,
	nf = length(GpSlot) - 1;
    double
	*D = REAL(GET_SLOT(x, Matrix_DSym)),
	*DIsqrt = REAL(GET_SLOT(x, Matrix_DIsqrtSym)),
	*Y = Calloc(n, double),
	*xcp = Calloc(Ap[n], double);

    Memcpy(xcp, REAL(GET_SLOT(x, Matrix_xSym)), Ap[n]);
    for (j = 0; j < nf; j++) {
	int  diag, i, ii, k, G2 = Gp[j + 1], ncj = nc[j];
	double *omgj = REAL(VECTOR_ELT(Omega, j));
	
	for (i = Gp[j]; i < G2; i += ncj) {
	    for (k = 0; k < ncj; k++) {
		diag = Ap[i + k + 1] - 1;
		if (Ai[diag] != i+k)
		    error("Expected Ai[%d] to be %d (i.e on diagonal) not %d",
			  diag, i+k, Ai[diag]);
		for (ii = 0; ii <= k; ii++) {
		    xcp[diag + ii - k] += omgj[k*ncj + ii];
		}
	    }
	}
    }
    j = ldl_numeric(n, Ap, Ai, xcp,
		    INTEGER(GET_SLOT(x, Matrix_LpSym)),
		    INTEGER(GET_SLOT(x, Matrix_ParentSym)),
		    Lnz, INTEGER(GET_SLOT(x, Matrix_LiSym)),
		    REAL(GET_SLOT(x, Matrix_LxSym)),
		    D, Y, Pattern, Flag,
		    (int *) NULL, (int *) NULL); /* P & Pinv */
    if (j != n)
	error("rank deficiency of ZtZ+W detected at column %d",
	      j + 1);
    for (j = 0; j < n; j++) DIsqrt[j] = 1./sqrt(D[j]);
    Free(Lnz); Free(Flag); Free(Pattern); Free(Y); Free(xcp);
    return R_NilValue;
}


/** 
 * If status[["factored"]] is FALSE, create and factor Z'Z+Omega, then
 * create RZX and RXX, the deviance components, and the value of the
 * deviance for both ML and REML.
 * 
 * @param x pointer to an ssclme object
 * 
 * @return NULL
 */
SEXP ssclme_factor(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    
    if (!status[0]) {
	SEXP
	    GpSlot = GET_SLOT(x, Matrix_GpSym),
	    Omega = GET_SLOT(x, Matrix_OmegaSym);
	int
	    *Gp = INTEGER(GpSlot),
	    *Li = INTEGER(GET_SLOT(x, Matrix_LiSym)),
	    *Lp = INTEGER(GET_SLOT(x, Matrix_LpSym)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    i,
	    n = INTEGER(GET_SLOT(x, Matrix_DimSym))[1],
	    nf = length(GpSlot) - 1,
	    nobs = nc[nf + 1],
	    nreml = nobs + 1 - nc[nf],
	    pp1 = nc[nf],
	    pp2 = pp1 + 1;
	double
	    *D = REAL(GET_SLOT(x, Matrix_DSym)),
	    *DIsqrt = REAL(GET_SLOT(x, Matrix_DIsqrtSym)),
	    *Lx = REAL(GET_SLOT(x, Matrix_LxSym)),
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(GET_SLOT(x, Matrix_RZXSym)),
	    *dcmp = REAL(getAttrib(x, Matrix_devCompSym)),
	    *deviance = REAL(getAttrib(x, Matrix_devianceSym)),
	    minus1 = -1.,
	    one = 1.;
	
	ssclme_inflate_and_factor(x);
				/* Accumulate logdet of ZtZ+W */
	dcmp[0] = dcmp[1] = dcmp[2] = dcmp[3] = 0.;
	for (i = 0; i < n; i++) dcmp[0] += log(D[i]);
				/* Accumulate logdet of W */
	for (i = 0; i < nf; i++) {
	    int nci = nc[i],
		mi = (Gp[i+1] - Gp[i])/nci;
	
	    if (nci < 2) {
		dcmp[1] += mi * log(REAL(VECTOR_ELT(Omega, i))[0]);
	    } else {
		int j;
		double
		    *tmp = Calloc(nci * nci, double),
		    accum = 0.;
		F77_CALL(dpotrf)("U", &nci,
				 Memcpy(tmp, REAL(VECTOR_ELT(Omega, i)),
					nci * nci),
				 &nci, &j);
		if (j) 
		    error("Omega[%d] is not positive definite", i + 1);
		for (j = 0; j < nci; j++) {
		    accum += 2 * log(tmp[j * (nci + 1)]);
		}
		dcmp[1] += mi * accum;
		Free(tmp);
	    }
	}
				/* ldl_lsolve on Z'X */
	Memcpy(RZX, REAL(GET_SLOT(x, Matrix_ZtXSym)), n * pp1);
	for (i = 0; i < pp1; i++) {
	    int j;
	    double *RZXi = RZX + i * n;
	    ldl_lsolve(n, RZXi, Lp, Li, Lx);
	    for (j = 0; j < n; j++) RZXi[j] *= DIsqrt[j];
	}
				/* downdate and factor X'X */
	Memcpy(RXX, REAL(GET_SLOT(x, Matrix_XtXSym)), pp1 * pp1);
	F77_CALL(dsyrk)("U", "T", &pp1, &n, &minus1,
			RZX, &n, &one, RXX, &pp1);
	F77_CALL(dpotrf)("U", &pp1, RXX, &pp1, &i);
	if (i) {
	    warning("Could not factor downdated X'X, code %d", i);
	    dcmp[2] = dcmp[3] = deviance[0] = deviance[1] = NA_REAL;
	} else {
				/* logdet of RXX */
	    for (i = 0; i < (pp1 - 1); i++)
		dcmp[2] += 2 * log(RXX[i*pp2]);
				/* logdet of Ryy */
	    dcmp[3] = 2. * log(RXX[pp1 * pp1 - 1]);
	    deviance[0] =	/* ML criterion */
		dcmp[0] - dcmp[1] + nobs*(1+dcmp[3]+log(2*PI/nobs));
	    deviance[1] = dcmp[0] - dcmp[1] + /* REML */
		dcmp[2] + nreml * (1. + dcmp[3] + log(2. * PI/nreml));
	}
	status[0] = 1;
	status[1] = 0;
    }
    return R_NilValue;
}

/** 
 * Return the position of probe in the sorted index vector ind.  It is
 * known that the position is greater than or equal to start so a linear
 * search from start is used.
 * 
 * @param probe value to be matched
 * @param start index at which to start
 * @param ind vector of indices
 * 
 * @return index of the entry matching probe
 */
static
int ldl_update_ind(int probe, int start, const int ind[])
{
    while (ind[start] < probe) start++;
    if (ind[start] > probe) error("logic error in ldl_inverse");
    return start;
}

/** 
 * Update the diagonal blocks of the inverse of LDL' (=Z'Z+W).  The
 * lower Cholesky factors of the updated blocks are stored in the bVar
 * slot.
 * 
 * @param x pointer to an ssclme object
 *
 * @return R_NilValue (x is updated in place)

 */
static
SEXP ldl_inverse(SEXP x)
{
    SEXP
	Gpsl = GET_SLOT(x, Matrix_GpSym),
	bVar = GET_SLOT(x, Matrix_bVarSym);
    int *Gp = INTEGER(Gpsl),
	*Parent = INTEGER(GET_SLOT(x, Matrix_ParentSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i,
	nf = length(Gpsl) - 1,
	nzc = INTEGER(GET_SLOT(x, Matrix_DimSym))[1];
    int maxod = Parent[nzc];
    double *DIsqrt = REAL(GET_SLOT(x, Matrix_DIsqrtSym));
	
    ssclme_factor(x);
    if (maxod == 0) {		/* L and L^{-1} are the identity */
	for (i = 0; i < nf; i++) {
	    Memcpy(REAL(VECTOR_ELT(bVar, i)), DIsqrt + Gp[i],
		   Gp[i+1] - Gp[i]);
	}
    } else {
	int *Lp = INTEGER(GET_SLOT(x, Matrix_LpSym)),
	    *Li = INTEGER(GET_SLOT(x, Matrix_LiSym));
	    
	double one = 1.0, zero = 0.,
	    *Lx = REAL(GET_SLOT(x, Matrix_LxSym));
	
	for (i = 0; i < nf; i++) {
	    int j, jj, k, kk, nci = nc[i], nr, p, p2, pj, pp,
		m = maxod + 1,
		*ind = Calloc(m, int), G1 = Gp[i], G2 = Gp[i+1];
	    double
		*tmp = Calloc(m * nci, double),
		*bVi = REAL(VECTOR_ELT(bVar, i));
	    
				/* initialize bVi to zero */
	    memset(bVi, 0, sizeof(double) * (G2 - G1) * nci);

	    for (j = G1; j < G2; j += nci) {
		kk = 0;		/* ind gets indices of non-zeros */
		jj = j;		/* in this block of columns */
		while (jj >= 0) {
		    ind[kk++] = jj;
		    jj = Parent[jj];
		}
		nr = kk;       /* number of non-zeros in this block */
		while (kk < m) ind[kk++] = nzc; /* placeholders */

		for (k = 0; k < nci; k++) {
		    double *ccol = tmp + k * nr;
		    
		    for (kk = 0; kk < nr; kk++) ccol[kk] = 0.;
		    ccol[k] = 1.; /* initialize from unit diagonal */
		    for (jj = j + k; jj >= 0; jj = Parent[jj]) {
			p2 = Lp[jj+1];
			pp = pj = ldl_update_ind(jj, 0, ind);
			for (p = Lp[jj]; p < p2; p++) {
			    pp = ldl_update_ind(Li[p], pp, ind);
			    ccol[pp] -= Lx[p] * ccol[pj];
			}
		    }
		}
				
		for (kk = 0; kk < nr; kk++) { /* scale rows */
		    for (k = 0; k < nci; k++) {
			tmp[k * nr + kk] *= DIsqrt[ind[kk]];
		    }
		}
		F77_CALL(dsyrk)("L", "T", &nci, &nr, &one, tmp, &nr,
				&zero, bVi + (j - G1)*nci, &nci);
		F77_CALL(dpotrf)("L", &nci, bVi + (j - G1)*nci,
				 &nci, &jj);
		if (jj)		/* should never happen */
		    error(
			"Rank deficient variance matrix at group %d, level %d, error code %d",
			i + 1, j + 1, jj);
	    }
	    Free(tmp); Free(ind);
	}
    }
    return R_NilValue;
}

/** 
 * If necessary, factor Z'Z+Omega, ZtX, and XtX then, if necessary,
 * form RZX, RXX, and bVar for the inverse of the Cholesky factor.
 * 
 * @param x pointer to an ssclme object
 * 
 * @return NULL (x is updated in place)
 */
SEXP ssclme_invert(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    if (!status[0]) ssclme_factor(x);
    if (!R_FINITE(REAL(GET_SLOT(x, Matrix_devianceSym))[0]))
	error("Unable to invert singular factor of downdated X'X");
    if (!status[1]) {
	SEXP
	    RZXsl = GET_SLOT(x, Matrix_RZXSym);
	int
	    *dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	    *Li = INTEGER(GET_SLOT(x, Matrix_LiSym)),
	    *Lp = INTEGER(GET_SLOT(x, Matrix_LpSym)),
	    i,
	    n = dims[0],
	    pp1 = dims[1];
	double
	    *DIsqrt = REAL(GET_SLOT(x, Matrix_DIsqrtSym)),
	    *Lx = REAL(GET_SLOT(x, Matrix_LxSym)),
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(RZXsl),
	    one = 1.;

	F77_CALL(dtrtri)("U", "N", &pp1, RXX, &pp1, &i);
	if (i)
	    error("DTRTRI returned error code %d", i);
	F77_CALL(dtrmm)("R", "U", "N", "N", &n, &pp1, &one,
			RXX, &pp1, RZX, &n);
	for (i = 0; i < pp1; i++) {
	    int j; double *RZXi = RZX + i * n;
	    for (j = 0; j < n; j++) RZXi[j] *= DIsqrt[j];
	    ldl_ltsolve(n, RZXi, Lp, Li, Lx);
	}
	ldl_inverse(x);
	status[1] = 1;
    }
    return R_NilValue;
}

/** 
 * Create and insert initial values for Omega_i.
 * 
 * @param x pointer to an ssclme object
 * 
 * @return NULL
 */
SEXP ssclme_initial(SEXP x)
{
    SEXP Gpsl = GET_SLOT(x, Matrix_GpSym),
	Omg = GET_SLOT(x, Matrix_OmegaSym);
    int *Ai = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(x, Matrix_pSym)),
	*Gp = INTEGER(Gpsl),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	i, nf = length(Gpsl) - 1;
	
    double *Ax = REAL(GET_SLOT(x, Matrix_xSym));

    for (i = 0; i < nf; i++) {
	int
	    Gpi = Gp[i],
	    j, k,
	    nci = nc[i],
	    ncip1 = nci + 1,
	    p2 = Gp[i+1];
	double
	    mi = 0.375 * ((double) nci)/((double) (p2 - Gpi)),
	    *mm = REAL(VECTOR_ELT(Omg, i));

	memset((void *) mm, 0, sizeof(double) * nci * nci);
	for (j = Gpi; j < p2; j += nci) {
	    for (k = 0; k < nci; k++) {
		int jk = j+k, jj = Ap[jk+1] - 1;
		if (Ai[jj] != jk) error("malformed ZtZ structure");
		mm[k * ncip1] += Ax[jj] * mi;
	    }
	}
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}

/** 
 * Extract the conditional estimates of the fixed effects
 * 
 * @param x Pointer to an ssclme object
 * 
 * @return a numeric vector containing the conditional estimates of
 * the fixed effects
 */
SEXP ssclme_fixef(SEXP x)
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym);
    int pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1];
    int j, p = pp1 - 1;
    SEXP val = PROTECT(allocVector(REALSXP, p));
    double
	*RXX = REAL(RXXsl),
	*beta = REAL(val),
	nryyinv;		/* negative ryy-inverse */

    ssclme_invert(x);
    Memcpy(beta, RXX + p * pp1, p);
    nryyinv = -RXX[pp1*pp1 - 1];
    for (j = 0; j < p; j++) beta[j] /= nryyinv;
    UNPROTECT(1);
    return val;
}

/** 
 * Extract the conditional modes of the random effects.
 * 
 * @param x Pointer to an ssclme object
 * 
 * @return a vector containing the conditional modes of the random effects
 */
SEXP ssclme_ranef(SEXP x)
{
    SEXP RZXsl = GET_SLOT(x, Matrix_RZXSym),
	GpSl = GET_SLOT(x, Matrix_GpSym);
    int *dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	*Gp = INTEGER(GpSl),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, j,
	n = dims[0],
	nf = length(GpSl) - 1,
	pp1 = dims[1],
	p = pp1 - 1;
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double
	*b = REAL(RZXsl) + n * p,
	ryyinv;		/* ryy-inverse */

    ssclme_invert(x);
    ryyinv = REAL(GET_SLOT(x, Matrix_RXXSym))[pp1*pp1 - 1];
    for (i = 0; i < nf; i++) {
	int nci = nc[i], Mi = Gp[i+1] - Gp[i];
	double *mm;
	
	SET_VECTOR_ELT(val, i, allocMatrix(REALSXP, nci, Mi/nci));
	mm = Memcpy(REAL(VECTOR_ELT(val, i)), b, Mi);
	b += Mi;
	for (j = 0; j < Mi; j++) mm[j] /= ryyinv;
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Extract the ML or REML conditional estimate of sigma
 * 
 * @param x pointer to an ssclme object
 * @param REML logical scalar - TRUE if REML estimates are requested
 * 
 * @return numeric scalar 
 */
SEXP ssclme_sigma(SEXP x, SEXP REML)
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym);
    int pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1],
	nobs = INTEGER(GET_SLOT(x, Matrix_ncSym))[
	    length(GET_SLOT(x, Matrix_GpSym))];
   
    ssclme_invert(x);
    return ScalarReal(1./(REAL(RXXsl)[pp1*pp1 - 1] *
			  sqrt((double)(asLogical(REML) ?
					nobs + 1 - pp1 : nobs))));
}

/** 
 * Calculate the length of the parameter vector, which is called coef
 * for historical reasons.
 * 
 * @param nf number of factors
 * @param nc number of columns in the model matrices for each factor
 * 
 * @return total length of the coefficient vector
 */
static
int coef_length(int nf, const int nc[])
{
    int i, ans = 0;
    for (i = 0; i < nf; i++) ans += (nc[i] * (nc[i] + 1))/2;
    return ans;
}

/** 
 * Extract the upper triangles of the Omega matrices.  These aren't
 * "coefficients" but the extractor is called coef for historical
 * reasons.  Within each group these values are in the order of the
 * diagonal entries first then the strict upper triangle in row
 * order.
 * 
 * @param x pointer to an ssclme object
 * 
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP ssclme_coef(SEXP x)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), vind;
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double *vv = REAL(val);

    vind = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	if (nci == 1) {
	    vv[vind++] = REAL(VECTOR_ELT(Omega, i))[0];
	} else {
	    int j, k, odind = vind + nci, ncip1 = nci + 1;
	    double *omgi = REAL(VECTOR_ELT(Omega, i));
	    
	    for (j = 0; j < nci; j++) {
		vv[vind++] = omgi[j * ncip1];
		for (k = j + 1; k < nci; k++) {
		    vv[odind++] = omgi[k*nci + j];
		}
	    }
	    vind = odind;
	}
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Extract the unconstrained parameters that determine the
 * Omega matrices. (Called coef for historical reasons.)  The
 * unconstrained parameters are derived from the LDL' decomposition of
 * Omega_i.  The first nc[i] entries in each group are the diagonals
 * of log(D) followed by the strict lower triangle of L in column
 * order.
 * 
 * @param x pointer to an ssclme object
 * 
 * @return numeric vector of unconstrained parameters that determine the
 * Omega matrices
 */
SEXP ssclme_coefUnc(SEXP x)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), vind;
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double *vv = REAL(val);

    vind = 0;
    for (i = 0; i < nf; i++) {
	int  nci = nc[i];
	if (nci == 1) {
	    vv[vind++] = log(REAL(VECTOR_ELT(Omega, i))[0]);
	} else {
	    int j, k, ncip1 = nci + 1, ncisq = nci * nci;
	    double *tmp = Memcpy(Calloc(ncisq, double), 
				 REAL(VECTOR_ELT(Omega, i)), ncisq);
	    F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
	    if (j)		/* should never happen */
		error("DPOTRF returned error code %d on Omega[[%d]]",
		      j, i+1);
	    for (j = 0; j < nci; j++) {
		double diagj = tmp[j * ncip1];
		vv[vind++] = 2. * log(diagj);
		for (k = j + 1; k < nci; k++) {
		    tmp[j + k * nci] /= diagj;
		}
	    }
	    for (j = 0; j < nci; j++) {
		for (k = j + 1; k < nci; k++) {
		    vv[vind++] = tmp[j + k * nci];
		}
	    }
	    Free(tmp);
	}
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Assign the Omega matrices from the unconstrained parameterization.
 * 
 * @param x pointer to an ssclme object
 * @param coef pointer to an numeric vector of appropriate length
 * 
 * @return R_NilValue
 */
SEXP ssclme_coefGetsUnc(SEXP x, SEXP coef)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	cind, i, nf = length(Omega),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    double *cc = REAL(coef);

    if (length(coef) != coef_length(nf, nc) || !isReal(coef))
	error("coef must be a numeric vector of length %d",
	      coef_length(nf, nc));
    cind = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	if (nci == 1) {
	    REAL(VECTOR_ELT(Omega, i))[0] = exp(cc[cind++]);
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		j, k,
		ncip1 = nci + 1,
		ncisq = nci * nci;
	    double
		*omgi = REAL(VECTOR_ELT(Omega, i)),
		*tmp = Calloc(ncisq, double),
		diagj, one = 1., zero = 0.;

	    memset(omgi, 0, sizeof(double) * ncisq);
	    for (j = 0; j < nci; j++) {
		tmp[j * ncip1] = diagj = exp(cc[cind++]/2.);
		for (k = j + 1; k < nci; k++) {
		    tmp[k*nci + j] = cc[odind++] * diagj;
		}
	    }
	    F77_CALL(dsyrk)("U", "T", &nci, &nci, &one,
			    tmp, &nci, &zero, omgi, &nci);
	    Free(tmp);
	    cind = odind;
	}
    }
    status[0] = status[1] = 0;
    return x;
}

/** 
 * Assign the upper triangles of the Omega matrices.
 * (Called coef for historical reasons.)
 * 
 * @param x pointer to an ssclme object
 * @param coef pointer to an numeric vector of appropriate length
 * 
 * @return R_NilValue
 */
SEXP ssclme_coefGets(SEXP x, SEXP coef)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	cind, i, nf = length(Omega),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    double *cc = REAL(coef);

    if (length(coef) != coef_length(nf, nc) || !isReal(coef))
	error("coef must be a numeric vector of length %d",
	      coef_length(nf, nc));
    cind = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	if (nci == 1) {
	    REAL(VECTOR_ELT(Omega, i))[0] = cc[cind++];
	} else {
	    int j, k, odind = cind + nci, ncip1 = nci + 1;
	    double *omgi = REAL(VECTOR_ELT(Omega, i));
	
	    for (j = 0; j < nci; j++) {
		omgi[j * ncip1] = cc[cind++];
		for (k = j + 1; k < nci; k++) {
		    omgi[k*nci + j] = cc[odind++];
		}
	    }
	    cind = odind;
	}
    }
    status[0] = status[1] = 0;
    return x;
}

/** 
 * Perform a number of ECME steps for the REML or ML criterion.
 * 
 * @param x pointer to an ssclme object
 * @param nsteps pointer to an integer scalar giving the number of ECME steps to perform
 * @param REMLp pointer to a logical scalar indicating if REML is to be used
 * @param verb pointer to a logical scalar indicating verbose mode
 * 
 * @return NULL
 */
SEXP ssclme_EMsteps(SEXP x, SEXP nsteps, SEXP REMLp, SEXP verb)
{
    SEXP
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	RZXsl = GET_SLOT(x, Matrix_RZXSym),
	ncsl = GET_SLOT(x, Matrix_ncSym),
	bVar = GET_SLOT(x, Matrix_bVarSym);
    int
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	*nc = INTEGER(ncsl),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	REML = asLogical(REMLp),
	i, info, iter,
	n = dims[0],
	nEM = asInteger(nsteps),
	nf = length(ncsl) - 2,
	nobs = nc[nf + 1],
	p,
	pp1 = dims[1],
	verbose = asLogical(verb);
    double
	*RZX = REAL(RZXsl),
	*dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	*b,
        alpha,
	one = 1.,
	zero = 0.;

    p = pp1 - 1;
    b = RZX + p * n;
    if (verbose) {
	SEXP coef = PROTECT(ssclme_coef(x));
	int lc = length(coef); double *cc = REAL(coef);

	ssclme_factor(x);
	Rprintf("  EM iterations\n");
	Rprintf("%3d %.3f", 0, dev[REML ? 1 : 0]);
	for (i = 0; i < lc; i++) Rprintf(" %#8g", cc[i]);
	Rprintf("\n");
	UNPROTECT(1);
    }
    for (iter = 0; iter < nEM; iter++) {
	ssclme_invert(x);
	for (i = 0; i < nf; i++) {
	    int ki = Gp[i+1] - Gp[i],
		nci = nc[i],
		mi = ki/nci;
	    double
	        *vali = REAL(VECTOR_ELT(Omega, i));
	    
	    alpha = ((double)(REML?(nobs-p):nobs))/((double)mi);
	    F77_CALL(dsyrk)("U", "N", &nci, &mi,
			    &alpha, b + Gp[i], &nci,
			    &zero, vali, &nci);
	    alpha = 1./((double) mi);
	    F77_CALL(dsyrk)("U", "N", &nci, &ki,
			    &alpha, REAL(VECTOR_ELT(bVar, i)), &nci,
			    &one, vali, &nci);
	    if (REML) {
		int j;
		for (j = 0; j < p; j++) { 
		    F77_CALL(dsyrk)("U", "N", &nci, &mi,
				&alpha, RZX + Gp[i] + j*n, &nci,
				&one, vali, &nci);
		}
	    }
	    F77_CALL(dpotrf)("U", &nci, vali, &nci, &info);
	    if (info)
		error("DPOTRF returned error code %d in Omega[[%d]] update",
		      info, i + 1);
	    F77_CALL(dpotri)("U", &nci, vali, &nci, &info);
	    if (info)
		error("DPOTRI returned error code %d in Omega[[%d]] update",
		      info, i + 1);
	}
	status[0] = status[1] = 0;
	if (verbose) {
	    SEXP coef = PROTECT(ssclme_coef(x));
	    int lc = length(coef); double *cc = REAL(coef);

	    ssclme_factor(x);
	    Rprintf("%3d %.3f", iter + 1, dev[REML ? 1 : 0]);
	    for (i = 0; i < lc; i++) Rprintf(" %#8g", cc[i]);
	    Rprintf("\n");
	    UNPROTECT(1);
	}
    }
    ssclme_factor(x);
    return R_NilValue;
}

/** 
 * Return the gradient of the ML or REML deviance.
 * 
 * @param x pointer to an ssclme object
 * @param REMLp pointer to a logical scalar indicating if REML is to be used
 * @param Uncp pointer to a logical scalar indicating if the unconstrained parameterization is to be used
 * 
 * @return pointer to a numeric vector of the gradient.
 */
SEXP ssclme_gradient(SEXP x, SEXP REMLp, SEXP Uncp)
{
    SEXP
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	RZXsl = GET_SLOT(x, Matrix_RZXSym),
	bVar = GET_SLOT(x, Matrix_bVarSym);
    int
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)), 
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	REML = asLogical(REMLp),
	cind, i, n = dims[0],
	nf = length(Omega),
	nobs, p, pp1 = dims[1],
	uncst = asLogical(Uncp);
    double
	*RZX = REAL(RZXsl),
	*b,
        alpha,
	one = 1.;
    SEXP ans = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));

    ssclme_factor(x);
    if (!R_FINITE(REAL(GET_SLOT(x, Matrix_devianceSym))[0])) {
	int ncoef = coef_length(nf, nc);
	for (i = 0; i < ncoef; i++) REAL(ans)[i] = NA_REAL;
	UNPROTECT(1);
	return ans;
    }
    nobs = nc[nf + 1];
    p = pp1 - 1;
    b = RZX + p * n;
    ssclme_invert(x);
    cind = 0;
    for (i = 0; i < nf; i++) {
	int j, ki = Gp[i+1] - Gp[i],
	    nci = nc[i], ncip1 = nci + 1, ncisq = nci * nci,
	    mi = ki/nci;
	double
	    *chol = Memcpy(Calloc(ncisq, double),
			   REAL(VECTOR_ELT(Omega, i)), ncisq),
	    *tmp = Calloc(ncisq, double);
	
	    
	F77_CALL(dpotrf)("U", &nci, chol, &nci, &j);
	if (j)
	    error("DPOTRF gave error code %d on Omega[[%d]]", j, i + 1);
	Memcpy(tmp, chol, ncisq);
	F77_CALL(dpotri)("U", &nci, tmp, &nci, &j);
	if (j)
	    error("DPOTRI gave error code %d on Omega[[%d]]", j, i + 1);
	alpha = (double) -mi;
	F77_CALL(dsyrk)("U", "N", &nci, &ki,
			&one, REAL(VECTOR_ELT(bVar, i)), &nci,
			&alpha, tmp, &nci);
	alpha = (double)(REML ? (nobs-p) : nobs);
	F77_CALL(dsyrk)("U", "N", &nci, &mi,
			&alpha, b + Gp[i], &nci,
			&one, tmp, &nci);
	if (REML) {
	    for (j = 0; j < p; j++) { 
		F77_CALL(dsyrk)("U", "N", &nci, &mi,
				&one, RZX + Gp[i] + j*n, &nci,
				&one, tmp, &nci);
	    }
	}
	if (nci == 1) {
	    REAL(ans)[cind++] = *tmp *
		(uncst ? *REAL(VECTOR_ELT(Omega, i)) : 1.);
	} else {
	    int k, odind = cind + nci;
	    if (uncst) {
		int ione = 1, kk;
		double *rr = Calloc(nci, double); /* j'th row of R, the Cholesky factor */
		nlme_symmetrize(tmp, nci);
		for (j = 0; j < nci; j++, cind++) {
		    for (k = 0; k < j; k++) rr[k] = 0.;
		    for (k = j; k < nci; k++) rr[k] = chol[j + k*nci];
		    REAL(ans)[cind] = 0.;
		    for (k = j; k < nci; k++) {
			for (kk = j; kk < nci; kk++) {
			    REAL(ans)[cind] += rr[k] * rr[kk] *
				tmp[kk * nci + k];
			}
		    }
		    for (k = j + 1; k < nci; k++) {
			REAL(ans)[odind++] = 2. * rr[j] *
			    F77_CALL(ddot)(&nci, rr, &ione, tmp + k*nci, &ione);
		    }			
		}
		Free(rr);
	    } else {
		for (j = 0; j < nci; j++) {
		    REAL(ans)[cind++] = tmp[j * ncip1];
		    for (k = j + 1; k < nci; k++) {
			REAL(ans)[odind++] = tmp[k*nci + j] * 2.;
		    }
		}
	    }
	    cind = odind;
	}
	Free(tmp); Free(chol);
    }
    UNPROTECT(1);
    return ans;
}

/** 
 * Return the Hessian of the ML or REML deviance.  This is a
 * placeholder until I work out the evaluation of the analytic
 * Hessian, which probably will involve several helper functions.
 * 
 * @param x pointer to an ssclme object
 * @param REMLp pointer to a logical scalar indicating if REML is to be used
 * @param Uncp pointer to a logical scalar indicating if the
 * unconstrained parameterization is to be used
 * 
 * @return pointer to an approximate Hessian matrix
 */
SEXP ssclme_Hessian(SEXP x, SEXP REMLp, SEXP Uncp)
{
    int j, ncoef = coef_length(length(GET_SLOT(x, Matrix_OmegaSym)),
			       INTEGER(GET_SLOT(x, Matrix_ncSym))),
	unc = asLogical(Uncp);
    SEXP ans = PROTECT(allocMatrix(REALSXP, ncoef, ncoef)),
	base = PROTECT(unc ? ssclme_coefUnc(x) : ssclme_coef(x)),
	current = PROTECT(duplicate(base)),
	gradient;

    for (j = 0; j < ncoef; j++) {
	double delta = (REAL(base)[j] ? 1.e-7 * REAL(base)[j] : 1.e-7);
	int i;

	for (i = 0; i < ncoef; i++) REAL(current)[i] = REAL(base)[i];
	REAL(current)[j] += delta/2.;
	if (unc) {
	    ssclme_coefGetsUnc(x, current);
	} else {
	    ssclme_coefGets(x, current);
	}
	PROTECT(gradient = ssclme_gradient(x, REMLp, Uncp));
	for (i = 0; i < ncoef; i++) REAL(ans)[j * ncoef + i] = REAL(gradient)[i];
	UNPROTECT(1);
	REAL(current)[j] -= delta;
	if (unc) {
	    ssclme_coefGetsUnc(x, current);
	} else {
	    ssclme_coefGets(x, current);
	}
	PROTECT(gradient = ssclme_gradient(x, REMLp, Uncp));
	for (i = 0; i < ncoef; i++)
	    REAL(ans)[j * ncoef + i] = (REAL(ans)[j * ncoef + i] - REAL(gradient)[i])/
		delta;
	UNPROTECT(1);
	/* symmetrize */
	for (i = 0; i < j; i++) {
	    REAL(ans)[j * ncoef + i] = REAL(ans)[i * ncoef + j] =
		(REAL(ans)[j * ncoef + i] + REAL(ans)[i * ncoef + j])/2.;
	}
    }
    UNPROTECT(3);
    return ans;
}

/** 
 * Calculate and return the fitted values.
 * 
 * @param x pointer to an ssclme object
 * @param facs list of grouping factors
 * @param mmats list of model matrices
 * @param useRf pointer to a logical scalar indicating if the random effects should be used
 * 
 * @return pointer to a numeric array of fitted values
 */
SEXP ssclme_fitted(SEXP x, SEXP facs, SEXP mmats, SEXP useRf)
{
    SEXP val, b;
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ione = 1, nf = length(facs), nobs, p;
    double *vv, one = 1.0, zero = 0.0;

    if (nf < 1)
	error("null factor list passed to ssclme_fitted");
    nobs = length(VECTOR_ELT(facs, 0));
    val = PROTECT(allocVector(REALSXP, nobs));
    vv = REAL(val);
    p = nc[nf] - 1;
    if (p > 0) {
	F77_CALL(dgemm)("N", "N", &nobs, &ione, &p, &one,
			REAL(VECTOR_ELT(mmats, nf)), &nobs,
			REAL(PROTECT(ssclme_fixef(x))), &p,
			&zero, vv, &nobs);
	UNPROTECT(1);
    } else {
	memset(vv, 0, sizeof(double) * nobs);
    }
    if (asLogical(useRf)) {
	b = PROTECT(ssclme_ranef(x));
	for (i = 0; i < nf; i++) {
	    int *ff = INTEGER(VECTOR_ELT(facs, i)), j, nci = nc[i];
	    double *bb = REAL(VECTOR_ELT(b, i)),
		*mm = REAL(VECTOR_ELT(mmats, i));
	    for (j = 0; j < nobs; ) {
		int nn = 1, lev = ff[j];
		/* check for adjacent rows with same factor level */
		while (ff[j + nn] == lev) nn++; 
		F77_CALL(dgemm)("N", "N", &nn, &ione, &nci,
				&one, mm + j, &nobs,
				bb + (lev - 1) * nci, &nci,
				&one, vv + j, &nobs);
		j += nn;
	    }
	}
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Return the unscaled variances
 * 
 * @param x pointer to an ssclme object
 * 
 * @return a list similar to the Omega list with the unscaled variances
 */
SEXP ssclme_variances(SEXP x)
{
    SEXP Omg = PROTECT(duplicate(GET_SLOT(x, Matrix_OmegaSym)));
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omg);
    
    for (i = 0; i < nf; i++) {
	double *mm = REAL(VECTOR_ELT(Omg, i));
	int j, nci = nc[i];

	F77_CALL(dpotrf)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error("DPOTRF returned error code %d on Omega[%d]",
		  j, i + 1);
	F77_CALL(dpotri)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error("DTRTRI returned error code %d on Omega[%d]",
		  j, i + 1);
	nlme_symmetrize(mm, nci);
    }
    UNPROTECT(1);
    return Omg;
}

/** 
 * Copy an ssclme object collapsing the fixed effects slots to the response only.
 * 
 * @param x pointer to an ssclme object
 * 
 * @return a duplicate of x with the fixed effects slots collapsed to the response only
 */
SEXP ssclme_collapse(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("ssclme"))),
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	Dim = GET_SLOT(x, Matrix_DimSym);
    int nf = length(Omega), nz = INTEGER(Dim)[1];

    slot_dup(ans, x, Matrix_DSym);
    slot_dup(ans, x, Matrix_DIsqrtSym);
    slot_dup(ans, x, Matrix_DimSym);
    slot_dup(ans, x, Matrix_GpSym);
    slot_dup(ans, x, Matrix_LiSym);
    slot_dup(ans, x, Matrix_LpSym);
    slot_dup(ans, x, Matrix_LxSym);
    slot_dup(ans, x, Matrix_OmegaSym);
    slot_dup(ans, x, Matrix_ParentSym);
    slot_dup(ans, x, Matrix_bVarSym);
    slot_dup(ans, x, Matrix_devianceSym);
    slot_dup(ans, x, Matrix_devCompSym);
    slot_dup(ans, x, Matrix_iSym);
    slot_dup(ans, x, Matrix_ncSym);
    slot_dup(ans, x, Matrix_statusSym);
    slot_dup(ans, x, Matrix_pSym);
    slot_dup(ans, x, Matrix_xSym);
    INTEGER(GET_SLOT(ans, Matrix_ncSym))[nf] = 1;
    SET_SLOT(ans, Matrix_XtXSym, allocMatrix(REALSXP, 1, 1));
    REAL(GET_SLOT(ans, Matrix_XtXSym))[0] = NA_REAL;
    SET_SLOT(ans, Matrix_RXXSym, allocMatrix(REALSXP, 1, 1));
    REAL(GET_SLOT(ans, Matrix_RXXSym))[0] = NA_REAL;
    SET_SLOT(ans, Matrix_ZtXSym, allocMatrix(REALSXP, nz, 1));
    SET_SLOT(ans, Matrix_RZXSym, allocMatrix(REALSXP, nz, 1));
    LOGICAL(GET_SLOT(ans, Matrix_statusSym))[0] = 0;
    UNPROTECT(1);
    return ans;
}


/** 
 * Create an lme object from its components.  This is not done by
 * new("lme", ...) at the R level because of the possibility of
 * causing the copying of very large objects.
 * 
 * @param call Pointer to the original call
 * @param facs pointer to the list of grouping factors
 * @param x pointer to the model matrices (may be of length zero)
 * @param model pointer to the model frame
 * @param REML pointer to a logical scalar indicating if REML is used
 * @param rep pointer to the converged ssclme object
 * @param fitted pointer to the fitted values
 * @param residuals pointer to the residuals
 * 
 * @return an lme object
 */
SEXP ssclme_to_lme(SEXP call, SEXP facs, SEXP x, SEXP model, SEXP REML,
		   SEXP rep, SEXP fitted, SEXP residuals)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lme")));

    SET_SLOT(ans, install("call"), call);
    SET_SLOT(ans, install("facs"), facs);
    SET_SLOT(ans, Matrix_xSym, x);
    SET_SLOT(ans, install("model"), model);
    SET_SLOT(ans, install("REML"), REML);
    SET_SLOT(ans, install("rep"), rep);
    SET_SLOT(ans, install("fitted"), fitted);
    SET_SLOT(ans, install("residuals"), residuals);
    UNPROTECT(1);
    return ans;
}
