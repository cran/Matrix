#include "ssclme.h"

/** 
 * Check for a nested series of grouping factors in the sparse,
 *  symmetric representation of the pairwise cross-tabulations.
 * 
 * @param n size of pairwise cross-tabulation matrix
 * @param nf number of groups of columns in pairwise cross-tabulation
 * @param upper non-zero if the upper triangle is stored
 * @param Ap array of pointers to columns
 * @param Ai row indices
 * @param Gp array of pointers to groups
 * 
 * @return 0 for non-nested groups, 1 for nested groups
 */
static
int ctab_isNested(int n, int nf, int upper,
		  const int Ap[], const int Ai[], const int Gp[])
{
    if (nf > 1) {  /* single factor always nested */
	int  i;
	if (upper) {
	    int *nnz = (int *) R_alloc(n, sizeof(int)), nz = Ap[n];
				/* count number of nonzeros in each row */
	    for (i = 0; i < n; i++) nnz[i] = 0;
	    for (i = 0; i < nz; i++) nnz[Ai[i]]++;
	    for (i = 0; i < nf; i++) {
		int j, p2 = Gp[i+1], target = nf - i;
		for (j = Gp[i]; j < p2; j++) {
		    if (nnz[j] != target) return 0;
		}
	    }
	} else {		/* lower triangle - the easy case */
	    for (i = 0; i < nf; i++) {
		int j, p2 = Gp[i+1], target = nf - i;
		for (j = Gp[i]; j < p2; j++) {
		    if ((Ap[j+1] - Ap[j]) != target)
			return 0;
		}
	    }
	}
    }
    return 1;
}

/** 
 * Determine if a fill-reducing permutation is needed for the pairwise
 * cross-tabulation matrix.  If so, determine such a permutation
 * (using Metis) then separate the groups.
 * 
 * @param ctab pointer to a pairwise cross-tabulation object
 * 
 * @return pointer to an integer R vector.
 */
static
SEXP ctab_permute(SEXP ctab)
{
    SEXP val, GpSl = GET_SLOT(ctab, Matrix_GpSym);
    int *Ai = INTEGER(GET_SLOT(ctab, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(ctab, Matrix_pSym)),
	*Gp = INTEGER(GpSl),
	*perm,
	*work,
	i,
	j,
	n = INTEGER(GET_SLOT(ctab, Matrix_DimSym))[1],
	nf = length(GpSl) - 1,
	nz = Ap[n],		/* number of non-zeros */
	pos;

    if (ctab_isNested(n, nf, 1, Ap, Ai, Gp))
	return allocVector(INTSXP, 0);
    val =  allocVector(INTSXP, n);
    perm = INTEGER(val);
    work = (int *) R_alloc(n, sizeof(int));
    ssc_metis_order(n, nz, Ap, Ai, work, perm);	/* perm gets inverse perm */
    /* work now contains desired permutation but with groups scrambled */

    /* copy work into perm preserving the order of the groups */
    pos = 0;		/* position in new permutation */
    for (i = 0; i < nf; i++) {
	for (j = 0; j < n; j++) {
	    int jj = work[j];
	    if (Gp[i] <= jj && jj < Gp[i+1]) {
		perm[pos] = jj;
		pos++;
	    }
	}
    }
    return val;
}

static
void ssclme_copy_ctab(int nf, const int nc[], SEXP ctab, SEXP ssc)
{
    int *snc, i, copyonly = 1;

    for (i = 0; i < nf; i++) {
	if (nc[i] > 1) copyonly = 0;
    }
    if (copyonly) {
	SET_SLOT(ssc, Matrix_pSym, duplicate(GET_SLOT(ctab, Matrix_pSym)));
	SET_SLOT(ssc, Matrix_iSym, duplicate(GET_SLOT(ctab, Matrix_iSym)));
	SET_SLOT(ssc, Matrix_xSym, duplicate(GET_SLOT(ctab, Matrix_xSym)));
	SET_SLOT(ssc, Matrix_DimSym,
		 duplicate(GET_SLOT(ctab, Matrix_DimSym)));
	SET_SLOT(ssc, Matrix_GpSym, duplicate(GET_SLOT(ctab, Matrix_GpSym)));
    } else {
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
	SET_SLOT(ssc, Matrix_DimSym,
		 duplicate(GET_SLOT(ctab, Matrix_DimSym)));
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
		int ii = AiIn[j], mj = map[j], ncci = ncc[ii],
		    pos = ApOut[mj];
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
		    AiOut[ApOut[mjj + 1] - 1] = mjj; /* maybe mjj-1? */
		}
	    }
	}
	Free(map); Free(ncc);
    }
    SET_SLOT(ssc, Matrix_ncSym, allocVector(INTSXP, nf + 2));
    snc = INTEGER(GET_SLOT(ssc, Matrix_ncSym));
    for (i = 0; i <= nf; i++) {
	snc[i] = nc[i];
    }
}

static
void ssclme_fill_LIp(int n, const int Parent[], int LIp[])
{
    int *sz = Calloc(n, int), i;
    for (i = n - 1; i >= 0; i--) {
	sz[i] = (Parent[i] < 0) ? 0 : 1 + sz[Parent[i]];
    }
    LIp[0] = 0;
    for (i = 0; i < n; i++) LIp[i+1] = LIp[i] + sz[i];
    Free(sz);
}

static
void ssclme_fill_LIi(int n, const int Parent[], const int LIp[], int LIi[])
{
    int i;
    for (i = n; i > 0; i--) {
	int im1 = i - 1, Par = Parent[im1];
	if (Par >= 0) {
	    LIi[LIp[im1]] = Par;
	    Memcpy(LIi + LIp[im1] + 1, LIi + LIp[Par],
		   LIp[Par + 1] - LIp[Par]);
	}
    }
}

SEXP
ssclme_create(SEXP facs, SEXP ncv, SEXP threshold)
{
    SEXP ctab, nms, ssc, tmp,
	val = PROTECT(allocVector(VECSXP, 2)),
	dd = PROTECT(allocVector(INTSXP, 3));	/* dimensions of 3-D arrays */
    int *Ai, *Ap, *Gp, *LIp, *Lp, *Parent,
	*nc, Lnz, i, nf = length(facs), nzcol, pp1,
	*dims = INTEGER(dd);

    if (length(ncv) != (nf + 1))
	error("length of nc (%d) should be length of facs (%d) + 1",
	      length(ncv), nf);
    SET_VECTOR_ELT(val, 0, NEW_OBJECT(MAKE_CLASS("ssclme")));
    ssc = VECTOR_ELT(val, 0);
				/* Pairwise cross-tabulation */
    ctab = PROTECT(sscCrosstab(facs, ScalarLogical(1)));
    SET_VECTOR_ELT(val, 1, ctab_permute(ctab));
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
        /* Zero the symmetric matrices (for cosmetic reasons only). */
    memset(REAL(GET_SLOT(ssc, Matrix_XtXSym)), 0,
	   sizeof(double) * pp1 * pp1); 
    memset(REAL(GET_SLOT(ssc, Matrix_RXXSym)), 0,
	   sizeof(double) * pp1 * pp1);
    SET_SLOT(ssc, Matrix_LpSym, allocVector(INTSXP, nzcol + 1));
    Lp = INTEGER(GET_SLOT(ssc, Matrix_LpSym));
    SET_SLOT(ssc, Matrix_ParentSym, allocVector(INTSXP, nzcol));
    Parent = INTEGER(GET_SLOT(ssc, Matrix_ParentSym));
    SET_SLOT(ssc, Matrix_DSym, allocVector(REALSXP, nzcol));
    SET_SLOT(ssc, Matrix_DIsqrtSym, allocVector(REALSXP, nzcol));
    ldl_symbolic(nzcol, Ap, Ai, Lp, Parent,
		 (int *) R_alloc(nzcol, sizeof(int)), /* Lnz */
		 (int *) R_alloc(nzcol, sizeof(int)), /* Flag */
		 (int *) NULL, (int *) NULL); /* P & Pinv */
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
    SET_SLOT(ssc, Matrix_LIpSym, allocVector(INTSXP, nzcol + 1));
    LIp = INTEGER(GET_SLOT(ssc, Matrix_LIpSym));
    ssclme_fill_LIp(nzcol, Parent, LIp);
    if (asInteger(threshold) > (Lnz = LIp[nzcol])) {
	SET_SLOT(ssc, Matrix_LIiSym, allocVector(INTSXP, Lnz));
	ssclme_fill_LIi(nzcol, Parent, LIp,
			INTEGER(GET_SLOT(ssc, Matrix_LIiSym)));
    	SET_SLOT(ssc, Matrix_LIxSym, allocVector(REALSXP, Lnz));
	memset(REAL(GET_SLOT(ssc, Matrix_LIxSym)), 0,
	       sizeof(double) * Lnz);
    }
    UNPROTECT(2);
    return val;
}

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

SEXP
ssclme_update_mm(SEXP x, SEXP facs, SEXP mmats)
{
    SEXP bVar = GET_SLOT(x, Matrix_bVarSym);
    int
	*Ai = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(x, Matrix_pSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
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
		nck = nc[k];
	    double
		*Zk = Z[k];

	    for (i = 0; i < nobs; i++) {
		int ii, ind = -1, fpji = fpj[i] - 1,
		    row = Gp[j] + fpji * ncj,
		    fpki = fpk[i] - 1,
		    lastind = Apk[fpki + 1];
		for (ii = Apk[fpki]; ii < lastind; ii++) {
		    if (Ai[ii] == row) {
			ind = ii;
			break;
		    }
		}
		if (ind < 0) error("logic error in ssclme_update_mm");
		if (Ncj || nck > 1) {
				/* FIXME: run a loop to update */
		    error("code not yet written");
		} else {	/* update scalars directly */
		    Ax[ind] += Zj[fpji] * Zk[fpki];
		}
	    }
	}
    }
    Free(Z);
    ssclme_transfer_dimnames(x, facs, mmats);
    return R_NilValue;
}

SEXP ssclme_inflate_and_factor(SEXP lme)
{
    SEXP
	GpSlot = GET_SLOT(lme, Matrix_GpSym),
	Omega = GET_SLOT(lme, Matrix_OmegaSym);
    int n = INTEGER(GET_SLOT(lme, Matrix_DimSym))[1];
    int
	*Ai = INTEGER(GET_SLOT(lme, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(lme, Matrix_pSym)),
	*Flag = Calloc(n, int),
	*Gp = INTEGER(GpSlot),
	*Lnz = Calloc(n, int),
	*Pattern = Calloc(n, int),
	*nc = INTEGER(GET_SLOT(lme, Matrix_ncSym)),
	j,
	nf = length(GpSlot) - 1;
    double
	*D = REAL(GET_SLOT(lme, Matrix_DSym)),
	*DIsqrt = REAL(GET_SLOT(lme, Matrix_DIsqrtSym)),
	*Y = Calloc(n, double),
	*xcp = Calloc(Ap[n], double);

    Memcpy(xcp, REAL(GET_SLOT(lme, Matrix_xSym)), Ap[n]);
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
		    INTEGER(GET_SLOT(lme, Matrix_LpSym)),
		    INTEGER(GET_SLOT(lme, Matrix_ParentSym)),
		    Lnz, INTEGER(GET_SLOT(lme, Matrix_LiSym)),
		    REAL(GET_SLOT(lme, Matrix_LxSym)),
		    D, Y, Pattern, Flag,
		    (int *) NULL, (int *) NULL); /* P & Pinv */
    if (j != n)
	error("rank deficiency of ZtZ+W detected at column %d",
	      j + 1);
    for (j = 0; j < n; j++) DIsqrt[j] = 1./sqrt(D[j]);
    Free(Lnz); Free(Flag); Free(Pattern); Free(Y); Free(xcp);
    return R_NilValue;
}

SEXP ssclme_factor(SEXP lme)
{
    int *status = LOGICAL(GET_SLOT(lme, Matrix_statusSym));
    
    if (!status[0]) {
	SEXP
	    GpSlot = GET_SLOT(lme, Matrix_GpSym),
	    Omega = GET_SLOT(lme, Matrix_OmegaSym);
	int
	    *Gp = INTEGER(GpSlot),
	    *Li = INTEGER(GET_SLOT(lme, Matrix_LiSym)),
	    *Lp = INTEGER(GET_SLOT(lme, Matrix_LpSym)),
	    *nc = INTEGER(GET_SLOT(lme, Matrix_ncSym)),
	    i,
	    n = INTEGER(GET_SLOT(lme, Matrix_DimSym))[1],
	    nf = length(GpSlot) - 1,
	    nobs = nc[nf + 1],
	    nreml = nobs + 1 - nc[nf],
	    pp1 = nc[nf],
	    pp2 = pp1 + 1;
	double
	    *D = REAL(GET_SLOT(lme, Matrix_DSym)),
	    *DIsqrt = REAL(GET_SLOT(lme, Matrix_DIsqrtSym)),
	    *Lx = REAL(GET_SLOT(lme, Matrix_LxSym)),
	    *RXX = REAL(GET_SLOT(lme, Matrix_RXXSym)),
	    *RZX = REAL(GET_SLOT(lme, Matrix_RZXSym)),
	    *dcmp = REAL(getAttrib(lme, Matrix_devCompSym)),
	    *deviance = REAL(getAttrib(lme, Matrix_devianceSym)),
	    minus1 = -1.,
	    one = 1.;
	
	ssclme_inflate_and_factor(lme);
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
	Memcpy(RZX, REAL(GET_SLOT(lme, Matrix_ZtXSym)), n * pp1);
	for (i = 0; i < pp1; i++) {
	    int j;
	    double *RZXi = RZX + i * n;
	    ldl_lsolve(n, RZXi, Lp, Li, Lx);
	    for (j = 0; j < n; j++) RZXi[j] *= DIsqrt[j];
	}
				/* downdate and factor X'X */
	Memcpy(RXX, REAL(GET_SLOT(lme, Matrix_XtXSym)), pp1 * pp1);
	F77_CALL(dsyrk)("U", "T", &pp1, &n, &minus1,
			RZX, &n, &one, RXX, &pp1);
	F77_CALL(dpotrf)("U", &pp1, RXX, &pp1, &i);
	if (i)
	    error("DPOTRF returned error code %d", i);
				/* logdet of RXX */
	for (i = 0; i < (pp1 - 1); i++)
	    dcmp[2] += 2 * log(RXX[i*pp2]);
				/* logdet of Ryy */
	dcmp[3] = 2. * log(RXX[pp1 * pp1 - 1]);
	deviance[0] =		/* ML criterion */
	    dcmp[0] - dcmp[1] + nobs*(1+dcmp[3]+log(2*PI/nobs));
	deviance[1] = dcmp[0] - dcmp[1] + /* REML */
	    dcmp[2] + nreml * (1. + dcmp[3] + log(2. * PI/nreml));
	status[0] = 1;
	status[1] = 0;
    }
    return R_NilValue;
}

static
int ldl_update_ind(int probe, int start, const int ind[])
{
    while (ind[start] < probe) start++;
    if (ind[start] > probe) error("logic error in ldl_inverse");
    return start;
}

/** 
 * Create the inverse of L and update the diagonal blocks of the inverse
 * of LDL' (=Z'Z+W)
 * 
 * @param x pointer to an ssclme object
 *
 * @return R_NilValue (x is updated in place)

 */
SEXP ldl_inverse(SEXP x)
{
    SEXP
	Gpsl = GET_SLOT(x, Matrix_GpSym),
	LIisl = GET_SLOT(x, Matrix_LIiSym),
	LIpsl = GET_SLOT(x, Matrix_LIpSym),
	bVar = GET_SLOT(x, Matrix_bVarSym);
    int *Gp = INTEGER(Gpsl),
	*Li,
	*LIp = INTEGER(LIpsl), *Lp,
	i,
	nf = length(Gpsl) - 1,
	nzc = length(LIpsl) - 1;
    double
	*DIsqrt = REAL(GET_SLOT(x, Matrix_DIsqrtSym)),
	*Lx;
	
    ssclme_factor(x);
    if (LIp[nzc] == 0) {	/* L and LI are the identity */
	for (i = 0; i < nf; i++) {
	    Memcpy(REAL(VECTOR_ELT(bVar, i)), DIsqrt + Gp[i],
		   Gp[i+1] - Gp[i]);
	}
	return R_NilValue;
    }
    Lp = INTEGER(GET_SLOT(x, Matrix_LpSym));
    Li = INTEGER(GET_SLOT(x, Matrix_LiSym));
    Lx = REAL(GET_SLOT(x, Matrix_LxSym));
    if (length(LIisl) == LIp[nzc]) { /* LIi is filled */
	int *LIi = INTEGER(LIisl),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    j, jj, k, kk, p1, p2, pi1, pi2;

	double *LIx = REAL(GET_SLOT(x, Matrix_LIxSym)),
	    one = 1., zero = 0.;

	memset(LIx, 0, sizeof(double) * LIp[nzc]);
				/* calculate inverse */
	for (i = 0; i < nzc; i++) {
	    p1 = Lp[i]; p2 = Lp[i+1]; pi1 = LIp[i]; pi2 = LIp[i+1];
				/* initialize from unit diagonal term */
	    kk = pi1;
	    for (j = p1; j < p2; j++) {
		k = Li[j];
		while (LIi[kk] < k && kk < pi2) kk++;
		if (LIi[kk] != k) error("logic error in ldl_inverse");
		LIx[kk] = -Lx[j];
	    }
	    for (j = pi1; j < pi2; j++) {
		jj = LIi[j];
		p1 = Lp[jj]; p2 = Lp[jj+1];
		kk = j;
		for (jj = p1; jj < p2; jj++) {
		    k = Li[jj];
		    while (LIi[kk] < k && kk < pi2) kk++;
		    if (LIi[kk] != k) error("logic error in ldl_inverse");
		    LIx[kk] -= Lx[jj]*LIx[j];
		}
	    }
	}
	for (i = 0; i < nf; i++) { /* accumulate bVar */
	    int G1 = Gp[i], G2 = Gp[i+1], j, k, kk,
		nci = nc[i], nr, nr1, rr;
	    double *bVi = REAL(VECTOR_ELT(bVar, i)), *tmp;

	    nr = -1;
	    for (j = G1; j < G2; j += nci) {
		rr = 1 + LIp[j + 1] - LIp[j];
		if (rr > nr) nr = rr;
	    }
	    tmp = Calloc(nr * nci, double); /* scratch storage */
	    nr1 = nr + 1;
				/* initialize bVi to zero (cosmetic) */
	    memset(bVi, 0, sizeof(double) * (G2 - G1) * nci);
	    for (j = G1; j < G2; j += nci) {
		memset(tmp, 0, sizeof(double) * nr * nci);
		rr = 1 + LIp[j + 1] - LIp[j];
		for (k = 0; k < nci; k++) { /* copy columns */
		    tmp[k * nr1] = 1.; /* (unstored) diagonal elt  */
		    Memcpy(tmp + k*nr1 + 1, LIx + LIp[j + k], rr - k - 1);
		}
				/* scale the rows */
		tmp[0] = DIsqrt[j]; /* first row only has one non-zero */
		for (kk = 1; kk < rr; kk++) {
		    for (k = 0; k < nci; k++) {
			tmp[k * nr + kk] *= DIsqrt[LIi[LIp[j] + kk - 1]];
		    }
		}
		F77_CALL(dsyrk)("U", "T", &nci, &rr, &one, tmp, &nr,
				&zero, bVi + (j - G1) * nci, &nci);
		F77_CALL(dpotrf)("U", &nci, bVi + (j - G1) * nci,
				 &nci, &kk);
		if (kk)		/* should never happen */
		    error(
			"Rank deficient variance matrix at group %d, level %d",
			i + 1, j + 1);
	    }
	}
	return R_NilValue;
    }
    if (length(LIisl)) error("logic error in ssclme_ldl_inverse");
    else {		    /* LIi and LIx are too big and not used */
	int *counts = Calloc(nzc, int), info, maxod = -1;
	int *Parent = INTEGER(GET_SLOT(x, Matrix_ParentSym));
	int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym));
	double one = 1.0, zero = 0.;
				/* determine maximum # of off-diagonals */
	for (i = nzc - 1; i >= 0; i--) { /* in a column of L^{-1} */
	    counts[i] = (Parent[i] < 0) ? 0 : 1 + counts[Parent[i]];
	    if (counts[i] > maxod) maxod = counts[i];
	}
	Free(counts);

	for (i = 0; i < nf; i++) {
	    int j, jj, k, kk, nci = nc[i], nr, p, p2, pp,
		m = maxod + nci,
		*ind = Calloc(m, int);
	    double
		*tmp = Calloc(m * nci, double),
		*mpt = REAL(VECTOR_ELT(bVar, i));
	    
	    for (j = Gp[i]; j < Gp[i+1]; j += nci) {
		memset((void *) tmp, 0, sizeof(double) * m * nci);
		
		kk = 0;		/* ind holds indices of non-zeros */
		jj = j;		/* in this block of columns */
		while (jj >= 0) {
		    ind[kk++] = jj;
		    jj = Parent[jj];
		}
		nr = kk;       /* number of non-zeros in this block */
		while (kk < m) ind[kk++] = nzc; /* placeholders */

		for (k = 0; k < nci; k++) {
		    double *ccol = tmp + k * nr;
		    
		    ccol[k] = 1.;
		    kk = k;
		    for (jj = j + k; jj >= 0; jj = Parent[jj]) {
			p2 = Lp[jj+1];
			pp = kk;
			for (p = Lp[jj]; p < p2; p++) {
			    pp = ldl_update_ind(Li[p], pp, ind);
			    ccol[pp] -= Lx[p] * ccol[kk];
			}
		    }
		}
				/* scale rows */
		for (kk = 0; kk < nr; kk++) {
		    for (k = 0; k < nci; k++) {
			tmp[k * nr + kk] *= DIsqrt[ind[kk]];
		    }
		}
		F77_CALL(dsyrk)("U", "T", &nci, &nr, &one, tmp, &nr,
				&zero, mpt + (j - Gp[i])*nci, &nci);
		F77_CALL(dpotrf)("U", &nci, mpt + (j - Gp[i])*nci,
				 &nci, &info);
		if (info)	/* should never happen */
		    error(
			"Rank deficient variance matrix at group %d, level %d",
			i + 1, j + 1);
	    }
	    Free(tmp); Free(ind);
	}
    }
    return R_NilValue;
}

SEXP ssclme_invert(SEXP lme)
{
    int *status = LOGICAL(GET_SLOT(lme, Matrix_statusSym));
    if (!status[0]) ssclme_factor(lme);
    if (!status[1]) {
	SEXP
	    RZXsl = GET_SLOT(lme, Matrix_RZXSym);
	int
	    *dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	    *Li = INTEGER(GET_SLOT(lme, Matrix_LiSym)),
	    *Lp = INTEGER(GET_SLOT(lme, Matrix_LpSym)),
	    i,
	    n = dims[0],
	    pp1 = dims[1];
	double
	    *DIsqrt = REAL(GET_SLOT(lme, Matrix_DIsqrtSym)),
	    *Lx = REAL(GET_SLOT(lme, Matrix_LxSym)),
	    *RXX = REAL(GET_SLOT(lme, Matrix_RXXSym)),
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
	ldl_inverse(lme);
	status[1] = 1;
    }
    return R_NilValue;
}

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
 * FIXME: Add names
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
 * FIXME: Change the returned value to be a named list of matrices
 *        with dimnames.
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
	nobs = INTEGER(GET_SLOT(x, Matrix_ncSym))
	[length(GET_SLOT(x, Matrix_GpSym))];
   
    ssclme_invert(x);
    return ScalarReal(1./(REAL(RXXsl)[pp1*pp1 - 1] *
			  sqrt((double)(asLogical(REML) ?
					nobs + 1 - pp1 : nobs))));
}

static
int coef_length(int nf, const int nc[])
{
    int i, ans = 0;
    for (i = 0; i < nf; i++) ans += (nc[i] * (nc[i] + 1))/2;
    return ans;
}

/** 
 * Extract the upper triangles of the Omega matrices.
 * (These are not in any sense "coefficients" but the extractor is
 * called coef for historical reasons.)
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
	int j, k, nci = nc[i];
	double *omgi = REAL(VECTOR_ELT(Omega, i));
	for (j = 0; j < nci; j++) {
	    for (k = 0; k <= j; k++) {
		vv[vind++] = omgi[j*nci + k];
	    }
	}
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Extract the upper triangles of the Omega matrices in the unconstrained
 * parameterization.
 * (These are not in any sense "coefficients" but the extractor is
 * called coef for historical reasons.)
 * 
 * @param x pointer to an ssclme object
 * 
 * @return numeric vector of the values in the upper triangles of the
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
		error("DPOTRF returned error code %d", j);
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
 * Assign the upper triangles of the Omega matrices in the unconstrained
 * parameterization.
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
		diagj, one = 1.;
	    /* FIXEME: Change this to use a factor and dsyrk */
				/* LD in omgi and L' in tmp */
	    memset(omgi, 0, sizeof(double) * ncisq);
	    for (j = 0; j < nci; j++) {
		omgi[j * ncip1] = diagj = exp(cc[cind++]);
		for (k = j + 1; k < nci; k++) {
		    omgi[j*nci + k] = diagj * (tmp[k*nci + j] = cc[odind++]);
		}
	    }
	    F77_CALL(dtrmm)("R", "U", "N", "U", &nci, &nci, &one,
			    tmp, &nci, omgi, &nci);
	    Free(tmp);
	    cind = odind;
	}
    }
    status[0] = status[1] = 0;
    return x;
}

/** 
 * Assign the upper triangles of the Omega matrices.
 * (These are not in any sense "coefficients" but are
 * called coef for historical reasons.)
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
	int j, k, nci = nc[i];
	double *omgi = REAL(VECTOR_ELT(Omega, i));
	for (j = 0; j < nci; j++) {
	    for (k = 0; k <= j; k++) {
		omgi[j*nci + k] = cc[cind++];
	    }
	}
    }
    status[0] = status[1] = 0;
    return x;
}

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
    if (verbose) Rprintf("  EM iterations\n");
    for (iter = 0; iter <= nEM; iter++) {
	ssclme_invert(x);
	if (verbose) {
	    SEXP coef = PROTECT(ssclme_coef(x));
	    int lc = length(coef); double *cc = REAL(coef);
	    Rprintf("%3d %.3f", iter, dev[REML ? 1 : 0]);
	    for (i = 0; i < lc; i++) Rprintf(" %#8g", cc[i]);
	    Rprintf("\n");
	    UNPROTECT(1);
	}
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
		int mp = mi * p;
		F77_CALL(dsyrk)("U", "N", &nci, &mp,
				&alpha, RZX + Gp[i], &nci,
				&one, vali, &nci);
	    }
	    F77_CALL(dpotrf)("U", &nci, vali, &nci, &info);
	    if (info) error("DPOTRF returned error code %d", info);
	    F77_CALL(dpotri)("U", &nci, vali, &nci, &info);
	    if (info) error("DPOTRF returned error code %d", info);
	}
	status[0] = status[1] = 0;
    }
    ssclme_factor(x);
    return R_NilValue;
}

SEXP ssclme_fitted(SEXP x, SEXP facs, SEXP mmats)
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
    UNPROTECT(2);
    return val;
}

SEXP ssclme_variances(SEXP x, SEXP REML)
{
    SEXP Omg = PROTECT(duplicate(GET_SLOT(x, Matrix_OmegaSym))),
	val = PROTECT(allocVector(VECSXP, 2));
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omg);
    double sigmasq;
    

    SET_VECTOR_ELT(val, 0, Omg);
    SET_VECTOR_ELT(val, 1, ssclme_sigma(x, REML));
    sigmasq = REAL(VECTOR_ELT(val, 1))[0];
    sigmasq = (sigmasq) * (sigmasq);
    for (i = 0; i < nf; i++) {
	double *mm = REAL(VECTOR_ELT(Omg, i));
	int j, k, nci = nc[i], ncip1 = nci+1;

	F77_CALL(dpotrf)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error("DPOTRF returned error code %d on Omega[%d]",
		  j, i + 1);
	F77_CALL(dpotri)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error("DTRTRI returned error code %d on Omega[%d]",
		  j, i + 1);
	for (j = 0; j < nci; j++) {
	    mm[j * ncip1] *= sigmasq;
	    for (k = 0; k < j; k++) {
		mm[j + k * nci] = (mm[k + j * nci] *= sigmasq);
	    }
	}
    }
    UNPROTECT(2);
    return val;
}
