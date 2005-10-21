#include "lmer.h"
#include <float.h>

/**
 * Check validity of an lmer object.
 *
 * @param x Pointer to an lmer object
 *
 * @return TRUE if the object is a valid lmer object, else a string
 * describing the nature of the violation.
 */
SEXP lmer_validate(SEXP x)
{
    SEXP
	/* ZZxP = GET_SLOT(x, Matrix_ZZxSym), */
	ZtXP = GET_SLOT(x, Matrix_ZtXSym),
	XtXP = GET_SLOT(x, Matrix_XtXSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	RXXP = GET_SLOT(x, Matrix_RXXSym)
	/* , cnames = GET_SLOT(x, Matrix_cnamesSym) */
	;
    int *ZtXd = INTEGER(getAttrib(ZtXP, R_DimSymbol)),
	*XtXd = INTEGER(getAttrib(XtXP, R_DimSymbol));

    if (!(isReal(ZtXP) && isReal(XtXP) && isReal(RZXP) && isReal(RXXP) ))
	return mkString(_("Slots ZtX, XtX, RZX, and RXX must be real matrices"));
    if (!match_mat_dims(ZtXd, INTEGER(getAttrib(RZXP, R_DimSymbol))))
	return mkString(_("Dimensions of slots ZtX and RZX must match"));
    if (!match_mat_dims(XtXd, INTEGER(getAttrib(RXXP, R_DimSymbol))))
	return mkString(_("Dimensions of slots XtX and RXX must match"));
    if (ZtXd[1] != XtXd[0] || XtXd[0] != XtXd[1])
	return mkString(_("Slot XtX must be a square matrix with same ncol as ZtX"));
    return ScalarLogical(1);
}

/**
 * Create the pairwise crosstabulation of the elements of flist.
 *
 * @param flist pointer to the factor list.
 * @param nobs number of observations.
 * @param nc number of columns in the model matrices.
 *
 * @return the pairwise crosstabulation in the form of the ZtZ array.
 * This version does not fill in the counts as they are not needed.
 */
static SEXP
internal_crosstab(SEXP flist, int nobs, const int nc[])
{
    int i, nf = length(flist);
    int npairs = (nf * (nf + 1))/2;
    SEXP val = PROTECT(allocVector(VECSXP, npairs));
    SEXP cscbCl = PROTECT(MAKE_CLASS("dgBCMatrix"));
    int *Ti = Calloc(nobs, int),
	*nlevs = Calloc(nf, int),
	**zb = Calloc(nf, int*); /* zero-based indices */

    for (i = 0; i < nf; i++) {	/* populate the zb vectors */
	SEXP fi = VECTOR_ELT(flist, i);
	int j;

	zb[i] = Calloc(nobs, int);
	nlevs[i] = length(getAttrib(fi, R_LevelsSymbol));
	for (j = 0; j < nobs; j++) zb[i][j] = INTEGER(fi)[j] - 1;
	for (j = 0; j <= i; j++) {
	    int *ijp, ind = Lind(i, j), nnz;
	    SEXP ZZij;

	    SET_VECTOR_ELT(val, ind, ZZij = NEW_OBJECT(cscbCl));
	    ijp = INTEGER(ALLOC_SLOT(ZZij, Matrix_pSym,
				     INTSXP, nlevs[j] + 1));
	    triplet_to_col(nlevs[i], nlevs[j], nobs,
			   zb[i], zb[j], (double *) NULL,
			   ijp, Ti, (double *) NULL);
	    nnz = ijp[nlevs[j]];
	    Memcpy(INTEGER(ALLOC_SLOT(ZZij, Matrix_iSym, INTSXP, nnz)),
		   Ti, nnz);
	}
    }

    for (i = 0; i < nf; i++) Free(zb[i]);
    Free(zb); Free(nlevs); Free(Ti);
    UNPROTECT(2);
    return val;
}

SEXP lmer_Crosstab(SEXP flist)
{
    SEXP val;
    int i, nf = length(flist), nobs;
    int *nc = Calloc(nf, int);

    if (!(nf > 0 && isNewList(flist)))
	error(_("flist must be a non-empty list"));
    nobs = length(VECTOR_ELT(flist, 0));
    if (nobs < 1) error(_("flist[[1]] must be a non-null factor"));
    for (i = 0; i < nf; i++) {
	SEXP fi = VECTOR_ELT(flist, i);
	if (!(isFactor(fi) && length(fi) == nobs))
	    error(_("flist[[%d]] must be a factor of length %d"),
		  i + 1, nobs);
	nc[i] = 1;
    }
    val = internal_crosstab(flist, nobs, nc);
    Free(nc);
    return val;
}

/** 
 * Allocate the x slot in an dgBCMatrix object
 * 
 * mm Pointer to a dgBCMatrix object
 * nr number of rows per block
 * nc number of columns per block
 */
#define ALLOC_X_SLOT(mm, nr, nc) \
    SET_SLOT(mm, Matrix_xSym, alloc3Darray(REALSXP, nr, nc, \
					   length(GET_SLOT(mm, Matrix_iSym))))

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
	int *mj = expand_cmprPt(ncol, cp, Calloc(nnz, int));
	int *mi = Memcpy(Calloc(nnz, int), INTEGER(cscbi), nnz);

	if (j <= jj) int_permute(mi, nnz, iperm);
	if (j >= jj) int_permute(mj, nnz, iperm);
	if (j == jj) make_upper_triangular(mi, mj, nnz);
	triplet_to_col(nrow, ncol, nnz, mi, mj, (double *) NULL,
		       cp, INTEGER(cscbi), (double *) NULL);
	Free(mi); Free(mj);
    }
}

/** 
 * Apply a permutation of the rows and columns to a sparse symmetric
 * matrix object.
 * 
 * @param A A sparse, symmetric matrix object stored in the upper
 * triangle
 * @param nlev order of A
 * @param iperm A 0-based permutation of length nlev
 */
static void
symmetric_permute(int Ap[], int Ai[], int n, const int iperm[])
{
    int nnz = Ap[n];
    int *mj = expand_cmprPt(n, Ap, Calloc(nnz, int));
    int *mi = Memcpy(Calloc(nnz, int), Ai, nnz);

    int_permute(mi, nnz, iperm);
    int_permute(mj, nnz, iperm);
    make_upper_triangular(mi, mj, nnz);
    triplet_to_col(n, n, nnz, mi, mj, (double *) NULL,
		   Ap, Ai, (double *) NULL);
    Free(mi); Free(mj);
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
    SEXP D, L, Parent, ZZpO, 
	flist = GET_SLOT(val, Matrix_flistSym),
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
    L = ALLOC_SLOT(val, Matrix_LSym, VECSXP, npairs);
    Gp = INTEGER(ALLOC_SLOT(val, Matrix_GpSym, INTSXP, nf + 1));
    Gp[0] = 0;
    for (j = 0; j < nf; j++) {
	nlev[j] = length(getAttrib(VECTOR_ELT(flist, j), R_LevelsSymbol));
	Gp[j + 1] = Gp[j] + nc[j] * nlev[j];
	SET_VECTOR_ELT(D, j, alloc3Darray(REALSXP, nc[j], nc[j], nlev[j]));
	AZERO(REAL(VECTOR_ELT(D, j)), nc[j] * nc[j] * nlev[j]);
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
	    *ci = INTEGER(cip),
	    ncj = length(cpp) - 1,
	    nnz = length(cip);
				
	SET_VECTOR_ELT(Parent, j, Matrix_make_named(VECSXP, pnms));
	parent = VECTOR_ELT(Parent, j);
	SET_VECTOR_ELT(parent, 0, allocVector(INTSXP, ncj));
	SET_VECTOR_ELT(parent, 1, allocVector(INTSXP, ncj));
	SET_VECTOR_ELT(perm, j, allocVector(INTSXP, ncj));
	Perm = INTEGER(VECTOR_ELT(perm, j));
	if (nnz > ncj) {	/* calculate fill-reducing permutation */
	    SEXP fac = VECTOR_ELT(flist, j);
	    SEXP fcp = PROTECT(duplicate(fac));
	    int *iPerm = Calloc(ncj, int);

	    ssc_metis_order(ncj, cp, ci, Perm, iPerm);
				/* apply to the crosstabulation, L, and ZZpO */
	    bCrosstab_permute(ZtZ, nf, j, nlev, iPerm);
	    bCrosstab_permute(L, nf, j, nlev, iPerm);
	    symmetric_permute(cp, ci, nlev[j], iPerm);
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
	}
	for (k = j+1; k < nf; k++) { /* Update other blocks in this column */
	    SEXP Lkj = VECTOR_ELT(L, Lind(k,j));
	    SET_SLOT(Lkj, Matrix_iSym,
		     lCholClgCsm(RGT, TRN, nlev[k], nlev[j],
				 INTEGER(VECTOR_ELT(parent, 0)),
				 GET_SLOT(Lkj, Matrix_iSym),
				 INTEGER(GET_SLOT(Lkj, Matrix_pSym))));
	}
	for (k = j + 1; k < nf; k++) { /* Update remaining columns */
	    SEXP db = VECTOR_ELT(ZZpO, k), Lkj = VECTOR_ELT(L, Lind(k, j));
	    int *Lkji = INTEGER(GET_SLOT(Lkj, Matrix_iSym)),
		*Lkjp = INTEGER(GET_SLOT(Lkj, Matrix_pSym));
	    SET_SLOT(db, Matrix_iSym,
		     Matrix_lgCsyrk(1, 0, nlev[k], nlev[j], Lkji, Lkjp,
				    1, GET_SLOT(db, Matrix_iSym),
				    INTEGER(GET_SLOT(db, Matrix_pSym))));
	    for (i = k + 1; i < nf; i++) {
		SEXP Lij = VECTOR_ELT(L, Lind(i, j)),
		    Lik = VECTOR_ELT(L, Lind(i, k));
		SET_SLOT(Lik, Matrix_iSym,
			 Matrix_lgClgCmm(0, 1, nlev[i], nlev[k], nlev[j],
					 INTEGER(GET_SLOT(Lij, Matrix_iSym)),
					 INTEGER(GET_SLOT(Lij, Matrix_pSym)),
					 Lkji, Lkjp,
					 1, GET_SLOT(Lik, Matrix_iSym),
					 INTEGER(GET_SLOT(Lik, Matrix_pSym))));
	    }
	}
    }
				
    for (j = 0; j < nf; j++) {	/* allocate x slots in dgBCMatrix objects */
	ALLOC_X_SLOT(VECTOR_ELT(ZZpO, j), nc[j], nc[j]);
	for (k = j; k < nf; k++) {
	    int indkj = Lind(k,j);
	    ALLOC_X_SLOT(VECTOR_ELT(L, indkj), nc[k], nc[j]);
	    ALLOC_X_SLOT(VECTOR_ELT(ZtZ, indkj), nc[k], nc[j]);
	}
    }
/* FIXME: Use these macros from Tim Davis instead */
#define EMPTY -1
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) (((i) < EMPTY) ? FLIP(i) : (i))
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

/**
 * Update the arrays ZtZ, ZtX, and XtX in an lme object
 * according to a list of model matrices.
 *
 * @param x pointer to an lmer object
 * @param mmats pointer to a list of model matrices
 *
 * @return NULL
 */
SEXP lmer_update_mm(SEXP x, SEXP mmats)
{
    SEXP
	ZtZP = GET_SLOT(x, Matrix_ZtZSym),
	ZtXP = GET_SLOT(x, Matrix_ZtXSym),
	flist = GET_SLOT(x, Matrix_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*dims = INTEGER(getAttrib(ZtXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	nf = length(flist), nfp1 = nf + 1,
	i, ione = 1,
	nobs = nc[nfp1],
	pp1 = nc[nf];
    double
	*X,
	*XtX = REAL(GET_SLOT(x, Matrix_XtXSym)),
	*ZtX = REAL(ZtXP),
	one = 1.0, zero = 0.0;

    if (pp1 < 0) pp1 = -pp1;
    if (!isNewList(mmats) || length(mmats) != nfp1)
	error(_("mmats must be a list of %d model matrices"), nfp1);
    for (i = 0; i <= nf; i++) {
	SEXP mmat = VECTOR_ELT(mmats, i);
	int *mdims = INTEGER(getAttrib(mmat, R_DimSymbol));

	if (!isMatrix(mmat) || !isReal(mmat))
	    error(_("mmats[[%d]] is not a numeric matrix"), i + 1);
	if (nobs != mdims[0])
	    error(_("Expected %d rows in mmats[[%d]]. Got %d"),
		  nobs, i+1, mdims[0]);
	if ((nc[i] < 0 ? -nc[i] : nc[i]) != mdims[1])
	    error(_("Expected %d columns in mmats[[%d]]. Got %d"),
		  nc[i], i+1, mdims[1]);
    }
				/* Create XtX */
    X = REAL(VECTOR_ELT(mmats, nf));
    F77_CALL(dsyrk)("U", "T", &pp1, &nobs, &one, X, &nobs, &zero, XtX, &pp1);
				/* Zero an accumulator */
    AZERO(ZtX, pp1 * Gp[nf]);
    for (i = 0; i < nf; i++) {
	int *fac = INTEGER(VECTOR_ELT(flist, i)),
	    j, k, nci = nc[i], ZtXrows = Gp[i+1] - Gp[i];
	int ncisqr = nci * nci, nlev = ZtXrows/nci;
	double *Z = REAL(VECTOR_ELT(mmats, i)), *ZZx;

	for (k = 0; k < i; k++) {
	    SEXP ZZxM = VECTOR_ELT(ZtZP, Lind(i, k));
	    int *rowind = INTEGER(GET_SLOT(ZZxM, Matrix_iSym)),
		*colptr = INTEGER(GET_SLOT(ZZxM, Matrix_pSym));
	    int *f2 = INTEGER(VECTOR_ELT(flist, k)), nck = nc[k];
	    double *Zk = REAL(VECTOR_ELT(mmats, k));

	    ZZx = REAL(GET_SLOT(ZZxM, Matrix_xSym));
	    AZERO(ZZx, length(GET_SLOT(ZZxM, Matrix_xSym)));
	    for (j = 0; j < nobs; j++) {
		F77_CALL(dgemm)("T", "N", nc + i, nc + k, &ione, &one,
				Z + j, &nobs, Zk + j, &nobs, &one,
				ZZx + check_csc_index(colptr, rowind,
						      fac[j] - 1, f2[j] - 1, 0)
				* (nci * nck), &nci);
	    }
	}
	ZZx = REAL(GET_SLOT(VECTOR_ELT(ZtZP, Lind(i, i)), Matrix_xSym));
	AZERO(ZZx, nci * nci * nlev);
	if (nci == 1) {		/* single column in Z */
	    for (j = 0; j < nobs; j++) {
		int fj = fac[j] - 1; /* factor indices are 1-based */
		ZZx[fj] += Z[j] * Z[j];
		F77_CALL(daxpy)(&pp1, Z + j, X + j, &nobs, ZtX + fj, dims);
	    }
	} else {
	    for (j = 0; j < nobs; j++) {
		int fj = fac[j] - 1; /* factor indices are 1-based */

		F77_CALL(dsyr)("U", nc + i, &one, Z + j, &nobs,
			       ZZx + fj * ncisqr, nc + i);
		F77_CALL(dgemm)("T", "N", nc + i, &pp1, &ione,
				&one, Z + j, &nobs,
				X + j, &nobs, &one,
				ZtX + fj * nci, dims);
	    }
	}
	ZtX += ZtXrows;
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}

/**
 * Create an lmer object from a list of grouping factors and a list of model
 * matrices.  There is one more model matrix than grouping factor.  The last
 * model matrix is the fixed effects and the response.
 *
 * @param flist pointer to a list of grouping factors
 * @param mmats pointer to a list of model matrices
 *
 * @return pointer to an lmer object
 */
SEXP lmer_create(SEXP flist, SEXP mmats, SEXP method)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("mer")));
    SEXP ZtZ, cnames, fnms, nms;
    int *nc, i, nf = length(flist), nobs;

				/* Check validity of flist */
    if (!(nf > 0 && isNewList(flist)))
	error(_("flist must be a non-empty list"));
    nobs = length(VECTOR_ELT(flist, 0));
    if (nobs < 1) error(_("flist[[0]] must be a non-null factor"));
    for (i = 0; i < nf; i++) {
	SEXP fi = VECTOR_ELT(flist, i);
	if (!(isFactor(fi) && length(fi) == nobs))
	    error(_("flist[[%d]] must be a factor of length %d"),
		  i + 1, nobs);
    }
    SET_SLOT(val, Matrix_flistSym, duplicate(flist));
				/* Check mmats; allocate and populate nc */
    if (!(isNewList(mmats) && length(mmats) == (nf + 1)))
	error(_("mmats must be a list of length %d"), nf + 1);
    nc = INTEGER(ALLOC_SLOT(val, Matrix_ncSym, INTSXP, nf + 2));
    nc[nf + 1] = nobs;
    for (i = 0; i <= nf; i++) {
	SEXP mi = VECTOR_ELT(mmats, i);
	int *dims;

	if (!(isMatrix(mi) && isReal(mi)))
	    error(_("mmats[[%d]] must be a numeric matrix"), i + 1);
	dims = INTEGER(getAttrib(mi, R_DimSymbol));
	if (dims[0] != nobs)
	    error(_("mmats[[%d]] must have %d rows"), i + 1, nobs);
	if (dims[1] < 1)
	    error(_("mmats[[%d]] must have at least 1 column"), i + 1);
	nc[i] = dims[1];
    }   /* Arguments have now been checked for type, dimension, etc. */
				/* Create pairwise crosstabulation in ZtZ */
    SET_SLOT(val, Matrix_ZtZSym, internal_crosstab(flist, nobs, nc));
    SET_SLOT(val, Matrix_methodSym, duplicate(method));
    lmer_populate(val);
    ZtZ = GET_SLOT(val, Matrix_ZtZSym);
    /* FIXME: Check for possible reordering of the factors to maximize the
     * number of levels (columns?) in the leading nested sequence. */
    fnms = getAttrib(flist, R_NamesSymbol);
				/* Allocate and populate cnames */
    cnames = ALLOC_SLOT(val, Matrix_cnamesSym, VECSXP, nf + 1);
    setAttrib(cnames, R_NamesSymbol, allocVector(STRSXP, nf + 1));
    nms = getAttrib(cnames, R_NamesSymbol);
    for (i = 0; i <= nf; i++) {
	SEXP mi = VECTOR_ELT(mmats, i);
	SET_VECTOR_ELT(cnames, i,
		       duplicate(VECTOR_ELT(getAttrib(mi, R_DimNamesSymbol),
					    1)));
	SET_STRING_ELT(nms, i, (i < nf) ? duplicate(STRING_ELT(fnms, i)) :
		       mkChar(".fixed"));
    }
    lmer_update_mm(val, mmats);
    SET_SLOT(val, Matrix_bVarSym, duplicate(GET_SLOT(val, Matrix_DSym)));
    UNPROTECT(1);
    return val;
}

/**
 * Create and insert initial values for Omega.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP lmer_initial(SEXP x)
{
    SEXP Omg = GET_SLOT(x, Matrix_OmegaSym);
    int	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)), i, nf = length(Omg);

    for (i = 0; i < nf; i++) {
	SEXP ZZxP = GET_SLOT(VECTOR_ELT(GET_SLOT(x, Matrix_ZtZSym), Lind(i, i)),
			     Matrix_xSym);
	int *dims = INTEGER(getAttrib(ZZxP, R_DimSymbol));
	int j, k, nzc = dims[0], nlev = dims[2];
	int nzcsqr = nzc * nzc, nzcp1 = nzc + 1;
	double *Omega = REAL(VECTOR_ELT(Omg, i)),
	    mi = 0.375 / ((double) nlev);

	AZERO(Omega, nzc * nzc);
	for (j = 0; j < nlev; j ++) {
	    for (k = 0; k < nzc; k++) {
		Omega[k * nzcp1] += REAL(ZZxP)[k * nzcp1 + j * nzcsqr] * mi;
	    }
	}
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}

/**
 * Create and insert initial values for Omega.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP lmer_set_initial(SEXP x, SEXP iv)
{
    SEXP Omg = GET_SLOT(x, Matrix_OmegaSym);
    int	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)), i, nf = LENGTH(Omg);

    if (!isNewList(iv) || LENGTH(iv) != nf)
	error(_("Initial values must be a list of length %d"), nf);
    for (i = 0; i < nf; i++) {
	SEXP Omgi = VECTOR_ELT(Omg, i), ivi = VECTOR_ELT(iv, i);
	int *dims = INTEGER(getAttrib(ivi, R_DimSymbol)), j, k,
	    nci = *INTEGER(getAttrib(Omgi, R_DimSymbol));  
	double *Omega = REAL(Omgi), *Ivi = REAL(ivi);

	if (!isMatrix(ivi) || !isReal(ivi) || dims[0] != nci ||
	    dims[1] != nci)
	    error(_("iv[[%d]] must be a real %d by %d matrix"),
		  i + 1, nci, nci);
	AZERO(Omega, nci * nci);
	for (j = 0; j < nci; j++) {
	    for (k = 0; k <= j; k++) {
		int ind = k + j * nci;
		Omega[ind] += Ivi[ind];
	    }
	}
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}


/**
 * Copy ZtZ to ZZpO and L.  Inflate diagonal blocks of ZZpO by Omega.
 * Update devComp[1].
 *
 * @param x pointer to an lmer object
 */
SEXP
lmer_inflate(SEXP x)
{
    SEXP Omg = GET_SLOT(x, Matrix_OmegaSym),
	ZZpO = GET_SLOT(x, Matrix_ZZpOSym),
	ZtZ = GET_SLOT(x, Matrix_ZtZSym),
	LP = GET_SLOT(x, Matrix_LSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, k, nf = length(Omg);
    double *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym));

    for (i = 0; i < nf; i++) {
	SEXP ZZOel = VECTOR_ELT(ZZpO, i);
	SEXP ZZOm = GET_SLOT(ZZOel, Matrix_xSym);
	SEXP ZZel = VECTOR_ELT(ZtZ, Lind(i, i));
	int *Di = INTEGER(GET_SLOT(ZZOel, Matrix_iSym)),
	    *Dp = INTEGER(GET_SLOT(ZZOel, Matrix_pSym)),
	    *Si = INTEGER(GET_SLOT(ZZel, Matrix_iSym)),
	    *Sp = INTEGER(GET_SLOT(ZZel, Matrix_pSym)),
	    *dims = INTEGER(getAttrib(ZZOm, R_DimSymbol));
	int sz = dims[0] * dims[1];
	int ii, j, nci = nc[i], ncisqr = nci * nci;
	int nlev = (Gp[i + 1] - Gp[i])/nci;
	double *Omega = REAL(VECTOR_ELT(Omg, i)),
	    *ZZ = REAL(GET_SLOT(ZZel, Matrix_xSym)),
	    *tmp = Memcpy(Calloc(ncisqr, double), Omega, ncisqr);

	F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
	if (j)
	    error(_("Leading %d minor of Omega[[%d]] not positive definite"),
		  j, i + 1);
				/* update dcmp[1] */
	for (j = 0; j < nci; j++) { /* nlev * logDet(Omega_i) */
	    dcmp[1] += nlev * 2. * log(tmp[j * (nci + 1)]);
	}
	Free(tmp);
	AZERO(REAL(ZZOm), dims[0] * dims[1] * dims[2]);
	for (j = 0; j < nlev; j++) { /* copy diagonal block and inflate */
	    double *ZZOkk = REAL(ZZOm) + check_csc_index(Dp, Di, j, j, 0) * sz;
	    int kk, k2 = Sp[j + 1];
	    for (kk = Sp[j]; kk < k2; kk++) {
		Memcpy(REAL(ZZOm) + check_csc_index(Dp, Di, Si[kk], j, 0) * sz,
		       ZZ + kk * sz, sz);
	    }
	    for (kk = 0; kk < nci; kk++) {
		for (ii = 0; ii <= kk; ii++) {
		    int ind = ii + kk * nci;
		    ZZOkk[ind] += Omega[ind];
		}
	    }
	}
	for (k = i + 1; k < nf; k++) {
	    int ind = Lind(k, i);
	    SEXP Lel = VECTOR_ELT(LP, ind),
		Lm = GET_SLOT(Lel, Matrix_xSym);
	    double *L = REAL(Lm);

	    dims = INTEGER(getAttrib(Lm, R_DimSymbol));
	    ZZel = VECTOR_ELT(ZtZ, ind);
	    ZZ = REAL(GET_SLOT(ZZel, Matrix_xSym));
	    Di = INTEGER(GET_SLOT(Lel, Matrix_iSym));
	    Dp = INTEGER(GET_SLOT(Lel, Matrix_pSym));
	    Si = INTEGER(GET_SLOT(ZZel, Matrix_iSym));
	    Sp = INTEGER(GET_SLOT(ZZel, Matrix_pSym));
	    sz = dims[0] * dims[1];

	    AZERO(L, sz * dims[2]); /* zero L  */
	    for (j = 0; j < nlev; j++) { /* copy src blocks to dest */
		int kk, k2 = Sp[j + 1];
		for (kk = Sp[j]; kk < k2; kk++) {
		    Memcpy(L + check_csc_index(Dp, Di, Si[kk], j, 0) * sz,
			   ZZ + kk * sz, sz);
		}
	    }
	}
    }
    return R_NilValue;
}

/**
 * Convert the extended parent pair (Parent, Block) to a parent array
 * for the jth diagonal block of size n.
 *
 * @param j index (0-based) of the diagonal outer block
 * @param n number of inner column blocks in the outer block
 * @param par array of length n to be filled with the parent array
 * @param ParP pointer to the extended parent structure
 *
 * @return par
 */
static R_INLINE
int *block_parent(int j, int n, int par[], SEXP ParP)
{
    SEXP Parj = VECTOR_ELT(ParP, j);
    int *Parent = INTEGER(VECTOR_ELT(Parj, 0)),
	*Block = INTEGER(VECTOR_ELT(Parj, 1)), i;
    for (i = 0; i < n; i++) par[i] = (Block[i] == j) ? Parent[i] : -1;
    return par;
}

/**
 * If status[["factored"]] is FALSE, create and factor Z'Z+Omega.  Also
 * create RZX and RXX, the deviance components, and the value of the
 * deviance for both ML and REML.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP lmer_factor(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));

    if (!status[0]) {
	SEXP DP = GET_SLOT(x, Matrix_DSym),
	    LP = GET_SLOT(x, Matrix_LSym),
	    RZXP = GET_SLOT(x, Matrix_RZXSym),
	    ZZOP = GET_SLOT(x, Matrix_ZZpOSym),
	    Parent = GET_SLOT(x, Matrix_ParentSym);
	int *dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	    i, j, nf = length(DP);
	int nml = nc[nf + 1], nreml = nml + 1 - nc[nf], ncX = dims[1];
	double
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(RZXP),
	    *ZtX = REAL(GET_SLOT(x, Matrix_ZtXSym)),
	    *XtX = REAL(GET_SLOT(x, Matrix_XtXSym)),
	    *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	    *deviance = REAL(GET_SLOT(x, Matrix_devianceSym)),
	    minus1 = -1., one = 1.;

	if (nc[nf] < 0) {	/* skip the fixed effects part */
	    int ZXshift = dims[0] * (dims[1] - 1),
		XXshift = dims[1]*dims[1] - 1;
	    
	    RZX += ZXshift; ZtX += ZXshift; RXX += XXshift; XtX += XXshift;
	    ncX = 1;
	    nreml = nml + 1 + nc[nf];
	}
	dcmp[0] = dcmp[1] = dcmp[2] = dcmp[3] = 0.;
	Memcpy(RZX, ZtX, dims[0] * ncX);
	lmer_inflate(x);	/* initialize ZZpO and L */
	for (i = 0; i < nf; i++) {
	    SEXP ZZOiP = VECTOR_ELT(ZZOP, i);
	    SEXP DiP = VECTOR_ELT(DP, i);
	    SEXP LiP = VECTOR_ELT(LP, Lind(i, i));
	    int nlev = INTEGER(getAttrib(DiP, R_DimSymbol))[2];
	    int jj, nci = nc[i], ncisqr = nci * nci;
	    int *Pari = block_parent(i, nlev, Calloc(nlev, int), Parent);
	    double *D = REAL(DiP);

	    jj = cscb_ldl(ZZOiP, Pari, LiP, DiP);
	    if (jj != nlev) error(_("cscb_ldl returned %d < nlev = %d"), jj, nlev);
	    for (j = 0; j < nlev; j++) { /* accumulate dcmp[0] */
		double *Dj = D + j * ncisqr;
		for (jj = 0; jj < nci; jj++) /* accumulate determinant */
		    dcmp[0] += 2. * log(Dj[jj * (nci + 1)]);
	    }
	    /* Solve L_{i,i} %*% RZX_i := RZX_i */
	    cscb_trsm(LOW, NTR, UNT, 1., LiP,
		      Gp[i+1] - Gp[i], ncX, RZX + Gp[i], dims[0]);
	    /* Solve D_i^{T/2} %*% RZX_i := RZX_i */
	    for (jj = 0; jj < nlev; jj++) {
		F77_CALL(dtrsm)("L", "U", "T", "N", &nci, &ncX,
				&one, D + jj * ncisqr, &nci,
				RZX + Gp[i] + jj * nci, &dims[0]);
	    }
	    for (j = i + 1; j < nf; j++) { /*  further blocks */
		SEXP Lji = VECTOR_ELT(LP, Lind(j, i));
		SEXP Lx = GET_SLOT(Lji, Matrix_xSym);
		double *L = REAL(Lx);
		int *xdims = INTEGER(getAttrib(Lx, R_DimSymbol)),
		    *Lp = INTEGER(GET_SLOT(Lji, Matrix_pSym));
		int ntot = xdims[0] * xdims[1];

		/* L_{j,i} := L_{j,i} %*% L_{i,i}^{-T} %*% D_i^{-1/2} */
		cscb_trcbsm(RGT, LOW, TRN, UNT, 1.0, LiP, Pari, Lji);
		for (jj = 0; jj < nlev; jj++) {
		    int k, k2 = Lp[jj + 1];
		    for (k = Lp[jj]; k < k2; k++)
			F77_CALL(dtrsm)("R", "U", "N", "N", xdims, xdims + 1,
					&one, D + jj * ncisqr, &nci,
					L + k * ntot, xdims);
		}
		/* RZX_j := RZX_j - (L_{j,i} %*% D_i^{T/2}) %*% RZX_i */
		/* At this point Lji contains L_{j,i} %*% D_i^{T/2} */
		cscb_mm(LFT, NTR, Gp[j + 1] - Gp[j], ncX, Gp[i+1] - Gp[i],
			-1.0, Lji, RZX + Gp[i], dims[0],
			1.0, RZX + Gp[j], dims[0]);
	    }
	    for (j = i + 1; j < nf; j++) { /* block pairs and final update */
		SEXP Lji = VECTOR_ELT(LP, Lind(j, i));
		SEXP Lx = GET_SLOT(Lji, Matrix_xSym);
		double *L = REAL(Lx);
		int *xdims = INTEGER(getAttrib(Lx, R_DimSymbol)),
		    *Lp = INTEGER(GET_SLOT(Lji, Matrix_pSym));
		int ntot = xdims[0] * xdims[1];


		/* ZZpO_{j,j} := ZZpO_{j,j} - L{j,i} %*% L_{j,i}^T */
		cscb_syrk(UPP, NTR, -1.0, Lji, 1.0, VECTOR_ELT(ZZOP, j));
		for (jj = j+1; jj < nf; jj++) {
		    /* L_{jj,j} := L_{jj,j} - L{jj,i} %*% L_{j,i}^T */
		    cscb_cscbm(NTR, TRN, -1.0, VECTOR_ELT(LP, Lind(jj, i)),
			    Lji, 1.0, VECTOR_ELT(LP, Lind(jj, j)));
		}
		/* L_{j,i} := L_{j,i} %*% D_i^{-T/2} */
		for (jj = 0; jj < nlev; jj++) {
		    int k, k2 = Lp[jj + 1];
		    for (k = Lp[jj]; k < k2; k++)
			F77_CALL(dtrsm)("R", "U", "T", "N", xdims, xdims + 1,
					&one, D + jj * ncisqr, &nci,
					L + k * ntot, xdims);
		}
	    }
	    Free(Pari);
	}
				/* downdate and factor XtX */
	Memcpy(RXX, XtX, ncX * ncX);
	F77_CALL(dsyrk)("U", "T", &ncX, &dims[0],
			&minus1, RZX, &dims[0], &one, RXX, &ncX);
	F77_CALL(dpotrf)("U", &ncX, RXX, &ncX, &j);
	if (j) {
	    warning(_("Leading minor of size %d of downdated X'X is indefinite"),
		    j);
	    dcmp[2] = dcmp[3] = deviance[0] = deviance[1] = NA_REAL;
	} else {
	    for (j = 0; j < (ncX - 1); j++) /* 2 logDet(RXX) */
		dcmp[2] += 2 * log(RXX[j * (dims[1] + 1)]);
	    dcmp[3] = 2. * log(RXX[ncX * ncX - 1]); /* 2 log(ryy) */
	    deviance[0] =	/* ML criterion */
		dcmp[0] - dcmp[1] + nml*(1.+dcmp[3]+log(2.*PI/nml));
	    deviance[1] = dcmp[0] - dcmp[1] + /* REML */
		dcmp[2] + nreml*(1.+dcmp[3]+log(2.*PI/nreml));
	}
	status[0] = 1; status[1] = 0; /* factored but not inverted */
    }
    return R_NilValue;
}

/**
 * Solve one of the matrix equations op(L)*X=alpha*B or
 * X*op(L)=alpha*B where L is a sparse, blocked, unit lower triangular matrix.
 *
 * @param side LFT or RGT for left or right
 * @param trans TRN or NTR for transpose or no transpose
 * @param nf number of grouping factors
 * @param Gp group pointers for the rows
 * @param n number of columns
 * @param alpha multiplier
 * @param L pointer to the L cscb object
 * @param B pointer to the matrix of right-hand sides
 * @param ldb leading dimension of array B as declared in the caller
 */
static void
internal_sm(enum CBLAS_SIDE side, enum CBLAS_TRANSPOSE trans, int nf,
	    const int Gp[], int n, double alpha, SEXP L, double B[], int ldb)
{
    int j, k;

    if (side == LFT) {
	if (trans == TRN) {
	    for (j = nf - 1; j >= 0; j--) {
		int nrj = Gp[j + 1] - Gp[j];

		cscb_trsm(LOW, TRN, UNT, alpha, VECTOR_ELT(L, Lind(j, j)),
			  nrj, n, B + Gp[j], ldb);
		for (k = 0; k < j; k++) {
		    cscb_mm(LFT, TRN, Gp[k + 1] - Gp[k], n, nrj,
			    -1., VECTOR_ELT(L, Lind(j, k)),
			    B + Gp[j], ldb, alpha, B + Gp[k], ldb);
		}
	    }
	} else error(_("Code for non-transpose case not yet written"));
    } else error(_("Code for right-side solutions not yet written"));
}

/** 
 * Determine the maximum number of nonzero elements in a column and
 * allocate storage for the tmp and ind arrays.
 * 
 * @param j level
 * @param Parent Parent list
 * 
 * @return Maximum number of nonzero elements in a column
 */
static void
alloc_tmp_ind(int nf, const int nc[], const int nlevs[], SEXP Parent,
	      double *tmp[], int *ind[])
{
    int j, maxnc;
    for (maxnc = -1, j = 0; j < nf; j++) {
	SEXP lst = VECTOR_ELT(Parent, j);
	SEXP blk = VECTOR_ELT(lst, 1), par = VECTOR_ELT(lst, 0);
	int *nfj = Calloc(nlevs[j], int), i, val;

	
	if (nc[j] > maxnc) maxnc = nc[j];
	for (val = -1, i = nlevs[j] - 1; i >= 0; i--) {
	    int thisnnz = (INTEGER(blk)[i] != j) ? 1 : nfj[INTEGER(par)[i]] + 1;
	    if (thisnnz > val) val = thisnnz;
	    nfj[i] = thisnnz;
	}
	ind[j] = Calloc(val, int);
	tmp[j] = Calloc(val * nc[j] * maxnc, double);
	Free(nfj);
    }
}

#define BLK(i,j) INTEGER(VECTOR_ELT(VECTOR_ELT(Parent, i), 1))[j]
#define PAR(i,j) INTEGER(VECTOR_ELT(VECTOR_ELT(Parent, i), 0))[j]

/**
 * Fill the nnz array with the number of nonzero inner blocks in each
 * outer block of the jth inner column block of the ith outer block of
 * L^{-1}.  Also fill the ind array.
 *
 * @param i outer block index
 * @param j inner block index within the ith outer block
 * @param nf number of factors
 * @param Parent pointer to the extended parent pairs
 * @param nc
 * @param nnz array of length nf
 * @param tmp array of length nf of pointers to doubles
 * @param ind array of length nf of pointers to ints
 *
 */
static
void fill_ind(int i, int j, int nf, SEXP Parent, int nnz[], int *ind[])
{
    int blk, k, par;

    AZERO(nnz, nf);
    for (blk = BLK(i,j), par = PAR(i,j); blk >= 0;
	 k = BLK(blk,par), par = PAR(blk,par), blk = k) {
	ind[blk][nnz[blk]++] = par;
    }
}

static R_INLINE
int fsrch(int target, const int vals[], int nvals)
{
    int i;
    for (i = 0; i < nvals; i++) if (vals[i] == target) return i;
    error(_("fsrch: unable to find target %d in nvals %d "), target, nvals);
    return -1;			/* -Wall */
}

/**
 * If necessary, factor Z'Z+Omega, ZtX, and XtX then, if necessary,
 * replace the RZX and RXX slots by the corresponding parts of the
 * inverse of the Cholesky factor.  Replace the elements of the D slot
 * by the blockwise inverses and evaluate bVar.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL (x is updated in place)
 */
SEXP lmer_invert(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    if (!status[0]) lmer_factor(x);
    if (!R_FINITE(REAL(GET_SLOT(x, Matrix_devianceSym))[0]))
	error(_("Unable to invert singular factor of downdated X'X"));
    if (!status[1]) {
	SEXP DP = GET_SLOT(x, Matrix_DSym),
	    LP = GET_SLOT(x, Matrix_LSym),
	    ParP = GET_SLOT(x, Matrix_ParentSym),
	    RZXP = GET_SLOT(x, Matrix_RZXSym),
	    bVarP = GET_SLOT(x, Matrix_bVarSym);
	int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	    *dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    i, nf = length(DP);
	int **ind = Calloc(nf, int *),
	    *nlevs = Calloc(nf, int),
	    *nnz = Calloc(nf, int), ncX = dims[1];
	double **tmp = Calloc(nf, double *),
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(RZXP),
	    minus1 = -1., one = 1., zero = 0.;

	if (nc[nf] < 0) {	/* skip the fixed effects part */
	    RZX += dims[0] * (dims[1] - 1);
	    RXX += dims[1]*dims[1] - 1;
	    ncX = 1;
	}
	/* RXX := RXX^{-1} */
	F77_CALL(dtrtri)("U", "N", &ncX, RXX, &ncX, &i);
	if (i)
	    error(_("Leading minor of size %d of downdated X'X,is indefinite"),
		  i + 1);

				/* RZX := - RZX %*% RXX */
	F77_CALL(dtrmm)("R", "U", "N", "N", &dims[0], &ncX, &minus1,
			RXX, &ncX, RZX, &dims[0]);
	for(i = 0; i < nf; i++) {
	    int info, j, jj, nci = nc[i];
	    int ncisqr = nci * nci;
	    double *Di = REAL(VECTOR_ELT(DP, i)),
		*RZXi = RZX + Gp[i];

	    nlevs[i] = (Gp[i+1] - Gp[i])/nci;
				/* D_i := D_i^{-1}; RZX_i := D_i %*% RZX_i */
	    if (nci == 1) {
		for (j = 0; j < nlevs[i]; j++) {
		    Di[j] = 1./Di[j];
		    for (jj = 0; jj < ncX; jj++)
			RZXi[j + jj * dims[0]] *= Di[j];
		}
	    } else {
		for (j = 0; j < nlevs[i]; j++) {
		    F77_CALL(dtrtri)("U", "N", &nci, Di + j * ncisqr,
				     &nci, &info);
		    if (info) error(_("D[,,%d] for factor %d is singular"),
				    j + 1, i + 1);
		    F77_CALL(dtrmm)("L", "U", "N", "N", &nci, &ncX, &one,
				    Di + j * ncisqr, &nci, RZXi + j * nci,
				    &dims[0]);
		}
	    }
	}

				/* RZX := L^{-T} %*% RZX */
	internal_sm(LFT, TRN, nf, Gp, ncX, 1.0, LP, RZX, dims[0]);

	alloc_tmp_ind(nf, nc, nlevs, ParP, tmp, ind);
	/* Create bVar as crossprod of column blocks of D^{-T/2}%*%L^{-1} */
	for (i = 0; i < nf; i++) { /* ith column of outer blocks */
	    int j, k, kj, nci = nc[i];
	    int ncisqr = nci * nci;
	    double *Di = REAL(VECTOR_ELT(DP, i)),
		*bVi = REAL(VECTOR_ELT(bVarP, i));

	    AZERO(bVi, ncisqr * nlevs[i]);
	    for (j = 0; j < nlevs[i]; j++) {
		double *bVij = bVi + j * ncisqr, *Dij = Di + j * ncisqr;
		
		F77_CALL(dsyrk)("U", "N", &nci, &nci, &one, Dij,
				&nci, &zero, bVij, &nci);
		/* count non-zero blocks; allocate and zero storage */
		fill_ind(i, j, nf, ParP, nnz, ind);

		for (k = i; k < nf; k++) { /* kth row of outer blocks */
		    SEXP Lki = VECTOR_ELT(LP, Lind(k, i));
		    int *Lkii = INTEGER(GET_SLOT(Lki, Matrix_iSym)),
			*Lkip = INTEGER(GET_SLOT(Lki, Matrix_pSym));
		    double *Lkix = REAL(GET_SLOT(Lki, Matrix_xSym));
		    int kk, sz = nc[i] * nc[k];
		    
		    AZERO(tmp[k], sz * nnz[k]);
		    /* initialize tmp from jth column of (k,i)th block */
		    /* - sign in sol'n incorporated in dtrmm call below */
		    for (kk = Lkip[j]; kk < Lkip[j + 1]; kk++)
			Memcpy(tmp[k] + fsrch(Lkii[kk], ind[k], nnz[k]) * sz,
			       Lkix + kk * sz, sz);
		    /* columns in ind[kk] for (k,kk)th block */
		    for (kk = i; kk <= k; kk++) {
			int szk = nc[k] * nc[kk];
			/* skip getting slots if not using them */
			if (!nnz[kk]) continue;
			Lki = VECTOR_ELT(LP, Lind(k, kk));
			Lkii = INTEGER(GET_SLOT(Lki, Matrix_iSym));
			Lkip = INTEGER(GET_SLOT(Lki, Matrix_pSym));
			Lkix = REAL(GET_SLOT(Lki, Matrix_xSym));
			for (kj = 0; kj < nnz[kk]; kj++) {
			    int col = ind[kk][kj], k1, szkk = nc[i] * nc[kk];
			    
			    for (k1 = Lkip[col]; k1 < Lkip[col + 1]; k1++) {
				if ((kk == k) && col >= Lkii[k1]) break;
				F77_CALL(dgemm)("N", "N", &nc[k], &nci,
						&nc[kk], &minus1,
						Lkix + k1 * szk, &nc[k],
						tmp[kk] + kj * szkk, &nc[kk],
						&one, tmp[k] +
						fsrch(Lkii[k1],ind[k],nnz[k]) *
						sz, &nc[k]);
			    }
			}
		    }
		}
		for (k = 0; k < nf; k++) {
		    for (kj = 0; kj < nnz[k]; kj++) {
			F77_CALL(dtrmm)("L", "U", "T", "N", nc + k, &nci,
					&minus1,
					REAL(VECTOR_ELT(DP, k))+
					ind[k][kj]*nc[k]*nc[k],
					nc + k, tmp[k] + kj * nc[i] * nc[k],
					nc + k);
		    }
		    if (nnz[k] > 0) {
			kj = nc[k] * nnz[k];
			F77_CALL(dsyrk)("U", "T", &nci, &kj, &one, tmp[k], &kj,
					&one, bVij, &nci);
		    }
		}
	    }
	}
	for (i = 0; i < nf; i++) {
	    if (tmp[i]) Free(tmp[i]);
	    if (ind[i]) Free(ind[i]);
	}
	Free(tmp); Free(nlevs); Free(nnz); Free(ind);
	status[1] = 1;
    }
    return R_NilValue;
}

static double
internal_sigma(SEXP x, int REML)
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym);
    int pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1],
	nobs = INTEGER(GET_SLOT(x, Matrix_ncSym))
	[length(GET_SLOT(x, Matrix_OmegaSym)) + 1];

    lmer_invert(x);
    return 1./(REAL(RXXsl)[pp1*pp1 - 1] *
	       sqrt((double)(REML ? nobs + 1 - pp1 : nobs)));
}

/**
 * Extract the ML or REML conditional estimate of sigma
 *
 * @param x pointer to an lme object
 * @param REML logical scalar - TRUE if REML estimates are requested
 *
 * @return pointer to a numeric scalar
 */
SEXP lmer_sigma(SEXP x, SEXP REML)
{
    return ScalarReal(internal_sigma(x, asLogical(REML)));
}

/**
 * Calculate the length of the parameter vector (historically called "coef"
 * even though these are not coefficients).
 *
 * @param nf number of factors
 * @param nc number of columns in the model matrices for each factor
 *
 * @return total length of the coefficient vector
 */
static R_INLINE
int coef_length(int nf, const int nc[])
{
    int i, ans = 0;
    for (i = 0; i < nf; i++) ans += (nc[i] * (nc[i] + 1))/2;
    return ans;
}

static double *
internal_coef(SEXP x, int ptyp, double ans[])
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), vind;

    vind = 0;			/* index in ans */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1;
	if (nci == 1) {
	    double dd = REAL(VECTOR_ELT(Omega, i))[0];
	    ans[vind++] = ptyp ? ((ptyp == 1) ? log(dd) : 1./dd) : dd;
	} else {
	    if (ptyp) {	/* L log(D) L' factor of Omega[,,i] */
		int j, k, ncisq = nci * nci;
		double *tmp = Memcpy(Calloc(ncisq, double),
				     REAL(VECTOR_ELT(Omega, i)), ncisq);
		F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
		if (j)		/* should never happen */
		    error(_("DPOTRF returned error code %d on Omega[[%d]]"),
			  j, i+1);
		for (j = 0; j < nci; j++) {
		    double diagj = tmp[j * ncip1];
		    ans[vind++] = (ptyp == 1) ? (2. * log(diagj)) :
			1./(diagj * diagj);
		    for (k = j + 1; k < nci; k++) {
			tmp[j + k * nci] /= diagj;
		    }
		}
		for (j = 0; j < nci; j++) {
		    for (k = j + 1; k < nci; k++) {
			ans[vind++] = tmp[j + k * nci];
		    }
		}
		Free(tmp);
	    } else {		/* upper triangle of Omega[,,i] */
		int j, k, odind = vind + nci;
		double *omgi = REAL(VECTOR_ELT(Omega, i));

		for (j = 0; j < nci; j++) {
		    ans[vind++] = omgi[j * ncip1];
		    for (k = j + 1; k < nci; k++) {
			ans[odind++] = omgi[k*nci + j];
		    }
		}
		vind = odind;
	    }
	}
    }
    return ans;
}

/**
 * Extract parameters from the Omega matrices.  These aren't
 * "coefficients" but the extractor is called coef for historical
 * reasons.  Within each group these values are in the order of the
 * diagonal entries first then the strict upper triangle in row
 * order.
 * 
 * The parameters can be returned in three forms:
 *   0 - nonlinearly constrained - elements of the relative precision matrix
 *   1 - unconstrained - from the LDL' decomposition - logarithms of
 *       the diagonal elements of D
 *   2 - box constrained - also from the LDL' decomposition - inverses
 *       of the diagonal elements of D
 *
 * @param x pointer to an lme object
 * @param pType pointer to an integer scalar indicating the form of the 
 *        parameters to be returned.
 *
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP lmer_coef(SEXP x, SEXP pType)
{
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	nf = LENGTH(GET_SLOT(x, Matrix_OmegaSym));
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));

    internal_coef(x, asInteger(pType), REAL(val));
    UNPROTECT(1);
    return val;
}

static
void internal_coefGets(SEXP x, const double cc[], int ptyp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	cind, i, nf = length(Omega);

    cind = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	if (nci == 1) {
	    double dd = cc[cind++];
	    REAL(VECTOR_ELT(Omega, i))[0] =
		ptyp ? ((ptyp == 1) ? exp(dd) : 1./dd) : dd;
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		j, k,
		ncip1 = nci + 1,
		ncisq = nci * nci;
	    double
		*omgi = REAL(VECTOR_ELT(Omega, i));
	    if (ptyp) {
		double *tmp = Calloc(ncisq, double),
		    diagj, one = 1., zero = 0.;

		AZERO(omgi, ncisq);
		for (j = 0; j < nci; j++) {
		    double dd = cc[cind++];
		    tmp[j * ncip1] = diagj =
			(ptyp == 1) ? exp(dd/2.) : sqrt(1./dd);
		    for (k = j + 1; k < nci; k++) {
			tmp[k*nci + j] = cc[odind++] * diagj;
		    }
		}
		F77_CALL(dsyrk)("U", "T", &nci, &nci, &one,
				tmp, &nci, &zero, omgi, &nci);
		Free(tmp);
	    } else {
		for (j = 0; j < nci; j++) {
		    omgi[j * ncip1] = cc[cind++];
		    for (k = j + 1; k < nci; k++) {
			omgi[k*nci + j] = cc[odind++];
		    }
		}
	    }
	    cind = odind;
	}
    }
    status[0] = status[1] = 0;
}

/**
 * Assign the upper triangles of the Omega matrices according to a
 * vector of parameters.
 *
 * @param x pointer to an lme object
 * @param coef pointer to an numeric vector of appropriate length
 * @param pType pointer to an integer scalar 
 *
 * @return R_NilValue
 */
SEXP lmer_coefGets(SEXP x, SEXP coef, SEXP pType)
{
    int clen = coef_length(LENGTH(GET_SLOT(x, Matrix_flistSym)),
			   INTEGER(GET_SLOT(x, Matrix_ncSym)));   
    if (LENGTH(coef) != clen || !isReal(coef))
	error(_("coef must be a numeric vector of length %d"), clen);
    internal_coefGets(x, REAL(coef), asInteger(pType));
    return x;
}

static double*
internal_fixef(SEXP x, double beta[])
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym);
    int j, pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1];
    int p = pp1 - 1;
    double nryyinv;		/* negative ryy-inverse */

    lmer_invert(x);
    Memcpy(beta, REAL(RXXsl) + pp1 * (pp1 - 1), p);
    nryyinv = -REAL(RXXsl)[pp1*pp1 - 1];
    for (j = 0; j < p; j++) beta[j] /= nryyinv;
    return beta;
}

/**
 * Extract the conditional estimates of the fixed effects
 *
 * @param x Pointer to an mer object
 *
 * @return a numeric vector containing the conditional estimates of
 * the fixed effects
 */
SEXP lmer_fixef(SEXP x)
{
    SEXP cnames = GET_SLOT(x, Matrix_cnamesSym);
    int *dims = INTEGER(getAttrib(GET_SLOT(x, Matrix_RZXSym), R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	nf = LENGTH(cnames) - 1;
    int i, p = dims[1] - 1;
    SEXP ans, nms, cnms;
		
    if (nc[nf] < 0 || !p) return(allocVector(REALSXP, 0));
    ans = PROTECT(allocVector(REALSXP, p));
    internal_fixef(x, REAL(ans));
    nms = PROTECT(allocVector(STRSXP, p));
    cnms = VECTOR_ELT(cnames, nf);
    for (i = 0; i < p; i++) SET_STRING_ELT(nms, i, STRING_ELT(cnms, i));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/**
 * Extract the conditional modes of the random effects.
 *
 * @param x Pointer to an lme object
 *
 * @return a list of matrices containing the conditional modes of the
 * random effects
 */
SEXP lmer_ranef(SEXP x)
{
    SEXP RZXP = GET_SLOT(x, Matrix_RZXSym),
	cnames = GET_SLOT(x, Matrix_cnamesSym),
	flist = GET_SLOT(x, Matrix_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ii, jj,
	nf = length(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double
	*b = REAL(RZXP) + dims[0] * (dims[1] - 1),
	nryyinv;		/* negative ryy-inverse */

    lmer_invert(x);
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(flist, R_NamesSymbol)));
    nryyinv = -REAL(GET_SLOT(x, Matrix_RXXSym))[dims[1] * dims[1] - 1];
    for (i = 0; i < nf; i++) {
	SEXP nms, rnms = getAttrib(VECTOR_ELT(flist, i), R_LevelsSymbol);
	int nci = nc[i], mi = length(rnms);
	double *bi = b + Gp[i], *mm;

	SET_VECTOR_ELT(val, i, allocMatrix(REALSXP, mi, nci));
	setAttrib(VECTOR_ELT(val, i), R_DimNamesSymbol, allocVector(VECSXP, 2));
	nms = getAttrib(VECTOR_ELT(val, i), R_DimNamesSymbol);
	SET_VECTOR_ELT(nms, 0, duplicate(rnms));
	SET_VECTOR_ELT(nms, 1, duplicate(VECTOR_ELT(cnames, i)));
	mm = REAL(VECTOR_ELT(val, i));
	for (jj = 0; jj < nci; jj++)
	    for(ii = 0; ii < mi; ii++)
		mm[ii + jj * mi] = bi[jj + ii * nci]/nryyinv;
    }
    UNPROTECT(1);
    return val;
}

/**
 * Fill in four symmetric matrices for each level, providing the
 * information to generate the gradient or the ECME step.  The four
 * matrices are
 *  1) -m_i\bOmega_i^{-1}
 *  2) \bB_i\bB_i\trans
 *  3) \tr\left[\der_{\bOmega_i}\bOmega\left(\bZ\trans\bZ+\bOmega\right)\inv\right]
 *  4) The term added to 3) to get \tr\left[\der_{\bOmega_i}\bOmega\vb\right]
 *
 * @param x pointer to an mer object
 * @param val pointer to a list of matrices of the correct sizes
 *
 * @return val
 */
/* static */
SEXP lmer_firstDer(SEXP x, SEXP val)
{
    SEXP bVarP = GET_SLOT(x, Matrix_bVarSym),
	OmegaP = GET_SLOT(x, Matrix_OmegaSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym);
    int *dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	i, nf = length(OmegaP), p = dims[1] - 1;
    double *RZX = REAL(RZXP),
	*b = REAL(RZXP) + dims[0] * p;

    lmer_invert(x);
    /* FIXME: Why is this loop run backwards?
       It appears it could run forwards. */
    for (i = nf - 1; i >= 0; i--) {
	SEXP bVPi = VECTOR_ELT(bVarP, i);
	int *ddims = INTEGER(getAttrib(bVPi, R_DimSymbol)), j, k;
	int nci = ddims[0];
	int ncisqr = nci * nci, RZXrows = Gp[i + 1] - Gp[i];
	int nlev = RZXrows/nci;
	double *RZXi = RZX + Gp[i], *bVi = REAL(bVPi),
	    *bi = b + Gp[i], *mm = REAL(VECTOR_ELT(val, i)),
	    *tmp = Memcpy(Calloc(ncisqr, double),
			  REAL(VECTOR_ELT(OmegaP, i)), ncisqr),
	    dlev = (double) nlev,
	    one = 1., zero = 0.;

 	if (nci == 1) {
	    int ione = 1;
 	    mm[0] = ((double) nlev)/tmp[0];
 	    mm[1] = F77_CALL(ddot)(&nlev, bi, &ione, bi, &ione);
	    mm[2] = 0.;
	    for (k = 0; k < nlev; k++) mm[2] += bVi[k];
	    mm[3] = 0.;
  	    for (j = 0; j < p; j++) {
  		mm[3] += F77_CALL(ddot)(&RZXrows, RZXi + j * dims[0], &ione,
					RZXi + j * dims[0], &ione);
  	    }
 	} else {
	    AZERO(mm, 4 * ncisqr);
	    F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
	    if (j)
		error(_("Omega[[%d]] is not positive definite"), i + 1);
	    F77_CALL(dtrtri)("U", "N", &nci, tmp, &nci, &j);
	    if (j)
		error(_("Omega[[%d]] is not positive definite"), i + 1);
	    F77_CALL(dsyrk)("U", "N", &nci, &nci, &dlev, tmp, &nci,
			    &zero, mm, &nci);
	    mm += ncisqr;	/* \bB_i term */
	    F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one, bi, &nci,
			    &zero, mm, &nci);
	    mm += ncisqr;     /* Sum of diagonal blocks of the inverse
			       * (Z'Z+Omega)^{-1} */
	    for (j = 0; j < ncisqr; j++) {
		for (k = 0; k < nlev; k++) mm[j] += bVi[j + k*ncisqr];
	    }
	    mm += ncisqr;	/* Extra term for \vb */
	    for (j = 0; j < p; j++) {
		F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one,
				RZXi + j * dims[0], &nci,
				&one, mm, &nci);
	    }
	}
	Free(tmp);
    }
    return val;
}

/**
 * Return a length nf list of arrays of dimension (nci, nci, 4).  The
 * values of these arrays are assigned in lmer_firstDer.
 *
 * @param nf number of factors
 * @param nc vector of number of columns per factor
 *
 * @return pointer to a list of REAL arrays
 */
static
SEXP EM_grad_array(int nf, const int nc[])
{
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    int i;

    for (i = 0; i < nf; i++) {
	SET_VECTOR_ELT(val, i, alloc3Darray(REALSXP, nc[i], nc[i], 4));
    }
    UNPROTECT(1);
    return val;
}

/**
 * Fill in the 4-dimensional vector of linear combinations of the
 * firstDer array according to whether ECME steps or the gradient are
 * needed and to whether or not REML is being used.
 *
 * @param cc coefficient vector to be filled in
 * @param EM non-zero for ECME steps, zero for gradient
 * @param REML non-zero for REML, zero for ML
 * @param ns ns[0] is p+1, ns[1] is n
 *
 * @return cc with the coefficients filled in
 */
static R_INLINE
double *EM_grad_lc(double *cc, int EM, int REML, int ns[])
{
    cc[0] = EM ? 0. : -1.;
    cc[1] = (double)(ns[1] - (REML ? ns[0] - 1 : 0));
    cc[2] = 1.;
    cc[3] = REML ? 1. : 0.;
    return cc;
}

/**
 * Print the verbose output in the ECME iterations
 *
 * @param x pointer to an ssclme object
 * @param iter iteration number
 * @param REML non-zero for REML, zero for ML
 * @param firstDer arrays for calculating ECME steps and the first
 * derivative
 * @param val Pointer to a list of arrays to receive the calculated
 * values
 */
static void
EMsteps_verbose_print(SEXP x, int iter, int REML, SEXP firstDer)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ifour = 4, ii, ione = 1, jj, nf = LENGTH(Omega);
    double
	*cc = EM_grad_lc(Calloc(4, double), 0, REML, nc + nf),
	*dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	one = 1., zero = 0.;

    lmer_factor(x);
    if (iter == 0) Rprintf("  EM iterations\n");
    Rprintf("%3d %.3f", iter, dev[REML ? 1 : 0]);
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1, ncisqr = nci * nci;
	double
	    *Omgi = REAL(VECTOR_ELT(Omega, i)),
	    *Grad = Calloc(ncisqr, double);

				/* diagonals */
	Rprintf(" (%#8g", Omgi[0]);
	for (jj = 1; jj < nci; jj++) {
	    Rprintf(" %#8g", Omgi[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++)
		Rprintf(" %#8g", Omgi[ii + jj * nci]);
				/* Evaluate and print the gradient */
	F77_CALL(dgemv)("N", &ncisqr, &ifour, &one,
			REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			cc, &ione, &zero, Grad, &ione);
	Rprintf(":%#8.3g", Grad[0]);
				
	for (jj = 1; jj < nci; jj++) { /* diagonals */
	    Rprintf(" %#8.3g", Grad[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++)
		Rprintf(" %#8.3g", Grad[ii + jj * nci]);
	Rprintf(")");
	Free(Grad);
    }
    Rprintf("\n");
    Free(cc);
}

static void
internal_ECMEsteps(SEXP x, int nEM, int verb)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	flist = GET_SLOT(x, Matrix_flistSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),
		       "REML"),
	i, ifour = 4, info, ione = 1, iter,
	nf = LENGTH(Omega);
    double
	*cc = EM_grad_lc(Calloc(4, double), 1, REML, nc + nf),
	zero = 0.0;
    SEXP firstDer = PROTECT(EM_grad_array(nf, nc));

    lmer_firstDer(x, firstDer);
    if (verb)
	EMsteps_verbose_print(x, 0, REML, firstDer);
    for (iter = 0; iter < nEM; iter++) {
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], ncisqr = nci * nci;
	    double *Omgi = REAL(VECTOR_ELT(Omega, i)),
		mult = 1./
		((double) length(getAttrib(VECTOR_ELT(flist, i),
				 R_LevelsSymbol)));

	    F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &mult,
			    REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			    cc, &ifour, &zero, Omgi, &ncisqr);
	    F77_CALL(dpotrf)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error(_("DPOTRF in ECME update gave code %d"),
		      info);
	    F77_CALL(dpotri)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error(_("Matrix inverse in ECME update gave code %d"), info);
	}
	status[0] = status[1] = 0;
	lmer_firstDer(x, firstDer);
	if (verb)
	    EMsteps_verbose_print(x, iter + 1, REML, firstDer);
    }
    lmer_factor(x);
    Free(cc);
    UNPROTECT(1);
}

/**
 * Perform ECME steps for the REML or ML criterion.
 *
 * @param x pointer to an mer object
 * @param nsteps pointer to an integer scalar - the number of ECME
 * steps to perform
 * @param Verbp pointer to a logical scalar indicating verbose output
 *
 * @return R_NilValue
 */
SEXP lmer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp)
{
    internal_ECMEsteps(x, asInteger(nsteps), asLogical(Verbp));
    return R_NilValue;
}

/** 
 * Evaluate the gradient vector
 * 
 * @param x Pointer to an lmer object
 * @param pType Pointer to an integer indicator of the
 * parameterization being used
 * 
 * @return pointer to a gradient vector
 */
SEXP lmer_gradient(SEXP x, SEXP pType)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	dind, i, ifour = 4, info, ione = 1, nf = length(Omega),
	odind, ptyp = asInteger(pType);
    SEXP
	firstDer = lmer_firstDer(x, PROTECT(EM_grad_array(nf, nc))),
	val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double
	*cc = EM_grad_lc(Calloc(4, double), 0,
			 !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),
				 "REML"), nc + nf),
	one = 1.0, zero = 0.0;

    dind = 0;			/* index into val for diagonals */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncisqr = nci * nci;
	double
	    *Omgi = REAL(VECTOR_ELT(Omega, i)),
	    *tmp = Calloc(ncisqr, double);

	F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &one,
			REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			cc, &ifour, &zero, tmp, &ncisqr);
	if (nci == 1) {
	    REAL(val)[dind++] =
		(ptyp?((ptyp == 1)?Omgi[0]: -Omgi[0] * Omgi[0]) : 1) * tmp[0];
	} else {
	    int ii, j, ncip1 = nci + 1;

	    odind = dind + nci; /* index into val for off-diagonals */
	    if (ptyp) {
		double *chol = Memcpy(Calloc(ncisqr, double),
				      REAL(VECTOR_ELT(Omega, i)), ncisqr),
		    *tmp2 = Calloc(ncisqr, double);

		/* Overwrite the gradient with respect to positions in
		 * Omega[[i]] by the gradient with respect to the
		 * unconstrained parameters.*/

		F77_CALL(dpotrf)("U", &nci, chol, &nci, &info);
		if (info)
		    error(_("Omega[[%d]] is not positive definite"), i + 1);
		/* tmp2 := chol %*% tmp using only upper triangle of tmp */
		F77_CALL(dsymm)("R", "U", &nci, &nci, &one, tmp, &nci,
				chol, &nci, &zero, tmp2, &nci);
		/* full symmetric product gives diagonals */
		F77_CALL(dtrmm)("R", "U", "T", "N", &nci, &nci, &one, chol,
				&nci, Memcpy(tmp, tmp2, ncisqr), &nci);
		/* overwrite upper triangle with gradients for L' */
		for (ii = 1; ii < nci; ii++) {
		    for (j = 0; j < ii; j++) {
			tmp[j + ii*nci] = chol[j*ncip1] * tmp2[j + ii*nci];
			tmp[ii + j*nci] = 0.;
		    }
		}
		if (ptyp > 1)
		    for (ii = 0; ii < nci; ii++) {
			int ind = ii * ncip1;
			double sqrtd = chol[ind];
			tmp[ind] *= -(sqrtd*sqrtd);
		    }
	    }
	    for (j = 0; j < nci; j++) {
		REAL(val)[dind + j] = tmp[j * ncip1];
		for (ii = 0; ii < j; ii++) /* offdiagonals count twice */
		    REAL(val)[odind++] = 2. * tmp[ii + j * nci];
	    }
	    dind = odind;
	}
	Free(tmp);
    }
    UNPROTECT(2);
    Free(cc);
    return val;
}

/**
 * Fill in five symmetric matrices, providing the information to
 * generate the Hessian.
 *
 * @param x pointer to an mer object
 *
 * @return an array consisting of five symmetric faces
 */
/* static */
SEXP lmer_secondDer(SEXP x)
{
    SEXP
	D = GET_SLOT(x, Matrix_DSym),
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	val;
    int *dRZX = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	Q, Qsqr, RZXpos, facepos,
	i, ione = 1, j, nf = length(Omega), p = dRZX[1] - 1, pos;
    SEXP
	firstDer = lmer_firstDer(x, PROTECT(EM_grad_array(nf, nc)));
    double
	*RZX = REAL(RZXP),
	*b = REAL(RZXP) + dRZX[0] * p,
	*bbface,		/* vec of second faces of firstDer elts */
	one = 1.,
	zero = 0.;

    Q = 0;			/* number of rows and columns in the result */
    for (i = 0; i < nf; i++) Q += nc[i] * nc[i];
    Qsqr = Q * Q;
    bbface = Calloc(Q, double);
    val = PROTECT(alloc3Darray(REALSXP, Q, Q, 5));
    AZERO(REAL(val), Qsqr * 5);

    pos = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int ncisqr = nci * nci;
	double *fDi = REAL(VECTOR_ELT(firstDer, i)),
	    mult = 1./((double)(Gp[i + 1] - Gp[i])/nci);

	Memcpy(bbface + pos, fDi + ncisqr, ncisqr);
	/* outer product of the third face of firstDer on the diagonal
	 * of the third face of val */
	F77_CALL(dsyr)("U", &ncisqr, &mult, fDi + 2 * ncisqr, &ione,
		       REAL(val) + 2 * Qsqr + pos * Q, &Q);
	pos += ncisqr;
    }
				/* fifth face is outer product of bbface */
    F77_CALL(dsyr)("U", &Q, &one, bbface, &ione, REAL(val) + 4 * Qsqr, &Q);
				/* fourth face from \bb\trans\der\vb\der\bb */
    AZERO(REAL(val) + 3 * Qsqr, Qsqr); /* zero accumulator */
    RZXpos = 0;
    facepos = 0;
    for (i = 0; i < nf; i++) {
	int ii, jj, nci = nc[i];
	int ncisqr = nci * nci, nctp = nci * p, 
	    nlev = (Gp[i + 1] - Gp[i])/nci;
	int maxpq = (p > nci) ? p : nci;
	double
	    *Di = REAL(VECTOR_ELT(D, i)),
	    *arr = Calloc(ncisqr * maxpq, double), /* tmp 3Darray */
	    *face = REAL(val) + 3 * Qsqr,
	    *mat = Calloc(nci * maxpq, double); /* tmp matrix */

	for (j = 0; j < nlev; j++) {
	    F77_CALL(dgemm)("T", "T", &p, &nci, &nci,
			    &one, RZX + j * nci, dRZX, Di + j * ncisqr, &nci,
			    &zero, mat, &p);
	    F77_CALL(dgemm)("N", "N", &nctp, &nci, &ione,
			    &one, mat, &nctp, b + j * nci, &ione,
			    &zero, arr, &nctp);
	    F77_CALL(dsyrk)("U", "T", &ncisqr, &p, &one, arr, &p,
			    &one, face + facepos, &Q);
				/* Add the D_{i,j}^{-T/2} term */
	    Memcpy(mat, Di + j * ncisqr, ncisqr);
	    for (jj = 1; jj < nci; jj++) { /* transpose mat */
		for (ii = 0; ii < jj; ii++) {
		    mat[jj + ii * nci] = mat[ii + jj * nci];
		    mat[ii + jj * nci] = 0.;
		}
	    }
	    F77_CALL(dgemm)("N", "N", &ncisqr, &nci, &ione,
			    &one, mat, &ncisqr, b + j * nci, &ione,
			    &zero, arr, &ncisqr);
	    /* FIXME: Next call could be dsyr (it's rank one). */
	    F77_CALL(dsyrk)("U", "T", &ncisqr, &nci, &one, arr, &nci,
			    &one, face + facepos, &Q);

	}
	RZXpos += nci * nlev;
	facepos += ncisqr;
	Free(arr); Free(mat);
    }
    UNPROTECT(2);
    Free(bbface);
    return val;
}

/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static double*
internal_symmetrize(double *a, int nc)
{
    int i, j;

    for (i = 1; i < nc; i++)
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

/**
 * Return the unscaled variances
 *
 * @param x pointer to an lmer object
 *
 * @return a list similar to the Omega list with the unscaled variances
 */
SEXP lmer_variances(SEXP x)
{
    SEXP Omg = PROTECT(duplicate(GET_SLOT(x, Matrix_OmegaSym)));
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omg);

    for (i = 0; i < nf; i++) {
	double *mm = REAL(VECTOR_ELT(Omg, i));
	int j, nci = nc[i];

	F77_CALL(dpotrf)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error(_("DPOTRF returned error code %d on Omega[%d]"),
		  j, i + 1);
	F77_CALL(dpotri)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error(_("DTRTRI returned error code %d on Omega[%d]"),
		  j, i + 1);
	internal_symmetrize(mm, nci);
    }
    UNPROTECT(1);
    return Omg;
}

/** 
 * Add the contribution from the random effects to the fitted values
 * 
 * @param flist pointer to a list of factors
 * @param mmats pointer to a list of model matrices
 * @param b pointer to a list of random effects matrices
 * @param nc number of columns per model matrix
 * @param val fitted values that are incremented
 * 
 * @return val (updated)
 */
static double*
fitted_ranef(SEXP flist, SEXP mmats, SEXP b,
	     const int nc[], double val[])
{
    int i, ione = 1, nf = LENGTH(flist),
	nobs = LENGTH(VECTOR_ELT(flist, 0));
    double one = 1;
    
    for (i = 0; i < nf; i++) {
	SEXP bi = VECTOR_ELT(b, i);
	int mi = INTEGER(getAttrib(bi, R_DimSymbol))[0];
	int *ff = INTEGER(VECTOR_ELT(flist, i)), j, nci = nc[i];
	double *mm = REAL(VECTOR_ELT(mmats, i));

	for (j = 0; j < nobs; ) {
	    int nn = 1, lev = ff[j];
	    /* check for adjacent rows with same factor level */
	    while ((j + nn) < nobs && ff[j + nn] == lev) nn++; 
	    F77_CALL(dgemm)("N", "T", &nn, &ione, &nci,
			    &one, mm + j, &nobs,
			    REAL(bi) + (lev - 1), &mi,
			    &one, &val[j], &nobs);
	    j += nn;
	}
    }
    return val;
}
    
/** 
 * Calculate the fitted values.
 * 
 * @param x pointer to an lmer object
 * @param mmats list of model matrices
 * @param useRand logical scalar indicating if the random
 * effects should be used
 * @param val array to hold the fitted values
 * 
 * @return pointer to a numeric array of fitted values
 */
static double*
internal_fitted(SEXP x, SEXP mmats, int useRand, double val[])
{
    SEXP flist = GET_SLOT(x, Matrix_flistSym),
	fixef = PROTECT(lmer_fixef(x));
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)), ione = 1,
	nf = LENGTH(flist), nobs = length(VECTOR_ELT(flist, 0)),
	p = LENGTH(fixef);
    double one = 1.0, zero = 0.0;

    if (p) F77_CALL(dgemm)("N", "N", &nobs, &ione, &p, &one,
			   REAL(VECTOR_ELT(mmats, nf)), &nobs,
			   REAL(fixef), &p, &zero, val, &nobs);
    else AZERO(val, nobs);
    if (useRand) {
	fitted_ranef(flist, mmats, PROTECT(lmer_ranef(x)),
		     nc, val);
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Return the fitted values as an SEXP
 * 
 * @param x pointer to an lmer object
 * @param mmats list of model matrices
 * @param useRf pointer to a logical scalar indicating if the random
 * effects should be used
 * 
 * @return pointer to a numeric array of fitted values
 */

SEXP lmer_fitted(SEXP x, SEXP mmats, SEXP useRf)
{
    int nobs = LENGTH(VECTOR_ELT(GET_SLOT(x, Matrix_flistSym), 0));
    SEXP ans = PROTECT(allocVector(REALSXP, nobs));

    internal_fitted(x, mmats, asLogical(useRf), REAL(ans));
    UNPROTECT(1);
    return ans;
}
				 
/** 
 * Apply weights to a list of model matrices
 * 
 * @param MLin Input list of model matrices
 * @param wts Weights to be applied
 * @param adjst Working residual to be appended to the last model matrix
 * @param n number of observations
 * @param MLout Output list of model matrices
 */
static void
internal_weight_list(SEXP MLin, double *wts, double *adjst, int n,
		     SEXP MLout)
{
    int i, j, nf = LENGTH(MLin);
    SEXP lastM;

    for (i = 0; i < nf; i++) {
	SEXP Min = VECTOR_ELT(MLin, i),
	    Mout = VECTOR_ELT(MLout, i);
	int *din, *dout, k, nc;

	if (!(isMatrix(Min) && isReal(Min)))
	    error(_("component %d of MLin is not a numeric matrix"), i + 1);
	din = INTEGER(getAttrib(Min, R_DimSymbol));
	nc = din[1];
	if (din[0] != n)
	    error(_("component %d of MLin has %d rows, expected %d"), i + 1,
		  din[0], n);
	if (!(isMatrix(Mout) && isReal(Mout)))
	    error(_("component %d of MLout is not a numeric matrix"), i + 1);
	dout = INTEGER(getAttrib(Mout, R_DimSymbol));
	if (dout[0] != n)
	    error(_("component %d of MLout has %d rows, expected %d"), i + 1,
		  dout[0], n);
	if (dout[1] != nc)
	    error(_("component %d of MLout has %d columns, expected %d"), i + 1,
		  dout[1], nc);
	for (k = 0; k < nc; k++) {
	    for (j = 0; j < n; j++) {
		REAL(Mout)[j + k * n] = REAL(Min)[j + k * n] * wts[j];
	    }
	}
    }
    lastM = VECTOR_ELT(MLout, nf - 1);
    j = INTEGER(getAttrib(lastM, R_DimSymbol))[1] - 1;
    for (i = 0; i < n; i++)
	REAL(lastM)[j*n + i] = adjst[i] * wts[i];
}    

/** 
 * Find a variable of a given name in a given environment and check
 * that its length and mode are correct.
 * 
 * @param rho Environment in which to find the variable
 * @param nm Name of the variable to find
 * @param mode Desired mode
 * @param len Desired length
 * 
 * @return 
 */
static
SEXP find_and_check(SEXP rho, SEXP nm, SEXPTYPE mode, int len)
{    
    SEXP ans;
    if (R_NilValue == PROTECT(ans = findVarInFrame(rho, nm)))
	error(_("environment `rho' must contain an object `%s'"),
	      CHAR(PRINTNAME(nm)));
    if (TYPEOF(ans) != mode)
	error(_("object `%s' of incorrect type"),
	      CHAR(PRINTNAME(nm)));
    if (len && LENGTH(ans) != len)
	error(_("object `%s' must be of length `%d'"),
	      CHAR(PRINTNAME(nm)), len);
    UNPROTECT(1);
    return ans;
}

/** 
 * Evaluate an expression in an environment, check that the length and
 * mode are as expected and store the result.
 * 
 * @param fcn expression to evaluate
 * @param rho environment in which to evaluate it
 * @param vv position to store the result
 * 
 * @return vv with new contents
 */
static
SEXP eval_check_store(SEXP fcn, SEXP rho, SEXP vv)
{
    SEXP v = PROTECT(eval(fcn, rho));
    if (TYPEOF(v) != TYPEOF(vv) || LENGTH(v) != LENGTH(vv))
	error(_("fcn produced mode %d, length %d - wanted mode %d, length %d"),
	      TYPEOF(v), LENGTH(v), TYPEOF(vv), LENGTH(vv));
    switch (TYPEOF(v)) {
    case LGLSXP:
	Memcpy(LOGICAL(vv), LOGICAL(v), LENGTH(vv));
	break;
    case INTSXP:
	Memcpy(INTEGER(vv), INTEGER(v), LENGTH(vv));
	break;
    case REALSXP:
	Memcpy(REAL(vv), REAL(v), LENGTH(vv));
	break;
    default:
	error(_("invalid type for eval_check_store"));
    }
    UNPROTECT(1);
    return vv;
}

/** 
 * Evaluate an expression in an environment, check that the length and
 * mode are as expected and return the result.
 * 
 * @param fcn expression to evaluate
 * @param rho environment in which to evaluate it
 * @param mode desired mode
 * @param len desired length
 * 
 * @return evaluated expression
 */
static SEXP
eval_check(SEXP fcn, SEXP rho, SEXPTYPE mode, int len) {
    SEXP v = PROTECT(eval(fcn, rho));
    if (TYPEOF(v) != mode || LENGTH(v) != len)
	error(_("fcn produced mode %d, length %d - wanted mode %d, length %d"),
	      TYPEOF(v), LENGTH(v), mode, len);
    UNPROTECT(1);
    return v;
}

typedef struct glmer_struct
{
    SEXP cv;         /* control values */
    SEXP mer;	     /* mixed-effects representation */
    SEXP rho;        /* where to evaluate the calls */
    SEXP y;          /* response vector (REALSXP) */
    SEXP eta;        /* linear predictor */
    SEXP mu;         /* mean vector */
    SEXP wts;        /* prior weights */
    SEXP offset;     /* offset for original glm problem */
    SEXP off;        /* offset for random effects only */
    SEXP unwtd;      /* unweighted model matrices */
    SEXP wtd;        /* weighted model matrices */
    SEXP linkinv;    /* expression for inverse link evaluation */
    SEXP mu_eta;     /* expression for dmu/deta evaluation */
    SEXP LMEopt;     /* expression for LME optimization */
    SEXP dev_resids; /* expression for deviance residuals */
    SEXP var;        /* expression for variance evaluation */
    SEXP x;          /* fixed-effects model matrix */
    int n;	     /* length of the response vector */
    int p;	     /* length of fixed effects vector */
    int nf;	     /* number of grouping factors */
    int npar;        /* total length of the parameter */
    int niterEM;     /* default number of ECME iterations */
    int EMverbose;   /* logical indicator */
    int maxiter;     /* maximum number of IRLS iterations */
    double tol;      /* convergence tolerance for IRLS iterations */
} glmer_struct, *GlmerStruct;

/** 
 * Return an external pointer object to a GlmerStruct created in
 * environment rho
 * 
 * @param rho An environment
 * 
 * @return An external pointer to a GlmerStruct
 */
SEXP glmer_init(SEXP rho) {
    GlmerStruct GS;
    int *dims;
    SEXP wtdSym = install("wtd"), offSym = install("off");
    
    GS = (GlmerStruct) Calloc(1, glmer_struct);
    if (!isEnvironment(rho))
	error(_("`rho' must be an environment"));
    GS->rho = rho;
    GS->y = find_and_check(rho, install("y"), REALSXP, 0);
    if ((GS->n = LENGTH(GS->y)) < 1)
	error(_("`%s' must be a nonempty, numeric vector"),
	      "y");
    GS->mer = find_and_check(rho, install("mer"), VECSXP, 0);
    GS->nf = LENGTH(GET_SLOT(GS->mer, Matrix_flistSym));
    GS->cv = find_and_check(rho, install("cv"), VECSXP, 0);
    GS->niterEM = asInteger(Matrix_getElement(GS->cv, "niterEM"));
    GS->EMverbose = asLogical(Matrix_getElement(GS->cv, "EMverbose"));
    GS->tol = asReal(Matrix_getElement(GS->cv, "tolerance"));
    GS->maxiter = asInteger(Matrix_getElement(GS->cv, "maxIter"));

    GS->mu = find_and_check(rho, install("mu"), REALSXP, GS->n);
    GS->offset = find_and_check(rho, install("offset"),
				REALSXP, GS->n);
    defineVar(offSym, duplicate(GS->offset), rho);
    GS->off = find_and_check(rho, offSym, REALSXP, GS->n);
    GS->x = find_and_check(rho, install("x"), REALSXP, 0);
    if (!isMatrix(GS->x))
	error(_("%s must be a model matrix"), "x");
    dims = INTEGER(getAttrib(GS->x, R_DimSymbol));
    if (dims[0] != GS->n)
	error(_("nrow(%s) must be %d"), "x", GS->n);
    GS->p = dims[1];
    GS->npar = dims[1] +
	coef_length(GS->nf, INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)));
    GS->eta = find_and_check(rho, install("eta"),
			     REALSXP, GS->n);
    GS->unwtd = find_and_check(rho, install("mmats"),
			       VECSXP, GS->nf + 1);
    defineVar(wtdSym, duplicate(GS->unwtd), rho);
    GS->wtd = find_and_check(rho, wtdSym, VECSXP, GS->nf + 1);
    GS->wts = find_and_check(rho, install("wts"), REALSXP, GS->n);
    GS->linkinv = find_and_check(rho, install("linkinv"),
				 LANGSXP, 0);
    GS->mu_eta = find_and_check(rho, install("mu.eta"),
				LANGSXP, 0);
    GS->var = find_and_check(rho, install("variance"),
			     LANGSXP, 0);
    GS->LMEopt = find_and_check(rho, install("doLMEopt"),
				LANGSXP, 0);
    GS->dev_resids = find_and_check(rho, install("dev.resids"),
				    LANGSXP, 0);
    return R_MakeExternalPtr(GS, R_NilValue, GS->mer);
}

/** 
 * Release the storage for a GlmerStruct
 * 
 * @param GSp External pointer to a  GlmerStruct
 * 
 * @return R_NilValue
 */
SEXP glmer_finalize(SEXP GSp) {
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    
    Free(GS);
    return R_NilValue;
}

/** 
 * Evaluate new weights and working residuals, apply the weights to a
 * copy of the model matrices and update the mer object.
 * 
 * @param GS a GlmerStruct object
 */
static void
reweight_update(GlmerStruct GS) {
    SEXP dmu_deta, var;
    int i;
    double *w = Calloc(GS->n, double),
	*z = Calloc(GS->n, double);
				/* reweight mer */
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				  REALSXP, GS->n));
    var = PROTECT(eval_check(GS->var, GS->rho,
			     REALSXP, GS->n));
    for (i = 0; i < GS->n; i++) {
	w[i] = REAL(GS->wts)[i] *
	    REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	z[i] = REAL(GS->eta)[i] - REAL(GS->off)[i] +
	    (REAL(GS->y)[i] - REAL(GS->mu)[i]) /
	    REAL(dmu_deta)[i];
    }
    internal_weight_list(GS->unwtd, w, z, GS->n, GS->wtd);
    lmer_update_mm(GS->mer, GS->wtd);
    UNPROTECT(2);
    Free(w); Free(z);
}

/** 
 * Update eta, evaluate the convergence criterion, then copy eta to
 * etaold
 * 
 * @param GS a GlmerStruct object
 * @param etaold previous values of the linear predictors
 * 
 * @return convergence criterion
 */
static double
conv_crit(double etaold[], double eta[], int n) {
    double max_abs_eta = -1, max_abs_diff = -1;
    int i;

    for (i = 0; i < n; i++) {
	double abs_eta, abs_diff;

	abs_eta = fabs(eta[i]);
	if (abs_eta > max_abs_eta) max_abs_eta = abs_eta;
	abs_diff = fabs(eta[i] - etaold[i]);
	if (abs_diff > max_abs_diff) max_abs_diff = abs_diff;
	etaold[i] = eta[i];
    }
    return max_abs_diff / (0.1 + max_abs_eta);
}

/** 
 * Perform the PQL optimization
 * 
 * @param GSp pointer to a GlmerStruct object
 * 
 * @return R_NilValue
 */
SEXP glmer_PQL(SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    int i;
    double *etaold = Memcpy(Calloc(GS->n, double),
			    REAL(GS->eta), GS->n),
	*fitted = Calloc(GS->n, double), crit;

    for (i = 0, crit = GS->tol + 1;
	 i < GS->maxiter && crit > GS->tol; i++) {
	reweight_update(GS);
	if (!i) lmer_initial(GS->mer); /* initialize first fit */
	internal_ECMEsteps(GS->mer, i ? 2 : GS->niterEM,
			   GS->EMverbose);
	eval(GS->LMEopt, GS->rho);
	vecSum(REAL(GS->eta), REAL(GS->off), 
	       internal_fitted(GS->mer, GS->unwtd, 1, fitted),
	       GS->n);
	crit = conv_crit(etaold, REAL(GS->eta), GS->n);
    }
    if (crit > GS->tol)
	warning(_("IRLS iterations for PQL did not converge"));

    Free(etaold); Free(fitted);
    return R_NilValue;
}
    
/** 
 * Establish off, the effective offset for the fixed effects, and
 * iterate to determine the conditional modes.
 * 
 * @param pars parameter vector
 * @param GS a GlmerStruct object
 * 
 * @return An indicator of whether the iterations converged
 */
static int
internal_bhat(GlmerStruct GS, const double fixed[], const double varc[])
{
    SEXP flist = GET_SLOT(GS->mer, Matrix_flistSym);
    int *nc = INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)),
	i, ione = 1;
    double *etaold = Calloc(GS->n, double),
	*fitted = Calloc(GS->n, double),
	crit, one = 1, zero = 0;

    /* skip this step if varc == (double*) NULL */	
    if (varc)			
	internal_coefGets(GS->mer, varc, 2);

    F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
		    REAL(GS->x), &(GS->n), fixed,
		    &ione, &zero, REAL(GS->off), &ione);
    vecIncrement(REAL(GS->off), REAL(GS->offset), GS->n);
    Memcpy(REAL(GS->eta), REAL(GS->off), GS->n);
    Memcpy(etaold, REAL(GS->eta), GS->n);

    for (i = 0, crit = GS->tol + 1;
	 i < GS->maxiter && crit > GS->tol; i++) {
	reweight_update(GS);
	AZERO(fitted, GS->n);
	fitted_ranef(flist, GS->unwtd, lmer_ranef(GS->mer),
		     nc, fitted);
	vecSum(REAL(GS->eta), REAL(GS->off), fitted, GS->n);
	crit = conv_crit(etaold, REAL(GS->eta), GS->n);
    }
    Free(etaold); Free(fitted);
    return (crit > GS->tol) ? 0 : 1;
}

/** 
 * Determine the deviance components associated with each of the
 * levels of a grouping factor at the conditional modes or a value
 * offset from the conditional modes by delb.
 * 
 * @param GS pointer to a GlmerStruct
 * @param b pointer to the conditional modes of the random effects 
 * @param nlev number of levels of the kth grouping factor
 * @param nc number of columns in the model matrix for the kth
 * grouping factor
 * @param k index (0-based) of the grouping factor
 * @param delb vector of length nc giving the changes in the
 * orthonormalized random effects 
 * @param OmgFac Cholesky factor of the inverse of the penalty matrix
 * for this grouping factor 
 * @param bVfac 3-dimensional array holding the factors of the
 * conditional variance-covariance matrix of the random effects 
 * @param devcmp array to hold the deviance components
 * 
 * @return devcmp
 */
static double*
rel_dev_1(GlmerStruct GS, SEXP b, int nlev, int nc, int k,
	  const double delb[], const double OmgFac[],
	  const double bVfac[], double devcmp[])
{
    SEXP devs, flist = GET_SLOT(GS->mer, Matrix_flistSym);
    int *fv = INTEGER(VECTOR_ELT(flist, k)), i, ione = 1,
	j, ntot = nlev * nc;
    double *bb = REAL(VECTOR_ELT(b, k)), *bcp = (double *) NULL;

    AZERO(devcmp, nlev);
    if (delb) {
	double sumsq = 0;
				/* copy the contents of b */
	bcp = Memcpy(Calloc(ntot, double), bb, ntot);
	if (nc == 1) {
	    sumsq = delb[0] * delb[0];
	    for (i = 0; i < nlev; i++) bb[i] += delb[0] * bVfac[i];
	} else {
	    int ncsq = nc * nc;
	    double *tmp = Calloc(nc, double);
	    for (i = 0; i < nlev; i++) {
		Memcpy(tmp, delb, nc);
		F77_CALL(dtrmv)("U", "N", "N", &nc, &(bVfac[i * ncsq]),
				&nc, tmp, &ione);
		for (j = 0; j < nc; j++) bb[i + j * nc] = tmp[j];
	    }
				/* sum of squares of delb */
	    for (j = 0; j < nc; j++) sumsq += delb[j] * delb[j];
	}
	for (i = 0; i < nlev; i++) devcmp[i] = -sumsq;
    }
    Memcpy(REAL(GS->eta), REAL(GS->off), GS->n);
    fitted_ranef(flist, GS->unwtd, b, &nc, REAL(GS->eta));
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0; i < GS->n; i++)
	devcmp[fv[i] - 1] += REAL(devs)[i];
    UNPROTECT(1);
    if (nc == 1) {
	for (i = 0; i < nlev; i++) {
	    double tmp = *OmgFac * bb[i];
	    devcmp[i] += tmp * tmp;
	}
    } else {
	double *tmp = Calloc(nc, double);
	int ione = 1;
	
	for (i = 0; i < nlev; i++) {
	    for (j = 0; j < nc; j++) tmp[j] = bb[i + j * nlev];
	    F77_CALL(dtrmv)("U", "N", "N", &nc, OmgFac, &nc,
			    tmp, &ione);
	    for (j = 0; j < nc; j++) 
		devcmp[i] += tmp[j] * tmp[j];
	}
    }
    if (delb) {
	Memcpy(bb, bcp, ntot);
	Free(bcp);
    }
    return devcmp;
}


/* Gauss-Hermite Quadrature x positions and weights */
static const double
    GHQ_x1[1] = {0},
    GHQ_w1[1] = {1},
    GHQ_x2[1] = {1},
    GHQ_w2[1] = {0.5},
    GHQ_x3[2] = {1.7320507779261, 0},
    GHQ_w3[2] = {0.166666666666667, 0.666666666666667},
    GHQ_x4[2] = {2.3344141783872, 0.74196377160456},
    GHQ_w4[2] = {0.0458758533899086, 0.454124131589555},
    GHQ_x5[3] = {2.85696996497785, 1.35562615677371, 0},
    GHQ_w5[3] = {0.0112574109895360, 0.222075915334214,
		 0.533333317311434},
    GHQ_x6[3] = {3.32425737665988, 1.88917584542184,
		 0.61670657963811},
    GHQ_w6[3] = {0.00255578432527774, 0.0886157433798025,
		 0.408828457274383},
    GHQ_x7[4] = {3.7504396535397, 2.36675937022918,
		 1.15440537498316, 0},
    GHQ_w7[4] = {0.000548268839501628, 0.0307571230436095,
		 0.240123171391455, 0.457142843409801},
    GHQ_x8[4] = {4.14454711519499, 2.80248581332504,
		 1.63651901442728, 0.539079802125417},
    GHQ_w8[4] = {0.000112614534992306, 0.0096352198313359,
		 0.117239904139746, 0.373012246473389},
    GHQ_x9[5] = {4.51274578616743, 3.20542894799789,
		 2.07684794313409, 1.02325564627686, 0},
    GHQ_w9[5] = {2.23458433364535e-05, 0.00278914123744297,
		 0.0499164052656755, 0.244097495561989,
		 0.406349194142045},
    GHQ_x10[5] = {4.85946274516615, 3.58182342225163,
		  2.48432579912153, 1.46598906930182,
		  0.484935699216176},
    GHQ_w10[5] = {4.31065250122166e-06, 0.000758070911538954,
		  0.0191115799266379, 0.135483698910192,
		  0.344642324578594},
    GHQ_x11[6] = {5.18800113558601, 3.93616653976536,
		  2.86512311160915, 1.87603498804787,
		  0.928868981484148, 0},
    GHQ_w11[6] = {8.12184954622583e-07, 0.000195671924393029,
		  0.0067202850336527, 0.066138744084179,
		  0.242240292596812, 0.36940835831095};

static const double
    *GHQ_x[12] = {(double *) NULL, GHQ_x1, GHQ_x2, GHQ_x3, GHQ_x4,
		  GHQ_x5, GHQ_x6, GHQ_x7, GHQ_x8, GHQ_x9, GHQ_x10,
		  GHQ_x11},
    *GHQ_w[12] = {(double *) NULL, GHQ_w1, GHQ_w2, GHQ_w3, GHQ_w4,
		  GHQ_w5, GHQ_w6, GHQ_w7, GHQ_w8, GHQ_w9, GHQ_w10,
		  GHQ_w11};

/** 
 * Compute the approximation to the deviance using adaptive
 * Gauss-Hermite quadrature (AGQ).  When nAGQ == 1 this is the Laplace
 * approximation.
 * 
 * @param pars pointer to a numeric vector of parameters
 * @param GSp pointer to a GlmerStruct object
 * @param nAGQp pointer to a scalar integer representing the number of
 * points in AGQ to use
 * 
 * @return the approximation to the deviance as computed using AGQ
 */
SEXP glmer_devAGQ(SEXP pars, SEXP GSp, SEXP nAGQp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    SEXP bhat, Omega, bVar;
    int i, j, k, nAGQ = asInteger(nAGQp);
    int n2 = (nAGQ + 1)/2;
    double *f0, LaplaceDev = 0, AGQadjst = 0;
	
    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (GS->nf > 1 && nAGQ > 1) {
	warning(_("AGQ not available for multiple grouping factors - using Laplace"));
	nAGQ = 1;
    }
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	return ScalarReal(DBL_MAX);
    bhat = PROTECT(lmer_ranef(GS->mer)); /* forces an inversion of GS->mer */
    bVar = GET_SLOT(GS->mer, Matrix_bVarSym);
    Omega = GET_SLOT(GS->mer, Matrix_OmegaSym);

    for (i = 0; i < GS->nf; i++) {
	int *dims = INTEGER(getAttrib(VECTOR_ELT(bVar, i),
				      R_DimSymbol));
	int nci = dims[0], nlev = dims[2];
	int ncip1 = nci + 1, ncisqr = nci * nci;
	double *omg = Memcpy(Calloc(ncisqr, double),
			     REAL(VECTOR_ELT(Omega, i)), ncisqr),
	    *bvar = Memcpy(Calloc(ncisqr * nlev, double),
			   REAL(VECTOR_ELT(bVar, i)), ncisqr * nlev);

	LaplaceDev = 0;
				/* Calculate difference of determinants */
        F77_CALL(dpotrf)("U", &nci, omg, &nci, &j);
        if (j)
            error(_("Leading %d minor of Omega[[%d]] not positive definite"),
                  j, i + 1);
        for (j = 0; j < nci; j++) { /* nlev * logDet(Omega_i) */
            /* Note: we subtract because Omega has been inverted */
            LaplaceDev -= 2 * nlev * log(omg[j * ncip1]);
        }
        for (k = 0; k < nlev; k++) {
	    double *bVik = bvar + k * ncisqr;
            F77_CALL(dpotrf)("U", &nci, bVik, &nci, &j);
            if (j)
                error(_("Leading %d minor of bVar[[%d]][,,%d] not positive definite"),
                      j, i + 1, k + 1);
            for (j = 0; j < nci; j++) LaplaceDev -= 2 * log(bVik[j * ncip1]);
        }

	f0 = Calloc(nlev, double);
	rel_dev_1(GS, bhat, nlev, nci, i, (double *) NULL,
		  omg, bvar, f0);
	for (k = 0; k < nlev; k++) LaplaceDev += f0[k];
	if (nAGQ > 1) {
	    double *fx = Calloc(nlev, double),
		*rellik = Calloc(nlev, double),
		*delb = Calloc(nci, double);

	    if (nci > 1) error(_("code not yet written"));
	    AZERO(rellik, nlev);	/* zero accumulator */
	    for (k = 0; k < n2; k++) {	
		delb[0] = GHQ_x[nAGQ][k];
		if (delb[0]) {
		    rel_dev_1(GS, bhat, nlev, nci, i, delb,
			      omg, bvar, fx);
		    for (j = 0; j < nlev; j++) {
			rellik[j] += GHQ_w[nAGQ][k] *
			    exp(-(fx[j] - f0[j])/2);
		    }
		    delb[0] *= -1;
		    rel_dev_1(GS, bhat, nlev, nci, i, delb,
			      omg, bvar, fx);
		    for (j = 0; j < nlev; j++) {
			rellik[j] += GHQ_w[nAGQ][k] *
			    exp(-(fx[j] - f0[j])/2);
		    }
		} else {
		    for (j = 0; j < nlev; j++)
			rellik[j] += GHQ_w[nAGQ][k];
		}
	    }
	    for (j = 0; j < nlev; j++)
		AGQadjst -= 2 * log(rellik[j]);
	    Free(fx); Free(rellik);
	}
	Free(f0); Free(omg); Free(bvar);
    }
    UNPROTECT(1);
    return ScalarReal(LaplaceDev + AGQadjst);
}

/** 
 * Evaluate the conditional deviance for a given set of fixed effects.
 * 
 * @param GS Pointer to a GlmerStruct
 * @param fixed value of the fixed effects
 * 
 * @return conditional deviance
 */
static double
random_effects_deviance(GlmerStruct GS, SEXP b)
{
    SEXP devs;
    int i;
    double ans;

    Memcpy(REAL(GS->eta), REAL(GS->off), GS->n);
    fitted_ranef(GET_SLOT(GS->mer, Matrix_flistSym), GS->unwtd, b,
		 INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)), REAL(GS->eta));
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0, ans = 0; i < GS->n; i++) ans += REAL(devs)[i];
    UNPROTECT(1);
    return ans;
}

static void
internal_glmer_ranef_update(GlmerStruct GS, SEXP b)
{
    SEXP bhat, bprop = PROTECT(duplicate(b)), 
	bVar = GET_SLOT(GS->mer, Matrix_bVarSym),
	Omega = GET_SLOT(GS->mer, Matrix_OmegaSym);
    int i, ione = 1, j, k;
    double devr, one = 1;

    bhat = PROTECT(lmer_ranef(GS->mer));
				/* subtract deviance at b */
    devr = -random_effects_deviance(GS, b);
    for (i = 0; i < GS->nf; i++) {
	SEXP Bi = VECTOR_ELT(b, i);
	int *dims = INTEGER(getAttrib(Bi, R_DimSymbol));
	int nlev = dims[0], nci = dims[1];
	int ncisqr = nci * nci, ntot = nlev * nci;
	double *bcopy = Calloc(ntot, double),
	    *bi = REAL(Bi),
	    *bhati = REAL(VECTOR_ELT(bhat, i)),
	    *bpropi = REAL(VECTOR_ELT(bprop, i)),
	    *bVari = REAL(VECTOR_ELT(bVar, i)),
	    *chol = Calloc(ncisqr, double),
	    *delta = Calloc(nci, double),
	    *omgfac = Memcpy(Calloc(ncisqr, double),
			     REAL(VECTOR_ELT(Omega, i)),
			     ncisqr);

				/* subtract quadratic form in
				 * Omega[[i]] at b  */
	F77_CALL(dpotrf)("U", &nci, omgfac, &nci, &j);
	if (j)
	    error(_("Leading %d minor of Omega[[%d]] not positive definite"),
                      j, i + 1);
	Memcpy(bcopy, bi, ntot);
	F77_CALL(dtrmm)("R", "U", "T", "N", &nlev, &nci, &one,
			omgfac, &nci, bcopy, &nlev);
	for (k = 0; k < ntot; k++) devr -= bcopy[k] * bcopy[k];
				/* form bprop and proposal density */
	for (k = 0; k < nlev; k++) {
				/* proposal density at b */
	    for (j = 0; j < nci; j++)
		delta[j] = bi[k + j * nlev] - bhati[k + j * nlev];
	    Memcpy(chol, &(bVari[k * ncisqr]), ncisqr);
	    F77_CALL(dpotrf)("U", &nci, chol, &nci, &j);
	    if (j)
		error(_("Leading %d minor of bVar[[%d]][,,%d] not positive definite"),
                      j, i + 1, k + 1);
	    F77_CALL(dtrsv)("U", "T", "N", &nci, chol, &nci,
			    delta, &ione);
	    for (j = 0; j < nci; j++) {
		double nrm = norm_rand(); /* proposal deviate */
		devr += delta[j] * delta[j] - nrm * nrm;
		delta[j] = nrm;
	    }
				/* scale by Cholesky inverse */
	    F77_CALL(dtrmv)("U", "T", "N", &nci, chol, &nci,
			    delta, &ione);
				/* add mean */
	    for (j = 0; j < nci; j++)
		bpropi[k + j * nlev] = bhati[k + j * nlev] + delta[j];
	}
				/* add quadratic form in
				 * Omega[[i]] at bprop  */
	Memcpy(bcopy, bpropi, ntot);
	F77_CALL(dtrmm)("R", "U", "T", "N", &nlev, &nci, &one,
			omgfac, &nci, bcopy, &nlev);
	for (k = 0; k < ntot; k++) devr += bcopy[k] * bcopy[k];

	Free(bcopy); Free(chol); Free(delta); Free(omgfac);
    }
				/* add deviance at bprop */
    devr += random_effects_deviance(GS, bprop);

    if (unif_rand() < exp(-0.5 * devr))
	for (i = 0; i < GS->nf; i++) { /* copy each face of b */
	    SEXP Bi = VECTOR_ELT(b, i);
	    int *dims = INTEGER(getAttrib(Bi, R_DimSymbol));

	    Memcpy(REAL(Bi), REAL(VECTOR_ELT(bprop, i)),
		   dims[0] * dims[1]);
	}
    
    if (asLogical(Matrix_getElement(GS->cv, "msVerbose"))) {
	double *b0 = REAL(VECTOR_ELT(bprop, 0));
	Rprintf("%5.3f:", exp(-0.5 * devr));
	for (k = 0; k < 5; k++) Rprintf("%#10g ", b0[j]);
	Rprintf("\n");
    }
	
    UNPROTECT(2);
}

/** 
 * Determine the conditional modes and the conditional variance of the
 * random effects given the data and the current fixed effects and
 * variance components.
 *
 * Create a Metropolis-Hasting proposal step from a multivariate
 * normal density centered at bhat with variance-covariance matrix
 * from bVar, determine the acceptance probability and return the
 * current value or the proposed value.
 * 
 * @param GSp pointer to a GlmerStruct
 * @param fixed pointer to a numeric vector of the fixed effects
 * @param varc pointer to a numeric vector of the variance components
 * @param varc pointer to current values of b
 * 
 * @return updated b from the Metropolis-Hastings step
 */
SEXP glmer_ranef_update(SEXP GSp, SEXP fixed, SEXP varc, SEXP b)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    int nvarc = GS->npar - GS->p;

    if (!isReal(fixed) || LENGTH(fixed) != GS->p)
	error(_("`%s' must be a numeric vector of length %d"),
	      "fixed", GS->p);
    if (INTEGER(GET_SLOT(GS->mer, Matrix_ncSym))[GS->nf] > 0)
	error(_("the mer object must be set to skip fixed effects"));
    if (!isReal(varc) || LENGTH(varc) != nvarc)
	error(_("`%s' must be a numeric vector of length %d"),
	      "varc", nvarc);

    GetRNGstate();
    /* Don't check for convergence failure after internal_bhat.
     * It is determining the mean of the proposal density and
     * does not need exact convergence. */
    internal_bhat(GS, REAL(fixed), REAL(varc));
    internal_glmer_ranef_update(GS, b);
    PutRNGstate();

    UNPROTECT(1);
    return b;
}

/** 
 * Evaluate the quadratic form (x-mn)'A'A(x-mn) from the multivariate
 * normal kernel.
 * 
 * @param n dimension of random variate
 * @param mn mean vector
 * @param a upper Cholesky factor of the precision matrix
 * @param lda leading dimension of A
 * @param x vector at which to evaluate the kernel
 * 
 * @return value of the normal kernel
 */
static double
normal_kernel(int n, const double mn[],
	      const double a[], int lda, const double x[])
{
    int i, ione = 1;
    double *tmp = Calloc(n, double), ans;
    
    for (i = 0; i < n; i++) tmp[i] = x[i] - mn[i];
    F77_CALL(dtrmv)("U", "N", "N", &n, a, &lda, tmp, &ione);
    for (i = 0, ans = 0; i < n; i++) ans += tmp[i] * tmp[i];
    Free(tmp);
    return ans;
}  

/** 
 * Evaluate the conditional deviance for a given set of fixed effects.
 * 
 * @param GS Pointer to a GlmerStruct
 * @param fixed value of the fixed effects
 * 
 * @return conditional deviance
 */
static double
fixed_effects_deviance(GlmerStruct GS, const double fixed[])
{
    SEXP devs;
    int i, ione = 1;
    double ans, one = 1, zero = 0;

    F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
		    REAL(GS->x), &(GS->n), fixed,
		    &ione, &zero, REAL(GS->eta), &ione);
				/* add in random effects and offset */
    vecIncrement(REAL(GS->eta), REAL(GS->off), GS->n);
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0, ans = 0; i < GS->n; i++) ans += REAL(devs)[i];
    UNPROTECT(1);
    return ans;
}
    
/** 
 * Determine the conditional modes and the conditional variance of the
 * fixed effects given the data and the current random effects.
 * Create a Metropolis-Hasting proposal step from the multivariate
 * normal density, determine the acceptance probability and, if the
 * step is to be accepted, overwrite the contents of fixed with the
 * new contents.
 * 
 * @param GS a GlmerStruct
 * @param b list of random effects
 * @param fixed current value of the fixed effects
 * 
 * @return updated value of the fixed effects
 */
static double *
internal_glmer_fixef_update(GlmerStruct GS, SEXP b,
			    double fixed[])
{
    SEXP dmu_deta, var;
    int i, ione = 1, it, j, lwork = -1;
    double *ans = Calloc(GS->p, double), /* proposal point */
	*etaold = Calloc(GS->n, double),
	*md = Calloc(GS->p, double), /* conditional modes */
	*w = Calloc(GS->n, double), *work,
	*wtd = Calloc(GS->n * GS->p, double),
	*z = Calloc(GS->n, double),
	crit, devr, one = 1, tmp, zero = 0;
    
    if (!isNewList(b) || LENGTH(b) != GS->nf)
	error(_("%s must be a %s of length %d"), "b", "list", GS->nf);
    for (i = 0; i < GS->nf; i++) {
	SEXP bi = VECTOR_ELT(b, i);
	if (!isReal(bi) || !isMatrix(bi))
	    error(_("b[[%d]] must be a numeric matrix"), i);
    }
    AZERO(z, GS->n);		/* -Wall */
    Memcpy(md, fixed, GS->p);
				/* calculate optimal size of work array */
    F77_CALL(dgels)("N", &(GS->n), &(GS->p), &ione, wtd, &(GS->n),
		    z,  &(GS->n), &tmp, &lwork, &j);
    if (j)			/* shouldn't happen */
	error(_("%s returned error code %d"), "dgels", j);
    lwork = (int) tmp;
    work = Calloc(lwork, double);
				
    AZERO(REAL(GS->off), GS->n); /* fitted values from random effects */
    fitted_ranef(GET_SLOT(GS->mer, Matrix_flistSym), GS->unwtd, b,
		 INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)), REAL(GS->off));
    for (i = 0; i < GS->n; i++)
	etaold[i] = (REAL(GS->off)[i] += REAL(GS->offset)[i]);
    
    for (it = 0, crit = GS->tol + 1;
	 it < GS->maxiter && crit > GS->tol; it++) {
				/* fitted values from current beta */
	F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
			REAL(GS->x), &(GS->n), md,
			&ione, &zero, REAL(GS->eta), &ione);
				/* add in random effects and offset */
	vecIncrement(REAL(GS->eta), REAL(GS->off), GS->n);
				/* check for convergence */
	crit = conv_crit(etaold, REAL(GS->eta), GS->n);
				/* obtain mu, dmu_deta, var */
	eval_check_store(GS->linkinv, GS->rho, GS->mu);
	dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				      REALSXP, GS->n));
	var = PROTECT(eval_check(GS->var, GS->rho, REALSXP, GS->n));
				/* calculate weights and working residual */
	for (i = 0; i < GS->n; i++) {
	    w[i] = REAL(GS->wts)[i] *
		REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	    z[i] = w[i] * (REAL(GS->eta)[i] - REAL(GS->off)[i] +
			   (REAL(GS->y)[i] - REAL(GS->mu)[i]) /
			   REAL(dmu_deta)[i]);
	}
	UNPROTECT(2);
				/* weighted copy of the model matrix */
	for (j = 0; j < GS->p; j++)
	    for (i = 0; i < GS->n; i++)
		wtd[i + j * GS->n] = REAL(GS->x)[i + j * GS->n] * w[i];
				/* weighted least squares solution */
	F77_CALL(dgels)("N", &(GS->n), &(GS->p), &ione, wtd, &(GS->n),
			z, &(GS->n), work, &lwork, &j);
	if (j) error(_("%s returned error code %d"), "dgels", j);
	Memcpy(md, z, GS->p);
    }
				/* wtd contains the Cholesky factor of
				 * the precision matrix */
    devr = normal_kernel(GS->p, md, wtd, GS->n, fixed);
    devr -= fixed_effects_deviance(GS, fixed);
    for (i = 0; i < GS->p; i++) {
	double var = norm_rand();
	ans[i] = var;
	devr -= var * var;
    }
    F77_CALL(dtrsv)("U", "N", "N", &(GS->p), wtd, &(GS->n), ans, &ione);
    for (i = 0; i < GS->p; i++) ans[i] += md[i];
    devr += fixed_effects_deviance(GS, ans);
    crit = exp(-0.5 * devr);	/* acceptance probability */
    tmp = unif_rand();
    if (asLogical(Matrix_getElement(GS->cv, "msVerbose"))) {
	Rprintf("%5.3f: ", crit);
	for (j = 0; j < GS->p; j++) Rprintf("%#10g ", ans[j]);
	Rprintf("\n");
    }
    if (tmp < crit) Memcpy(fixed, ans, GS->p);
    Free(ans); Free(etaold); Free(md); Free(w);
    Free(work); Free(wtd); Free(z);
    return fixed;
}

/** 
 * Determine the conditional modes and the conditional variance of the
 * fixed effects given the data and the current random effects.
 * Create a Metropolis-Hasting proposal step from the multivariate
 * normal density, determine the acceptance probability and return the
 * current value or the proposed value.
 * 
 * @param GSp pointer to a GlmerStruct
 * @param b list of random effects
 * @param fixed current value of the fixed effects
 * 
 * @return updated value of the fixed effects
 */
SEXP glmer_fixed_update(SEXP GSp, SEXP b, SEXP fixed)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
 
    if (!isReal(fixed) || LENGTH(fixed) != GS->p)
	error(_("%s must be a %s of length %d"), "fixed",
		"numeric vector", GS->p);
    GetRNGstate();
    internal_glmer_fixef_update(GS, b, REAL(fixed));
    PutRNGstate();
    return fixed;
}

/** 
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and df degrees of freedom.
 * 
 * @param df degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param ans array of size p * p to hold the result
 * 
 * @return ans
 */
static double*
std_rWishart_factor(double df, int p, double ans[])
{
    int i, j, pp1 = p + 1;

    if (df < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");
    for (j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(df - (double) j));
	for (i = 0; i < j; i++) ans[i + j * p] = norm_rand();
    }
    return ans;
}

/* FIXME: Combine internal_Omega_update and lmer_Omega_update.
 * Consider changing the internal storage of b in glmer_MCMCsamp to a
 * single vector.
 */

/** 
 * Simulate Omega in an mer object from a Wishart distribution with
 * scale given by the new values of the random effects.
 * 
 * @param x Pointer to an mer object
 * @param bnew Pointer to the new values of the random effects in the
 * chain
 */
static void
internal_Omega_update(SEXP x, SEXP bnew)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int i, info, nf = LENGTH(bnew);
    double one = 1, zero = 0;

    for (i = 0; i < nf; i++) {
	SEXP bobj = VECTOR_ELT(bnew, i);
	int *dims = INTEGER(getAttrib(bobj, R_DimSymbol));
	int ncisq = dims[1] * dims[1];
	double *tmp = Calloc(ncisq, double),
	    *omgi = REAL(VECTOR_ELT(Omega, i));
	
	AZERO(tmp, ncisq);
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[0]),
			&one, REAL(bobj), &(dims[0]),
			&zero, omgi, &(dims[1]));
	F77_CALL(dpotrf)("U", &(dims[1]), omgi, &(dims[1]), 
			 &info);
	if (info)
	    error(_("Singular random effects varcov at level %d"), i + 1);
	std_rWishart_factor((double) (dims[0] - dims[1] + 1),
			    dims[1], tmp);
	F77_CALL(dtrsm)("R", "U", "T", "N", &(dims[1]), &(dims[1]),
			&one, omgi, &(dims[1]), tmp, &(dims[1]));
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, omgi, &(dims[1]));
	Free(tmp);
    }
}

/** 
 * Create a Markov Chain Monte Carlo sample from a fitted generalized
 * linear mixed model
 * 
 * @param GSpt External pointer to a GlmerStruct
 * @param b Conditional modes of the random effects at the parameter
 * estimates
 * @param fixed Estimates of the fixed effects
 * @param varc Estimates of the variance components
 * @param savebp Logical indicator of whether or not to save the
 * random effects in the MCMC sample
 * @param nsampp Integer value of the number of samples to generate
 * 
 * @return 
 */
SEXP 
glmer_MCMCsamp(SEXP GSpt, SEXP b, SEXP fixed, SEXP varc,
	       SEXP savebp, SEXP nsampp) 
{ 
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSpt);
    int i, j, nf = LENGTH(b), nsamp = asInteger(nsampp),
	p = LENGTH(fixed), q = LENGTH(varc),
	saveb = asLogical(savebp);
    int *mcpar, nc = p + q;
    SEXP ans, mcparSym = install("mcpar");
    
    if (nsamp <= 0) nsamp = 1;
    nc = p + q;
    if (saveb)
	for (i = 0; i < nf; i++) {
	    int *dims = INTEGER(getAttrib(VECTOR_ELT(b, i),
					  R_DimSymbol));
	    nc += dims[0] * dims[1];
	}
    ans = PROTECT(allocMatrix(REALSXP, nsamp, nc));
    GetRNGstate();
    for (i = 0; i < nsamp; i++) {
	internal_glmer_fixef_update(GS, b, REAL(fixed));
	internal_bhat(GS, REAL(fixed), REAL(varc));
	internal_glmer_ranef_update(GS, b);
	internal_Omega_update(GS->mer, b);
	internal_coef(GS->mer, 2, REAL(varc));
	for (j = 0; j < p; j++)
	    REAL(ans)[i + j * nsamp] = REAL(fixed)[j];
	for (j = 0; j < q; j++)
	    REAL(ans)[i + (j + p) * nsamp] = REAL(varc)[j];
	if (saveb) {
	    int base = p + q, k;
	    for (k = 0; k < nf; k++) {
		SEXP bk = VECTOR_ELT(b, k);
		int *dims = INTEGER(getAttrib(bk, R_DimSymbol));
		int klen = dims[0] * dims[1];

		for (j = 0; j < klen; j++)
		    REAL(ans)[i + (j + base) * nsamp] = REAL(bk)[j];
		base += klen;
	    }
	}
    }
    PutRNGstate();
    UNPROTECT(1);
				/* set (S3) class and mcpar attribute */
    setAttrib(ans, R_ClassSymbol, mkString("mcmc"));
    setAttrib(ans, mcparSym, allocVector(INTSXP, 3));
    mcpar = INTEGER(getAttrib(ans, mcparSym));
    mcpar[0] = mcpar[2] = 1;
    mcpar[1] = nsamp;

    return ans;
} 

/** 
 * Simulate a sample of random matrices from a Wishart distribution
 * 
 * @param ns Number of samples to generate
 * @param dfp Degrees of freedom
 * @param scal Positive-definite scale matrix
 * 
 * @return 
 */
SEXP
Matrix_rWishart(SEXP ns, SEXP dfp, SEXP scal)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), j,
	n = asInteger(ns), psqr;
    double *scCp, *tmp, df = asReal(dfp), one = 1, zero = 0;

    if (!isMatrix(scal) || !isReal(scal) || dims[0] != dims[1])
	error("scal must be a square, real matrix");
    if (n <= 0) n = 1;
    psqr = dims[0] * dims[0];
    tmp = Calloc(psqr, double);
    AZERO(tmp, psqr);
    scCp = Memcpy(Calloc(psqr, double), REAL(scal), psqr);
    F77_CALL(dpotrf)("U", &(dims[0]), scCp, &(dims[0]), &j);
    if (j)
	error("scal matrix is not positive-definite");
    PROTECT(ans = alloc3Darray(REALSXP, dims[0], dims[0], n));
  
    GetRNGstate();
    for (j = 0; j < n; j++) {
	double *ansj = REAL(ans) + j * psqr;
	std_rWishart_factor(df, dims[0], tmp);
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims,
			&one, scCp, dims, tmp, dims);
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, ansj, &(dims[1]));
	internal_symmetrize(ansj, dims[0]);
	}

    PutRNGstate();
    Free(scCp); Free(tmp);
    UNPROTECT(1);
    return ans;
}

/** 
 * Update the relative precision matrices by sampling from a Wishart
 * distribution with scale factor determined by the current sample of
 * random effects.
 * 
 * @param x pointer to an mer object
 * @param b current sample from the random effects
 * @param sigma current value of sigma
 */
static void
lmer_Omega_update(SEXP Omega, const double b[], double sigma, int nf,
		  const int nc[], const int Gp[], double *bi, int nsamp,
		  int trans)
{
    int i, j, k, info;
    double one = 1, zero = 0;

    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int nlev = (Gp[i + 1] - Gp[i])/nci, ncip1 = nci + 1,
	    ncisqr = nci * nci;
	double
	    *scal = Calloc(ncisqr, double), /* factor of scale matrix */
	    *tmp = Calloc(ncisqr, double),
	    *var = Calloc(ncisqr, double), /* simulated variance-covariance */
	    *wfac = Calloc(ncisqr, double); /* factor of Wishart variate */

	/* generate and factor the scale matrix */
	AZERO(scal, ncisqr);
	F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one, b + Gp[i], &nci,
			&zero, scal, &nci);
	F77_CALL(dpotrf)("U", &nci, scal, &nci, &info);
	if (info)
	    error(_("Singular random effects varcov at level %d"), i + 1);

	/* generate a random factor from a standard Wishart distribution */
	AZERO(wfac, ncisqr);
	std_rWishart_factor((double) (nlev - nci + 1), nci, wfac);

	/* form the variance-covariance matrix and store elements */
	Memcpy(tmp, scal, ncisqr);
	F77_CALL(dtrsm)("L", "U", "T", "N", &nci, &nci,
			&one, wfac, &nci, tmp, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, tmp, &nci,
			&zero, var, &nci);
	for (j = 0; j < nci; j++) {
	    *bi = (trans ? log(var[j * ncip1]) : var[j * ncip1]);
	    bi += nsamp;
	}
	for (j = 1; j < nci; j++) {
	    for (k = 0; k < j; k++) {
		*bi = (trans ? atanh(var[k + j * nci]/
				     sqrt(var[j * ncip1] * var[k * ncip1]))
		       : var[k + j * nci]);
		bi += nsamp;
	    }
	}

	/* calculate and store the relative precision matrix */
	Memcpy(tmp, wfac, ncisqr);
	F77_CALL(dtrsm)("R", "U", "T", "N", &nci, &nci,
			&sigma, scal, &nci, tmp, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, tmp, &nci,
			&zero, REAL(VECTOR_ELT(Omega, i)), &nci);
	Free(scal); Free(tmp); Free(wfac); Free(var); 
    }
}

/** 
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
 * 
 * @param mer pointer to a mixed-effects representation object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * 
 * @return a matrix
 */
SEXP
lmer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp)
{
    SEXP ans, DP = GET_SLOT(x, Matrix_DSym),
	LP = GET_SLOT(x, Matrix_LSym),
	Omega = GET_SLOT(x, Matrix_OmegaSym), 
	RXXsl = GET_SLOT(x, Matrix_RXXSym),
	RZXsl = GET_SLOT(x, Matrix_RZXSym),
	mcparSym = install("mcpar");
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	*mcpar,
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),
		       "REML"),
	i, ione = 1, j, nf = LENGTH(Omega),
	nsamp = asInteger(nsampp),
	saveb = asLogical(savebp),
	trans = asLogical(transp);
    int ncol = coef_length(nf, nc) + dims[1] + (saveb ? dims[0] : 0),
	nobs = nc[nf + 1],
	p = dims[1] - 1, pp1 = dims[1];
    double 
	*bhat = Calloc(dims[0], double),
	*betahat = Calloc(p, double), one = 1,
	*bnew = Calloc(dims[0], double),
	*betanew = Calloc(p, double),
	df = (double)(REML ? nobs + 1 - pp1 : nobs);
    
    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nsamp, ncol));
    GetRNGstate();
    for (i = 0; i < nsamp; i++) {
	double sigma, ryyinv; /* ryy-inverse */

	/* decompose and invert the mer object */
	status[0] = status[1] = 0;
	lmer_invert(x);

	/* simulate and store new value of sigma */
	ryyinv = REAL(RXXsl)[pp1 * pp1 - 1];
	sigma = 1/(ryyinv * sqrt(rchisq(df)));
	REAL(ans)[i + p * nsamp] =
	    (trans ? 2 * log(sigma) : sigma * sigma);

	/* determine bhat and betahat */
	Memcpy(bhat, REAL(RZXsl) + dims[0] * p, dims[0]);
	for (j = 0; j < dims[0]; j++) bhat[j] /= -ryyinv;
	Memcpy(betahat, REAL(RXXsl) + pp1 * (pp1 - 1), p);
	for (j = 0; j < p; j++) betahat[j] /= -ryyinv;

	/* simulate scaled, independent normals */
	for (j = 0; j < p; j++) betanew[j] = sigma * norm_rand();
	for (j = 0; j < dims[0]; j++) bnew[j] = sigma * norm_rand();

	/* bnew := D^{-1/2} %*% bnew */
	for (j = 0; j < nf; j++) {
	    int k, ncj = nc[j];
	    int ncjsqr = ncj * ncj,
		nlev = (Gp[j+1] - Gp[j])/ncj;
	    double *Dj = REAL(VECTOR_ELT(DP, j)), *bj = bnew + Gp[j];

	    for (k = 0; k < nlev; k++)
		F77_CALL(dtrmv)("U", "N", "N", &ncj, Dj + k * ncjsqr,
				&ncj, bj + k * ncj, &ione);
	}

	/* bnew := L^{-T} %*% bnew */
	internal_sm(LFT, TRN, nf, Gp, 1, 1.0, LP, bnew, dims[0]);

	/* bnew := bnew + RZX %*% betanew*/
	F77_CALL(dgemv)("N", dims, &p, &one, REAL(RZXsl), dims,
			betanew, &ione, &one, bnew, &ione);

	/* betanew := RXX^{-1} %*% betanew */
	F77_CALL(dtrmv)("U", "N", "N", &p, REAL(RXXsl), &pp1,
			betanew, &ione);

	/* Add conditional modes and store */
	for (j = 0; j < p; j++) {
	    REAL(ans)[i + j * nsamp] = (betanew[j] += betahat[j]);
	}
	for (j = 0; j < dims[0]; j++) bnew[j] += bhat[j];

	/* Update and store variance-covariance of the random effects */
	lmer_Omega_update(Omega, bnew, sigma, nf, nc, Gp,
			  REAL(ans) + i + (p + 1) * nsamp, nsamp, trans);
    }
    PutRNGstate();
    Free(betanew); Free(bnew);
				/* set (S3) class and mcpar attribute */
    setAttrib(ans, R_ClassSymbol, mkString("mcmc"));
    setAttrib(ans, mcparSym, allocVector(INTSXP, 3));
    mcpar = INTEGER(getAttrib(ans, mcparSym));
    mcpar[0] = mcpar[2] = 1;
    mcpar[1] = nsamp;

    UNPROTECT(1);
    return ans;
}

/** 
 * Simulate a set of linear predictors from an mer object
 * 
 * @param x Pointer to an mer object
 * @param np Pointer to an integer giving the number of values to simulate
 * @param fxd Pointer to the fixed-effects prediction
 * @param mmats Pointer to a list of model matrices
 * @param useScale Logical indicator of whether to use the scale
 * 
 * @return a matrix of simulated linear predictors
 */
SEXP lmer_simulate(SEXP x, SEXP np, SEXP fxd, SEXP mmats, SEXP useScale)
{
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	i, ii, ione = 1, n = asInteger(np), nobs = LENGTH(fxd);
    SEXP ans = PROTECT(allocMatrix(REALSXP, nobs, n)),
	flist = GET_SLOT(x, Matrix_flistSym),
	Omega = GET_SLOT(x, Matrix_OmegaSym);
    double *aa = REAL(ans),
	one = 1.,
	scale = (asLogical(useScale)
		 ? internal_sigma(x,
				  !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),
					  "REML"))
		 : 1.);

    GetRNGstate();

    /* initialize the columns to the fixed effects prediction */
    if (!isReal(fxd))
	error(_("fxd must be a numeric vector"));
    for (ii = 0; ii < n; ii++) Memcpy(aa + ii * nobs, REAL(fxd), nobs);

				/* simulate and add random effects */
    for (i = 0; i < LENGTH(flist); i++) {
	SEXP mmi = VECTOR_ELT(mmats, i);
	int *fv = INTEGER(VECTOR_ELT(flist, i)), 
	    *dims = INTEGER(getAttrib(mmi, R_DimSymbol)),
	    nci = nc[i], relen = Gp[i + 1] - Gp[i];
	int ncisqr = nci * nci, ntot = relen * n;
	int ncol = ntot/nci, bsz = n * nci;
	double *mm = REAL(mmi),
	    *rev = Calloc(ntot, double),
	    *Rmat = Memcpy(Calloc(ncisqr, double),
			   REAL(VECTOR_ELT(Omega, i)), ncisqr);
	
	if (!isMatrix(mmi) || !isReal(mmi) || dims[0] != nobs || dims[1] != nci)
	    error(_("mmats[[%d]] must be a numeric matrix of dimension %d by %d"),
		    i + 1, nobs, nci);
	F77_CALL(dpotrf)("U", &nci, Rmat, &nci, &ii);
	if (ii)
	    error(_("Leading %d minor of Omega[[%d]] not positive definite"),
		  ii, i + 1);
	for (ii = 0; ii < ntot; ii++) rev[ii] = norm_rand();
	F77_CALL(dtrsm)("L", "U", "N", "N", &nci, &ncol, &scale, Rmat, &nci,
			rev, &nci);
	for (ii = 0; ii < nobs; ii++) /* loop over rows */
	    F77_CALL(dgemm)("N", "N", &ione, &n, &nci, &one, mm + ii, &nobs,
			    rev + (fv[ii] - 1) * bsz, &nci, &one, aa + ii, &nobs);
	Free(rev); Free(Rmat); 
    }

    PutRNGstate();
    UNPROTECT(1);
    return ans;
}

/** 
 * Update only the response in an mer object.
 * 
 * @param x Pointer to an mer object
 * @param y Numeric vector of responses
 * @param mmats List of model matrices.  Only the first p columns from
 * the fixed-effects model matrix are used.
 * 
 * @return x, updated for the new response vector
 */
SEXP lmer_update_y(SEXP x, SEXP y, SEXP mmats)
{
    SEXP ZtXP = GET_SLOT(x, Matrix_ZtXSym),
	XtXP = GET_SLOT(x, Matrix_XtXSym),
	flist = GET_SLOT(x, Matrix_flistSym);
    int *dims = INTEGER(getAttrib(ZtXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	i, ione = 1, nf = LENGTH(flist);
    int nfp1 = nf + 1, nobs = nc[nf + 1], p = dims[1] - 1;
    double *X = REAL(VECTOR_ELT(mmats, nf)),
	*ycol = REAL(ZtXP) + dims[0] * p,
	one = 1, zero = 0;
    
    if (!isReal(y) || LENGTH(y) != nobs)
	error(_("Argument `y' must be a numeric vector of length %d"),
		nobs);
    if (!isNewList(mmats) || LENGTH(mmats) != nfp1)
	error(_("mmats must be a list of %d model matrices"), nfp1);
    if (!isMatrix(VECTOR_ELT(mmats, nf)) || !isReal(VECTOR_ELT(mmats, nf)))
	error(_("mmats[[%d]] is not a numeric matrix"), nf + 1);
				/* Update XtX */
    F77_CALL(dgemv)("T", &nobs, &p, &one, X, &nobs, REAL(y), &ione,
		    &zero, REAL(XtXP) + p * dims[1], &ione);
    REAL(XtXP)[dims[1] * dims[1] - 1] =
	F77_CALL(ddot)(&nobs, REAL(y), &ione, REAL(y), &ione);
				/* Update ZtX */
    AZERO(ycol, Gp[nf]);	/* Zero the accumulator */
    for (i = 0; i < nf; i++) {
	SEXP mmat = VECTOR_ELT(mmats, i);
	int *mdims = INTEGER(getAttrib(mmat, R_DimSymbol)),
	    *fac = INTEGER(VECTOR_ELT(flist, i)),
	    j, k, nci = nc[i];
	double *yy = ycol + Gp[i], *Z = REAL(mmat);

	if (!isMatrix(mmat) || !isReal(mmat))
	    error(_("mmats[[%d]] is not a numeric matrix"), i + 1);
	if (nobs != mdims[0])
	    error(_("Expected %d rows in mmats[[%d]]. Got %d"),
		  nobs, i+1, mdims[0]);
	if (nc[i] != mdims[1])
	    error(_("Expected %d columns in mmats[[%d]]. Got %d"),
		  nc[i], i+1, mdims[1]);

	for (j = 0; j < nobs; j++) {
	    double *yp = yy + (fac[j] - 1) * nci, yv = REAL(y)[j];
	    for (k = 0; k < nci; k++) yp[k] += Z[j + nobs * k] * yv;
	}
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}

/* FIXME: Set this up so that there is the ability to update without reanalyzing */
/* Change the class definition to include Zt and X slots.  Allocate Zt as a dgCMatrix
   and save it. */
/* Do the checking of the factor list separately. */
static cholmod_sparse*
list_to_sparse(SEXP rand, int nobs, int *nc, int *Gp, SEXP flist)
{
    int *dims, *p, *ii, i, nctot = 0, nf = LENGTH(rand);
    int *offset = Calloc(nf + 1, int);
    double *x;
    SEXP fact, tmmat, randi;
    cholmod_sparse *A;

    if (!isNewList(rand) || nf < 1)
	error(_("random must be a non-empty list"));
    offset[0] = Gp[0] = 0;
    for (i = 0; i < nf; i++) {	/* check consistency and establish dimensions */
	randi = VECTOR_ELT(rand, i);
	if (!isNewList(randi) || LENGTH(randi) != 2)
	    error(_("random[[%d]] must be a list of length 2"),
		i + 1);
				/* transpose of model matrix */
	tmmat = VECTOR_ELT(randi, 0);
	if (!isMatrix(tmmat) || !isReal(tmmat))
	    error(_("random[[%d]][[1]] must be real matrix"),
		  i + 1);
	dims = INTEGER(getAttrib(tmmat, R_DimSymbol));
	nctot += (nc[i] = dims[0]); /* matrix has been transposed */
	if (dims[1] != nobs)
	    error(_("random[[%d]][[1]] must have %d columns"),
		  i + 1, nobs);
				/* grouping factor */
	SET_VECTOR_ELT(flist, i, fact = VECTOR_ELT(randi, 1));
	if (!isFactor(fact) || LENGTH(fact) != nobs)
	    error(_("random[[%d]][[2]] must be a factor of length %d"),
		  i + 1, nobs);
	Gp[i + 1] = Gp[i] + nc[i] * LENGTH(getAttrib(fact, R_LevelsSymbol));
	offset[i + 1] = offset[i] + nc[i];
    }
				/* allocate A, extract pointers, fill A->p */
    A = cholmod_allocate_sparse((size_t) Gp[nf], (size_t) nobs,
				(size_t) nctot * nobs, TRUE/* sorted */,
				TRUE/* packed */, 0/* stype */, CHOLMOD_REAL, &c);
    p = (int *) A->p; ii = (int *) A->i; x = (double *) A->x;
    p[0] = 0; for(i = 0; i < nobs; i++) p[i + 1] = p[i] + nctot;

    for (i = 0; i < nf; i++) {	/* fill A */
	int *vals = INTEGER(VECTOR_ELT(flist, i)), j;
	double *xvals = REAL(VECTOR_ELT(VECTOR_ELT(rand, i), 0));

	for (j = 0; j < nobs; j++) {
	    int jbase = Gp[i] + nc[i] * (vals[j] - 1), k;
	    for (k = 0; k < nc[i]; k++) {
		int ind = j * nctot + offset[i] + k;
		ii[ind] = jbase + k;
		x[ind] = xvals[j * nc[i] + k];
	    }
	}
    }

    Free(offset);
    return A;
}

static cholmod_dense
*numeric_as_chm_dense(double *v, int n)
{
    cholmod_dense *ans = Calloc(1, cholmod_dense);
    
    ans->d = ans->nzmax = ans->nrow = n;
    ans->ncol = 1;
    ans->x = (void *) v;
    ans->xtype = CHOLMOD_REAL;
    ans->dtype = CHOLMOD_DOUBLE;
    return ans;
}

#if 0
static
SEXP alloc_dgeMatrix(int m, int n)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = m; dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, m * n);
    UNPROTECT(1);
    return ans;
}
#endif

static
SEXP alloc_dsyMatrix(int n, char *uplo)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsyMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * n);
    SET_SLOT(ans, Matrix_uploSym, mkString(uplo));
    UNPROTECT(1);
    return ans;
}

static
SEXP alloc_dtrMatrix(int n, char *uplo, char *diag)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * n);
    SET_SLOT(ans, Matrix_uploSym, mkString(uplo));
    SET_SLOT(ans, Matrix_diagSym, mkString(diag));
    UNPROTECT(1);
    return ans;
}

/**
 * Create and insert initial values for Omega.
 *
 * @param x pointer to an mer2 object
 *
 * @return NULL
 */
SEXP mer2_initial(SEXP x)
{
    SEXP Omg = GET_SLOT(x, Matrix_OmegaSym),
	ZtZ = GET_SLOT(x, Matrix_ZtZSym);
    int	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*p = INTEGER(GET_SLOT(ZtZ, Matrix_pSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	i, nf = length(Omg);
    double *xx = REAL(GET_SLOT(ZtZ, Matrix_xSym));

    for (i = 0; i < nf; i++) {
	double *omgi = REAL(GET_SLOT(VECTOR_ELT(Omg, i),
				     Matrix_xSym));
	int bb = Gp[i], j, k, nci = nc[i];
	int ncip1 = nci + 1, nlev = (Gp[i + 1] - bb)/nci;

	AZERO(omgi, nci * nci);
	for (j = 0; j < nlev; j++) {
	    int base = bb + j * nci;
	    for (k = 0; k < nci; k++)
				/* add the last element in the column */
		omgi[k * ncip1] += xx[p[base + k + 1] - 1];
	}
	for (k = 0; k < nci; k++) omgi[k * ncip1] *= 0.375/nlev;
    }
    status[0] = status[1] = FALSE;
    return R_NilValue;
}

/** 
 * Create a copy of ZtZ with the diagonal blocks inflated according to Omega
 * 
 * @param zz cholmod_sparse version of ZtZ
 * @param nf number of factors
 * @param Omega Omega list
 * @param nc number of columns in model matrices
 * @param Gp group pointers
 * 
 * @return a freshly allocated cholmod_sparse version of the sum
 */
static cholmod_sparse *
ZZ_inflate(cholmod_sparse *zz, int nf, SEXP Omega, int *nc, int *Gp)
{
    cholmod_sparse *Omg, *ans;
    int *omp, *nnz = Calloc(nf + 1, int), i;
    double one = 1;

    for (nnz[0] = 0, i = 0; i < nf; i++)
	nnz[i + 1] = nnz[i] + (Gp[i + 1] - Gp[i])*(nc[i] + 1)/2;
    Omg = cholmod_allocate_sparse(zz->nrow, zz->ncol, (size_t) nnz[nf],
				  TRUE, TRUE, 1, CHOLMOD_REAL, &c);
    omp = (int *) Omg->p;
    for (i = 0; i < nf; i++) {
	int bb = Gp[i], j, jj, k, nci = nc[i];
	int nlev = (Gp[i + 1] - bb)/nci;
	double *Omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym));

	for (j = 0; j < nlev; j++) { /* column of result */
	    int col0 = bb + j * nci; /* absolute column number */

	    for (jj = 0; jj < nci; jj++) { /* column of Omega_i */
		int coljj = col0 + jj;

		omp[coljj + 1] = omp[coljj] + jj + 1;
		for (k = 0; k <= jj; k++) { /* row of Omega_i */
		    int ind = omp[coljj];
		    ((int *)Omg->i)[ind + k] = col0 + k;
		    ((double *)Omg->x)[ind + k] = Omgi[jj * nci + k];
		}
	    }
	}
    }
    ans = cholmod_add(zz, Omg, &one, &one, TRUE, TRUE, &c);

    Free(nnz); cholmod_free_sparse(&Omg, &c);
    return ans;
}

/**
 * Create an mer2 object from a list of grouping factors and a list of model
 * matrices.  There is one more model matrix than grouping factor.  The last
 * model matrix is the fixed effects and the response.
 *
 * @param flist pointer to a list of grouping factors
 * @param mmats pointer to a list of model matrices
 *
 * @return pointer to an mer2 object
 */
SEXP mer2_create(SEXP random, SEXP Xp, SEXP yp, SEXP method)
{
    SEXP LL, Omega, XtXp, fl,
	fnms = getAttrib(random, R_NamesSymbol), tt,
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("mer2")));
    cholmod_sparse *Zt, *ZtZ, *A;
    cholmod_dense *X = as_cholmod_dense(Xp), *tmp1, *tmp2;
    cholmod_factor *F;
    int *nc, *Gp, *xdims, ione = 1, j, nf = LENGTH(random), nobs = LENGTH(yp), p, q;
    double *XtX, one = 1, zero = 0;
    char *statnms[] = {"factored", "inverted", ""},
	*devnms[] = {"ML", "REML", ""};

    if (!isReal(yp)) error(_("yp must be a real vector"));
    if (!isMatrix(Xp) || !isReal(Xp))
	error(_("Xp must be a real matrix"));
    xdims = INTEGER(getAttrib(Xp, R_DimSymbol));
    if (xdims[0] != nobs) error(_("Xp must have %d rows"), nobs);
    p = xdims[1];
				/* allocate slots in val */
    fl = ALLOC_SLOT(val, Matrix_flistSym, VECSXP, nf);
    nc = INTEGER(ALLOC_SLOT(val, Matrix_ncSym, INTSXP, nf + 1));
    nc[nf] = nobs;
    Gp = INTEGER(ALLOC_SLOT(val, Matrix_GpSym, INTSXP, nf + 1));
				/* check random, create Zt and ZtZ */
    Zt = list_to_sparse(random, nobs, nc, Gp, fl);
				/* analyze ZtZ */
    q = Zt->nrow;
    j = c.supernodal;
    c.supernodal = CHOLMOD_SUPERNODAL;
    F = cholmod_analyze(Zt, &c);
    c.supernodal = j;
    LL = ALLOC_SLOT(val, Matrix_LSym, VECSXP, 1);
/* FIXME: Need to set up a finalizer for F.  Right now that storage is not being released. */
    SET_VECTOR_ELT(LL, 0, R_MakeExternalPtr(F, R_NilValue, val));
    A = cholmod_aat(Zt, (int *) NULL, (size_t) 0, 1/* mode */, &c);
    ZtZ = cholmod_copy(A, 1/* stype */, 1/* mode */, &c);
    cholmod_free_sparse(&A, &c);
    SET_SLOT(val, Matrix_ZtZSym, chm_sparse_to_SEXP(ZtZ, 0)); /* do not free ZtZ */
				/* allocate other slots */
    SET_SLOT(val, Matrix_methodSym, duplicate(method));
    SET_SLOT(val, Matrix_devianceSym, Matrix_make_named(REALSXP, devnms));
    SET_SLOT(val, Matrix_statusSym, Matrix_make_named(LGLSXP, statnms));
    LOGICAL(GET_SLOT(val, Matrix_statusSym))[0] = 0;
    tt = ALLOC_SLOT(val, Matrix_devCompSym, REALSXP, 5);
    AZERO(REAL(tt), 5);		/* cosmetic */
				/* save y'y */
    REAL(tt)[0] = F77_CALL(ddot)(&nobs, REAL(yp), &ione, REAL(yp), &ione);
				/* allocate and populate Omega */
    Omega = ALLOC_SLOT(val, Matrix_OmegaSym, VECSXP, nf);
    setAttrib(Omega, R_NamesSymbol, duplicate(fnms));
    for (j = 0; j < nf; j++)
	SET_VECTOR_ELT(Omega, j, alloc_dsyMatrix(nc[j], "U"));
				/* initialize Omega and factorize Z'Z+Omega */
    mer2_initial(val);
    A = ZZ_inflate(ZtZ, nf, Omega, nc, Gp);
    if (!cholmod_factorize(A, F, &c)) /* should never happen */
	error(_("rank_deficient Z'Z+Omega"));
    cholmod_free_sparse(&A, &c);
    cholmod_free_sparse(&ZtZ, &c);
				/* create ZtX, RZX, XtX, RXX */
    tmp1 = cholmod_allocate_dense(q, p, q, CHOLMOD_REAL, &c);
    if (!cholmod_sdmult(Zt, 0, &one, &zero, X, tmp1, &c))
	error(_("cholmod_sdmult failed"));
    tmp2 = cholmod_solve(CHOLMOD_P, F, tmp1, &c);
    cholmod_free_dense(&tmp1, &c);
    SET_SLOT(val, Matrix_ZtXSym, chm_dense_to_SEXP(tmp2, 1)); /* frees tmp2 */
    SET_SLOT(val, Matrix_RZXSym,
	     duplicate(GET_SLOT(val, Matrix_ZtXSym)));
    XtXp = alloc_dsyMatrix(p, "U");
    SET_SLOT(val, Matrix_XtXSym, XtXp);
    XtX = REAL(GET_SLOT(XtXp, Matrix_xSym));
    AZERO(XtX, p * p);
    F77_CALL(dsyrk)("U", "T", &p, &nobs, &one, REAL(Xp), &nobs,
		    &zero, XtX, &p);
    XtXp = alloc_dtrMatrix(p, "U", "N");
    SET_SLOT(val, Matrix_RXXSym, XtXp);
    AZERO(REAL(GET_SLOT(XtXp, Matrix_xSym)), p * p);
				/* create Zty, rZy, Xty, rXy */
    tmp1 = numeric_as_chm_dense(REAL(yp), nobs);
    tmp2 = cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!cholmod_sdmult(Zt, 0, &one, &zero, tmp1, tmp2, &c))
	error(_("cholmod_sdmult failed"));
    Free(tmp1); 
    tmp1 = cholmod_solve(CHOLMOD_P, F, tmp2, &c);
    Memcpy(REAL(ALLOC_SLOT(val, Matrix_ZtySym, REALSXP, q)),
	   (double *)tmp1->x, q);
    cholmod_free_dense(&tmp2, &c);
    cholmod_free_dense(&tmp1, &c);
    SET_SLOT(val, Matrix_rZySym,
	     duplicate(GET_SLOT(val, Matrix_ZtySym)));
    tt = ALLOC_SLOT(val, Matrix_XtySym, REALSXP, p);
    F77_CALL(dgemm)("T", "N", &p, &ione, &nobs, &one, REAL(Xp), &nobs,
		    REAL(yp), &nobs, &zero, REAL(tt), &p);
    SET_SLOT(val, Matrix_rXySym, duplicate(tt));
    
    Free(X); cholmod_free_sparse(&Zt, &c);
    UNPROTECT(1);
    return val;
}


static double chm_log_abs_det(cholmod_factor *F)
{
    double ans = 0;

    if (F->is_super) {
	int i;
	for (i = 0; i < F->nsuper; i++) {
	    int j, nrp1 = 1 + ((int *)(F->pi))[i + 1] - ((int *)(F->pi))[i],
		nc = ((int *)(F->super))[i + 1] - ((int *)(F->super))[i];
	    double *x = (double *)(F->x) + ((int *)(F->px))[i];

	    for (j = 0; j < nc; j++) ans += log(fabs(x[j * nrp1]));
	}
    } else
	error(_("code for simplicial factorization not yet written"));
    return ans;
}

static double Omega_log_det(SEXP Omega, int nf, int *nc, int *Gp)
{
    double ans = 0;
    int i;

    for (i = 0; i < nf; i++) {
	int j, nci = nc[i], ncip1 = nc[i] + 1, ncisq = nc[i] * nc[i],
	    nlev = (Gp[i + 1] - Gp[i])/nc[i];
	double *omgi = Memcpy(Calloc(ncisq, double),
			      REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym)),
			      ncisq);
	F77_CALL(dpotrf)("U", &nci, omgi, &nci, &j);
	if (j)
	    error(_("Leading minor of size %d in Omega[[%d]] is not positive definite"),
		  j, i + 1);
	for (j = 0; j < nci; j++) ans += 2. * nlev * log(fabs(omgi[j * ncip1]));
	Free(omgi);
    }
    return ans;
}

/**
 * If status[["factored"]] is FALSE, create and factor Z'Z+Omega.  Also
 * create RZX and RXX, the deviance components, and the value of the
 * deviance for both ML and REML.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP mer2_factor(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));

    if (!status[0]) {
	SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
	cholmod_sparse *A,
	    *zz = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtZSym));
	cholmod_factor *L = (cholmod_factor *)
	    R_ExternalPtrAddr(VECTOR_ELT(GET_SLOT(x, Matrix_LSym), 0));
	cholmod_dense *ZtX = as_cholmod_dense(GET_SLOT(x, Matrix_ZtXSym)),
	    *Zty = numeric_as_chm_dense(REAL(GET_SLOT(x, Matrix_ZtySym)), L->n),
	    *rZy, *RZX;
	int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)), i, info, ione = 1,
	    nf = LENGTH(Omega), p = ZtX->ncol, q = L->n;
	int nml = nc[nf]; int nreml = nml - p;
	double *RXX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym), Matrix_xSym)),
	    *rXy = REAL(GET_SLOT(x, Matrix_rXySym)),
	    *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	    *dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	    one[2] = {1, 0}, m1[2] = {-1, 0};
	    
 	dcmp[1] = Omega_log_det(Omega, nf, nc, Gp); /* logDet(Omega) */
	A = ZZ_inflate(zz, nf, Omega, nc, Gp); Free(zz);
	if (!cholmod_factorize(A, L, &c))
	    error(_("rank_deficient Z'Z+Omega"));
	cholmod_free_sparse(&A, &c);
	dcmp[3] = 2 * chm_log_abs_det(L); /* 2 * logDet(L) */
				/* calculate and store RZX and rZy */
	RZX = cholmod_solve(CHOLMOD_L, L, ZtX, &c); Free(ZtX);
	rZy = cholmod_solve(CHOLMOD_L, L, Zty, &c); Free(Zty);
	Memcpy(REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXSym), Matrix_xSym)),
	       (double *) RZX->x, q * p);
	Memcpy(REAL(GET_SLOT(x, Matrix_rZySym)), (double *) rZy->x, q);
				/* downdate XtX and factor */
	Memcpy(RXX, REAL(GET_SLOT(GET_SLOT(x, Matrix_XtXSym), Matrix_xSym)), p * p);
	F77_CALL(dsyrk)("U", "T", &p, &q, m1, (double*)RZX->x, &q, one, RXX, &p);
	F77_CALL(dpotrf)("U", &p, RXX, &p, &info);
	if (info) {
	    error(_("Leading minor of order %d in downdated X'X is not positive definite"),
		  info);
	    dcmp[2] = dcmp[4] = dev[0] = dev[1] = NA_REAL;
	} else {
	    for (dcmp[4] = 0, i = 0; i < p; i++) /* 2 * logDet(RXX) */
		dcmp[4] += 2. * log(RXX[i * (p + 1)]);
				/* solve for rXy  and ryy^2 */
	    Memcpy(rXy, REAL(GET_SLOT(x, Matrix_XtySym)), p);
	    F77_CALL(dgemv)("T", &q, &p, m1, (double*) RZX->x, &q,
			    (double*) rZy->x, &ione, one, rXy, &ione);
	    F77_CALL(dtrsv)("U", "T", "N", &p, RXX, &p, rXy, &ione);
	    dcmp[2] = log(dcmp[0]	/* dcmp[2] = log(ryy^2); dcmp[0] = y'y; */
			  - F77_CALL(ddot)(&p, rXy, &ione, rXy, &ione)
			  - F77_CALL(ddot)(&q, (double*)rZy->x, &ione,
					   (double*)rZy->x, &ione));
				/* evaluate ML and REML deviance */
	    dev[0] = dcmp[3] - dcmp[1] +
		nml*(1.+dcmp[2]+log(2.*PI/nml));
	    dev[1] = dcmp[3] - dcmp[1] + dcmp[4] +
		nreml*(1.+dcmp[2]+log(2.*PI/nreml));
	}
	    
	cholmod_free_dense(&RZX, &c); cholmod_free_dense(&rZy, &c);
	status[0] = TRUE; status[1] = FALSE;
    }
    return R_NilValue;
}

/** 
 * Return the permutation of the columns of Z as a pMatrix object
 * 
 * @param x pointer to an mer2 object
 * 
 * @return the permutation as an pMatrix object
 */
SEXP mer2_pMatrix(SEXP x)
{
    cholmod_factor *L = (cholmod_factor *)
	R_ExternalPtrAddr(VECTOR_ELT(GET_SLOT(x, Matrix_LSym), 0));
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("pMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	*perm = INTEGER(ALLOC_SLOT(ans, Matrix_permSym, INTSXP, L->n)), i;

    dims[0] = dims[1] = (int) L->n;
    for (i = 0; i < (int) L->n; i++) perm[i] = ((int *)(L->Perm))[i] + 1;
    UNPROTECT(1);
    return ans;
}

/** 
 * Return L as a dtCMatrix object
 * 
 * @param x pointer to an mer2 object
 * 
 * @return L as an dtCMatrix object
 */
SEXP mer2_dtCMatrix(SEXP x)
{
    cholmod_factor *L;
    cholmod_sparse *Lm;
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	nz, q;

    mer2_factor(x);
    L = cholmod_copy_factor((cholmod_factor *)
			    R_ExternalPtrAddr(VECTOR_ELT(GET_SLOT(x, Matrix_LSym),
							 0)), &c);
    dims[0] = dims[1] = q = (int)(L->n);
    Lm = cholmod_factor_to_sparse(L, &c); cholmod_free_factor(&L, &c);
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, q + 1)),
	   (int *) Lm->p, q + 1);
    nz = ((int *)(Lm->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)),
	   (int *) Lm->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)),
	   (double *) Lm->x, nz);
    cholmod_free_sparse(&Lm, &c);
    UNPROTECT(1);
    return ans;
}

static double *
internal_mer2_coef(SEXP x, int ptyp, double ans[])
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), vind;

    vind = 0;			/* index in ans */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1;
	if (nci == 1) {
	    double dd = REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym))[0];
	    ans[vind++] = ptyp ? ((ptyp == 1) ? log(dd) : 1./dd) : dd;
	} else {
	    if (ptyp) {	/* L log(D) L' factor of Omega[,,i] */
		int j, k, ncisq = nci * nci;
		double *tmp = Memcpy(Calloc(ncisq, double),
				     REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym)),
				     ncisq);
		F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
		if (j)		/* should never happen */
		    error(_("DPOTRF returned error code %d on Omega[[%d]]"),
			  j, i+1);
		for (j = 0; j < nci; j++) {
		    double diagj = tmp[j * ncip1];
		    ans[vind++] = (ptyp == 1) ? (2. * log(diagj)) :
			1./(diagj * diagj);
		    for (k = j + 1; k < nci; k++) {
			tmp[j + k * nci] /= diagj;
		    }
		}
		for (j = 0; j < nci; j++) {
		    for (k = j + 1; k < nci; k++) {
			ans[vind++] = tmp[j + k * nci];
		    }
		}
		Free(tmp);
	    } else {		/* upper triangle of Omega[,,i] */
		int j, k, odind = vind + nci;
		double *omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym));

		for (j = 0; j < nci; j++) {
		    ans[vind++] = omgi[j * ncip1];
		    for (k = j + 1; k < nci; k++) {
			ans[odind++] = omgi[k*nci + j];
		    }
		}
		vind = odind;
	    }
	}
    }
    return ans;
}

/**
 * Extract parameters from the Omega matrices.  These aren't
 * "coefficients" but the extractor is called coef for historical
 * reasons.  Within each group these values are in the order of the
 * diagonal entries first then the strict upper triangle in row
 * order.
 * 
 * The parameters can be returned in three forms:
 *   0 - nonlinearly constrained - elements of the relative precision matrix
 *   1 - unconstrained - from the LDL' decomposition - logarithms of
 *       the diagonal elements of D
 *   2 - box constrained - also from the LDL' decomposition - inverses
 *       of the diagonal elements of D
 *
 * @param x pointer to an mer2 object
 * @param pType pointer to an integer scalar indicating the form of the 
 *        parameters to be returned.
 *
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP mer2_coef(SEXP x, SEXP pType)
{
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	nf = LENGTH(GET_SLOT(x, Matrix_OmegaSym));
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));

    internal_mer2_coef(x, asInteger(pType), REAL(val));
    UNPROTECT(1);
    return val;
}

static
void internal_mer2_coefGets(SEXP x, const double cc[], int ptyp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	cind, i, nf = length(Omega);

    cind = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	if (nci == 1) {
	    double dd = cc[cind++];
	    REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym))[0] =
		ptyp ? ((ptyp == 1) ? exp(dd) : 1./dd) : dd;
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		j, k,
		ncip1 = nci + 1,
		ncisq = nci * nci;
	    double
		*omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym));
	    if (ptyp) {
		double *tmp = Calloc(ncisq, double),
		    diagj, one = 1., zero = 0.;

		AZERO(omgi, ncisq);
		for (j = 0; j < nci; j++) {
		    double dd = cc[cind++];
		    tmp[j * ncip1] = diagj =
			(ptyp == 1) ? exp(dd/2.) : sqrt(1./dd);
		    for (k = j + 1; k < nci; k++) {
			tmp[k*nci + j] = cc[odind++] * diagj;
		    }
		}
		F77_CALL(dsyrk)("U", "T", &nci, &nci, &one,
				tmp, &nci, &zero, omgi, &nci);
		Free(tmp);
	    } else {
		for (j = 0; j < nci; j++) {
		    omgi[j * ncip1] = cc[cind++];
		    for (k = j + 1; k < nci; k++) {
			omgi[k*nci + j] = cc[odind++];
		    }
		}
	    }
	    cind = odind;
	}
    }
    status[0] = status[1] = 0;
}

/**
 * Assign the upper triangles of the Omega matrices according to a
 * vector of parameters.
 *
 * @param x pointer to an lme object
 * @param coef pointer to an numeric vector of appropriate length
 * @param pType pointer to an integer scalar 
 *
 * @return R_NilValue
 */
SEXP mer2_coefGets(SEXP x, SEXP coef, SEXP pType)
{
    int clen = coef_length(LENGTH(GET_SLOT(x, Matrix_flistSym)),
			   INTEGER(GET_SLOT(x, Matrix_ncSym)));   
    if (LENGTH(coef) != clen || !isReal(coef))
	error(_("coef must be a numeric vector of length %d"), clen);
    internal_mer2_coefGets(x, REAL(coef), asInteger(pType));
    return x;
}
