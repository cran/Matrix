#include "lmeRep.h"

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

static
SEXP lmeRep_alloc3Darray(int TYP, int nr, int nc, int nf)
{
    SEXP val, dd = PROTECT(allocVector(INTSXP, 3));
    
    INTEGER(dd)[0] = nr; INTEGER(dd)[1] = nc; INTEGER(dd)[2] = nf;
    val = allocArray(TYP, dd);
    UNPROTECT(1);
    return val;
}

SEXP lmeRep_validate(SEXP x)
{
    /* FIXME: add checks for correct dimensions, modes, etc. */
    return ScalarLogical(1);
}

SEXP
lmeRep_crosstab(SEXP facs)
{
    int I = length(facs), nobs, pos = 0;
    int Ic2 = (I * (I - 1))/2;	/* I choose 2 */
    SEXP fac, levs, val = PROTECT(allocVector(VECSXP, Ic2));
    
    if (!isNewList(facs))
	error("Argument facs must be a list");
    nobs = length(VECTOR_ELT(facs, 0));

    if (Ic2 > 0) {
	int i, j, k, *nlevs = Calloc(I, int), *itmp = Calloc(nobs, int),
	    *zb = Calloc(nobs * I, int); /* zero-based indices */
	double *xtmp = Calloc(nobs, double), *atmp = Calloc(nobs, double);

	for (i = 0; i < nobs; i++) xtmp[i] = 1.;
	for (i = 0; i < I; i++) {
	    fac = VECTOR_ELT(facs, i);
	    if (!isFactor(fac) || length(fac) <= 0)
		error("All elements of facs must be nonnull factors");
	    if (length(fac) != nobs)
		error("All elements of facs must have the same length");
	    nlevs[i] = length(getAttrib(fac, R_LevelsSymbol));
	    for (k = 0; k < nobs; k++) zb[i * nobs + k] = INTEGER(fac)[k] - 1;
	    for (j = 0; j < i; j++) {
		SEXP mm;
		int *Dims, *mp, *mi, nz;
		double *mx;

		SET_VECTOR_ELT(val, pos, NEW_OBJECT(MAKE_CLASS("cscMatrix")));
		mm = VECTOR_ELT(val, pos);
		SET_SLOT(mm, Matrix_DimSym, allocVector(INTSXP, 2));
		Dims = INTEGER(GET_SLOT(mm, Matrix_DimSym));
		Dims[0] = nlevs[i]; Dims[1] = nlevs[j];
		SET_SLOT(mm, Matrix_pSym, allocVector(INTSXP, nlevs[j] + 1));
		mp = INTEGER(GET_SLOT(mm, Matrix_pSym));
		triplet_to_col(nlevs[i], nlevs[j], nobs,
			       zb + i * nobs, zb + j * nobs,
			       xtmp, mp, itmp, atmp);
		nz = mp[nlevs[j]];
		SET_SLOT(mm, Matrix_iSym, allocVector(INTSXP, nz));
		mi = INTEGER(GET_SLOT(mm, Matrix_iSym));
		SET_SLOT(mm, Matrix_xSym, allocVector(REALSXP, nz));
		mx = REAL(GET_SLOT(mm, Matrix_xSym));
		for (k = 0; k < nz; k++) {
		    mx[k] = atmp[k];
		    mi[k] = itmp[k];
		}
		pos++;
	    }
	}
	Free(nlevs); Free(itmp); Free(xtmp); Free(atmp); Free(zb);
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Create an lmeRep object from a list of one grouping factor and
 * an integer vector of length two giving the number of columns in the
 * model matrices. The last element is the number of columns in the
 * model matrix for the fixed effects plus the response.  (i.e. p+1)
 * 
 * @param facs pointer to a list of grouping factors
 * @param ncv pointer to an integer vector of number of columns per model matrix
 * 
 * @return pointer to an lmeRep object
 */
SEXP
lmeRep_create(SEXP facs, SEXP ncv)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("lmeRep")));
    int *nc, I = length(facs), i, nzcol;
    SEXP fac, levs, nms, tmp;

    if (!isNewList(facs) || length(facs) != 1)
	error("Argument facs must be a list of length 1");
    if (!isInteger(ncv) || length(ncv) != 2)
	error("Argument ncv must be an integer vector of length 2");
    if (INTEGER(ncv)[0] <= 0 || INTEGER(ncv)[1] <= 0)
	error("Number of columns in model matrices must be positive");

    SET_SLOT(val, Matrix_ncSym, allocVector(INTSXP, I + 2));
    nc = INTEGER(GET_SLOT(val, Matrix_ncSym));
    for (i = 0; i <= I; i++) nc[i] = INTEGER(ncv)[i];
    fac = VECTOR_ELT(facs, 0);
    if (!isFactor(fac) || length(fac) <= 0)
	error("All elements of facs must be nonnull factors");
    nc[I + 1] = length(fac); /* number of observations */
				/* allocate slots that are lists */
    nms = getAttrib(facs, R_NamesSymbol);
    SET_SLOT(val, R_LevelsSymbol, allocVector(VECSXP, I));
    levs = GET_SLOT(val, R_LevelsSymbol);
    setAttrib(levs, R_NamesSymbol, nms);
    SET_SLOT(val, Matrix_cnamesSym, allocVector(VECSXP, I + 1));
    tmp = PROTECT(allocVector(STRSXP, I + 1));
    for (i = 0; i < I; i++)
	SET_VECTOR_ELT(tmp, i, VECTOR_ELT(nms, i));
    SET_VECTOR_ELT(tmp, I, mkChar(".fixed"));
    setAttrib(GET_SLOT(val, Matrix_cnamesSym), R_NamesSymbol, duplicate(tmp));
    UNPROTECT(1);
    SET_SLOT(val, Matrix_DSym, allocVector(VECSXP, I));
    setAttrib(GET_SLOT(val, Matrix_DSym), R_NamesSymbol, duplicate(nms));
    SET_SLOT(val, Matrix_DIsqrtSym, allocVector(VECSXP, 0));
/*     setAttrib(GET_SLOT(val, Matrix_DIsqrtSym), R_NamesSymbol, duplicate(nms)); */
    SET_SLOT(val, Matrix_OmegaSym, allocVector(VECSXP, I));
    setAttrib(GET_SLOT(val, Matrix_OmegaSym), R_NamesSymbol, duplicate(nms));
    SET_SLOT(val, Matrix_ZZxSym, allocVector(VECSXP, (I * (I + 1))/2));
    nzcol = 0;
    for (i = 0; i < I; i++) {	/* check factors and fill lists */
	int nci = nc[i], nlev;
	SEXP LL;

	fac = VECTOR_ELT(facs, i);
	if (!isFactor(fac) || length(fac) != nc[I + 1])
	    error("All elements of facs must be factors of length %d",
		  nc[I + 1]);
	LL = getAttrib(fac, R_LevelsSymbol);
	SET_VECTOR_ELT(levs, i, LL);
	nlev = length(LL);
	nzcol += nlev * nci;
	SET_VECTOR_ELT(GET_SLOT(val, Matrix_OmegaSym), i,
		   allocMatrix(REALSXP, nci, nci));
	SET_VECTOR_ELT(GET_SLOT(val, Matrix_DSym), i,
		       lmeRep_alloc3Darray(REALSXP, nci, nci, nlev));
/* 	SET_VECTOR_ELT(GET_SLOT(val, Matrix_DIsqrtSym), i,	 */
/* 		       lmeRep_alloc3Darray(REALSXP, nci, nci, nlev)); */
	SET_VECTOR_ELT(GET_SLOT(val, Matrix_ZZxSym), (i * (i + 1))/2,
		       lmeRep_alloc3Darray(REALSXP, nci, nci, nlev));
    }
				/* Create slots */
    SET_SLOT(val, Matrix_XtXSym, allocMatrix(REALSXP, nc[I], nc[I]));
    SET_SLOT(val, Matrix_RXXSym, allocMatrix(REALSXP, nc[I], nc[I]));
    SET_SLOT(val, Matrix_ZtXSym, allocMatrix(REALSXP, nzcol, nc[I]));
    SET_SLOT(val, Matrix_RZXSym, allocMatrix(REALSXP, nzcol , nc[I]));
				/* Zero symmetric matrices (cosmetic) */
    memset(REAL(GET_SLOT(val, Matrix_XtXSym)), 0,
	   sizeof(double) * nc[I] * nc[I]); 
    memset(REAL(GET_SLOT(val, Matrix_RXXSym)), 0,
	   sizeof(double) * nc[I] * nc[I]);
				/*  3-D arrays */
    SET_SLOT(val, Matrix_devianceSym, allocVector(REALSXP, 2));
    tmp = GET_SLOT(val, Matrix_devianceSym);
    REAL(tmp)[0] = REAL(tmp)[1] = NA_REAL;
    setAttrib(tmp, R_NamesSymbol, allocVector(STRSXP, 2));
    nms = getAttrib(tmp, R_NamesSymbol);
    SET_STRING_ELT(nms, 0, mkChar("ML"));
    SET_STRING_ELT(nms, 1, mkChar("REML"));
    SET_SLOT(val, Matrix_devCompSym, allocVector(REALSXP, 4));
    tmp = GET_SLOT(val, Matrix_devCompSym);
    REAL(tmp)[0] = REAL(tmp)[1] = REAL(tmp)[2] = REAL(tmp)[3] = NA_REAL;
    SET_SLOT(val, Matrix_statusSym, allocVector(LGLSXP, 2));
    tmp = GET_SLOT(val, Matrix_statusSym);
    LOGICAL(tmp)[0] = LOGICAL(tmp)[1] = 0;
    setAttrib(tmp, R_NamesSymbol, allocVector(STRSXP, 2));
    nms = getAttrib(tmp, R_NamesSymbol);
    SET_STRING_ELT(nms, 0, mkChar("factored"));
    SET_STRING_ELT(nms, 1, mkChar("inverted"));
    UNPROTECT(1);
    return val;
}

/** 
 * Update the arrays ZZx, ZtX, and XtX in an lme object
 * according to a list of factors and a list of model matrices.
 * 
 * @param x pointer to an object inheriting from lmeCommon
 * @param facs pointer to a list of grouping factors
 * @param mmats pointer to a list of model matrices
 * 
 * @return NULL
 */
SEXP lmeRep_update_mm(SEXP x, SEXP facs, SEXP mmats)
{
    SEXP ZZxP = GET_SLOT(x, Matrix_ZZxSym),
	ZtXP = GET_SLOT(x, Matrix_ZtXSym),
	levs = GET_SLOT(x, R_LevelsSymbol),
	cnames = GET_SLOT(x, Matrix_cnamesSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	I = length(ZZxP), Ip1 = I + 1,
	i, ione = 1,
	nobs = nc[Ip1],
	ncZ = 0,
	pp1 = nc[I];
    double
	*X,
	*XtX = REAL(GET_SLOT(x, Matrix_XtXSym)),
	*ZtX = REAL(ZtXP),
	one = 1.0, zero = 0.0;

    if (!isNewList(facs) || length(facs) != I)
	error("facs must be a list of %d factors", I);
    if (!isNewList(mmats) || length(mmats) != Ip1)
	error("mmats must be a list of %d model matrices", Ip1);
    for (i = 0; i <= I; i++) {
	SEXP mmat = VECTOR_ELT(mmats, i);
	int *dims = INTEGER(getAttrib(mmat, R_DimSymbol));
	
	if (!isMatrix(mmat) || !isReal(mmat))
	    error("element %d of mmats is not a numeric matrix", i + 1);
	if (nobs != dims[0])
	    error("Expected %d rows in the %d'th model matrix. Got %d",
		  nobs, i+1, dims[0]);
	if (nc[i] != dims[1])
	    error("Expected %d columns in the %d'th model matrix. Got %d",
		  nc[i], i+1, dims[1]);
	SET_VECTOR_ELT(cnames, i,
		       duplicate(VECTOR_ELT(getAttrib(mmat, R_DimNamesSymbol),
					    1)));
    }
    for (i = 0; i < I; i++) {
	SEXP fac = VECTOR_ELT(facs, i),
	    lev = getAttrib(fac, R_LevelsSymbol);

	if (!isFactor(fac))
	    error("element %i in list facs is not a factor", i + 1);
	SET_VECTOR_ELT(levs, i, duplicate(lev));
	ncZ += length(lev) * nc[i];
    }
    if (ncZ != INTEGER(getAttrib(ZtXP, R_DimSymbol))[0])
	error("# rows of ZtX slot, %d, != sum of # levels * # columns, %d",
	      INTEGER(getAttrib(ZtXP, R_DimSymbol))[0], ncZ);
				/* Create XtX */
    X = REAL(VECTOR_ELT(mmats, I));
    F77_CALL(dsyrk)("U", "T", &pp1, &nobs, &one, X, &nobs, &zero, XtX, nc + 1);
				/* Zero an accumulator */
    memset((void *) ZtX, 0, sizeof(double) * pp1 * ncZ);
    for (i = 0; i < I; i++) {
	int *fac = INTEGER(VECTOR_ELT(facs, i)),
	    j, nci = nc[i], ncisqr = nci * nci,
	    nlev = length(VECTOR_ELT(levs, i));
	int ZtXrows = nci * nlev;
	double *Z = REAL(VECTOR_ELT(mmats, i)),
	    *ZZx = REAL(VECTOR_ELT(ZZxP, i));

	memset((void *) ZZx, 0, sizeof(double) * nci * nci * nlev);
	if (nci == 1) {		/* single column in Z */
	    for (j = 0; j < nobs; j++) {
		int fj = fac[j] - 1; /* factor indices are 1-based */
		ZZx[fj] += Z[j] * Z[j];
		F77_CALL(daxpy)(&pp1, Z + j, X + j, &nobs, ZtX + fj, &nlev);
	    }
	} else {
	    for (j = 0; j < nobs; j++) {
		int fj = fac[j] - 1; /* factor indices are 1-based */

		F77_CALL(dsyr)("U", nc + i, &one, Z + j, &nobs,
			       ZZx + fj * ncisqr, nc + i);
		F77_CALL(dgemm)("T", "N", nc + i, &pp1, &ione,
				&one, Z + j, &nobs,
				X + j, &nobs, &one,
				ZtX + fj * nci, &ZtXrows);
	    }
	}
	ZtX += ZtXrows;
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}


/** 
 * Create and insert initial values for Omega.
 * 
 * @param x pointer to an lmeRep object
 * 
 * @return NULL
 */
SEXP lmeRep_initial(SEXP x)
{
    SEXP ZZxP = VECTOR_ELT(GET_SLOT(x, Matrix_ZZxSym), 0);
    int *dims = INTEGER(getAttrib(ZZxP, R_DimSymbol)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	j, k, nzc = dims[0], nlev = dims[2], nzcsqr = nzc * nzc,
	nzcp1 = nzc + 1;
    double
	*Omega = REAL(VECTOR_ELT(GET_SLOT(x, Matrix_OmegaSym), 0)),
	*ZZx = REAL(ZZxP),
	mi = 0.375 / ((double) nlev);
    
    memset((void *) Omega, 0, sizeof(double) * nzc * nzc);
    for (j = 0; j < nlev; j ++) {
	for (k = 0; k < nzc; k++) {
	    Omega[k * nzcp1] += ZZx[k * nzcp1 + j * nzcsqr] * mi;
	}
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}

/** 
 * If status[["factored"]] is FALSE, create and factor Z'Z+Omega, then
 * create RZX and RXX, the deviance components, and the value of the
 * deviance for both ML and REML.
 * 
 * @param x pointer to an lmeRep object
 * 
 * @return NULL
 */
SEXP lmeRep_factor(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    
    if (!status[0]) {
	SEXP
	    levs = GET_SLOT(x, R_LevelsSymbol),
	    RZXsl = GET_SLOT(x, Matrix_RZXSym);
	int *dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    i, j, nf = length(levs);
	int nml = nc[nf + 1], nreml = nml + 1 - nc[nf];
	double
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(GET_SLOT(x, Matrix_RZXSym)),
	    *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	    *deviance = REAL(GET_SLOT(x, Matrix_devianceSym)),
	    minus1 = -1., one = 1.;


	Memcpy(RZX, REAL(GET_SLOT(x, Matrix_ZtXSym)), dims[0] * dims[1]);
	dcmp[0] = dcmp[1] = dcmp[2] = dcmp[3] = 0.;
	for (i = 0; i < nf; i++) {
	    int ii, jj, nci = nc[i], ncisqr = nci * nci,
		nlev = length(VECTOR_ELT(levs, i)),
		RZXrows = nlev * nci;
	    double
		*ZZx = REAL(VECTOR_ELT(GET_SLOT(x, Matrix_ZZxSym), i)),
		*D = REAL(VECTOR_ELT(GET_SLOT(x, Matrix_DSym), i)),
		*Omega = REAL(VECTOR_ELT(GET_SLOT(x, Matrix_OmegaSym), i));
	
	    if (nci == 1) {
		dcmp[1] += nlev * log(Omega[0]);
		for (j = 0; j < nlev; j++) { /* inflate and factor */
		    D[j] = sqrt(ZZx[j] + Omega[0]);
		    dcmp[0] += 2. * log(D[j]);
		    for (jj = 0; jj < dims[1]; jj++) RZX[j + jj*dims[0]] /= D[j];
		}
	    } else {
		double *tmp = Memcpy(Calloc(ncisqr, double), Omega, ncisqr);
		
		F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
		if (j)
		    error("Leading minor of size %d of Omega[[%d]] is not positive definite",
			  j, i + 1);
		for (j = 0; j < nci; j++) { /* nlev * logDet(Omega_i) */
		    dcmp[1] += nlev * 2. * log(tmp[j * (nci + 1)]);
		}
		Free(tmp);
	    
		for (j = 0; j < nlev; j++) {
		    double
			*Dj = D + j * ncisqr,
			*ZZxj = ZZx + j * ncisqr;

		    for (jj = 0; jj < nci; jj++) { /* Copy ZZx to Dj and inflate */
			for (ii = 0; ii <= jj; ii++) {
			    Dj[ii + jj * nci]
				= ZZxj[ii + jj * nci] + Omega[ii + jj * nci];
			}
			for (ii = jj + 1; ii < nci; ii++)
			    Dj[ii + jj * nci] = 0.;
		    }
		    F77_CALL(dpotrf)("U", &nci, Dj, &nci, &jj);
		    if (jj) 
			error("D[ , , %d] is not positive definite", j + 1);
		    for (jj = 0; jj < nci; jj++) /* accumulate determinant */
			dcmp[0] += 2. * log(Dj[jj * (nci + 1)]);
				/* Update RZX */
		    F77_CALL(dtrsm)("L", "U", "T", "N", &nci, dims + 1,
				    &one, Dj, &nci, RZX + j * nci, dims);
		    /* FIXME: Need another step here for L in the general case */
		}
	    }
	    RZX += RZXrows;
	}
				/* downdate and factor X'X */
	Memcpy(RXX, REAL(GET_SLOT(x, Matrix_XtXSym)), dims[1] * dims[1]);
	F77_CALL(dsyrk)("U", "T", dims + 1, dims,
			&minus1, REAL(RZXsl), dims,
			&one, RXX, dims + 1);
	F77_CALL(dpotrf)("U", dims + 1, RXX, dims + 1, &j);
	if (j) {
	    warning("Leading minor of size %d of downdated X'X is indefinite",
		    j + 1);
	    dcmp[2] = dcmp[3] = deviance[0] = deviance[1] = NA_REAL;
	} else {
	    for (j = 0; j < (dims[1] - 1); j++) /* 2 logDet(RXX) */
		dcmp[2] += 2 * log(RXX[j * (dims[1] + 1)]);
	    dcmp[3] = 2. * log(RXX[dims[1] * dims[1] - 1]); /* 2 log(ryy) */
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
 * If necessary, factor Z'Z+Omega, ZtX, and XtX then, if necessary,
 * replace the RZX and RXX slots by the corresponding parts of the
 * inverse of the Cholesky factor.  Replace the elements of the D slot
 * by the blockwise inverses.
 * 
 * @param x pointer to an lme object
 * 
 * @return NULL (x is updated in place)
 */
SEXP lmeRep_invert(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    if (!status[0]) lmeRep_factor(x);
    if (!R_FINITE(REAL(GET_SLOT(x, Matrix_devianceSym))[0]))
	error("Unable to invert singular factor of downdated X'X");
    if (!status[1]) {
	SEXP
	    RZXsl = GET_SLOT(x, Matrix_RZXSym),
	    Dsl = GET_SLOT(x, Matrix_DSym),
	    levs = GET_SLOT(x, R_LevelsSymbol);
	int
	    *dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    i, nf = length(levs);
	double
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(RZXsl),
	     minus1 = -1., one = 1.;

	F77_CALL(dtrtri)("U", "N", dims + 1, RXX, dims + 1, &i);
	if (i)
	    error("Leading minor of size %d of downdated X'X,is indefinite",
		    i + 1);
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims + 1, &minus1,
			RXX, dims + 1, RZX, dims);
	for(i = 0; i < nf; i++) {
	    int info, j, jj, nci = nc[i], ncisqr = nci * nci,
		nlev = length(VECTOR_ELT(levs, i));
	    double *Di = REAL(VECTOR_ELT(Dsl, i));
	    
	    if (nci == 1) {
		for (j = 0; j < nlev; j++) {
		    Di[j] = 1./Di[j];
		    for (jj = 0; jj < dims[1]; jj++)
			RZX[j + jj * dims[0]] *= Di[j];
		}
	    } else {
		for (j = 0; j < nlev; j++) {
		    F77_CALL(dtrtri)("U", "N", &nci, Di + j * ncisqr, &nci, &info);
		    if (info)
			error("D[,,%d] for factor %d is singular", j + 1, i + 1);
		    F77_CALL(dtrmm)("L", "U", "N", "N", &nci, dims + 1, &one,
				    Di + j * ncisqr, &nci, RZX + j * nci, dims);
		}
	    }
	    RZX += nci * nlev;
	}
	status[1] = 1;
    }
    return R_NilValue;
}

/** 
 * Extract the ML or REML conditional estimate of sigma
 * 
 * @param x pointer to an lme object
 * @param REML logical scalar - TRUE if REML estimates are requested
 * @param nf - number of grouping factors
 * 
 * @return pointer to a numeric scalar 
 */
SEXP lmeRep_sigma(SEXP x, SEXP REML)
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym);
    int pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1],
	nobs = INTEGER(GET_SLOT(x, Matrix_ncSym))
	[length(GET_SLOT(x, Matrix_OmegaSym)) + 1];

    if (!LOGICAL(GET_SLOT(x, Matrix_statusSym))[1])
	error("Call the 'solve' method before calling lmeRep_sigma");
    return ScalarReal(1./(REAL(RXXsl)[pp1*pp1 - 1] *
			  sqrt((double)(asLogical(REML) ?
					nobs + 1 - pp1 : nobs))));
}

/** 
 * Extract the upper triangles of the Omega matrices.  These aren't
 * "coefficients" but the extractor is called coef for historical
 * reasons.  Within each group these values are in the order of the
 * diagonal entries first then the strict upper triangle in row
 * order.
 * 
 * @param x pointer to an lme object
 * @param Unc pointer to a logical scalar indicating if the parameters
 * are in the unconstrained form.
 * 
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP lmeRep_coef(SEXP x, SEXP Unc)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), unc = asLogical(Unc), vind;
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double *vv = REAL(val);

    vind = 0;			/* index in vv */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1;
	if (nci == 1) {
	    vv[vind++] = (unc ?
			  log(REAL(VECTOR_ELT(Omega, i))[0]) :
			  REAL(VECTOR_ELT(Omega, i))[0]);
	} else {
	    if (unc) {		/* L log(D) L' factor of Omega[,,i] */
		int j, k, ncisq = nci * nci;
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
	    } else {		/* upper triangle of Omega[,,i] */
		int j, k, odind = vind + nci;
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
    }
    UNPROTECT(1);
    return val;
}


/** 
 * Assign the upper triangles of the Omega matrices.
 * (Called coef for historical reasons.)
 * 
 * @param x pointer to an lme object
 * @param coef pointer to an numeric vector of appropriate length
 * @param Unc pointer to a logical scalar indicating if the parameters
 * are in the unconstrained form.
 * 
 * @return R_NilValue
 */
SEXP lmeRep_coefGets(SEXP x, SEXP coef, SEXP Unc)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	cind, i, nf = length(Omega),
	unc = asLogical(Unc);
    double *cc = REAL(coef);

    if (length(coef) != coef_length(nf, nc) || !isReal(coef))
	error("coef must be a numeric vector of length %d",
	      coef_length(nf, nc));
    cind = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	if (nci == 1) {
	    REAL(VECTOR_ELT(Omega, i))[0] = (unc ?
					     exp(cc[cind++]) :
					     cc[cind++]);
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		j, k,
		ncip1 = nci + 1,
		ncisq = nci * nci;
	    double
		*omgi = REAL(VECTOR_ELT(Omega, i));
	    if (unc) {
		double
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
    return x;
}

/** 
 * Extract the conditional estimates of the fixed effects
 * 
 * @param x Pointer to an lme object
 * 
 * @return a numeric vector containing the conditional estimates of
 * the fixed effects
 */
SEXP lmeRep_fixef(SEXP x)
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym),
	cnames = GET_SLOT(x, Matrix_cnamesSym);
    int pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1];
    int j;
    SEXP val = PROTECT(allocVector(REALSXP, pp1));
    double
	*beta = REAL(val),
	nryyinv;		/* negative ryy-inverse */

    if (!LOGICAL(GET_SLOT(x, Matrix_statusSym))[1])
	error("Call the 'solve' method before calling lmeRep_fixef");
    Memcpy(beta, REAL(RXXsl) + pp1 * (pp1 - 1), pp1);
    nryyinv = -REAL(RXXsl)[pp1*pp1 - 1];
    for (j = 0; j < pp1; j++) beta[j] /= nryyinv;
    setAttrib(val, R_NamesSymbol, VECTOR_ELT(cnames, length(cnames) - 1));
    UNPROTECT(1);
    return val;
}

/** 
 * Extract the conditional modes of the random effects.
 * 
 * @param x Pointer to an lme object
 * 
 * @return a vector containing the conditional modes of the random effects
 */
SEXP lmeRep_ranef(SEXP x)
{
    SEXP RZXsl = GET_SLOT(x, Matrix_RZXSym),
	cnames = GET_SLOT(x, Matrix_cnamesSym),
	levs = GET_SLOT(x, R_LevelsSymbol);
    int *dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ii, jj,
	nf = length(levs);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double
	*b = REAL(RZXsl) + dims[0] * (dims[1] - 1),
	nryyinv;		/* negative ryy-inverse */
    /* FIXME: Check if this really should be negative ryy-inverse now
     * that the negative sign is used in the inversion of the RZX slot */

    if (!LOGICAL(GET_SLOT(x, Matrix_statusSym))[1])
	error("Call the 'solve' method before calling lmeRep_ranef");
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(GET_SLOT(x, Matrix_OmegaSym), R_NamesSymbol)));
    nryyinv = -REAL(GET_SLOT(x, Matrix_RXXSym))[dims[1] * dims[1] - 1];
    for (i = 0; i < nf; i++) {
	SEXP nms, rnms = VECTOR_ELT(levs, i);
	int nci = nc[i], mi = length(rnms), Mi = mi * nci;
	double *mm;
	
	SET_VECTOR_ELT(val, i, allocMatrix(REALSXP, mi, nci));
	setAttrib(VECTOR_ELT(val, i), R_DimNamesSymbol, allocVector(VECSXP, 2));
	nms = getAttrib(VECTOR_ELT(val, i), R_DimNamesSymbol);
	SET_VECTOR_ELT(nms, 0, duplicate(rnms));
	SET_VECTOR_ELT(nms, 1, duplicate(VECTOR_ELT(cnames, i)));
	mm = REAL(VECTOR_ELT(val, i));
	for (jj = 0; jj < nci; jj++)
	    for(ii = 0; ii < mi; ii++)
		mm[ii + jj * mi] = b[jj + ii * nci]/nryyinv;
	b += Mi;
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
 * @param x pointer to an lme object
 * @param val pointer to a list of matrices of the correct sizes
 * 
 * @return val
 */
static
SEXP lmeRep_firstDer(SEXP x, SEXP val)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	D = GET_SLOT(x, Matrix_DSym),
	RZXsl = GET_SLOT(x, Matrix_RZXSym),
	levels = GET_SLOT(x, R_LevelsSymbol);
    int *dims = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), p = dims[1] - 1;
    double *RZX = REAL(RZXsl),
	*b = REAL(RZXsl) + dims[0] * p;

    lmeRep_invert(x);
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(Omega, R_NamesSymbol)));
    for (i = 0; i < nf; i++) {
	int j, nlev = length(VECTOR_ELT(levels, i)),
	    nci = nc[i], ncisqr = nci * nci, RZXrows = nlev * nci;
	    
	double *Di = REAL(VECTOR_ELT(D, i)), *mm, dlev = (double) nlev;
	
	mm = memset(REAL(VECTOR_ELT(val, i)), 0, sizeof(double) * ncisqr * 4);
	
 	if (nci == 1) {
	    int ione = 1;
 	    mm[0] = ((double) nlev)/REAL(VECTOR_ELT(Omega, i))[0]; 
 	    mm[1] = F77_CALL(ddot)(&nlev, b, &ione, b, &ione); 
	    mm[2] = F77_CALL(ddot)(&nlev, Di, &ione, Di, &ione);
  	    for (j = 0; j < p; j++) {
  		mm[3] += F77_CALL(ddot)(&RZXrows, RZX + j * dims[0], &ione,
					RZX + j * dims[0], &ione);
  	    }
 	} else {
	    double *tmp = Memcpy(Calloc(ncisqr, double),
				 REAL(VECTOR_ELT(Omega, i)), ncisqr),
		one = 1., zero = 0.;

	    F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
	    if (j)
		error("Omega[[%d]] is not positive definite", i + 1);
	    F77_CALL(dtrtri)("U", "N", &nci, tmp, &nci, &j);
	    if (j)
		error("Omega[[%d]] is not positive definite", i + 1);
	    F77_CALL(dsyrk)("U", "N", &nci, &nci, &dlev, tmp, &nci,
			    &zero, mm, &nci);
	    Free(tmp);
	    mm += ncisqr;	/* \bB_i term */
	    F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one, b, &nci,
			    &zero, mm, &nci);
	    mm += ncisqr;     /* Sum of inverses of diagonal blocks */
	    F77_CALL(dsyrk)("U", "N", &nci, &RZXrows, &one, Di, &nci,
			    &zero, mm, &nci);
	    mm += ncisqr;	/* Extra term for \vb */
	    for (j = 0; j < p; j++) {
		F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one,
				RZX + j * dims[0], &nci,
				&one, mm, &nci);
	    }
	}
	RZX += RZXrows;
	b += RZXrows;
    }
    return val;
}

/** 
 * Generate and zero the contents of a list of length nf of arrays of
 * dimension (nci, nci, 4).  The values of these arrays are assigned
 * in lmeRep_firstDer.
 * 
 * @param nf number of factors
 * @param nc vector of number of columns per factor
 * 
 * @return pointer to a list of REAL arrays
 */
static
SEXP EM_grad_array(int nf, const int nc[])
{
    SEXP val = PROTECT(allocVector(VECSXP, nf)),
	dd = PROTECT(allocVector(INTSXP, 3));
    int *dims = INTEGER(dd), i;

    dims[2] = 4;
    for (i = 0; i < nf; i++) {
	dims[0] = dims[1] = nc[i];
	SET_VECTOR_ELT(val, i, allocArray(REALSXP, dd));
	memset(REAL(VECTOR_ELT(val, i)), 0, sizeof(double) * nc[i] * nc[i]);
    }
    UNPROTECT(2);
    return val;
}
	
/** 
 * Fill in the 4-dimensional vector of linear combinations of the
 * firstDer array according to whether ECME steps of the gradient are
 * needed and to whether or not REML is being used.
 * 
 * @param cc coefficient vector to be filled in
 * @param EM non-zero for ECME steps, zero for gradient
 * @param REML non-zero for REML, zero for ML
 * @param ns ns[0] is p+1, ns[1] is n
 * 
 * @return cc with the coefficients filled in
 */
static
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
 */
static
void EMsteps_verbose_print(SEXP x, int iter, int REML, SEXP firstDer, SEXP val)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	pMat = VECTOR_ELT(val, 2);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Its = INTEGER(VECTOR_ELT(val, 0)),
	i, ifour = 4, ii, ione = 1, jj, nf = length(Omega),
	niter = INTEGER(getAttrib(pMat, R_DimSymbol))[0];
    double
	*dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	*cc = EM_grad_lc(Calloc(4, double), 0, REML, nc + nf),
	*Devs = REAL(VECTOR_ELT(val, 1)),
	*pars = REAL(pMat) + iter,
	*grds = REAL(VECTOR_ELT(val, 3)) + iter,
	one = 1., zero = 0.;
			       
    /* FIXME: Check with MM for format. */
    lmeRep_factor(x);
    if (iter == 0) Rprintf("  EM iterations\n");
    Rprintf("%3d %.3f", Its[iter] = iter, Devs[iter] = dev[REML ? 1 : 0]);
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1, ncisqr = nci * nci;
	double
	    *Omgi = REAL(VECTOR_ELT(Omega, i)),
	    *Grad = Calloc(ncisqr, double);
	
				/* diagonals */
	for (jj = 0; jj < nci; jj++, pars += niter) {
	    Rprintf(" %#8g", *pars = Omgi[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++, pars += niter)
		Rprintf(" %#8g", *pars = Omgi[ii + jj * nci]);
				/* Evaluate and print the gradient */
	F77_CALL(dgemv)("N", &ncisqr, &ifour, &one,
			REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			cc, &ione, &zero, Grad, &ione);
	Rprintf(":");
				/* diagonals */
	for (jj = 0; jj < nci; jj++, grds += niter) {
	    Rprintf(" %#8g", *grds = Grad[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++, grds += niter)
		Rprintf(" %#8g", *grds = Grad[ii + jj * nci]);
	Free(Grad);
    }
    Rprintf("\n");
    Free(cc);
}

/** 
 * Perform ECME steps for the REML or ML criterion.
 * 
 * @param x pointer to an ssclme object
 * @param nsteps pointer to an integer scalar - the number of ECME steps to perform
 * @param REMLp pointer to a logical scalar indicating if REML is to be used
 * @param verb pointer to a logical scalar indicating verbose output
 * 
 * @return R_NilValue if verb == FALSE, otherwise a list of iteration
 *numbers, deviances, parameters, and gradients.
 */
SEXP lmeRep_ECMEsteps(SEXP x, SEXP nsteps, SEXP REMLp, SEXP Verbp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	levs = GET_SLOT(x, R_LevelsSymbol),
	val = R_NilValue;
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	REML = asLogical(REMLp),
	i, ifour = 4, info, ione = 1, iter,
	nEM = asInteger(nsteps),
	nf = length(Omega),
	verb = asLogical(Verbp);
    double
	*cc = EM_grad_lc(Calloc(4, double), 1, REML, nc + nf),
	zero = 0.0;
    SEXP firstDer = PROTECT(EM_grad_array(nf, nc));

    if (verb) {
	int nEMp1 = nEM + 1, npar = coef_length(nf, nc);
	val = PROTECT(allocVector(VECSXP, 4));
	SET_VECTOR_ELT(val, 0, allocVector(INTSXP, nEMp1));
	SET_VECTOR_ELT(val, 1, allocVector(REALSXP, nEMp1));
	SET_VECTOR_ELT(val, 2, allocMatrix(REALSXP, nEMp1, npar));
	SET_VECTOR_ELT(val, 3, allocMatrix(REALSXP, nEMp1, npar));
	EMsteps_verbose_print(x, 0, REML,
			      lmeRep_firstDer(x, firstDer), val);
    }
    for (iter = 0; iter < nEM; iter++) {
	lmeRep_firstDer(x, firstDer);
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], ncisqr = nci * nci;
	    double *Omgi = REAL(VECTOR_ELT(Omega, i)),
		mult = 1./((double) length(VECTOR_ELT(levs, i)));
	
	    F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &mult,
			    REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			    cc, &ifour, &zero, Omgi, &ncisqr);
	    F77_CALL(dpotrf)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error("DPOTRF in ECME update gave code %d", info);
	    F77_CALL(dpotri)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error("Matrix inverse in ECME update gave code %d", info);
	}
	status[0] = status[1] = 0;
	if (verb) EMsteps_verbose_print(x, iter + 1, REML, firstDer, val);
    }
    lmeRep_factor(x);
    if (verb) UNPROTECT(1);
    UNPROTECT(1);
    return val;
}	

SEXP lmeRep_gradient(SEXP x, SEXP REMLp, SEXP Uncp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	dind, i, ifour = 4, info, ione = 1, nf = length(Omega),
	odind, unc = asLogical(Uncp);
    SEXP
	firstDer = lmeRep_firstDer(x, PROTECT(EM_grad_array(nf, nc))),
	val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double
	*cc = EM_grad_lc(Calloc(4, double), 0,
			 asInteger(REMLp), nc + nf),
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
	    REAL(val)[dind++] = (unc ? Omgi[0] : 1.) * tmp[0];
	} else {
	    int ii, j, ncip1 = nci + 1;
	    
	    odind = dind + nci; /* index into val for off-diagonals */
	    if (unc) {
		double *chol = Memcpy(Calloc(ncisqr, double),
				      REAL(VECTOR_ELT(Omega, i)), ncisqr),
		    *tmp2 = Calloc(ncisqr, double);

		/* Overwrite the gradient with respect to positions in
		 * Omega[[i]] by the gradient with respect to the
		 * unconstrained parameters.*/

		F77_CALL(dpotrf)("U", &nci, chol, &nci, &info);
		if (info)
		    error("Omega[[%d]] is not positive definite", i + 1);
		/* tmp2 := chol %*% tmp using only upper triangle of tmp */
		F77_CALL(dsymm)("R", "U", &nci, &nci, &one, tmp, &nci,
				chol, &nci, &zero, tmp2, &nci);
		/* full symmetric product gives diagonals */
		F77_CALL(dtrmm)("R", "U", "T", "N", &nci, &nci, &one, chol, &nci,
				Memcpy(tmp, tmp2, ncisqr), &nci);
		/* overwrite upper triangle with gradients for positions in L' */
		for (ii = 1; ii < nci; ii++) { 
		    for (j = 0; j < ii; j++) {
			tmp[j + ii*nci] = chol[j*ncip1] * tmp2[j + ii*nci];
			tmp[ii + j*nci] = 0.;
		    }
		}
		Free(chol); Free(tmp2);
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
 * Fill in five symmetric matrices, providing the
 * information to generate the Hessian.

 * @param x pointer to an lme object
 * @param val ignored at present
 * 
 * @return val an array consisting of five symmetric faces
 */
SEXP lmeRep_secondDer(SEXP x, SEXP Valp)
{
    SEXP
	D = GET_SLOT(x, Matrix_DSym),
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	RZXsl = GET_SLOT(x, Matrix_RZXSym),
	levels = GET_SLOT(x, R_LevelsSymbol),
	val;
    int *dRZX = INTEGER(getAttrib(RZXsl, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	Q, Qsqr, RZXpos, facepos,
	i, ione = 1, j, nf = length(Omega), p = dRZX[1] - 1, pos;
    SEXP
	firstDer = lmeRep_firstDer(x, PROTECT(EM_grad_array(nf, nc)));
    double
	*RZX = REAL(RZXsl),
	*b = REAL(RZXsl) + dRZX[0] * p,
	*bbface,		/* vec of second faces of firstDer elts */
	one = 1.,
	zero = 0.;
    
    Q = 0;			/* number of rows and columns in the result */
    for (i = 0; i < nf; i++) Q += nc[i] * nc[i];
    Qsqr = Q * Q;
    bbface = Calloc(Q, double);
    val = PROTECT(lmeRep_alloc3Darray(REALSXP, Q, Q, 5));
    memset(REAL(val), 0, sizeof(double) * Qsqr * 5);
    
    pos = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncisqr = nci * nci;
	double *fDi = REAL(VECTOR_ELT(firstDer, i)),
	    mult = 1./((double) length(VECTOR_ELT(levels, i)));

	Memcpy(bbface + pos, fDi + ncisqr, ncisqr);
	/* outer product of the third face of firstDer on the diagonal
	 * of the third face of val */
	F77_CALL(dsyr)("U", &ncisqr, &mult, fDi + 2 * ncisqr, &ione,
		       REAL(val) + 2 * Qsqr + pos * Q, &Q);
	pos += ncisqr;
    }
				/* fifth face of val is outer product of bbface */
    F77_CALL(dsyr)("U", &Q, &one, bbface, &ione, REAL(val) + 4 * Qsqr, &Q);
				/* fourth face from \bb\trans\der\vb\der\bb */
    memset(REAL(val) + 3 * Qsqr, 0, sizeof(double) * Qsqr); /* zero accumulator */
    RZXpos = 0;
    facepos = 0;
    for (i = 0; i < nf; i++) {
	int ii, jj, nci = nc[i], ncisqr = nci * nci, nctp = nci * p,
	    nlev = length(VECTOR_ELT(levels, i));
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

