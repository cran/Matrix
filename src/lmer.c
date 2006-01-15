#include "lmer.h"
#include <float.h>

/* These should be moved to the R sources */

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


/* Internally used utilities */

#define flag_not_factored(x) *LOGICAL(GET_SLOT(x, Matrix_statusSym)) = 0

/** 
 * Create the diagonal blocks of the variance-covariance matrix of the
 * random effects
 * 
 * @param Linv - cholmod_sparse representation of L^{-1}
 * @param nf - number of grouping factors
 * @param Gp - group pointers
 * @param nc - number of columns per factor
 * @param bVar - list of 3-d arrays to be filled in
 * @param uplo - "U" or "L" for upper or lower triangle
 */
static void
Linv_to_bVar(cholmod_sparse *Linv, const int Gp[], const int nc[],
	     SEXP bVar, const char uplo[])
{
    int *Lii = (int*)(Linv->i), *Lip = (int*)(Linv->p), i, nf = LENGTH(bVar);
    double *Lix = (double*)(Linv->x), one[] = {1,0}, zero[] = {0,0};
    
    for (i = 0; i < nf; i++) {
	int *ind, j, nci = nc[i], maxnnz = 0;
	int ncisqr = nci * nci, nlev = (Gp[i + 1] - Gp[i])/nci;
	double *bVi = REAL(VECTOR_ELT(bVar, i)), *tmp;

	AZERO(bVi, nlev * ncisqr);
	for (j = 0; j < nlev; j++) {
	    int nzm = Lip[Gp[i] + (j + 1) * nci] - Lip[Gp[i] + j * nci];
	    if (nzm > maxnnz) maxnnz = nzm;
	}
	ind = Calloc(maxnnz, int);
	tmp = Calloc(maxnnz * nci, double);
	for (j = 0; j < nlev; j++) {
	    int jj, k, kk;
	    int *ap = Lip + Gp[i] + j * nci;
	    int nr = ap[1] - ap[0];

	    AZERO(tmp, maxnnz * nci);
	    Memcpy(ind, Lii + ap[0], nr);
	    Memcpy(tmp, Lix + ap[0], nr);
	    for (jj = 1; jj < nci; jj++) {
		for (k = ap[jj]; k < ap[jj + 1]; k++) {
		    int aik = Lii[k];
		    for (kk = 0; kk < nr; kk++) {
			if (aik == ind[kk]) {
			    tmp[kk + jj * maxnnz] = Lix[k];
			    aik = -1;
			    break;
			}
		    }
		    if (aik >= 0) {	/* did not find the row index */
			ind[nr] = aik;
			tmp[nr + jj * maxnnz] = Lix[k];
			nr++;
		    }
		}
	    }
	    F77_CALL(dsyrk)(uplo, "T", &nci, &nr, one, tmp, &maxnnz,
			    zero, bVi + j * ncisqr, &nci);
	}
	Free(ind); Free(tmp);
    }
}

static void
internal_mer_bVar(SEXP x)
{
    int q = LENGTH(GET_SLOT(x, Matrix_rZySym));
    cholmod_factor *L = as_cholmod_factor(GET_SLOT(x, Matrix_LSym));
    cholmod_sparse *eye = cholmod_speye(q, q, CHOLMOD_REAL, &c), *Linv;
    int *Perm = (int *)(L->Perm), *iperm = Calloc(q, int), i;
				/* create the inverse permutation */
    for (i = 0; i < q; i++) iperm[Perm[i]] = i;
				/* apply iperm to the identity matrix */
    for (i = 0; i < q; i++) ((int*)(eye->i))[i] = iperm[i];
				/* Create Linv */
    Linv = cholmod_spsolve(CHOLMOD_L, L, eye, &c);
    cholmod_free_sparse(&eye, &c);
    Linv_to_bVar(Linv, INTEGER(GET_SLOT(x, Matrix_GpSym)),
		 INTEGER(GET_SLOT(x, Matrix_ncSym)), 
		 GET_SLOT(x, Matrix_bVarSym), "U");
    cholmod_free_sparse(&Linv, &c);
    Free(L);
}

/** 
 * Evaluate the quadratic form in b defined by Omega
 * 
 * @param b vector of random effects
 * @param Omega - list of dpoMatrix objects defining the pattern for Omega
 * @param nf - number of grouping factors
 * @param Gp - group pointers
 * @param nc - number of columns per factor
 * 
 * @return 
 */
static double
b_quadratic(const double b[], SEXP Omega, const int Gp[], const int nc[])
{
    int i, ione = 1, nf = LENGTH(Omega);
    double ans = 0., one[] = {1.,0.};

    for (i = 0; i < nf; i++) {
	int nci = nc[i], ntot = Gp[i + 1] - Gp[i];
	int nlev = ntot/nci;
	double *bcp = Memcpy(Calloc(ntot, double), b + Gp[i], ntot),
	    *omgf = REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(Omega, i)), Matrix_xSym));

	F77_CALL(dtrmm)("L", "U", "N", "N", &nci, &nlev, one, omgf, &nci, bcp, &nci);
	ans += F77_CALL(ddot)(&ntot, bcp, &ione, bcp, &ione);
	Free(bcp);
    }
    return ans;
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
    SEXP rho;        /* environment in which to evaluate the calls */
    SEXP eta;        /* linear predictor */
    SEXP mu;         /* mean vector */
    SEXP linkinv;    /* expression for inverse link evaluation */
    SEXP mu_eta;     /* expression for dmu/deta evaluation */
    SEXP LMEopt;     /* expression for LME optimization */
    SEXP dev_resids; /* expression for deviance residuals */
    SEXP var;        /* expression for variance evaluation */
    double *offset;  /* offset for GLM */
    double *wts;     /* prior weights for GLM */
    double *y;       /* copy of response vector */
    double *etaold;  /* previous value of eta for evaluating  */
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
 * Calculate fitted values for the current fixed and random effects.
 * 
 * @param x pointer to an mer object
 * @param initial initial values (i.e. an offset) or (double *) NULL
 * @param val array to hold the fitted values
 * 
 * @return pointer to a numeric array of fitted values
 */
static double *
internal_mer_fitted(SEXP x, const double initial[], double val[])
{
    SEXP fixef = GET_SLOT(x, Matrix_fixefSym),
	ranef = GET_SLOT(x, Matrix_ranefSym);
    int ione = 1, n = LENGTH(GET_SLOT(x, Matrix_ySym)), p = LENGTH(fixef);
    double *X = REAL(GET_SLOT(x, Matrix_XSym)), one[] = {1,0};
    cholmod_sparse *Zt = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtSym));
    cholmod_dense *chv = numeric_as_chm_dense(val, n),
	*chb = as_cholmod_dense(ranef);
    
    mer_secondary(x);
    if (initial) Memcpy(val, initial, n);
    else AZERO(val, n);
    F77_CALL(dgemv)("N", &n, &p, one, X, &n, REAL(fixef), &ione, one, val, &ione);
    if (!cholmod_sdmult(Zt, 1, one, one, chb, chv, &c))
	error(_("Error return from sdmult"));
    Free(chv); Free(chb); Free(Zt);
    return val;
}
    
/** 
 * Extract the coefficients
 * 
 * @param x pointer to an mer object
 * @param ptyp parameter type to extract
 * @param ans vector to hold the extracted values
 * 
 * @return ans
 */
static double *
internal_mer_coef(SEXP x, int ptyp, double ans[])
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
				     REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(Omega, i)),
						   Matrix_xSym)), ncisq);
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
 * Set the coefficient vector and perform a factorization
 * 
 * @param x pointer to an mer object
 * @param cc vector of coefficients to assign
 * @param ptyp indicator of the parametrization being used
 */
static
void internal_mer_coefGets(SEXP x, const double cc[], int ptyp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	cind, i, nf = length(Omega);

    cind = 0;
    for (i = 0; i < nf; i++) {
	SEXP Omegai = VECTOR_ELT(Omega, i);
	int j, k, nci = nc[i], ncisq = nc[i] * nc[i];
	double *choli, *omgi = REAL(GET_SLOT(Omegai, Matrix_xSym));

	if (nci == 1) {
	    double dd = cc[cind++];
	    *omgi = ptyp ? ((ptyp == 1) ? exp(dd) : 1./dd) : dd;
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		ncip1 = nci + 1;

	    if (ptyp) {
		/* FIXME: Write this as an LDL decomposition */
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
	choli = REAL(GET_SLOT(dpoMatrix_chol(Omegai), Matrix_xSym));
	Memcpy(choli, omgi, ncisq);
	F77_CALL(dpotrf)("U", &nci, choli, &nci, &j);
	/* Yes, you really do need to do that decomposition.
	   The contents of choli before the decomposition are
	   from the previous value of Omegai. */
	if (j)
	    error(_("Omega[[%d]] is not positive definite"), i + 1);
    }
    flag_not_factored(x);
}

/** 
 * Evaluate current estimate of sigma from an mer object
 * 
 * @param x pointer to an mer object
 * @param REML indicator of whether to use REML.
 *           < 0  -> determine REML or ML from x@method
 *           == 0 -> use ML unconditionally
 *           > 0  -> use REML unconditionally
 * 
 * @return 
 */
static double
internal_mer_sigma(SEXP x, int REML)
{
    double *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym));

    if (REML < 0)		/* get REML from x */
	REML = !strcmp(CHAR(asChar(GET_SLOT(x,
					    Matrix_methodSym))),
		       "REML");
    mer_factor(x);
    return exp(dcmp[3]/2)/sqrt(dcmp[0] - (REML ? dcmp[1] : 0));
}

/** 
 * Update the derived quantities (ZtZ, ZtX, XtX, Zty, Xty 
 * and dcmp[2] = y'y when Z, X, y, wts or wrkres has been changed.
 * 
 * @param x pointer to an mer object
 * @param perm permutation from the Cholesky factor
 */

static void
internal_mer_update_ZXy(SEXP x, int *perm)
{
    SEXP Xp = GET_SLOT(x, Matrix_XSym), ZtZ = GET_SLOT(x, Matrix_ZtZSym),
	Ztyp = GET_SLOT(x, Matrix_ZtySym);
    SEXP ZtZx = GET_SLOT(ZtZ, Matrix_xSym),
	ZtZp = GET_SLOT(ZtZ, Matrix_pSym), ZtZi = GET_SLOT(ZtZ, Matrix_iSym);
    int *dims = INTEGER(getAttrib(Xp, R_DimSymbol)), i, ione = 1, j;
    int n = dims[0], nnz, p = dims[1], q = LENGTH(Ztyp);
    cholmod_sparse *ts1, *ts2,
	*Zt = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtSym));
    cholmod_sparse *Ztcp = cholmod_copy_sparse(Zt, &c);
    int *Zp = (int*)Ztcp->p;
    double *XtX = REAL(GET_SLOT(GET_SLOT(x, Matrix_XtXSym), Matrix_xSym)),
	*Xty = REAL(GET_SLOT(x, Matrix_XtySym)),
	*ZtX = REAL(GET_SLOT(GET_SLOT(x, Matrix_ZtXSym), Matrix_xSym)),
	*Zty = REAL(Ztyp),
	*wts = REAL(GET_SLOT(x, Matrix_wtsSym)),
	one[] = {1, 0}, zero[] = {0,0};
    cholmod_dense *td1, *Xd = as_cholmod_dense(Xp),
	*wkrd = as_cholmod_dense(GET_SLOT(x, Matrix_wrkresSym));
    cholmod_dense *Xcp = cholmod_copy_dense(Xd, &c),
	*wkrcp = cholmod_copy_dense(wkrd, &c);
    double *X = (double*)(Xcp->x), *Ztx = (double*)(Ztcp->x),
	*wtres = (double*)(wkrcp->x);
				/* Apply weights */
    for (i = 0; i < n; i++)
	wtres[i] = ((double*)(wkrd->x))[i] * wts[i];
    for (j = 0; j < p; j++)
	for (i = 0; i < n; i++)
	    X[i + j * n] = ((double*)(Xd->x))[i + j * n] * wts[i];
    for (j = 0; j < n; j++)
	for (i = Zp[j]; i < Zp[j + 1]; i++)
	    Ztx[i] = ((double*)(Zt->x))[i] * wts[j];
    Free(Zt); Free(Xd); Free(wkrd);
    
				/* y'y */
    REAL(GET_SLOT(x, Matrix_devCompSym))[2] =
	F77_CALL(ddot)(&n, wtres, &ione, wtres, &ione); 
				/* ZtZ */
    ts1 = cholmod_aat(Ztcp, (int *) NULL, (size_t) 0, 1/* mode */, &c);
    /* cholmod_aat returns stype == 0; copy to set stype == 1 */ 
    ts2 = cholmod_copy(ts1, 1/* stype */, 1/* mode */, &c);
    nnz = cholmod_nnz(ts2, &c);
    if (((int)(ts2->ncol) + 1) != LENGTH(ZtZp))
	error(_("Order of Z'Z has changed - was %d, now %d"),
	      LENGTH(ZtZp) - 1, (int)(ts2->ncol));
    Memcpy(INTEGER(ZtZp), (int*)(ts2->p), LENGTH(ZtZp));
    if (nnz != LENGTH(ZtZx))
	error(_("Number of nonzeros in Z'Z has changed - was %d, now %d"),
	      LENGTH(ZtZx), nnz);
    Memcpy(INTEGER(ZtZi), (int*)(ts2->i), nnz);
    Memcpy(REAL(ZtZx), (double*)(ts2->x), nnz);
    cholmod_free_sparse(&ts1, &c); cholmod_free_sparse(&ts2, &c);
				/* PZ'X into ZtX */
    td1 = cholmod_allocate_dense(q, p, q, CHOLMOD_REAL, &c);
    if (!cholmod_sdmult(Ztcp, 0, one, zero, Xcp, td1, &c))
	error(_("cholmod_sdmult failed"));
    for (j = 0; j < p; j++) { 	/* apply the permutation to each column */
	double *dcol = ZtX + j * q,
	    *scol = (double*)(td1->x) + j * q;
	for (i = 0; i < q; i++) dcol[i] = scol[perm[i]];
    }
    cholmod_free_dense(&td1, &c); 
				/* PZ'y into Zty */
    td1 = cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!cholmod_sdmult(Ztcp, 0, one, zero, wkrcp, td1, &c))
	error(_("cholmod_sdmult failed"));
    for (i = 0; i < q; i++) Zty[i] = ((double *)(td1->x))[perm[i]];
    cholmod_free_dense(&td1, &c); cholmod_free_sparse(&Ztcp, &c);
				/* XtX and Xty */
    AZERO(XtX, p * p);		
    F77_CALL(dsyrk)("U", "T", &p, &n, one, X, &n, zero, XtX, &p);
    F77_CALL(dgemv)("T", &n, &p, one, X, &n, wtres, &ione, zero, Xty, &ione);
    cholmod_free_dense(&Xcp, &c); cholmod_free_dense(&wkrcp, &c);
    flag_not_factored(x);
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

static double
Omega_log_det(SEXP Omega, int *nc, int *Gp)
{
    double ans = 0;
    int i;

    for (i = 0; i < LENGTH(Omega); i++) {
	int j, nci = nc[i], ncip1 = nc[i] + 1, nlev = (Gp[i + 1] - Gp[i])/nc[i];
	double *omgi = REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				     Matrix_xSym));

	for (j = 0; j < nci; j++) ans += 2. * nlev * log(fabs(omgi[j * ncip1]));
    }
    return ans;
}


/** 
 * Inflate Z'Z to Z'Z+Omega and factor. Form RZX and rZy and update
 * the status flags.
 * 
 * @param x pointer to an mer object.
 */
static void
internal_mer_Zfactor(SEXP x, cholmod_factor *L)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	rZyP = GET_SLOT(x, Matrix_rZySym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	nf = LENGTH(Omega), q = LENGTH(rZyP),
	p = LENGTH(GET_SLOT(x, Matrix_rXySym));
    cholmod_sparse *A, *Omg,
	*zz = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtZSym));
    cholmod_dense *ZtX = as_cholmod_dense(GET_SLOT(x, Matrix_ZtXSym)),
	*Zty = as_cholmod_dense(GET_SLOT(x, Matrix_ZtySym)),
	*rZy, *RZX;
    int *omp, *nnz = Calloc(nf + 1, int), i,
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    double *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)), one[] = {1, 0};
    

    dcmp[5] = Omega_log_det(Omega, nc, Gp); /* logDet(Omega) */
    for (nnz[0] = 0, i = 0; i < nf; i++)
	nnz[i + 1] = nnz[i] + ((Gp[i + 1] - Gp[i])*(nc[i] + 1))/2;
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
    Free(nnz);
    A = cholmod_add(zz, Omg, one, one, TRUE, TRUE, &c);
    Free(zz); cholmod_free_sparse(&Omg, &c);
    if (!cholmod_factorize(A, L, &c))
	error(_("rank_deficient Z'Z+Omega"));
    cholmod_free_sparse(&A, &c);
    dcmp[4] = 2 * chm_log_abs_det(L); /* 2 * logDet(L) */

				/* calculate and store RZX and rZy */
    RZX = cholmod_solve(CHOLMOD_L, L, ZtX, &c); Free(ZtX);
    rZy = cholmod_solve(CHOLMOD_L, L, Zty, &c); Free(Zty);
    Memcpy(REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXSym), Matrix_xSym)),
	   (double *) RZX->x, q * p);
    cholmod_free_dense(&RZX, &c);
    Memcpy(REAL(rZyP), (double *) rZy->x, q);
    cholmod_free_dense(&rZy, &c);
				/* signal that secondary slots are not valid */
    status[0] = 1;
    status[1] = status[2] = status[3] = 0;
}

/** 
 * Update the relative precision matrices by sampling from a Wishart
 * distribution with scale factor determined by the current sample of
 * random effects.
 * 
 * @param Omega pointer to the list of relative precision matrices
 * @param b current sample from the random effects
 * @param sigma current value of sigma
 * @param nf number of grouping factors
 * @param nc number columns per grouping factor
 * @param Gp integer vector pointing to the beginning of each outer
 * group of columns
 * @param vals vector in which to store values
 * @param trans logical value indicating if variance components are
 * stored in the transformed scale.
 */
static void
internal_Omega_update(SEXP Omega, const double b[], double sigma, int nf,
		  const int nc[], const int Gp[], double *vals, int trans)
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
	    *vals++ = (trans ? log(var[j * ncip1]) : var[j * ncip1]);
	}
	for (j = 1; j < nci; j++) {
	    for (k = 0; k < j; k++) {
		*vals++ = (trans ? atanh(var[k + j * nci]/
				     sqrt(var[j * ncip1] * var[k * ncip1]))
		       : var[k + j * nci]);
	    }
	}
	/* calculate and store the relative precision matrix */
	Memcpy(tmp, wfac, ncisqr);
	F77_CALL(dtrsm)("R", "U", "T", "N", &nci, &nci,
			&sigma, scal, &nci, tmp, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, tmp, &nci, &zero,
			REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym)),
			&nci);
	dpoMatrix_chol(VECTOR_ELT(Omega, i));
	Free(scal); Free(tmp); Free(wfac); Free(var); 
    }
}

/** 
 * Evaluate new weights and working residuals.
 * 
 * @param GS a GlmerStruct object
 */
static void
internal_glmer_reweight(GlmerStruct GS) {
    SEXP dmu_deta, var;
    int i;
    double *w = REAL(GET_SLOT(GS->mer, Matrix_wtsSym)),
	*y = REAL(GET_SLOT(GS->mer, Matrix_ySym)),
	*z = REAL(GET_SLOT(GS->mer, Matrix_wrkresSym));

				/* reweight mer */
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				  REALSXP, GS->n));
    var = PROTECT(eval_check(GS->var, GS->rho,
			     REALSXP, GS->n));
    for (i = 0; i < GS->n; i++) {
	w[i] = GS->wts[i] *
	    REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	z[i] = REAL(GS->eta)[i] - GS->offset[i] +
	    (y[i] - REAL(GS->mu)[i])/REAL(dmu_deta)[i];
    }
    UNPROTECT(2);
    mer_update_ZXy(GS->mer);
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
 * Update the ranef slot assuming that the fixef, rZy, RZX and L slots
 * are up to date.
 * 
 * @param x Pointer to an mer object
 * 
 * @return 
 */
static double*
internal_mer_ranef(SEXP x)
{
    SEXP ranef = GET_SLOT(x, Matrix_ranefSym);
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    if (!status[0]) {
	error("Applying internal_mer_ranef without factoring");
	return (double*)NULL;	/* -Wall */
    }
    if (!status[1]) {
	SEXP fixef = GET_SLOT(x, Matrix_fixefSym),
	    ranef = GET_SLOT(x, Matrix_ranefSym);
	int ione = 1, p = LENGTH(fixef), q = LENGTH(ranef);
	cholmod_factor *L = as_cholmod_factor(GET_SLOT(x, Matrix_LSym));
	cholmod_dense *td1, *td2,
	    *chranef = as_cholmod_dense(ranef);
	double *RZX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXSym), Matrix_xSym)),
	    m1[] = {-1,0}, one[] = {1,0};
	
	Memcpy(REAL(ranef), REAL(GET_SLOT(x, Matrix_rZySym)), q);
	F77_CALL(dgemv)("N", &q, &p, m1, RZX, &q, REAL(fixef), &ione,
			one, REAL(ranef), &ione);
	td1 = cholmod_solve(CHOLMOD_Lt, L, chranef, &c);
	td2 = cholmod_solve(CHOLMOD_Pt, L, td1, &c);
	Free(chranef); cholmod_free_dense(&td1, &c);
	Memcpy(REAL(ranef), (double *)(td2->x), q);
	cholmod_free_dense(&td2, &c);
	status[1] = 1;
	status[2] = status[3] = 0;
	Free(L);
    }
    return REAL(ranef);
}

/** 
 * Update the fixef slot on a factored mer object.
 * 
 * @param x Pointer to an mer object
 * 
 * @return fixed effects vector
 */
static double*
internal_mer_fixef(SEXP x)
{
    SEXP fixef = GET_SLOT(x, Matrix_fixefSym);
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    if (!status[0]) {
	error("Applying internal_mer_fixef without factoring");
	return (double*)NULL;	/* -Wall */
    }
    if (!status[1]) {
	int ione = 1, p = LENGTH(fixef);
	Memcpy(REAL(fixef), REAL(GET_SLOT(x, Matrix_rXySym)), p);
	F77_CALL(dtrsv)("U", "N", "N", &p,
			REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym),
				      Matrix_xSym)),
			&p, REAL(fixef), &ione);
    }
    return REAL(fixef);
}

/** 
 * Iterate to determine the conditional modes of the random effects.
 * 
 * @param GS a GlmerStruct object
 * @param fixed vector of fixed effects
 * @param varc vector of parameters for the variance-covariance
 * 
 * @return An indicator of whether the iterations converged
 */
static int
internal_bhat(GlmerStruct GS, const double fixed[], const double varc[])
{
    SEXP fixef = GET_SLOT(GS->mer, Matrix_fixefSym);
    cholmod_factor *L = as_cholmod_factor(GET_SLOT(GS->mer, Matrix_LSym));
    int i;
    double crit = GS->tol + 1;
    
    if (varc)	  /* skip this step if varc == (double*) NULL */	
	internal_mer_coefGets(GS->mer, varc, 2);
    Memcpy(REAL(fixef), fixed, LENGTH(fixef));
    internal_mer_Zfactor(GS->mer, L);
    internal_mer_ranef(GS->mer);
    internal_mer_fitted(GS->mer, GS->offset, REAL(GS->eta));
    Memcpy(GS->etaold, REAL(GS->eta), GS->n);
    
    for (i = 0; i < GS->maxiter && crit > GS->tol; i++) {
	internal_glmer_reweight(GS);
	internal_mer_Zfactor(GS->mer, L);
	internal_mer_ranef(GS->mer);
	internal_mer_fitted(GS->mer, GS->offset, REAL(GS->eta));
	crit = conv_crit(GS->etaold, REAL(GS->eta), GS->n);
    }
    Free(L);
    return (crit > GS->tol) ? 0 : i;
}

/**
 * Print the verbose output in the ECME iterations
 *
 * @param x pointer to an ssclme object
 * @param iter iteration number
 * @param REML non-zero for REML, zero for ML
 */
static void
EMsteps_verbose_print(SEXP x, int iter, int REML)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	gradComp = GET_SLOT(x, Matrix_gradCompSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ifour = 4, ii, ione = 1, jj,
	nf = LENGTH(Omega);
    double
	cc[] = {-1, 1, 1, REML ? 1 : 0},
	*dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	one = 1., zero = 0.;

    if (iter == 0) Rprintf("  EM iterations\n");
    Rprintf("%3d %.3f", iter, dev[REML ? 1 : 0]);
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1, ncisqr = nci * nci;
	double
	    *Omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym)),
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
			REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
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
}

/** 
 * Perform a number of ECME steps
 * 
 * @param x pointer to an mer object
 * @param nEM number of iteration to perform
 * @param verb indicator of verbose output
 */
static void
internal_ECMEsteps(SEXP x, int nEM, int verb)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	gradComp = GET_SLOT(x, Matrix_gradCompSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),
		       "REML"),
	i, ifour = 4, info, ione = 1, iter,
	nf = LENGTH(Omega);
    double
	cc[] = {0, 1, 1, REML ? 1 : 0},
	zero = 0.0;

    mer_gradComp(x);
    if (verb)
	EMsteps_verbose_print(x, 0, REML);
    for (iter = 0; iter < nEM; iter++) {
	for (i = 0; i < nf; i++) {
	    SEXP Omgi = VECTOR_ELT(Omega, i);
	    int nci = nc[i], ncisqr = nci * nci,
		nlev = (Gp[i + 1] - Gp[i])/nci;
	    double *Omgx = REAL(GET_SLOT(Omgi, Matrix_xSym)),
		mult = 1./((double) nlev);

	    F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &mult,
			    REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
			    cc, &ifour, &zero, Omgx, &ncisqr);
	    F77_CALL(dpotrf)("U", &nci, Omgx, &nci, &info);
	    if (info)
		error(_("DPOTRF in ECME update gave code %d"),
		      info);
	    F77_CALL(dpotri)("U", &nci, Omgx, &nci, &info);
	    if (info)
		error(_("Matrix inverse in ECME update gave code %d"), info);
	    SET_SLOT(Omgi, Matrix_factorSym, allocVector(VECSXP, 0));
	}
	flag_not_factored(x);
	mer_gradComp(x);
	if (verb)
	    EMsteps_verbose_print(x, iter + 1, REML);
    }
}

/** 
 * Evaluate the conditional deviance for the stored random effects.
 * 
 * @param GS Pointer to a GlmerStruct
 * @param fixed value of the fixed effects
 * 
 * @return conditional deviance
 */
static double
random_effects_deviance(GlmerStruct GS)
{
    SEXP devs;
    int i;
    double ans;

    internal_mer_fitted(GS->mer, GS->offset, REAL(GS->eta));
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0, ans = 0; i < GS->n; i++) ans += REAL(devs)[i];
    UNPROTECT(1);
    return ans;
}

/* Externally accessible functions */

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
 * Perform the PQL optimization
 * 
 * @param GSp pointer to a GlmerStruct object
 * 
 * @return R_NilValue
 */
SEXP glmer_PQL(SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    int i; double crit;

    Memcpy(GS->etaold, REAL(GS->eta), GS->n);
    for (i = 0, crit = GS->tol + 1;
	 i < GS->maxiter && crit > GS->tol; i++) {
	internal_glmer_reweight(GS);
	if (!i) mer_initial(GS->mer); /* initialize first fit */
	internal_ECMEsteps(GS->mer, i ? 2 : GS->niterEM,
			   GS->EMverbose);
	eval(GS->LMEopt, GS->rho);
	internal_mer_fitted(GS->mer, GS->offset, REAL(GS->eta));
	crit = conv_crit(GS->etaold, REAL(GS->eta), GS->n);
    }
    if (crit > GS->tol)
	warning(_("IRLS iterations for PQL did not converge"));

    return R_NilValue;
}

/** 
 * Compute the Laplace approximation to the deviance.
 * 
 * @param pars pointer to a numeric vector of parameters
 * @param GSp pointer to a GlmerStruct object
 * 
 * @return the Laplace approximation to the deviance
 */
SEXP glmer_devLaplace(SEXP pars, SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    SEXP Omega = GET_SLOT(GS->mer, Matrix_OmegaSym);
    int *Gp = INTEGER(GET_SLOT(GS->mer, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(GS->mer, Matrix_ncSym));
    double *bhat = REAL(GET_SLOT(GS->mer, Matrix_ranefSym)),
	*dcmp = REAL(GET_SLOT(GS->mer, Matrix_devCompSym)),
	*dev = REAL(GET_SLOT(GS->mer, Matrix_devianceSym));
    
    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	return ScalarReal(DBL_MAX);
    dev[0] = dcmp[4] - dcmp[5] + random_effects_deviance(GS) +
	b_quadratic(bhat, Omega, Gp, nc);
    dev[1] = NA_REAL;
    return ScalarReal(dev[0]);
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
    
    Free(GS->offset); Free(GS->wts); Free(GS->etaold);
    Free(GS);
    return R_NilValue;
}

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
    SEXP tmp, y, Ztx;
    
    
    GS = (GlmerStruct) Calloc(1, glmer_struct);
    if (!isEnvironment(rho))
	error(_("`rho' must be an environment"));
    GS->rho = rho;
    GS->mer = find_and_check(rho, install("mer"), VECSXP, 0);
    y = GET_SLOT(GS->mer, Matrix_ySym);
    GS->n = LENGTH(y);
    GS->p = LENGTH(GET_SLOT(GS->mer, Matrix_rXySym));
    GS->y = Memcpy(Calloc(GS->n, double), REAL(y), GS->n);
    Ztx = GET_SLOT(GET_SLOT(GS->mer, Matrix_ZtSym), Matrix_xSym);
    GS->mu = find_and_check(rho, install("mu"), REALSXP, GS->n);
    tmp = find_and_check(rho, install("offset"), REALSXP, GS->n);
    GS->offset = Memcpy(Calloc(GS->n, double), REAL(tmp), GS->n);
    tmp = find_and_check(rho, install("wts"), REALSXP, GS->n);
    GS->wts = Memcpy(Calloc(GS->n, double), REAL(tmp), GS->n);
    GS->etaold = Calloc(GS->n, double);
    GS->cv = find_and_check(rho, install("cv"), VECSXP, 0);
    GS->niterEM = asInteger(Matrix_getElement(GS->cv, "niterEM"));
    GS->EMverbose = asLogical(Matrix_getElement(GS->cv, "EMverbose"));
    GS->tol = asReal(Matrix_getElement(GS->cv, "tolerance"));
    GS->maxiter = asInteger(Matrix_getElement(GS->cv, "maxIter"));
    GS->nf = LENGTH(GET_SLOT(GS->mer, Matrix_flistSym));
    GS->npar = GS->p +
	coef_length(GS->nf, INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)));
    GS->eta = find_and_check(rho, install("eta"), REALSXP, GS->n);

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
 * Perform ECME steps for the REML or ML criterion.
 *
 * @param x pointer to an mer object
 * @param nsteps pointer to an integer scalar - the number of ECME
 * steps to perform
 * @param Verbp pointer to a logical scalar indicating verbose output
 *
 * @return R_NilValue
 */
SEXP mer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp)
{
    int nstp = asInteger(nsteps);
    if (nstp > 0) internal_ECMEsteps(x, nstp, asLogical(Verbp));
    return R_NilValue;
}

/**
 * Fill in five symmetric matrices, providing the information to
 * generate the Hessian.
 *
 * @param x pointer to an mer object
 *
 * @return an array consisting of five symmetric faces
 */
SEXP mer_Hessian(SEXP x)
{
    SEXP
	D = GET_SLOT(x, Matrix_DSym),
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	gradComp = GET_SLOT(x, Matrix_gradCompSym),
	val;
    int *dRZX = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	Q, Qsqr, RZXpos, facepos,
	i, ione = 1, j, nf = length(Omega), p = dRZX[1] - 1, pos;
    double
	*RZX = REAL(RZXP),
	*b = REAL(RZXP) + dRZX[0] * p,
	*bbface,		/* vec of second faces of gradComp elts */
	one = 1.,
	zero = 0.;

    mer_gradComp(x);
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
	double *fDi = REAL(VECTOR_ELT(gradComp, i)),
	    mult = 1./((double)(Gp[i + 1] - Gp[i])/nci);

	Memcpy(bbface + pos, fDi + ncisqr, ncisqr);
	/* outer product of the third face of gradComp on the diagonal
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
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
 * 
 * @param mer pointer to an mer object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 * 
 * @return a matrix
 */
SEXP
mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp)
{
    SEXP ans, Omega = GET_SLOT(x, Matrix_OmegaSym),
	Omegacp = PROTECT(duplicate(Omega));
    cholmod_factor *L = as_cholmod_factor(GET_SLOT(x, Matrix_LSym));
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))), "REML"),
	i, ione = 1, j, n = LENGTH(GET_SLOT(x, Matrix_ySym)),
	nf = LENGTH(Omega), nsamp = asInteger(nsampp),
	p = LENGTH(GET_SLOT(x, Matrix_rXySym)),
	q = LENGTH(GET_SLOT(x, Matrix_rZySym)),
	saveb = asLogical(savebp),
	trans = asLogical(transp);
    double 
	*RXX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym), Matrix_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXSym), Matrix_xSym)),
	*bhat = REAL(GET_SLOT(x, Matrix_ranefSym)),
	*betahat = REAL(GET_SLOT(x, Matrix_fixefSym)), 
	*bnew = Calloc(q, double), *betanew = Calloc(p, double),
	*dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	df = n - (REML ? p : 0), m1[] = {-1,0}, one[] = {1,0};
    int nrbase = p + 1 + coef_length(nf, nc); /* rows always included */
    int nrtot = nrbase + (saveb ? q : 0);
    cholmod_dense *chb, *chbnew = numeric_as_chm_dense(bnew, q);
    
    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    GetRNGstate();
    for (i = 0; i < nsamp; i++) {
	double *col = REAL(ans) + i * nrtot, sigma;
				/* factor x and get secondary values */
	mer_factor(x);
	mer_secondary(x);
				/* simulate and store new value of sigma */
	sigma = exp(dcmp[3]/2)/sqrt(rchisq(df));
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
				/* simulate scaled, independent normals */
	for (j = 0; j < p; j++) betanew[j] = sigma * norm_rand();
	for (j = 0; j < q; j++) bnew[j] = sigma * norm_rand();
				/* betanew := RXX^{-1} %*% betanew */
	F77_CALL(dtrsv)("U", "N", "N", &p, RXX, &p, betanew, &ione);
				/* bnew := bnew - RZX %*% betanew */
	F77_CALL(dgemv)("N", &q, &p, m1, RZX, &q, betanew, &ione,
			one, bnew, &ione);
				/* chb := L^{-T} %*% bnew */
	chb = cholmod_solve(CHOLMOD_Lt, L, chbnew, &c);
				/* Copy chb to bnew and free chb */
	Memcpy(bnew, (double *)(chb->x), q);
 	cholmod_free_dense(&chb, &c);
				/* Add conditional modes and store beta */
	for (j = 0; j < p; j++) {
	    col[j] = (betanew[j] += betahat[j]);
	}
				/* Add conditional modes and
				 * optionally store b */
	for (j = 0; j < q; j++) {
	    bnew[j] += bhat[j];
	    if (saveb) col[nrbase + j] = bnew[j];
	}
	/* Update and store variance-covariance of the random effects */
	internal_Omega_update(Omega, bnew, sigma, nf, nc, Gp, col + p + 1, trans);
    }
    PutRNGstate();
    Free(betanew); Free(bnew); Free(chbnew);
				/* Restore original Omega */
    SET_SLOT(x, Matrix_OmegaSym, Omegacp);
    mer_factor(x);

    Free(L);
    UNPROTECT(2);
    return ans;
}

/**
 * Create an mer object from a list of grouping factors and a list of model
 * matrices. 
 *
 * @param fl named list of grouping factors
 * @param Ztl list of transposes of model matrices
 * @param Xp model matrix for the fixed effects
 * @param yp response vector
 * @param method character vector describing the estimation method
 *
 * @return pointer to an mer object
 */
SEXP mer_create(SEXP fl, SEXP ZZt, SEXP Xp, SEXP yp, SEXP method,
		 SEXP ncp, SEXP cnames, SEXP useS, SEXP call,
		 SEXP family)
{
    SEXP Omega, bVar, gradComp, fnms = getAttrib(fl, R_NamesSymbol),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("mer"))), xnms;
    cholmod_sparse *ts1, *ts2, *Zt;
    cholmod_factor *F;
    int *nc = INTEGER(ncp), *Gp, *xdims, i,
	nf = LENGTH(fl), nobs = LENGTH(yp), p, q;
    char *devCmpnms[] = {"n", "p", "yty", "logryy2", "logDetL2",
			 "logDetOmega", "logDetRXX", ""};
    char *devnms[] = {"ML", "REML", ""};
    char *statusnms[] =
	{"factored", "secondary", "gradComp", "Hesscomp", ""};
    double *dcmp, *wts, *wrkres, *y;
				/* Check arguments to be duplicated */
    if (!isReal(yp)) error(_("yp must be a real vector"));
    SET_SLOT(val, Matrix_ySym, duplicate(yp));
    if (!isMatrix(Xp) || !isReal(Xp))
	error(_("Xp must be a real matrix"));
    xdims = INTEGER(getAttrib(Xp, R_DimSymbol));
    if (xdims[0] != nobs) error(_("Xp must have %d rows"), nobs);
    p = xdims[1];
    xnms = VECTOR_ELT(getAttrib(Xp, R_DimNamesSymbol), 1);
    SET_SLOT(val, Matrix_XSym, duplicate(Xp));
    if (!isNewList(fl) || nf < 1) error(_("fl must be a nonempty list"));
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);
    }
    SET_SLOT(val, Matrix_flistSym, duplicate(fl));
    if (!isString(method) || LENGTH(method) != 1)
	error(_("method must be a character vector of length 1"));
    SET_SLOT(val, Matrix_methodSym, duplicate(method));
    if (!isLogical(useS) || LENGTH(useS) != 1)
	error(_("useS must be a logical vector of length 1"));
    SET_SLOT(val, Matrix_useScaleSym, duplicate(useS));
    if (!isNewList(cnames) || LENGTH(cnames) != nf + 1)
	error(_("cnames must be a list of length %d"), nf + 1);
    SET_SLOT(val, Matrix_cnamesSym, duplicate(cnames));
    if (!isInteger(ncp) || LENGTH(ncp) != nf)
	error(_("ncp must be an integer vector of length %d"), nf);
    SET_SLOT(val, Matrix_callSym, duplicate(call));
    SET_SLOT(val, Matrix_familySym, duplicate(family));
    SET_SLOT(val, Matrix_ncSym, duplicate(ncp));
    Gp = INTEGER(ALLOC_SLOT(val, Matrix_GpSym, INTSXP, nf + 1));
    Gp[0] = 0;
    if (!isNewList(fl) || nf < 1) error(_("fl must be a nonempty list"));
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);

    }
    SET_SLOT(val, Matrix_ZtSym, duplicate(ZZt));
    Zt = as_cholmod_sparse(GET_SLOT(val, Matrix_ZtSym));
				/* analyze Zt */
    q = Zt->nrow;
    i = c.supernodal;
    c.supernodal = CHOLMOD_SUPERNODAL; /* force a supernodal decomposition */
    F = cholmod_analyze(Zt, &c);
    c.supernodal = i;
    ts1 = cholmod_aat(Zt, (int*)NULL/* fset */,(size_t)0,
		      CHOLMOD_PATTERN, &c);
    ts2 = cholmod_copy(ts1, 1/*upper triangle*/, CHOLMOD_PATTERN, &c);
    SET_SLOT(val, Matrix_ZtZSym,
	     alloc_dsCMatrix(q, cholmod_nnz(ts2, &c), "U", R_NilValue,
			     R_NilValue));
    cholmod_free_sparse(&ts1, &c); cholmod_free_sparse(&ts2, &c);
				/* allocate other slots */
    SET_SLOT(val, Matrix_devianceSym, Matrix_make_named(REALSXP, devnms));
    SET_SLOT(val, Matrix_devCompSym, Matrix_make_named(REALSXP, devCmpnms));
    SET_SLOT(val, Matrix_statusSym, Matrix_make_named(LGLSXP, statusnms));
    AZERO(LOGICAL(GET_SLOT(val, Matrix_statusSym)), 4);
    dcmp = REAL(GET_SLOT(val, Matrix_devCompSym));
    AZERO(dcmp, 7);		/* cosmetic */
    dcmp[0] = (double) nobs;
    dcmp[1] = (double) p;
				/* allocate and populate list slots */
    Omega = ALLOC_SLOT(val, Matrix_OmegaSym, VECSXP, nf);
    bVar = ALLOC_SLOT(val, Matrix_bVarSym, VECSXP, nf);
    gradComp = ALLOC_SLOT(val, Matrix_gradCompSym, VECSXP, nf);
    setAttrib(Omega, R_NamesSymbol, duplicate(fnms));
    setAttrib(bVar, R_NamesSymbol, duplicate(fnms));
    setAttrib(gradComp, R_NamesSymbol, duplicate(fnms));
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int nlev = LENGTH(getAttrib(VECTOR_ELT(fl, i), R_LevelsSymbol));
	SET_VECTOR_ELT(Omega, i,
		       alloc_dpoMatrix(nci, "U",
				       VECTOR_ELT(cnames, i),
				       VECTOR_ELT(cnames, i)));
	SET_VECTOR_ELT(bVar, i, alloc3Darray(REALSXP, nci, nci, nlev));
	SET_VECTOR_ELT(gradComp, i, alloc3Darray(REALSXP, nci, nci, 4));
	Gp[i + 1] = Gp[i] + nlev * nci;
    }
				/* create ZtX, RZX, XtX, RXX */
    SET_SLOT(val, Matrix_ZtXSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
    SET_SLOT(val, Matrix_RZXSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
    SET_SLOT(val, Matrix_XtXSym, alloc_dpoMatrix(p, "U", xnms, xnms));
    SET_SLOT(val, Matrix_RXXSym, alloc_dtrMatrix(p, "U", "N", xnms, xnms));
    SET_SLOT(val, Matrix_ZtySym, allocVector(REALSXP, q));
    SET_SLOT(val, Matrix_rZySym, allocVector(REALSXP, q));
    SET_SLOT(val, Matrix_XtySym, allocVector(REALSXP, p));
    SET_SLOT(val, Matrix_rXySym, allocVector(REALSXP, p));
				/* create weights and working residuals */
    wts = REAL(ALLOC_SLOT(val, Matrix_wtsSym, REALSXP, nobs));
    wrkres = REAL(ALLOC_SLOT(val, Matrix_wrkresSym, REALSXP, nobs));
    y = REAL(GET_SLOT(val, Matrix_ySym));
    for (i = 0; i < nobs; i++) {wts[i] = 1.; wrkres[i] = y[i];}
    internal_mer_update_ZXy(val, (int*)(F->Perm));
    Free(Zt);
				/* secondary slots */
    SET_SLOT(val, Matrix_ranefSym, allocVector(REALSXP, q));
    SET_SLOT(val, Matrix_fixefSym, allocVector(REALSXP, p));
    SET_SLOT(val, Matrix_RZXinvSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
				/* initialize */
    mer_initial(val);
    /* The next calls are simply to set up the L slot.  At present the
     * factor F is a symbolic factor.  We need to force it to become
     * numeric before allocating the L slot in the object. */
    internal_mer_Zfactor(val, F); 
    /* One side-effect of this call is to set the status as
     * factored.  We need to reset it */
    LOGICAL(GET_SLOT(val, Matrix_statusSym))[0] = 0;
    /* Create the dCHMfactor object and store it in the L slot.  This
     * also frees the storage. */
    SET_SLOT(val, Matrix_LSym, chm_factor_to_SEXP(F, 1));
    /* OK, done now. */
    UNPROTECT(1);
    return val;
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
 * @param x pointer to an mer object
 * @param pType pointer to an integer scalar indicating the form of the 
 *        parameters to be returned.
 *
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP mer_coef(SEXP x, SEXP pType)
{
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	nf = LENGTH(GET_SLOT(x, Matrix_OmegaSym));
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));

    internal_mer_coef(x, asInteger(pType), REAL(val));
    UNPROTECT(1);
    return val;
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
SEXP mer_coefGets(SEXP x, SEXP coef, SEXP pType)
{
    int clen = coef_length(LENGTH(GET_SLOT(x, Matrix_flistSym)),
			   INTEGER(GET_SLOT(x, Matrix_ncSym)));   
    if (LENGTH(coef) != clen || !isReal(coef))
	error(_("coef must be a numeric vector of length %d"), clen);
    internal_mer_coefGets(x, REAL(coef), asInteger(pType));
    return x;
}

/** 
 * Return L as a dtCMatrix object
 * 
 * @param x pointer to an mer object
 * 
 * @return L as an dtCMatrix object
 */
SEXP mer_dtCMatrix(SEXP x)
{
    cholmod_factor *L = as_cholmod_factor(GET_SLOT(x, Matrix_LSym)), *Lcp;
    cholmod_sparse *Lm;
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	nz, q;

    dims[0] = dims[1] = q = (int)(L->n);
    Lcp = cholmod_copy_factor(L, &c); Free(L); /* next call changes Lcp */
    Lm = cholmod_factor_to_sparse(Lcp, &c); cholmod_free_factor(&Lcp, &c);
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, q + 1)),
	   (int*) Lm->p, q + 1);
    nz = ((int*)(Lm->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)),
	   (int*) Lm->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)),
	   (double*) Lm->x, nz);
    cholmod_free_sparse(&Lm, &c);
    UNPROTECT(1);
    return ans;
}

/** 
 * Return L^{-1} as a dtCMatrix object
 * 
 * @param x pointer to an mer object
 * 
 * @return L^{-1} as an dtCMatrix object
 */
SEXP mer_dtCMatrix_inv(SEXP x)
{
    cholmod_factor *L = as_cholmod_factor(GET_SLOT(x, Matrix_LSym));
    cholmod_sparse
	*b = cholmod_allocate_sparse(L->n, L->n, L->n, 1, 1,
				     0, CHOLMOD_REAL, &c),
	*Linv;
    double *bx = (double *)(b->x);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *bi = (int *) (b->i), *bp = (int *) (b->p),
	*dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	j, nz, q;

    dims[0] = dims[1] = q = (int)(L->n);
    for (j = 0; j < q; j++) {
	bp[j] = bi[j] = j;
	bx[j] = 1;
    }
    bp[q] = q;
    Linv = cholmod_spsolve(CHOLMOD_L, L, b, &c);
    cholmod_free_sparse(&b, &c);
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, q + 1)),
	   (int *) Linv->p, q + 1);
    nz = ((int *)(Linv->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)),
	   (int *) Linv->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)),
	   (double *) Linv->x, nz);
    cholmod_free_sparse(&Linv, &c);
    UNPROTECT(1);
    return ans;
}

/**
 * Create and factor Z'Z+Omega.  Also create RZX and RXX, the deviance
 * components, and the value of the deviance for both ML and REML.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP mer_factor(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    if (!status[0]) {
	SEXP rXyP = GET_SLOT(x, Matrix_rXySym),
	    rZyP = GET_SLOT(x, Matrix_rZySym);
	int i, info, ione = 1, p = LENGTH(rXyP), q = LENGTH(rZyP);
	cholmod_factor *L = as_cholmod_factor(GET_SLOT(x, Matrix_LSym));
	double *RXX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym), Matrix_xSym)),
	    *RZX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXSym), Matrix_xSym)),
	    *rXy = REAL(rXyP), *rZy = REAL(rZyP),
	    *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	    *dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	    one[2] = {1, 0}, m1[2] = {-1, 0};
	double nml = dcmp[0], nreml = dcmp[0] - dcmp[1];
	
	/* Inflate Z'Z to Z'Z+Omega and factor to form L. Form RZX and
	 * rZy. Update status flags, dcmp[4] and dcmp[5]. */
	internal_mer_Zfactor(x, L); 
				/* downdate XtX and factor */
	Memcpy(RXX, REAL(GET_SLOT(GET_SLOT(x, Matrix_XtXSym), Matrix_xSym)), p * p);
	F77_CALL(dsyrk)("U", "T", &p, &q, m1, RZX, &q, one, RXX, &p);
	F77_CALL(dpotrf)("U", &p, RXX, &p, &info);
	if (info) {
	    error(_("Leading minor of order %d in downdated X'X is not positive definite"),
		  info);
	    dcmp[3] = dcmp[6] = dev[0] = dev[1] = NA_REAL;
	} else {
	    for (dcmp[6] = 0, i = 0; i < p; i++) /* 2 * logDet(RXX) */
		dcmp[6] += 2. * log(RXX[i * (p + 1)]);
				/* solve for rXy  and ryy^2 */
	    Memcpy(rXy, REAL(GET_SLOT(x, Matrix_XtySym)), p);
	    F77_CALL(dgemv)("T", &q, &p, m1, RZX, &q, rZy, &ione, one, rXy, &ione);
	    F77_CALL(dtrsv)("U", "T", "N", &p, RXX, &p, rXy, &ione);
	    dcmp[3] = log(dcmp[2] /* dcmp[3] = log(ryy^2); dcmp[2] = y'y; */
			  - F77_CALL(ddot)(&p, rXy, &ione, rXy, &ione)
			  - F77_CALL(ddot)(&q, rZy, &ione, rZy, &ione));
				/* evaluate ML and REML deviance */
	    dev[0] = dcmp[4] - dcmp[5] +
		nml*(1.+dcmp[3]+log(2.*PI/nml));
	    dev[1] = dcmp[4] - dcmp[5] + dcmp[6] +
		nreml*(1.+dcmp[3]+log(2.*PI/nreml));
	}
	Free(L);
    }
    return R_NilValue;
}

/** 
 * Return the fitted values as an SEXP
 * 
 * @param x pointer to an mer object
 * @param useFe pointer to a logical scalar indicating if the fixed
 * effects should be used
 * @param useRe pointer to a logical scalar indicating if the random
 * effects should be used
 * 
 * @return pointer to a numeric array of fitted values
 */

SEXP mer_fitted(SEXP x)
{
    int n = LENGTH(GET_SLOT(x, Matrix_ySym));
    SEXP ans = PROTECT(allocVector(REALSXP, n));

    internal_mer_fitted(x, (double*) NULL, REAL(ans));
    UNPROTECT(1);
    return ans;
}

/**
 * Extract the conditional estimates of the fixed effects
 *
 * @param x Pointer to an mer object
 *
 * @return a numeric vector containing the conditional estimates of
 * the fixed effects
 */
SEXP mer_fixef(SEXP x)
{
    int nf = LENGTH(GET_SLOT(x, Matrix_OmegaSym));
    SEXP ans;
    
    mer_secondary(x);
    ans = PROTECT(duplicate(GET_SLOT(x, Matrix_fixefSym)));
    setAttrib(ans, R_NamesSymbol,
	      duplicate(VECTOR_ELT(GET_SLOT(x, Matrix_cnamesSym), nf)));
    UNPROTECT(1);
    return ans;
}

/**
 * Fill in the gradComp and bVar slots.  Each component in the gradComp slot
 * consists of four symmetric matrices used to generate the gradient or the ECME
 * step.  They are
 *  1) -m_i\bOmega_i^{-1}
 *  2) \bB_i\bB_i\trans
 *  3) \tr\left[\der_{\bOmega_i}\bOmega\left(\bZ\trans\bZ+\bOmega\right)\inv\right]
 *  4) The term added to 3) to get \tr\left[\der_{\bOmega_i}\bOmega\vb\right]
 *
 * @param x pointer to an mer object
 * @param val pointer to a list of matrices of the correct sizes
 *
 * @return NULL
 */
SEXP mer_gradComp(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));

    mer_secondary(x);
    if (!status[2]) {
	SEXP bVarP = GET_SLOT(x, Matrix_bVarSym),
	    OmegaP = GET_SLOT(x, Matrix_OmegaSym),
	    RZXP = GET_SLOT(x, Matrix_RZXSym),
	    gradComp = GET_SLOT(x, Matrix_gradCompSym),
	    ranefP = GET_SLOT(x, Matrix_ranefSym);
	int q = LENGTH(ranefP), p = LENGTH(GET_SLOT(x, Matrix_rXySym));
	cholmod_factor *L = as_cholmod_factor(GET_SLOT(x, Matrix_LSym));
	cholmod_dense *RZX = as_cholmod_dense(RZXP), *tmp1;
	int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	    *Perm = (int *)(L->Perm),
	    *iperm = Calloc(q, int),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    i, j, k, nf = length(OmegaP);
	double *b = REAL(GET_SLOT(x, Matrix_ranefSym)),
	    *RZXinv = REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXinvSym), Matrix_xSym)),
	    m1[] = {-1, 0}, alpha;
    
	alpha = 1./internal_mer_sigma(x, -1);
	alpha = alpha * alpha;
				
	for (j = 0; j < q; j++) iperm[Perm[j]] = j; /* create the inverse permutation */
				/* form RZXinv, the section of R^{-1} from RZX */
	tmp1 = cholmod_solve(CHOLMOD_Lt, L, RZX, &c); Free(RZX);
	/* copy columns of tmp1 to RZXinv applying the inverse permutation */
	for (j = 0; j < p; j++) {
	    double *dest = RZXinv + j * q, *src = ((double*)(tmp1->x)) + j * q;
	    for (i = 0; i < q; i++) dest[i] = src[iperm[i]];
	}
	cholmod_free_dense(&tmp1, &c);
	F77_CALL(dtrsm)("R", "U", "N", "N", &q, &p, m1,
			REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym), Matrix_xSym)),
			&p, RZXinv, &q);

	internal_mer_bVar(x);
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], RZXrows = Gp[i + 1] - Gp[i];
	    int ncisq = nci * nci, nlev = RZXrows/nci;
	    double *bVi = REAL(VECTOR_ELT(bVarP, i)),
		*bi = b + Gp[i], *mm = REAL(VECTOR_ELT(gradComp, i)),
		*tmp = Memcpy(Calloc(ncisq, double),
			      REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(OmegaP, i)),
					    Matrix_xSym)), ncisq),
		*RZXi = RZXinv + Gp[i],
		dlev = (double) nlev,
		one[] = {1,0}, zero[] = {0,0};

	    if (nci == 1) {
		int ione = 1;
		mm[0] = ((double) nlev)/(tmp[0] * tmp[0]);
		mm[1] = alpha * F77_CALL(ddot)(&nlev, bi, &ione, bi, &ione);
		mm[2] = 0.;
		for (k = 0; k < nlev; k++) mm[2] += bVi[k];
		mm[3] = 0.;
		for (j = 0; j < p; j++) {
		    mm[3] += F77_CALL(ddot)(&RZXrows, RZXi + j * q, &ione,
					    RZXi + j * q, &ione);
		}
	    } else {
		AZERO(mm, 4 * ncisq);
		F77_CALL(dtrtri)("U", "N", &nci, tmp, &nci, &j);
		if (j)
		    error(_("Omega[[%d]] is not positive definite"), i + 1);
		F77_CALL(dsyrk)("U", "N", &nci, &nci, &dlev, tmp, &nci,
				zero, mm, &nci);
		mm += ncisq;	/* \bB_i term */
		F77_CALL(dsyrk)("U", "N", &nci, &nlev, &alpha, bi, &nci,
				zero, mm, &nci);
		mm += ncisq;     /* Sum of diagonal blocks of the inverse
				  * (Z'Z+Omega)^{-1} */
		for (j = 0; j < ncisq; j++) {
		    for (k = 0; k < nlev; k++) mm[j] += bVi[j + k*ncisq];
		}
		mm += ncisq;	/* Extra term for \vb */
		for (j = 0; j < p; j++) {
		    F77_CALL(dsyrk)("U", "N", &nci, &nlev, one,
				    RZXi + j * q, &nci,
				    one, mm, &nci);
		}
	    }
	    Free(tmp);
	}
	Free(iperm); Free(L);
	status[2] = 1; status[3] = 0;
    }
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
SEXP mer_gradient(SEXP x, SEXP pType)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    SEXP gradComp = GET_SLOT(x, Matrix_gradCompSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	dind, i, ifour = 4, ione = 1, nf = length(Omega),
	odind, ptyp = asInteger(pType);
    int REML =
	!strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))), "REML");
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double cc[] = {-1, 1, 1, REML ? 1 : 0},
	one = 1.0, zero = 0.0;

    mer_gradComp(x);
    dind = 0;			/* index into val for diagonals */
    for (i = 0; i < nf; i++) {
	SEXP Omgi = VECTOR_ELT(Omega, i);
	int nci = nc[i], ncisqr = nci * nci;
	double *tmp = Calloc(ncisqr, double);

	F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &one,
			REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
			cc, &ifour, &zero, tmp, &ncisqr);
	if (nci == 1) {
	    double omg = *REAL(GET_SLOT(Omgi, Matrix_xSym));
	    REAL(val)[dind++] =
		(ptyp?((ptyp == 1)? omg : -omg * omg) : 1) * tmp[0];
	} else {
	    int ii, j, ncip1 = nci + 1;

	    odind = dind + nci; /* index into val for off-diagonals */
	    if (ptyp) {
		double *chol = REAL(GET_SLOT(dpoMatrix_chol(Omgi), Matrix_xSym)),
		    *tmp2 = Calloc(ncisqr, double);

		/* Overwrite the gradient with respect to positions in
		 * Omega[[i]] by the gradient with respect to the
		 * unconstrained parameters.*/

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
    UNPROTECT(1);
    return val;
}

/**
 * Create and insert initial values for Omega.
 *
 * @param x pointer to an mer object
 *
 * @return NULL
 */
SEXP mer_initial(SEXP x)
{
    SEXP Omg = GET_SLOT(x, Matrix_OmegaSym),
	ZtZ = GET_SLOT(x, Matrix_ZtZSym);
    int	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*p = INTEGER(GET_SLOT(ZtZ, Matrix_pSym)),
	i, nf = length(Omg);
    double *xx = REAL(GET_SLOT(ZtZ, Matrix_xSym));

    for (i = 0; i < nf; i++) {
	SEXP Omgi = VECTOR_ELT(Omg, i);
	double *omgi = REAL(GET_SLOT(Omgi, Matrix_xSym));
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
	SET_SLOT(Omgi, Matrix_factorSym, allocVector(VECSXP, 0));
	dpoMatrix_chol(Omgi);
    }
    flag_not_factored(x);
    return R_NilValue;
}

/**
 * Extract the conditional modes of the random effects.
 *
 * @param x Pointer to an mer object
 *
 * @return a list of matrices containing the conditional modes of the
 * random effects
 */
SEXP mer_ranef(SEXP x)
{
    SEXP cnames = GET_SLOT(x, Matrix_cnamesSym),
	flist = GET_SLOT(x, Matrix_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ii, jj,
	nf = length(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double *b = REAL(GET_SLOT(x, Matrix_ranefSym));

    mer_secondary(x);
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(flist, R_NamesSymbol)));
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
		mm[ii + jj * mi] = bi[jj + ii * nci];
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Update the secondary slots - fixef and ranef
 * 
 * @param x pointer to a mer object
 * 
 */
SEXP mer_secondary(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));

    mer_factor(x);
    if (!status[1]) {
	internal_mer_fixef(x);
	internal_mer_ranef(x);
    }
    return R_NilValue;
}

/**
 * Extract the ML or REML conditional estimate of sigma
 *
 * @param x pointer to an mer object
 * @param REML logical scalar - TRUE if REML estimates are requested
 *
 * @return pointer to a numeric scalar
 */
SEXP mer_sigma(SEXP x, SEXP REML)
{
    return ScalarReal(
	internal_mer_sigma(x,
			   (REML == R_NilValue) ? -1 :
			   (asLogical(REML))));
}

/** 
 * Simulate a set of linear predictors from the random effects part of
 * an mer object
 * 
 * @param x Pointer to an mer object
 * @param np Pointer to an integer giving the number of values to simulate
 * @param useScale Logical indicator of whether to use the scale
 * 
 * @return a matrix of simulated linear predictors
 */
SEXP mer_simulate(SEXP x, SEXP nsimP)
{
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),"REML"),
	i, ii, j, nsim = asInteger(nsimP),
	nf = LENGTH(GET_SLOT(x, Matrix_OmegaSym)),
	n = LENGTH(GET_SLOT(x, Matrix_ySym)),
	q = LENGTH(GET_SLOT(x, Matrix_ZtySym));
    SEXP ans = PROTECT(allocMatrix(REALSXP, n, nsim)),
	Omega = GET_SLOT(x, Matrix_OmegaSym);
    cholmod_dense *cha = as_cholmod_dense(ans),
	*chb = cholmod_allocate_dense(q, nsim, q, CHOLMOD_REAL, &c);
    double one[] = {1,0}, zero[] = {0,0},
	scale = (asLogical(GET_SLOT(x, Matrix_useScaleSym)) ?
		 internal_mer_sigma(x, REML) : 1);
    cholmod_sparse *Zt = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtSym));
	
    GetRNGstate();
    for (ii = 0; ii < nsim; ii++) {
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], relen = Gp[i + 1] - Gp[i];
	    int nlev = relen/nci;
	    double *bi = (double *)(chb->x) + ii * q + Gp[i],
		*Rmat = REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				      Matrix_xSym));

	    for (j = 0; j < relen; j++) bi[j] = norm_rand();
	    F77_CALL(dtrsm)("L", "U", "N", "N", &nci, &nlev, &scale,
			    Rmat, &nci, bi, &nci);
	}
    }
    PutRNGstate();

    if (!cholmod_sdmult(Zt, 1, one, zero, chb, cha, &c))
	error(_("cholmod_sdmult failed"));
    cholmod_free_dense(&chb, &c);
    Free(Zt); Free(cha);
    UNPROTECT(1);
    return ans;
}

/** 
 * Externally callable version of internal_mer_update_ZXy
 * 
 * @param x pointer to an mer object
 * 
 * @return NULL
 */
SEXP mer_update_ZXy(SEXP x)
{
    internal_mer_update_ZXy(x,
			    INTEGER(GET_SLOT(GET_SLOT(x, Matrix_LSym),
					     Matrix_permSym)));
    return R_NilValue;
}

/** 
 * Update the y slot (and slots derived from it) in an mer object
 * 
 * @param x pointer to an mer object
 * @param ynew pointer to a numeric vector of length n
 * 
 * @return NULL
 */
SEXP mer_update_y(SEXP x, SEXP ynew)
{
    SEXP y = GET_SLOT(x, Matrix_ySym),
	Xty = GET_SLOT(x, Matrix_XtySym),
	Zty = GET_SLOT(x, Matrix_ZtySym);
    cholmod_factor *L = as_cholmod_factor(GET_SLOT(x, Matrix_LSym));
    int *perm = (int*)(L->Perm), i, ione = 1,
	n = LENGTH(y), p = LENGTH(Xty), q = LENGTH(Zty);
    cholmod_sparse *Zt = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtSym));
    cholmod_dense *td1, *yd = as_cholmod_dense(GET_SLOT(x, Matrix_ySym));
    double one[] = {1,0}, zero[] = {0,0};

    if (!isReal(ynew) || LENGTH(ynew) != n)
	error(_("ynew must be a numeric vector of length %d"), n);
    Memcpy(REAL(y), REAL(ynew), n);
    				/* y'y */
    REAL(GET_SLOT(x, Matrix_devCompSym))[2] =
	F77_CALL(ddot)(&n, REAL(y), &ione, REAL(y), &ione); 
				/* PZ'y into Zty */
    td1 = cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!cholmod_sdmult(Zt, 0, one, zero, yd, td1, &c))
	error(_("cholmod_sdmult failed"));
    for (i = 0; i < q; i++) REAL(Zty)[i] = ((double *)(td1->x))[perm[i]];
    cholmod_free_dense(&td1, &c); Free(yd); Free(Zt);
    				/* Xty */
    F77_CALL(dgemv)("T", &n, &p, one, REAL(GET_SLOT(x, Matrix_XSym)),
		    &n, REAL(y), &ione, zero, REAL(Xty), &ione);
    flag_not_factored(x);
    Free(L);
    return R_NilValue;

}

/**
 * Check validity of an mer object.
 *
 * @param x Pointer to an mer object
 *
 * @return TRUE if the object is a valid lmer object, else a string
 * describing the nature of the violation.
 */
SEXP mer_validate(SEXP x)
{
    SEXP
	/* ZZxP = GET_SLOT(x, Matrix_ZZxSym), */
	ZtXP = GET_SLOT(x, Matrix_ZtXSym),
	XtXP = GET_SLOT(x, Matrix_XtXSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	RXXP = GET_SLOT(x, Matrix_RXXSym)
	/* , cnames = GET_SLOT(x, Matrix_cnamesSym) */
	;
    int *ZtXd = INTEGER(getAttrib(ZtXP, Matrix_DimSym)),
	*XtXd = INTEGER(getAttrib(XtXP, Matrix_DimSym));

    if (!match_mat_dims(ZtXd, INTEGER(getAttrib(RZXP, Matrix_DimSym))))
	return mkString(_("Dimensions of slots ZtX and RZX must match"));
    if (!match_mat_dims(XtXd, INTEGER(getAttrib(RXXP, Matrix_DimSym))))
	return mkString(_("Dimensions of slots XtX and RXX must match"));
    if (ZtXd[1] != XtXd[0] || XtXd[0] != XtXd[1])
	return mkString(_("Slot XtX must be a square matrix with same ncol as ZtX"));
    return ScalarLogical(1);
}

/** 
 * Create the sparse Zt matrix from a factor list and list of model matrices
 * 
 * @param fl list of factors
 * @param Ztl list of transposes of model matrices
 * 
 * @return a freshly created sparse Zt object
 */
SEXP Zt_create(SEXP fl, SEXP Ztl)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))), fi, tmmat;
    int *dims, *p, *ii, i, nrtot = 0, nf = LENGTH(fl), nobs;
    int *Gp = Calloc(nf + 1, int), *nr = Calloc(nf, int),
	*offset = Calloc(nf + 1, int);
    double *x;
    
    if (!isNewList(fl) || nf < 1)
	error(_("fl must be a non-null list"));
    if (!isNewList(Ztl) || LENGTH(Ztl) != nf)
	error(_("Ztl must be a list of length %d"), nf);
    fi = VECTOR_ELT(fl, 0);
    nobs = LENGTH(fi);
    if (!isFactor(fi) || nobs < 1)
	error(_("fl[[1]] must be a non-empty factor"));
    offset[0] = Gp[0] = 0;
    for (i = 0; i < nf; i++) {	/* check consistency and establish dimensions */
	fi = VECTOR_ELT(fl, i);	/* grouping factor */
	if (!isFactor(fi) || LENGTH(fi) != nobs)
	    error(_("fl[[%d]] must be a factor of length %d"), i + 1, nobs);
	tmmat = VECTOR_ELT(Ztl, i); /* transpose of model matrix */
	if (!isMatrix(tmmat) || !isReal(tmmat))
	    error(_("Ztl[[%d]] must be real matrix"), i + 1);
	dims = INTEGER(getAttrib(tmmat, R_DimSymbol));
	if (dims[1] != nobs)
	    error(_("Ztl[[%d]] must have %d columns"), i + 1, nobs);
	nrtot += (nr[i] = dims[0]);
	offset[i + 1] = offset[i] + nr[i];
	Gp[i + 1] = Gp[i] + dims[0] * LENGTH(getAttrib(fi, R_LevelsSymbol));
    }
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = Gp[nf]; dims[1] = nobs;
    p = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, nobs + 1));
    ii = INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nrtot * nobs));
    x = REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nrtot * nobs));
    p[0] = 0; for(i = 0; i < nobs; i++) p[i + 1] = p[i] + nrtot;

    for (i = 0; i < nf; i++) {	/* fill ans */
	int *vals = INTEGER(VECTOR_ELT(fl, i)), j;
	double *xvals = REAL(VECTOR_ELT(Ztl, i));

	for (j = 0; j < nobs; j++) {
	    int jbase = Gp[i] + nr[i] * (vals[j] - 1), k;
	    for (k = 0; k < nr[i]; k++) {
		int ind = j * nrtot + offset[i] + k;
		ii[ind] = jbase + k;
		x[ind] = xvals[j * nr[i] + k];
	    }
	}
    }

    Free(offset); Free(Gp); Free(nr);
    UNPROTECT(1);
    return ans;
}

#if 0

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
		    GS->X, &(GS->n), fixed,
		    &ione, &zero, REAL(GS->eta), &ione);
				/* add in random effects and offset */
    vecIncrement(REAL(GS->eta), GS->off, GS->n);
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0, ans = 0; i < GS->n; i++) ans += REAL(devs)[i];
    UNPROTECT(1);
    return ans;
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

/* FIXME: Separate the calculation of the offset random effects from
 * the calculation of the deviance contributions. */

/** 
 * Determine the deviance components associated with each of the
 * levels of a grouping factor at the conditional modes or a value
 * offset from the conditional modes by delb.
 * 
 * @param GS pointer to a GlmerStruct
 * @param b conditional modes of the random effects 
 * @param Gp group pointers
 * @param nc number of columns in the model matrix for the kth
 * grouping factor
 * @param k index (0-based) of the grouping factor
 * @param delb vector of length nc giving the changes in the
 * orthonormalized random effects 
 * @param OmgFac Cholesky factor of the inverse of the penalty matrix
 * for this grouping factor 
 * @param bVfac 3-dimensional array holding the factors of the
 * conditional variance-covariance matrix of the random effects 
FIXME: This is wrong.  It is bVar[[i]] not bVfac that is being passed.
This only affects the AGQ method.
 * @param devcmp array to hold the deviance components
 * 
 * @return devcmp
 */
static double*
rel_dev_1(GlmerStruct GS, const double b[], int nlev, int nc, int k,
	  const double delb[], const double OmgFac[],
	  const double bVfac[], double devcmp[])
{
    SEXP devs;
    int *fv = INTEGER(VECTOR_ELT(GET_SLOT(GS->mer, Matrix_flistSym), k)),
	i, j;
    double *bcp = (double *) NULL; 

    AZERO(devcmp, nlev);
    if (delb) {
	int ione = 1, ntot = nlev * nc;
	double sumsq = 0;
				/* copy the contents of b */
	bcp = Memcpy(Calloc(ntot, double), b, ntot);
	if (nc == 1) {
	    sumsq = delb[0] * delb[0];
	    for (i = 0; i < nlev; i++) b[i] += delb[0] * bVfac[i];
	} else {
	    int ncsq = nc * nc;
	    double *tmp = Calloc(nc, double);
	    for (i = 0; i < nlev; i++) {
		Memcpy(tmp, delb, nc);
		F77_CALL(dtrmv)("U", "N", "N", &nc, &(bVfac[i * ncsq]),
				&nc, tmp, &ione);
		for (j = 0; j < nc; j++) b[i + j * nc] = tmp[j];
	    }
				/* sum of squares of delb */
	    for (j = 0; j < nc; j++) sumsq += delb[j] * delb[j];
	}
	for (i = 0; i < nlev; i++) devcmp[i] = -sumsq;
    }
    internal_mer_fitted(GS->mer, GS->offset, REAL(GS->eta));
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0; i < GS->n; i++)
	devcmp[fv[i] - 1] += REAL(devs)[i];
    UNPROTECT(1);
    if (nc == 1) {
	for (i = 0; i < nlev; i++) {
	    double tmp = *OmgFac * b[i];
	    devcmp[i] += tmp * tmp;
	}
    } else {
	double *tmp = Calloc(nc, double);
	int ione = 1;
	
	for (i = 0; i < nlev; i++) {
	    for (j = 0; j < nc; j++) tmp[j] = b[i + j * nlev];
	    F77_CALL(dtrmv)("U", "N", "N", &nc, OmgFac, &nc,
			    tmp, &ione);
	    for (j = 0; j < nc; j++) 
		devcmp[i] += tmp[j] * tmp[j];
	}
    }
    if (delb) {
	Memcpy(b, bcp, ntot);
	Free(bcp);
    }
    return devcmp;
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
	internal_mer_coef(GS->mer, 2, REAL(varc));
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

SEXP glmer_bhat(SEXP pars, SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);

    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	warning(_("internal_bhat did not converge"));
    return R_NilValue;
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
				
    AZERO(GS->off, GS->n); /* fitted values from random effects */
/*     fitted_ranef(GET_SLOT(GS->mer, Matrix_flistSym), GS->unwtd, b, */
/* 		 INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)), GS->off); */
    for (i = 0; i < GS->n; i++)
	(GS->etaold)[i] = ((GS->off)[i] += (GS->offset)[i]);
    
    for (it = 0, crit = GS->tol + 1;
	 it < GS->maxiter && crit > GS->tol; it++) {
				/* fitted values from current beta */
	F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
			GS->X, &(GS->n), md,
			&ione, &zero, REAL(GS->eta), &ione);
				/* add in random effects and offset */
	vecIncrement(REAL(GS->eta), (GS->off), GS->n);
				/* check for convergence */
	crit = conv_crit(GS->etaold, REAL(GS->eta), GS->n);
				/* obtain mu, dmu_deta, var */
	eval_check_store(GS->linkinv, GS->rho, GS->mu);
	dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				      REALSXP, GS->n));
	var = PROTECT(eval_check(GS->var, GS->rho, REALSXP, GS->n));
				/* calculate weights and working residual */
	for (i = 0; i < GS->n; i++) {
	    w[i] = GS->wts[i] * REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	    z[i] = w[i] * (REAL(GS->eta)[i] - (GS->off)[i] +
			   ((GS->y)[i] - REAL(GS->mu)[i]) /
			   REAL(dmu_deta)[i]);
	}
	UNPROTECT(2);
				/* weighted copy of the model matrix */
	for (j = 0; j < GS->p; j++)
	    for (i = 0; i < GS->n; i++)
		wtd[i + j * GS->n] = GS->X[i + j * GS->n] * w[i];
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
    Free(ans); Free(md); Free(w);
    Free(work); Free(wtd); Free(z);
    return fixed;
}
    
static void
internal_glmer_ranef_update(GlmerStruct GS, SEXP b)
{
    SEXP bhat, bprop = PROTECT(duplicate(b)), 
	bVar = GET_SLOT(GS->mer, Matrix_bVarSym),
	Omega = GET_SLOT(GS->mer, Matrix_OmegaSym);
    int i, ione = 1, j, k;
    double *b = REAL(GET_SLOT(GS->mer, Matrix_ranefSym);
    double devr, one = 1;

    bhat = R_NilValue;
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

 				/* subtract quadratic form in */
 				/* Omega[[i]] at b  */
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
 				/* add quadratic form in */
 				/* Omega[[i]] at bprop  */
 	Memcpy(bcopy, bpropi, ntot);
 	F77_CALL(dtrmm)("R", "U", "T", "N", &nlev, &nci, &one,
 			omgfac, &nci, bcopy, &nlev);
 	for (k = 0; k < ntot; k++) devr += bcopy[k] * bcopy[k];

 	Free(bcopy); Free(chol); Free(delta); Free(omgfac);
     }
 				/* add deviance at bprop */
/*     devr += random_effects_deviance(GS, bprop); */

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
	for (k = 0; k < 5; k++) Rprintf("%#10g ", b0[k]);
	Rprintf("\n");
    }
	
    UNPROTECT(2);
}


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
    SEXP Omega = GET_SLOT(GS->mer, Matrix_OmegaSym),
	bVar = GET_SLOT(GS->mer, Matrix_bVarSym);
    int i, j, k, nAGQ = asInteger(nAGQp);
    int n2 = (nAGQ + 1)/2,
	*Gp = INTEGER(GET_SLOT(GS->mer, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(GS->mer, Matrix_ncSym));
    double *f0, LaplaceDev = 0, AGQadjst = 0,
	*bhat = REAL(GET_SLOT(GS->mer, Matrix_ranefSym));
	
    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (GS->nf > 1 && nAGQ > 1) {
	warning(_("AGQ not available for multiple grouping factors - using Laplace"));
	nAGQ = 1;
    }
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	return ScalarReal(DBL_MAX);

    for (i = 0; i < GS->nf; i++) {
	int nci = nc[i];
	int ncip1 = nci + 1, ncisqr = nci * nci,
	    nlev = (Gp[i + 1] - Gp[i])/nci;
	double *omgf = REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(Omega, i)), Matrix_xSym)),
	    *bVi = Memcpy(Calloc(ncisqr * nlev, double),
			   REAL(VECTOR_ELT(bVar, i)), ncisqr * nlev);

        for (j = 0; j < nci; j++) { /* nlev * logDet(Omega_i) */
            LaplaceDev += 2 * nlev * log(omgf[j * ncip1]);
        }
        for (k = 0; k < nlev; k++) {
	    double *bVik = bVi + k * ncisqr;
            F77_CALL(dpotrf)("U", &nci, bVik, &nci, &j);
            if (j)
                error(_("Leading %d minor of bVar[[%d]][,,%d] not positive definite"),
                      j, i + 1, k + 1);
            for (j = 0; j < nci; j++) LaplaceDev -= 2 * log(bVik[j * ncip1]);
        }

	f0 = Calloc(nlev, double);
	rel_dev_1(GS, bhat, nlev, nci, i, (double *) NULL,
		  omgf, bVi, f0);
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
			      omgf, bVi, fx);
		    for (j = 0; j < nlev; j++) {
			rellik[j] += GHQ_w[nAGQ][k] *
			    exp(-(fx[j] - f0[j])/2);
		    }
		    delb[0] *= -1;
		    rel_dev_1(GS, bhat, nlev, nci, i, delb,
			      omgf, bVi, fx);
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
	Free(f0); Free(bVi);
    }
    return ScalarReal(LaplaceDev + AGQadjst);
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


#endif
