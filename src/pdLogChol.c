#include "Mutils.h"
#include <R_ext/Lapack.h>

/** 
 * Populate the factor from the parameter vector and return the logarithm
 * the determinant of the factor.
 * 
 * @param par vector of parameters
 * @param factor pointer to matrix to be overwritten with the factor
 * @param nc number of columns
 * 
 * @return logarithm of the determinant of the factor
 */
static
double ld_factor_from_par(const double *par, double *factor, int nc)
{
    int i, j, k;
    double ld = 0.;

    for (i = 0; i < nc; i++) {
	double theta = par[i];
	ld += theta;
	factor[i * (nc + 1)] = exp(theta);
    }
    for (j = 1, k = nc; j < nc; j++) {
	for (i = 0; i < j; i++) {
	    factor[i + j * nc] = par[k];
	    factor[j + i * nc] = 0.;
	    k++;
	}
    }
    return ld;
}

/** 
 * An internal function that calculates the gradient of the
 * positive-definite matrix with respect to the parameters.
 * This function is used in both pdLogChol_LMEgradient and
 * pdLogChol_pdgradient.
 * 
 * @param nc number of columns (and rows) in the matrix
 * @param pars parameter vector of length nc*(nc+1)/2
 * @param value array into which the results are written
 *
 * @return the gradient in value
 */
static double*
gradient(const int nc, const double *factor, const double *pars,
	 double *value)
{
    int i, j, k, ncsq = nc * nc, ncp1 = nc + 1;
    double *facprime = Calloc(ncsq, double);
    const double one_d = 1.0;
    const double zero_d = 0.0;

    for (j = 0; j < ncsq; j++) facprime[j] = 0.0; /* is this necessary? */
    for (i = 0; i < nc; i++) { /* pars for diagonal elements of factor */
	facprime[i * ncp1] = exp(pars[i]);
	F77_CALL(dsyr2k)("U", "T", &nc, &nc, &one_d, factor,
			 &nc, facprime, &nc, &zero_d,
			 value + i * ncsq, &nc);
	nlme_symmetrize(value + i * ncsq, nc);
	facprime[i * ncp1] = 0.;
    }
    i = nc;			/* off-diagonals */
    for (j = 1; j < nc; j++) {
	for (k = 0; k < j; k++) {
	    facprime[k + j * nc] = 1.;
	    F77_CALL(dsyr2k)("U", "T", &nc, &nc, &one_d, factor,
			     &nc, facprime, &nc, &zero_d,
			     value + i * ncsq, &nc);
	    nlme_symmetrize(value + i * ncsq, nc);
	    facprime[k + j * nc] = 0.;
	    i++;
	}
    }
    Free(facprime);
    return value;

}

SEXP
pdLogChol_LMEhessian(SEXP x, SEXP Ain,
                     SEXP Hin, SEXP nlev)
{
    SEXP param = GET_SLOT((SEXP) x, install("param"));
    double* pars;
    int parlen = LENGTH(param);
    int ncol = asInteger(GET_SLOT((SEXP) x, install("Ncol")));
    int ncolsq = ncol*ncol;
    SEXP retval = PROTECT(allocMatrix(REALSXP, parlen, parlen));
    double* dretval = REAL(retval);
    double* factor = REAL(GET_SLOT((SEXP) x, install("factor")));
    int nlevVal = asInteger((SEXP) nlev);
    int* dims = INTEGER(getAttrib((SEXP)Ain, R_DimSymbol));
    int m = dims[0];
    int n = dims[1];
    double *grad = Calloc(ncol * ncol * parlen, double);
    double* Amat = REAL(PROTECT((TYPEOF((SEXP)Ain) == REALSXP) ?
                                duplicate((SEXP) Ain):
                                coerceVector((SEXP) Ain, REALSXP)));
    double* Hmat = REAL(PROTECT((TYPEOF((SEXP)Hin) == REALSXP) ?
                                duplicate((SEXP) Hin):
                                coerceVector((SEXP) Hin, REALSXP)));
    int j, k;

    if (parlen <= 0) {
	error("Uninitialized pdLogChol object");
    }
    if (m != n || m != ncol) {
	error("A must be a %d by %d matrix", ncol, ncol);
    }
    if (nlevVal <= 0) {
	error("nlev must be > 0");
    }

    gradient(ncol, factor, REAL(param), grad);

    /* The terms not involving D_iD_j */
    for (j = 0; j < parlen; j++) {
        double* col = dretval + j*parlen;
        double* grad_j = grad + ncolsq*j;
        int i;
        
        for (i = j; i < parlen; i++, col++) {
            double* grad_i = grad + ncolsq*i;
            double sum = 0.0;

            for (k = 0; k < ncolsq; k++) {
                double grad_ik = grad_i[k];
                double* H_k = Hmat+k*ncolsq;
                int l;
                for (l = 0; l < ncolsq; l++) {
                    sum += grad_ik*grad_j[l]*H_k[l];
                }
            }

            *col = sum;
        }
    }

    /* The terms involving D_iD_j */

    pars = REAL(param);

    /* both parameters are from diagonal and they are same */
    for (j = 0; j < ncol; j++) {
        double sum = 4*exp(2*pars[j]);
        double* grad_j = grad+j*(ncolsq+ncol);
        int i;
        for (i = j+1; i < ncol; i++) {
            sum += grad_j[i]*(Amat[i+j*ncol]+Amat[j+i*ncol]);
        }
        dretval[i*(parlen+1)] += sum;
    }
    /* one is from diagonal and the other from off-diagonal */
    for (j = 0; j < ncol; j++) {
        double tmp = exp(pars[j]);
        double* col = dretval + parlen*j + ncol + j*(j-3)/2 - 1;
        int i;
        for (i = j+1; j < ncol; j++) {
            col[i] += tmp*(Amat[i+j*ncol]+Amat[j+i*ncol]);
        }
    }
    /* both parameters are off diagonal and they are same */
    k = 0;
    for (j = 1; j < ncol; j++) {
        int i;
        for (i = 0; i < j; i++) {
            dretval[k*(parlen+1)] += 2.0*Amat[j*(ncol+1)];
            k++;
        }
    }
    /* two different off diagonals */
    for (j = 1; j < ncol; j++) {
        int i;
        for (i = 1; i < j; i++) {
            double tmp = Amat[i+j*ncol]+Amat[j+i*ncol];
            for (k = 0; k < i; k++) {
                int ki = ncol+i*(i-1)/2+k;
                int kj = ncol+j*(j-1)/2+k;
                dretval[kj + ki*parlen] += tmp;
            }
        }
    }
    nlme_symmetrize(dretval, parlen);
    UNPROTECT(3);
    return retval;
}

/** 
 * LMEgradient implementation for the pdLogChol class
 * 
 * @param x Pointer to a pdLogChol object
 * @param Ain Pointer to an upper-triangular double precision square matrix
 * @param nlev Pointer to an integer scalar giving the number of levels
 * 
 * @return Pointer to a REAL gradient vector
 */
SEXP
pdLogChol_LMEgradient(SEXP x, SEXP Ain,
		      SEXP nlev)
{
    SEXP param = GET_SLOT((SEXP) x, install("param"));
    int parlen = LENGTH(param);
    int ncol = asInteger(GET_SLOT((SEXP) x, install("Ncol")));
    SEXP retval = PROTECT(allocVector(REALSXP, parlen));
    double* factor = REAL(GET_SLOT((SEXP) x, install("factor")));
    int nlevVal = asInteger((SEXP) nlev);
    int* dims = INTEGER(getAttrib((SEXP)Ain, R_DimSymbol));
    int m = dims[0];
    int n = dims[1];
    double *grad = Calloc(ncol * ncol * parlen, double);
    double* Amat = REAL((TYPEOF((SEXP)Ain) == REALSXP) ?
                        duplicate((SEXP) Ain):
                        coerceVector((SEXP) Ain, REALSXP));
    
    if (parlen <= 0) {
	error("Uninitialized pdLogChol object");
    }
    if (m != n || m != ncol) {
	error("A must be a %d by %d matrix", ncol, ncol);
    }
    if (nlevVal <= 0) {
	error("nlev must be > 0");
    }
    gradient(ncol, factor, REAL(param), grad);
/*     LMEgradient(factor, Amat, nlevVal, ncol, grad, parlen, */
/* 		REAL(retval)); */
    Free(grad);
    UNPROTECT(1);
    return retval;
}

/** 
 * Implementation of the pdgradient method for pdLogChol objects.
 * 
 * @param x Pointer to a pdLogChol object
 * 
 * @return SEXP of a three-dimensional array with the gradient of the
 * pdgradient with respect to the parameters.
 */
SEXP
pdLogChol_pdgradient(SEXP x)
{
    SEXP param = GET_SLOT((SEXP) x, install("param"));
    int parlen = LENGTH(param);
    int ncol = asInteger(GET_SLOT((SEXP) x, install("Ncol")));
    SEXP factor = GET_SLOT((SEXP) x, install("factor"));
    SEXP retval, dims;

    if (parlen <= 0) {
	error("Uninitialized pdLogChol object");
    }
    dims = allocVector(INTSXP, 3);
    INTEGER(dims)[0] = INTEGER(dims)[1] = ncol;
    INTEGER(dims)[2] = parlen;
    retval = allocArray(REALSXP, dims);
    gradient(ncol, REAL(factor), REAL(param), REAL(retval));

    return retval;
}

/** 
 * Perform an EM update on a pdLogChol object.
 * 
 * @param x Pointer to a pdLogChol object
 * @param nlev An integer object - the number of levels in the grouping factor
 * @param Ain An upper triangular matrix object
 * 
 * @return The updated pdLogChol object x
 */
SEXP
pdLogChol_EMupdate(SEXP x, SEXP nlev, SEXP Ain)
{
    SEXP param = GET_SLOT(x, install("param"));
    SEXP factor = GET_SLOT(x, install("factor"));
    double scal = asReal((SEXP) nlev);
    int ncol = asInteger(GET_SLOT(x, install("Ncol")));
    double* apt = REAL((TYPEOF((SEXP)Ain) == REALSXP) ? duplicate((SEXP)Ain) :
                       coerceVector((SEXP)Ain, REALSXP));
    int ncsq = ncol * ncol;
    int i, info, j, k, one_i = 1;
    double ld;
				/* overwrite A with inverse of A'A */
    F77_CALL(dpotri)("U", &ncol, apt, &ncol, &info);
    nlme_check_Lapack_error(info, "dpotri");
				/* scale by nlev */
    F77_CALL(dscal)(&ncsq, &scal, apt, &one_i);
				/* decompose A'A-inverse */
    F77_CALL(dpotrf)("U", &ncol, apt, &ncol, &info);
    nlme_check_Lapack_error(info, "dpotrf");
				/* copy to the factor */
    F77_CALL(dlacpy)("U", &ncol, &ncol, apt, &ncol,
		     REAL(factor), &ncol);
    ld = 0.;			/* create the parameter and the logDet */
    for (j = 0, k = ncol; j < ncol; j++) {
	ld += REAL(param)[j] = log(apt[j * (ncol + 1)]);
	for (i = 0; i < j; i++) {
	    REAL(param)[k] = apt[i + j * ncol];
	    k++;
	}
    }
    REAL(GET_SLOT((SEXP) x, install("logDet")))[0] = ld;

    return x;
}

SEXP
pdLogChol_coefGets(SEXP x, SEXP value)
{
    SEXP val = PROTECT((TYPEOF((SEXP) value) == REALSXP) ? (SEXP) value :
			coerceVector((SEXP) value, REALSXP));
    SEXP param = GET_SLOT((SEXP) x, install("param"));
    int npar = LENGTH(param);

    if (npar == 0) {		/* uninitialized */
	int lv = LENGTH(val);
	int nc = (int)(0.5 + (sqrt(8.*(double)lv + 1.) - 1.)/2.);
	SEXP factor;
	
	if (((nc * (nc + 1))/2) != lv)
	    error("parameter vector cannot have length %d", lv);
	SET_SLOT(x, install("param"), duplicate(val));
	SET_SLOT(x, install("Ncol"), ScalarInteger(nc));
	SET_SLOT(x, install("factor"), allocMatrix(REALSXP, nc, nc));
	factor = GET_SLOT(x, install("factor"));
	SET_SLOT(x, install("logDet"), ScalarReal(
                 ld_factor_from_par(REAL(val), REAL(factor), nc)));
    } else {
	if (npar != LENGTH(val))
	    error("Cannot change length of parameter vector from %d to %d", npar,
		  LENGTH(val));
	memcpy(REAL(param), REAL(val), npar * sizeof(double));
	REAL(GET_SLOT(x, install("logDet")))[0] =
	    ld_factor_from_par(REAL(param), REAL(GET_SLOT(x, install("factor"))),
			       asInteger(GET_SLOT(x, install("Ncol"))));
    }
    UNPROTECT(1);
    return x;
}
