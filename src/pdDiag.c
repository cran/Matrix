#include "Mutils.h"

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
double pdDiag_ld_factor_from_par(const double *par, double *factor, int nc)
{
    int i, j;
    double ld = 0.;

    for (i = 0; i < nc; i++) {
	double theta = par[i];
	ld += theta;
	factor[i * (nc + 1)] = exp(theta);
	for (j = 0; j < i; j++) {
	    factor[i * nc + j] = factor[j * nc + i] = 0.;
	}
    }
    return ld;
}

SEXP
pdDiag_coefGets(SEXP x, SEXP value)
{
    SEXP val = PROTECT((TYPEOF((SEXP) value) == REALSXP) ? (SEXP) value :
			coerceVector((SEXP) value, REALSXP));
    SEXP param = GET_SLOT(x, install("param"));
    int npar = LENGTH(param);

    if (npar == 0) {		/* uninitialized */
	int lv = LENGTH(val);
	SEXP factor;
	
	if (lv < 1) error("Replacement coefficients must have length > 0");
	SET_SLOT(x, install("param"), duplicate(val));
	SET_SLOT(x, install("Ncol"), ScalarInteger(lv));
	SET_SLOT(x, install("factor"), allocMatrix(REALSXP, lv, lv));
	factor = GET_SLOT(x, install("factor"));
	
	SET_SLOT(x, install("logDet"), ScalarReal(
                 pdDiag_ld_factor_from_par(REAL(val), REAL(factor), lv)));
    } else {
	if (npar != LENGTH(val))
	    error("Cannot change length of parameter vector from %d to %d", npar,
		  LENGTH(val));
	memcpy(REAL(param), REAL(val), npar * sizeof(double));
	REAL(GET_SLOT(x, install("logDet")))[0] =
	    pdDiag_ld_factor_from_par(REAL(param),
				      REAL(GET_SLOT(x, install("factor"))),
				      npar);
    }
    UNPROTECT(1);
    return x;
}

SEXP
pdDiag_LMEgradient(SEXP x, SEXP Ain,
		SEXP nlev)
{
    int nlevVal = asInteger((SEXP)nlev);
    SEXP param = GET_SLOT((SEXP)x, install("param"));
    double* Amat = REAL(PROTECT((TYPEOF((SEXP)Ain) == REALSXP) ? (SEXP)Ain :
                                coerceVector((SEXP)Ain, REALSXP)));
    SEXP dims = getAttrib((SEXP)Ain, R_DimSymbol);
    int m = INTEGER(dims)[0];
    int n = INTEGER(dims)[1];
    int q = LENGTH(param);
    SEXP retval = allocVector(REALSXP, q);
    int i, j;
    double sum;
    
    if (q <= 0) {
	error("Uninitialized pdDiag object");
    }
    if (m != n || m != q) {
	error("A must be a %d by %d matrix", q, q);
    }
    if (nlevVal <= 0) {
	error("nlev must by > 0");
    }
    for (j = 0; j < q; j++) {
	sum = 0.;
	for (i = 0; i < q; i++) {
	    sum += Amat[i + j * q] * Amat[i + j * q];
	}
	REAL(retval)[j] =
	    ((double) nlevVal) - (sum * exp(REAL(param)[j] * 2.));
    }
    UNPROTECT(1);
    return retval;
}

SEXP pdDiag_EMupdate(SEXP x, SEXP nlev, SEXP Ain)
{
    int nlevVal = asInteger((SEXP)nlev);
    SEXP param = GET_SLOT(x, install("param"));
    double* factor = REAL(GET_SLOT(x, install("factor")));
    double* Amat = REAL((TYPEOF((SEXP)Ain) == REALSXP)?duplicate((SEXP)Ain):
                        coerceVector((SEXP)Ain, REALSXP));
    int* dims = INTEGER(getAttrib((SEXP)Ain, R_DimSymbol));
    int m = dims[0];
    int n = dims[1];
    int q = LENGTH(param);
    int qp1 = q+1;
    int i, j;
    double ld = 0.0;
    double sum;
    
    if (q <= 0) {
	error("Uninitialized pdDiag object");
    }
    if (m != n || m != q) {
	error("A must be a %d by %d matrix", q, q);
    }
    if (nlevVal <= 0) {
	error("nlev must by > 0");
    }
    for (j = 0; j < q; j++) {
	sum = 0.;
	for (i = 0; i < q; i++) {
	    sum += Amat[i + j * q] * Amat[i + j * q];
	}
    ld += sum = -0.5*log((double) sum/nlevVal);
	REAL(param)[j] = sum;
    factor[j*qp1] = exp(sum);
    }
	REAL(GET_SLOT(x, install("logDet")))[0] = ld;
    
    return x;
}
