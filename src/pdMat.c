#include "Mutils.h"
#include <R_ext/Lapack.h>

SEXP pdCompSymm_pdFactor(SEXP pd)
{
    int i, j, ncol = asInteger(GET_SLOT(pd, install("ncol")));
    double *param = REAL(GET_SLOT(pd, install("param"))), sd,
        lcorr, corr, corr2, *L, offdiag;
    SEXP retval;
    
    retval = PROTECT(allocMatrix(REALSXP, ncol, ncol));
    L = REAL(retval);
    sd = exp(param[0]);
    lcorr = exp(param[1]);	/* logistic of the correlation */
    corr = (lcorr - 1.)/((double) ncol - 1.)/(lcorr + 1.);
    corr2 = sd * sqrt(1. - corr);
    corr = sd * sqrt((1. + (ncol - 1.) * corr) / ((double) ncol));
    
    for (i = 0; i < ncol; i++) {
	L[i * ncol] = corr;
    }
    
    for (i = 1; i < ncol; i++) {
	offdiag = -corr2/sqrt((double) (i * (i + 1)));
	for (j = 0; j < i; j++) {
	    L[i + (j * ncol)] = offdiag;
	}
	L[i * (ncol + 1)] = -sd * i;
    }
    UNPROTECT(1);
    return retval;
}

SEXP nlme_Chol(SEXP A)
{
    SEXP ans = PROTECT((TYPEOF(A) == REALSXP)?duplicate(A):
		       coerceVector(A, REALSXP));
    SEXP adims = getAttrib(A, R_DimSymbol);
    int m = INTEGER(adims)[0];
    int n = INTEGER(adims)[1];
    int i, j;

    if (!(isMatrix(A)))
	error("A must be a numeric (real) matrix");
    if (m != n) error("A must be a square matrix");
    for (j = 0; j < n; j++) {	/* zero the lower triangle */
	for (i = j+1; i < n; i++) {
	    REAL(ans)[i + j * n] = 0.;
	}
    }

    F77_CALL(dpotrf)("Upper", &m, REAL(ans), &m, &i);
    nlme_check_Lapack_error(i, "dpotrf");
    UNPROTECT(1);
    return ans;
}


