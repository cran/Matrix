#include "Mutils.h"

SEXP pdIdent_gradient(SEXP x, SEXP Ain,
		      SEXP nlev)
{
    int nlevVal = asInteger((SEXP)nlev);
    SEXP param = GET_SLOT((SEXP)x, install("param"));
    int* dims = INTEGER(getAttrib((SEXP)Ain, R_DimSymbol));
    int m = dims[0];
    int n = dims[1];
    int q = asInteger(GET_SLOT((SEXP)x, install("Ncol")));
    double* Amat = REAL((TYPEOF((SEXP)Ain) == REALSXP) ? (SEXP)Ain :
                        coerceVector((SEXP)Ain, REALSXP));
    int i, j;
    double sum;
    
    if (q <= 0) {
	error("Uninitialized pdIdent object");
    }
    if (m != n || m != q) {
	error("A must be a %d by %d matrix", q, q);
    }
    if (nlevVal <= 0) {
	error("nlev must by > 0");
    }
    sum = 0.;
    for (j = 0; j < q; j++) {
	for (i = 0; i < q; i++) {
	    sum += Amat[i + j * q] * Amat[i + j * q];
	}
    }
    return ScalarReal((double)(q * nlevVal) -
		      sum * exp(REAL(param)[0] * 2.));
}

SEXP pdIdent_EMupdate(SEXP x, SEXP nlev, SEXP Ain)
{
    int nlevVal = asInteger(nlev);
    SEXP param = GET_SLOT(x, install("param"));
    int* dims = INTEGER(getAttrib(Ain, R_DimSymbol));
    int m = dims[0];
    int n = dims[1];
    int q = asInteger(GET_SLOT(x, install("Ncol")));
    double* Amat = REAL((TYPEOF(Ain) == REALSXP)?duplicate(Ain):
                        coerceVector(Ain, REALSXP));
    int i, j;
    double sum;
    
    if (q <= 0) {
	error("Uninitialized pdIdent object");
    }
    if (m != n || m != q) {
	error("A must be a %d by %d matrix", q, q);
    }
    if (nlevVal <= 0) {
	error("nlev must by > 0");
    }
    sum = 0.;
    for (j = 0; j < q; j++) {
	for (i = 0; i < q; i++) {
	    sum += Amat[i + j * q] * Amat[i + j * q];
	}
    }
    REAL(param)[0] = -log(sum/((double) nlevVal*q))/2;
    return x;
}
