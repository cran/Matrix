#include "Mutils.h"

static
void corr_from_par(const double *par, double *corr, int nc)
{
    int i, j, k = nc;
    
    for (i = 0; i < nc; i++) {
	corr[i * (nc + 1)] = 1.;
	for (j = 0; j < i; j++) {
		/* exponential of Fisher's z transformation */
	    double expz = exp(par[k]);

	    corr[i + nc * j] = corr[j + nc * i] = (expz - 1.)/(expz + 1.);
	    k++;
	}
    }
}

/** 
 * Evaluate the pdMatrix from a pdNatural object
 * 
 * @param x Pointer to a pdNatural object
 * 
 * @return A newly allocated matrix 
 */

SEXP
pdNatural_pdmatrix(SEXP x)
{
    SEXP param = GET_SLOT((SEXP) x, install("param"));
    int q = asInteger(GET_SLOT((SEXP) x, install("Ncol")));
    SEXP retval = PROTECT(allocMatrix(REALSXP, q, q));
    double *stdDev = Calloc(q, double);
    int i, j;

    for(i = 0; i < q; i++) stdDev[i] = exp(REAL(param)[i]);
    corr_from_par(REAL(param), REAL(retval), q);
    for(i = 0; i < q; i++) {
	for (j = 0; j < q; j++) {
	    REAL(retval)[i + j * q] *= stdDev[i] * stdDev[j];
	}
    }
    Free(stdDev);
    UNPROTECT(1);
    return(retval);
}

SEXP
pdNatural_corrmatrix(SEXP x)
{
    SEXP param = GET_SLOT((SEXP) x, install("param"));
    int q = asInteger(GET_SLOT((SEXP) x, install("Ncol")));
    SEXP corrmat = PROTECT(allocMatrix(REALSXP, q, q));

    corr_from_par(REAL(param), REAL(corrmat), q);
    UNPROTECT(1);
    return corrmat;
}

/** 
 * An internal function that calculates the gradient of the
 * positive-definite matrix with respect to the parameters.
 * This function is used in pdNatural_LMEgradient
 * 
 * @param nc number of columns (and rows) in the matrix
 * @param mat the positive definite matrix
 * @param value array into which the results are written
 *
 * @return the gradient in value
 */
static double*
gradient(int nc, const double *param, double *value)
{
    int i, ncsq = nc * nc, ncp1 = nc + 1;
    double* mat = Calloc(ncsq, double);
    double* stdDev = Calloc(nc, double);
    double* grad_i;

    memset(value, 0, sizeof(double)*ncsq*nc*(nc+1)/2);

    for(i = 0; i < nc; i++)
        stdDev[i] = exp(param[i]);
    corr_from_par(param, mat, nc);
    for(i = 0; i < nc; i++) {
        int j;
        for (j = 0; j < nc; j++) {
            mat[i + j * nc] *= stdDev[i] * stdDev[j];
        }
    }

    for (i = 0, grad_i = value; i < nc;         /* diagonals */
         i++, grad_i += ncsq) {
        int j;
        for (j = 0; j < nc; j++)
            grad_i[i+j*nc] = grad_i[j+i*nc] = mat[i+j*nc];
        grad_i[i*ncp1] = 2*mat[i*ncp1];
    }

    grad_i = value + nc*ncsq;                   /* off-diagonals */
    for (i = 1; i < nc; i++) {
        int j;
        for (j = 0; j < i; j++) {
            grad_i[i+j*nc] = grad_i[j+i*nc] =
                mat[i+j*nc] + stdDev[i]*stdDev[j];
            grad_i += ncsq;
        }
    }
    Free(mat);
    Free(stdDev);
    return value;
}

/** 
 * LMEgradient implementation for the pdNatural class
 * 
 * @param x Pointer to a pdNatural object
 * @param Ain Pointer to an upper-triangular double precision square matrix
 * @param nlev Pointer to an integer scalar giving the number of levels
 * 
 * @return Pointer to a REAL gradient vector
 */
SEXP
pdNatural_LMEgradient(SEXP x, SEXP Ain, SEXP nlev)
{
    SEXP param = GET_SLOT((SEXP) x, install("param"));
    int parlen = LENGTH(param);
    int ncol = asInteger(GET_SLOT((SEXP) x, install("Ncol")));
    SEXP retval = PROTECT(allocVector(REALSXP, parlen));
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
    gradient(ncol, REAL(param), grad);
    LMEgradient(REAL(GET_SLOT((SEXP)x, install("factor"))),
                Amat, nlevVal, ncol, grad, parlen,
                REAL(retval));
    Free(grad);
    UNPROTECT(1);
    return retval;
}


