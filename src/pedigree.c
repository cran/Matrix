				/* Methods for pedigrees */
#include "pedigree.h"

/**
 * Create the left Cholesky factor of the numerator relationship
 * matrix from a pedigree.
 *
 * @param x a pedigree object
 * @param tinv the inverse of T, a unit lower dtCMatrix
 * @param ttrans the transpose of T, a unit upper dtCMatrix
 * @param ans T stored as a non-unit, lower dtCMatrix
 * the pedigree
 *
 * @return ans with elements modified to incorporate D
 */
SEXP pedigree_chol(SEXP x, SEXP ans)
{
    SEXP Sire = GET_SLOT(x, install("sire"));
    int *ai = INTEGER(GET_SLOT(ans, Matrix_iSym)),
	*ap = INTEGER(GET_SLOT(ans, Matrix_pSym)),
	*dam = INTEGER(GET_SLOT(x, install("dam"))),
	*sire = INTEGER(Sire), 
	i, j, n = LENGTH(Sire);
    double *ax = REAL(GET_SLOT(ans, Matrix_xSym)), *F, Di, tmp;

    setAttrib(ans, install("F"), allocVector(REALSXP, n));
    F = REAL(getAttrib(ans, install("F")));
    for (i = 0; i < n; i++) {
	int k, p = sire[i] - 1, q = dam[i] - 1;
	if (sire[i] == NA_INTEGER) {
	    F[i] = 0;
	    Di = (dam[i] == NA_INTEGER) ? 1 : sqrt(0.75 - 0.25 * F[q]);
	} else {
	    if (dam[i] == NA_INTEGER) { /* sire only */
		F[i] = 0;
		Di = sqrt(0.75 - 0.25 * F[p]);
	    } else {		/* both parents in pedigree */
		Di = sqrt(0.5 - 0.25 * (F[p] + F[q]));
		F[i] = NA_REAL;
		if ((ap[i + 1] - ap[i]) > 1) {	  /* skip if no progeny */
		    if (p > q) {j = p; p = q; q = j;} /* ensure p <= q */
		    for (j = 0, F[i] = 0; j <= p; j++) {
			for (k = ap[j], tmp = 0;
			     k < ap[j + 1] && ai[k] <= q; k++) {
			    int ii = ai[k];
			    if (ii == p) tmp = ax[k];
			    if (ii == q) F[i] += tmp * ax[k]/2;
			}
		    }
		}
	    }
	}
	for (j = ap[i]; j < ap[i + 1]; j++) ax[j] *= Di;
    }

    return ans;
}

    
