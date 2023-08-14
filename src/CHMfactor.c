				/* CHOLMOD factors */
#include "CHMfactor.h"

/**
 * Evaluate the logarithm of the square of the determinant of L
 *
 * @param f pointer to a CHMfactor object
 *
 * @return log(det(L)^2)
 *
 */
double chm_factor_ldetL2(CHM_FR f)
{
    int i, j, p;
    double ans = 0;

    if (f->is_super) {
	int *lpi = (int*)(f->pi), *lsup = (int*)(f->super);
	for (i = 0; i < f->nsuper; i++) { /* supernodal block i */
	    int nrp1 = 1 + lpi[i + 1] - lpi[i],
		nc = lsup[i + 1] - lsup[i];
	    double *x = (double*)(f->x) + ((int*)(f->px))[i];

	    for (R_xlen_t jn = 0, j = 0; j < nc; j++, jn += nrp1) { // jn := j * nrp1
		ans += 2 * log(fabs(x[jn]));
	    }
	}
    } else {
	int *li = (int*)(f->i), *lp = (int*)(f->p);
	double *lx = (double *)(f->x);

	for (j = 0; j < f->n; j++) {
	    for (p = lp[j]; li[p] != j && p < lp[j + 1]; p++) {};
	    if (li[p] != j) {
		error(_("diagonal element %d of Cholesky factor is missing"), j);
		break;		/* -Wall */
	    }
	    ans += log(lx[p] * ((f->is_ll) ? lx[p] : 1.));
	}
    }
    return ans;
}

/**
 * Update the numerical values in the factor f as A + mult * I, if A is
 * symmetric, otherwise AA' + mult * I
 *
 * @param f pointer to a CHM_FR object.  f is updated upon return.
 * @param A pointer to a CHM_SP object, possibly symmetric
 * @param mult multiple of the identity to be added to A or AA' before
 * decomposing.
 *
 * @note: A and f must be compatible.  There is no check on this
 * here.  Incompatibility of A and f will cause the CHOLMOD functions
 * to take an error exit.
 *
 */
CHM_FR chm_factor_update(CHM_FR f, CHM_SP A, double mult)
{
    int ll = f->is_ll;
    double mm[2] = {0, 0};
    mm[0] = mult;
    // NB: Result depends if A is "dsC" or "dgC"; the latter case assumes we mean AA' !!!
    if (!cholmod_factorize_p(A, mm, (int*)NULL, 0 /*fsize*/, f, &c))
	/* -> ./CHOLMOD/Cholesky/cholmod_factorize.c */
	error(_("cholmod_factorize_p failed: status %d, minor %d of ncol %d"),
	      c.status, f->minor, f->n);
    if (f->is_ll != ll)
	if(!cholmod_change_factor(f->xtype, ll, f->is_super, 1 /*to_packed*/,
				  1 /*to_monotonic*/, f, &c))
	    error(_("cholmod_change_factor failed"));
    return f;
}
