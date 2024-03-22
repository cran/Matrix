#include "Lapack-etc.h"
#include "Mdefines.h"
#include "expm.h"

/* For matrix exponential calculation : */
const static double padec [] =
{
  5.0000000000000000e-1,
  1.1666666666666667e-1,
  1.6666666666666667e-2,
  1.6025641025641026e-3,
  1.0683760683760684e-4,
  4.8562548562548563e-6,
  1.3875013875013875e-7,
  1.9270852604185938e-9,
};

/* Based on _corrected_ code for Octave function 'expm' : */
SEXP dgeMatrix_expm(SEXP x)
{
    const double one = 1.0, zero = 0.0;
    const int i1 = 1;
    int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    const int n = Dims[1];
    const R_xlen_t n_ = n, np1 = n + 1, nsqr = n_ * n_; // nsqr = n^2

    SEXP val = PROTECT(duplicate(x));
    int i, ilo, ilos, ihi, ihis, j, sqpow;
    int *pivot = R_Calloc(n, int);
    double *dpp = R_Calloc(nsqr, double), /* denominator power Pade' */
	*npp = R_Calloc(nsqr, double), /* numerator power Pade' */
	*perm = R_Calloc(n, double),
	*scale = R_Calloc(n, double),
	*v = REAL(GET_SLOT(val, Matrix_xSym)),
	*work = R_Calloc(nsqr, double), inf_norm, m1_j/*= (-1)^j */, trshift;
    R_CheckStack();

    if (n < 1 || Dims[0] != n)
	error(_("Matrix exponential requires square, non-null matrix"));
    if(n == 1) {
	v[0] = exp(v[0]);
	UNPROTECT(1);
	return val;
    }

    /* Preconditioning 1.  Shift diagonal by average diagonal if positive. */
    trshift = 0;		/* determine average diagonal element */
    for (i = 0; i < n; i++) trshift += v[i * np1];
    trshift /= n;
    if (trshift > 0.) {		/* shift diagonal by -trshift */
	for (i = 0; i < n; i++) v[i * np1] -= trshift;
    }

    /* Preconditioning 2. Balancing with dgebal. */
    F77_CALL(dgebal)("P", &n, v, &n, &ilo, &ihi, perm, &j FCONE);
    if (j) error(_("dgeMatrix_exp: LAPACK routine dgebal returned %d"), j);
    F77_CALL(dgebal)("S", &n, v, &n, &ilos, &ihis, scale, &j FCONE);
    if (j) error(_("dgeMatrix_exp: LAPACK routine dgebal returned %d"), j);

    /* Preconditioning 3. Scaling according to infinity norm */
    inf_norm = F77_CALL(dlange)("I", &n, &n, v, &n, work FCONE);
    sqpow = (inf_norm > 0) ? (int) (1 + log(inf_norm)/log(2.)) : 0;
    if (sqpow < 0) sqpow = 0;
    if (sqpow > 0) {
	double scale_factor = 1.0;
	for (i = 0; i < sqpow; i++) scale_factor *= 2.;
	for (R_xlen_t i = 0; i < nsqr; i++) v[i] /= scale_factor;
    }

    /* Pade' approximation. Powers v^8, v^7, ..., v^1 */
    Matrix_memset(npp, 0, nsqr, sizeof(double));
    Matrix_memset(dpp, 0, nsqr, sizeof(double));
    m1_j = -1;
    for (j = 7; j >=0; j--) {
	double mult = padec[j];
	/* npp = m * npp + padec[j] *m */
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, npp, &n,
			&zero, work, &n FCONE FCONE);
	for (R_xlen_t i = 0; i < nsqr; i++) npp[i] = work[i] + mult * v[i];
	/* dpp = m * dpp + (m1_j * padec[j]) * m */
	mult *= m1_j;
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, dpp, &n,
			&zero, work, &n FCONE FCONE);
	for (R_xlen_t i = 0; i < nsqr; i++) dpp[i] = work[i] + mult * v[i];
	m1_j *= -1;
    }
    /* Zero power */
    for (R_xlen_t i = 0; i < nsqr; i++) dpp[i] *= -1.;
    for (j = 0; j < n; j++) {
	npp[j * np1] += 1.;
	dpp[j * np1] += 1.;
    }

    /* Pade' approximation is solve(dpp, npp) */
    F77_CALL(dgetrf)(&n, &n, dpp, &n, pivot, &j);
    if (j) error(_("dgeMatrix_exp: dgetrf returned error code %d"), j);
    F77_CALL(dgetrs)("N", &n, &n, dpp, &n, pivot, npp, &n, &j FCONE);
    if (j) error(_("dgeMatrix_exp: dgetrs returned error code %d"), j);
    Memcpy(v, npp, nsqr);

    /* Now undo all of the preconditioning */
    /* Preconditioning 3: square the result for every power of 2 */
    while (sqpow--) {
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, v, &n,
			&zero, work, &n FCONE FCONE);
	Memcpy(v, work, nsqr);
    }
    /* Preconditioning 2: apply inverse scaling */
    for (j = 0; j < n; j++) {
	R_xlen_t jn = j * n_;
	for (i = 0; i < n; i++)
	    v[i + jn] *= scale[i]/scale[j];
    }


    /* 2 b) Inverse permutation  (if not the identity permutation) */
    if (ilo != 1 || ihi != n) {
	/* Martin Maechler's code */

#define SWAP_ROW(I,J) F77_CALL(dswap)(&n, &v[(I)], &n, &v[(J)], &n)

#define SWAP_COL(I,J) F77_CALL(dswap)(&n, &v[(I)*n_], &i1, &v[(J)*n_], &i1)

#define RE_PERMUTE(I)				\
	int p_I = (int) (perm[I]) - 1;		\
	SWAP_COL(I, p_I);			\
	SWAP_ROW(I, p_I)

	/* reversion of "leading permutations" : in reverse order */
	for (i = (ilo - 1) - 1; i >= 0; i--) {
	    RE_PERMUTE(i);
	}

	/* reversion of "trailing permutations" : applied in forward order */
	for (i = (ihi + 1) - 1; i < n; i++) {
	    RE_PERMUTE(i);
	}
    }

    /* Preconditioning 1: Trace normalization */
    if (trshift > 0.) {
	double mult = exp(trshift);
	for (R_xlen_t i = 0; i < nsqr; i++) v[i] *= mult;
    }

    /* Clean up */
    R_Free(work); R_Free(scale); R_Free(perm); R_Free(npp); R_Free(dpp); R_Free(pivot);
    UNPROTECT(1);
    return val;
}
