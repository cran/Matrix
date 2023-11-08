#include <Rmath.h> /* math.h, logspace_add, logspace_sub */
#include "Mdefines.h"
#include "cholmod-etc.h"
#include "determinant.h"

static
SEXP mkDet(double modulus, int logarithm, int sign)
{
	SEXP nms = PROTECT(allocVector(STRSXP, 2)),
		cl = PROTECT(mkString("det")),
		det = PROTECT(allocVector(VECSXP, 2)),
		det0 = PROTECT(ScalarReal((logarithm) ? modulus : exp(modulus))),
		det1 = PROTECT(ScalarInteger(sign)),
		det0a = PROTECT(ScalarLogical(logarithm));
	SET_STRING_ELT(nms, 0, mkChar("modulus"));
	SET_STRING_ELT(nms, 1, mkChar("sign"));
	setAttrib(det, R_NamesSymbol, nms);
	setAttrib(det, R_ClassSymbol, cl);
	setAttrib(det0, install("logarithm"), det0a);
	SET_VECTOR_ELT(det, 0, det0);
	SET_VECTOR_ELT(det, 1, det1);
	UNPROTECT(6);
	return det;
}

SEXP denseLU_determinant(SEXP obj, SEXP logarithm)
{

#define DETERMINANT_START \
	SEXP dim = GET_SLOT(obj, Matrix_DimSym); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	if (m != n) \
		error(_("determinant of non-square matrix is undefined")); \
	int givelog = asLogical(logarithm) != 0; \
	double modulus = 0.0; /* result for n == 0 */

	DETERMINANT_START;

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int sign = (TYPEOF(x) == CPLXSXP) ? NA_INTEGER : 1;

	if (n > 0) {
	int j;
	R_xlen_t n1a = (R_xlen_t) n + 1;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x);
		for (j = 0; j < n; ++j) {
			modulus += log(hypot((*px).r, (*px).i));
			px += n1a;
		}
	} else {
		SEXP pivot = GET_SLOT(obj, Matrix_permSym);
		int *ppivot = INTEGER(pivot);
			double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			if (ISNAN(*px) || *px >= 0.0) {
				modulus += log(*px);
				if (*ppivot != j + 1)
					sign = -sign;
			} else {
				modulus += log(-(*px));
				if (*ppivot == j + 1)
					sign = -sign;
			}
			px += n1a;
			ppivot += 1;
		}
	}
	}

	UNPROTECT(1); /* x */
	return mkDet(modulus, givelog, sign);
}

SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm)
{
	DETERMINANT_START;

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int sign = (TYPEOF(x) == CPLXSXP) ? NA_INTEGER : 1;

	if (n > 0) {
	SEXP pivot = GET_SLOT(obj, Matrix_permSym);
	int *ppivot = INTEGER(pivot);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = *CHAR(STRING_ELT(uplo, 0));

	int j = 0, unpacked = (Matrix_int_fast64_t) n * n <= R_XLEN_T_MAX &&
		XLENGTH(x) == (R_xlen_t) m * m;
	R_xlen_t n1a = (R_xlen_t) n + 1;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x), a, b, c;
		while (j < n) {
			if (ppivot[j] > 0) {
				modulus += log(hypot((*px).r, (*px).i));
				px += (unpacked) ? n1a : ((ul == 'U') ? j + 2 : n - j);
				j += 1;
			} else {
				a = *px;
				if (ul == 'U') {
					px += (unpacked) ? n1a : j + 2;
					b = *px;
					c = *(px - 1);
					px += (unpacked) ? n1a : j + 3;
				} else {
					c = *(px + 1);
					px += (unpacked) ? n1a : n - j;
					b = *px;
					px += (unpacked) ? n1a : n - j - 1;
				}
				modulus += log(hypot(a.r * b.r - a.i * b.i -
									 c.r * c.r + c.i * c.i,
									 a.r * b.i + a.i * b.r -
									 2.0 * c.r * c.i));
				j += 2;
			}
		}
	} else {
		double *px = REAL(x), a, b, c, logab, logcc;
		while (j < n) {
			if (ppivot[j] > 0) {
				if (*px >= 0.0)
					modulus += log(*px);
				else {
					modulus += log(-(*px));
					sign = -sign;
				}
				px += (unpacked) ? n1a : ((ul == 'U') ? j + 2 : n - j);
				j += 1;
			} else {
				a = *px;
				if (ul == 'U') {
					px += (unpacked) ? n1a : j + 2;
					b = *px;
					c = *(px - 1);
					px += (unpacked) ? n1a : j + 3;
				} else {
					c = *(px + 1);
					px += (unpacked) ? n1a : n - j;
					b = *px;
					px += (unpacked) ? n1a : n - j - 1;
				}
				logab = log(fabs(a)) + log(fabs(b));
				logcc = 2.0 * log(fabs(c));
				if ((a < 0.0) != (b < 0.0)) {
					/* det = ab - cc = -(abs(ab) + cc) < 0 */
					modulus += logspace_add(logab, logcc);
					sign = -sign;
				} else if (logab < logcc) {
					/* det = ab - cc = -(cc - ab) < 0 */
					modulus += logspace_sub(logcc, logab);
					sign = -sign;
				} else {
					/* det = ab - cc > 0 */
					modulus += logspace_sub(logab, logcc);
				}
				j += 2;
			}
		}
	}
	}

	UNPROTECT(1); /* x */
	return mkDet(modulus, givelog, sign);
}

SEXP Cholesky_determinant(SEXP obj, SEXP logarithm)
{
	DETERMINANT_START;

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int sign = (TYPEOF(x) == CPLXSXP) ? NA_INTEGER : 1;

	if (n > 0) {
	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = *CHAR(STRING_ELT(uplo, 0));

	int j, unpacked = (Matrix_int_fast64_t) n * n <= R_XLEN_T_MAX &&
		XLENGTH(x) == (R_xlen_t) m * m;
	R_xlen_t n1a = (R_xlen_t) n + 1;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x);
		for (j = 0; j < n; ++j) {
			modulus += log(hypot((*px).r, (*px).i));
			px += (unpacked) ? n1a : ((ul == 'U') ? j + 2 : n - j);
		}
	} else {
		double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			if (ISNAN(*px) || *px >= 0.0)
				modulus += log(*px);
			else {
				modulus += log(-(*px));
				sign = -sign;
			}
			px += (unpacked) ? n1a : ((ul == 'U') ? j + 2 : n - j);
		}
	}
	modulus *= 2.0;
	}

	UNPROTECT(1); /* x */
	return mkDet(modulus, givelog, sign);
}

SEXP sparseLU_determinant(SEXP obj, SEXP logarithm)
{
	DETERMINANT_START;

	SEXP U = PROTECT(GET_SLOT(obj, Matrix_USym)),
		x = PROTECT(GET_SLOT(U, Matrix_xSym));
	int sign = (TYPEOF(x) == CPLXSXP) ? NA_INTEGER : 1;

	if (n > 0) {
	SEXP p = PROTECT(GET_SLOT(U, Matrix_pSym)),
		i = PROTECT(GET_SLOT(U, Matrix_iSym));
	int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k = 0, kend;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x);
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && pi[kend - 1] == j)
				modulus += log(hypot(px[kend - 1].r, px[kend - 1].i));
			else {
				UNPROTECT(4); /* i, p, x, U */
				return mkDet(R_NegInf, givelog, 1);
			}
			k = kend;
		}
	} else {
		double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && pi[kend - 1] == j) {
				if (ISNAN(px[kend - 1]) || px[kend - 1] >= 0.0)
					modulus += log(px[kend - 1]);
				else {
					modulus += log(-px[kend - 1]);
					sign = -sign;
				}
			} else {
				UNPROTECT(4); /* i, p, x, U */
				return mkDet(R_NegInf, givelog, 1);
			}
			k = kend;
		}

		/* defined in ./perm.c : */
		int signPerm(const int *, int, int);

		p = GET_SLOT(obj, Matrix_pSym);
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		p = GET_SLOT(obj, Matrix_qSym);
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
	}
	UNPROTECT(2); /* i, p */
	}

	UNPROTECT(2); /* x, U */
	return mkDet(modulus, givelog, sign);
}

SEXP sparseQR_determinant(SEXP obj, SEXP logarithm)
{
	DETERMINANT_START;

	SEXP R = PROTECT(GET_SLOT(obj, Matrix_RSym)),
		x = PROTECT(GET_SLOT(R, Matrix_xSym));
	int sign = (TYPEOF(x) == CPLXSXP) ? NA_INTEGER : 1;

	dim = GET_SLOT(R, Matrix_DimSym);
	if (INTEGER(dim)[0] > n)
		error(_("%s(<%s>) does not support structurally rank deficient case"),
			  "determinant", "sparseQR");

	if (n > 0) {
	SEXP p = PROTECT(GET_SLOT(R, Matrix_pSym)),
		i = PROTECT(GET_SLOT(R, Matrix_iSym));
	int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k = 0, kend;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x);
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && pi[kend - 1] == j)
				modulus += log(hypot(px[kend - 1].r, px[kend - 1].i));
			else {
				UNPROTECT(4); /* i, p, x, U */
				return mkDet(R_NegInf, givelog, 1);
			}
			k = kend;
		}
	} else {
		double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && pi[kend - 1] == j) {
				if (ISNAN(px[kend - 1]) || px[kend - 1] >= 0.0)
					modulus += log(px[kend - 1]);
				else {
					modulus += log(-px[kend - 1]);
					sign = -sign;
				}
			} else {
				UNPROTECT(4); /* i, p, x, R */
				return mkDet(R_NegInf, givelog, 1);
			}
			k = kend;
		}

		/* defined in ./perm.c : */
		int signPerm(const int *, int, int);

		p = GET_SLOT(obj, Matrix_pSym);
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		p = GET_SLOT(obj, Matrix_qSym);
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		if (n % 2)
			sign = -sign;
	}
	UNPROTECT(2); /* i, p */
	}

	UNPROTECT(2); /* x, R */
	return mkDet(modulus, givelog, sign);
}

SEXP CHMfactor_determinant(SEXP obj, SEXP logarithm, SEXP sqrt)
{
	DETERMINANT_START;

	cholmod_factor *L = M2CHF(obj, 1);
	int sign = (L->xtype == CHOLMOD_COMPLEX) ? NA_INTEGER : 1;

	if (n > 0) {
	int j, sqrt_ = asLogical(sqrt);
	if (L->is_super) {
		int k, nc,
			nsuper = (int) L->nsuper,
			*psuper = (int *) L->super,
			*ppi = (int *) L->pi,
			*ppx = (int *) L->px;
		R_xlen_t nr1a;
		if (L->xtype == CHOLMOD_COMPLEX) {
			Rcomplex *px = (Rcomplex *) L->x, *px_;
			for (k = 0; k < nsuper; ++k) {
				nc = psuper[k + 1] - psuper[k];
				nr1a = (R_xlen_t) (ppi[k + 1] - ppi[k]) + 1;
				px_ = px + ppx[k];
				for (j = 0; j < nc; ++j) {
					modulus += log(hypot((*px_).r, (*px_).i));
					px_ += nr1a;
				}
			}
		} else {
			double *px = (double *) L->x, *px_;
			for (k = 0; k < nsuper; ++k) {
				nc = psuper[k + 1] - psuper[k];
				nr1a = (R_xlen_t) (ppi[k + 1] - ppi[k]) + 1;
				px_ = px + ppx[k];
				for (j = 0; j < nc; ++j) {
					modulus += log(*px_);
					px_ += nr1a;
				}
			}
		}
		modulus *= 2.0;
	} else {
		int *pp = (int *) L->p;
		if (L->xtype == CHOLMOD_COMPLEX) {
			Rcomplex *px = (Rcomplex *) L->x;
			for (j = 0; j < n; ++j)
				modulus += log(hypot(px[pp[j]].r, px[pp[j]].i));
			if (L->is_ll)
				modulus *= 2.0;
		} else {
			double *px = (double *) L->x;
			if (L->is_ll) {
				for (j = 0; j < n; ++j)
					modulus += log(px[pp[j]]);
				modulus *= 2.0;
			} else {
				for (j = 0; j < n; ++j) {
					if (ISNAN(px[pp[j]]) || px[pp[j]] >= 0.0)
						modulus += log(px[pp[j]]);
					else {
						if (sqrt_)
							return mkDet(R_NaN, givelog, 1);
						modulus += log(-px[pp[j]]);
						sign = -sign;
					}
				}
			}
		}
	}
	if (sqrt_)
		modulus *= 0.5;
	}

	return mkDet(modulus, givelog, sign);
}
