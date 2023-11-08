#include "Mdefines.h"
#include "sparseVector.h"

SEXP v2spV(SEXP from)
{
	SEXP to = NULL, length = NULL, i = NULL, x = NULL;
	R_xlen_t n_ = XLENGTH(from);

#define V2SPV(_KIND_, _NZ_, \
	          _CTYPE1_, _SEXPTYPE1_, _PTR1_, \
	          _CTYPE2_, _SEXPTYPE2_, _PTR2_) \
	do { \
		PROTECT(to = newObject(#_KIND_ "sparseVector")); \
		_CTYPE1_ *py = _PTR1_(from); \
		for (k = 0; k < n; ++k) \
			if (_NZ_(py[k])) \
				++nnz; \
		PROTECT(i = allocVector(_SEXPTYPE2_, nnz)); \
		PROTECT(x = allocVector(_SEXPTYPE1_, nnz)); \
		_CTYPE2_ *pi = _PTR2_(i); \
		_CTYPE1_ *px = _PTR1_(x); \
		for (k = 0; k < n; ++k) { \
			if (_NZ_(py[k])) { \
				*(pi++) = (_CTYPE2_) (k + 1); \
				*(px++) = py[k]; \
			} \
		} \
	} while (0)

#define V2SPV_CASES(_CTYPE2_, _SEXPTYPE2_, _PTR2_) \
	do { \
		switch (TYPEOF(from)) { \
		case LGLSXP: \
			V2SPV(l, ISNZ_LOGICAL, int, LGLSXP, LOGICAL, \
			      _CTYPE2_, _SEXPTYPE2_, _PTR2_); \
			break; \
		case INTSXP: \
			V2SPV(i, ISNZ_INTEGER, int, INTSXP, INTEGER, \
			      _CTYPE2_, _SEXPTYPE2_, _PTR2_); \
			break; \
		case REALSXP: \
			V2SPV(d, ISNZ_REAL, double, REALSXP, REAL, \
			      _CTYPE2_, _SEXPTYPE2_, _PTR2_); \
			break; \
		case CPLXSXP: \
			V2SPV(z, ISNZ_COMPLEX, Rcomplex, CPLXSXP, COMPLEX, \
			      _CTYPE2_, _SEXPTYPE2_, _PTR2_); \
			break; \
		default: \
			ERROR_INVALID_TYPE(from, __func__); \
			break; \
		} \
	} while (0)

	if (n_ <= INT_MAX) {
		int k, n = (int) n_, nnz = 0;
		PROTECT(length = ScalarInteger(n));
		V2SPV_CASES(int, INTSXP, INTEGER);
	} else {
		R_xlen_t k, n = n_, nnz = 0;
		PROTECT(length = ScalarReal((double) n));
		V2SPV_CASES(double, REALSXP, REAL);
	}

#undef V2SPV_CASES
#undef V2SPV

	SET_SLOT(to, Matrix_lengthSym, length);
	SET_SLOT(to, Matrix_iSym, i);
	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(4); /* x, i, length, to */
	return to;
}

SEXP CR2spV(SEXP from)
{
	static const char *valid[] = { VALID_CSPARSE, VALID_RSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid];

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	Matrix_int_fast64_t mn = (Matrix_int_fast64_t) m * n;
	UNPROTECT(1); /* dim */

	if (mn > 0x1.0p+53)
		error(_("%s length cannot exceed %s"), "sparseVector", "2^53");

	/* defined in ./coerce.c : */
	SEXP sparse_as_general(SEXP, const char *);
	PROTECT(from = sparse_as_general(from, cl));

	char vcl[] = ".sparseVector";
	vcl[0] = cl[0];
	SEXP to = PROTECT(newObject(vcl));

	SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
	int *pp = INTEGER(p), nnz = (cl[2] == 'C') ? pp[n] : pp[m];

	SEXP vlength, vi;
	if (mn <= INT_MAX) {
		PROTECT(vlength = ScalarInteger(m * n));
		PROTECT(vi = allocVector(INTSXP, nnz));
	} else {
		PROTECT(vlength = ScalarReal((double) m * n));
		PROTECT(vi = allocVector(REALSXP, nnz));
	}
	SET_SLOT(to, Matrix_lengthSym, vlength);
	SET_SLOT(to, Matrix_iSym, vi);

	if (cl[2] == 'C') {
		SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
		int *pi = INTEGER(i), k, kend, j;
		if (TYPEOF(vi) == INTSXP) {
			int *pvi = INTEGER(vi), mj1a = 1;
			for (j = 0, k = 0; j < n; ++j) {
				kend = *(++pp);
				while (k < kend) {
					*(pvi++) = mj1a + *(pi++);
					++k;
				}
				mj1a += m;
			}
		} else {
			double *pvi = REAL(vi), mj1a = 1.0, m_ = (double) m;
			for (j = 0, k = 0; j < n; ++j) {
				kend = *(++pp);
				while (k < kend) {
					*(pvi++) = mj1a + (double) *(pi++);
					++k;
				}
				mj1a += m_;
			}
		}
		if (cl[0] != 'n') {
			SEXP vx = PROTECT(GET_SLOT(from, Matrix_xSym));
			SET_SLOT(to, Matrix_xSym, vx);
			UNPROTECT(1); /* vx */
		}
		UNPROTECT(1); /* i */
	} else {
		SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pj = INTEGER(j), k, kend, tmp, *work;
		Matrix_Calloc(work, n, int);
		for (k = 0; k < nnz; ++k)
			work[pj[k]]++;
		nnz = 0;
		for (k = 0; k < n; ++k) {
			tmp = work[k];
			work[k] = nnz;
			nnz += tmp;
		}

#define R2SPV(_CTYPE_, _PTR_, _MASK_) \
		do { \
			_MASK_(_CTYPE_ *px  = _PTR_(x )); \
			_MASK_(_CTYPE_ *pvx = _PTR_(vx)); \
			if (TYPEOF(vi) == INTSXP) { \
				int *pvi = INTEGER(vi), i; \
				k = 0; \
				for (i = 1; i <= m; i += 1) { \
					kend = *(++pp); \
					while (k < kend) { \
						pvi[work[pj[k]]] = m * pj[k] + i; \
						_MASK_(pvx[work[pj[k]]] = px[k]); \
						work[pj[k++]]++; \
					} \
				} \
			} else { \
				double *pvi = REAL(vi), i_, m_ = (double) m; \
				k = 0; \
				for (i_ = 1.0; i_ <= m_; i_ += 1.0) { \
					kend = *(++pp); \
					while (k < kend) { \
						pvi[work[pj[k]]] = m_ * pj[k] + i_; \
						_MASK_(pvx[work[pj[k]]] = px[k]); \
						work[pj[k++]]++; \
					} \
				} \
			} \
		} while (0)

		if (cl[0] == 'n')
			R2SPV(, , HIDE);
		else {
			SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym)),
				vx = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
			switch (TYPEOF(x)) {
			case LGLSXP:
				R2SPV(int, LOGICAL, SHOW);
				break;
			case INTSXP:
				R2SPV(int, INTEGER, SHOW);
				break;
			case REALSXP:
				R2SPV(double, REAL, SHOW);
				break;
			case CPLXSXP:
				R2SPV(Rcomplex, COMPLEX, SHOW);
				break;
			default:
				break;
			}
			SET_SLOT(to, Matrix_xSym, vx);
			UNPROTECT(2); /* vx, x */
		}

#undef R2SV

		UNPROTECT(1); /* j */
		Matrix_Free(work, n);
	}

	UNPROTECT(5); /* vi, vlength, p, to, from */
	return to;
}
