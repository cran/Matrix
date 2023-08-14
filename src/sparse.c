#include <math.h> /* fabs, hypot */
#include "sparse.h"

SEXP sparse_drop0(SEXP from, const char *class, double tol)
{
	if (class[0] == 'n')
		return from;

	SEXP to, x0 = PROTECT(GET_SLOT(from, Matrix_xSym));

#define TOLBASED_ISNZ_REAL(_X_) \
	(ISNA_REAL(_X_) || fabs(_X_) > tol)

#define TOLBASED_ISNZ_COMPLEX(_X_) \
	(ISNA_COMPLEX(_X_) || hypot((_X_).r, (_X_).i) > tol)

#define DROP0_CASES(_DO_) \
	do { \
		switch (class[0]) { \
		case 'l': \
			_DO_(int, LOGICAL, ISNZ_LOGICAL); \
			break; \
		case 'i': \
			_DO_(int, INTEGER, ISNZ_INTEGER); \
			break; \
		case 'd': \
			if (tol > 0.0) \
				_DO_(double, REAL, TOLBASED_ISNZ_REAL); \
			else \
				_DO_(double, REAL,          ISNZ_REAL); \
			break; \
		case 'z': \
			if (tol > 0.0) \
				_DO_(Rcomplex, COMPLEX, TOLBASED_ISNZ_COMPLEX); \
			else \
				_DO_(Rcomplex, COMPLEX,          ISNZ_COMPLEX); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	if (class[2] != 'T') {

		SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym));
		int *pp0 = INTEGER(p0), k, n = (int) (XLENGTH(p0) - 1),
			nnz0 = pp0[n], nnz1 = 0;

#undef DROP0_LOOP1
#define DROP0_LOOP1(_CTYPE_, _PTR_, _NZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0); \
			for (k = 0; k < nnz0; ++k) { \
				if (_NZ_(*px0)) \
					++nnz1; \
				++px0; \
			} \
		} while (0)

		DROP0_CASES(DROP0_LOOP1);
		if (nnz1 == nnz0) {
			UNPROTECT(2); /* p0, x0 */
			return from;
		}
		PROTECT(to = NEW_OBJECT_OF_CLASS(class));

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			i0 = PROTECT(GET_SLOT(from, iSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
			i1 = PROTECT(allocVector(INTSXP, nnz1)),
			x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
		int *pi0 = INTEGER(i0), *pp1 = INTEGER(p1), *pi1 = INTEGER(i1),
			j, kend;
		pp0++; *(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to,        iSym, i1);
		SET_SLOT(to, Matrix_xSym, x1);

#undef DROP0_LOOP2
#define DROP0_LOOP2(_CTYPE_, _PTR_, _NZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			for (j = 0, k = 0; j < n; ++j) { \
				pp1[j] = pp1[j - 1]; \
				kend = pp0[j]; \
				while (k < kend) { \
					if (_NZ_(*px0)) { \
						++pp1[j]; \
						*(pi1++) = *pi0; \
						*(px1++) = *px0; \
					} \
					++k; ++pi0; ++px0; \
				} \
			} \
		} while (0)

		DROP0_CASES(DROP0_LOOP2);
		UNPROTECT(7); /* x1, i1, p1, i0, to, p0, x0 */

	} else {

		R_xlen_t k, nnz0 = XLENGTH(x0), nnz1 = 0;

		DROP0_CASES(DROP0_LOOP1);
		if (nnz1 == nnz0) {
			UNPROTECT(1); /* x0 */
			return from;
		}
		PROTECT(to = NEW_OBJECT_OF_CLASS(class));

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1)),
			x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0),
			*pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		SET_SLOT(to, Matrix_xSym, x1);

#undef DROP0_LOOP2
#define DROP0_LOOP2(_CTYPE_, _PTR_, _NZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			for (k = 0; k < nnz0; ++k) { \
				if (_NZ_(*px0)) { \
					*(pi1++) = *pi0; \
					*(pj1++) = *pj0; \
					*(px1++) = *px0; \
				} \
				++pi0; ++pj0; ++px0; \
			} \
		} while (0)

		DROP0_CASES(DROP0_LOOP2);
		UNPROTECT(7); /* x1, j1, i1, j0, i0, to, x0 */

	}

	PROTECT(to);

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorSym, factors);
		UNPROTECT(1); /* factors */
	}

#undef TOLBASED_ISNZ_REAL
#undef TOLBASED_ISNZ_COMPLEX
#undef DROP0_CASES
#undef DROP0_LOOP1
#undef DROP0_LOOP2

	UNPROTECT(1); /* to */
	return to;
}

/* drop0(<[CRT]sparseMatrix>, tol) */
SEXP R_sparse_drop0(SEXP from, SEXP tol)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	double tol_;
	if (TYPEOF(tol) != REALSXP || LENGTH(tol) < 1 ||
	    ISNAN(tol_ = REAL(tol)[0]))
		error(_("'%s' is not a number"), "tol");

	return sparse_drop0(from, valid[ivalid], tol_);
}

SEXP sparse_band(SEXP from, const char *class, int a, int b)
{
	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */

	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
	if (a <= 1-m && b >= n-1 && (class[1] == 't' || m != n || m > 1 || n > 1))
		return from;

	int ge = 0, sy = 0, tr = 0;
	ge = m != n || !((tr = a >= 0 || b <= 0 || class[1] == 't') ||
	                 (sy = a == -b && class[1] == 's'));

	char ulf = 'U', ult = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ulf = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		if (class[1] == 't') {
			/* Be fast if band contains entire triangle */
			if ((ulf == 'U') ? (a <= 0 && b >= n-1) : (b >= 0 && a <= 1-m))
				return from;
			else if (a <= 0 && b >= 0) {
				SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
				di = *CHAR(STRING_ELT(diag, 0));
				UNPROTECT(1); /* diag */
			}
		}
	}

	/* band(<R>, a, b) is equivalent to t(band(t(<R>), -b, -a)) ! */

	if (class[2] == 'R') {
		int r;
		r = m; m =  n; n =  r;
		r = a; a = -b; b = -r;
		ulf = (ulf == 'U') ? 'L' : 'U';
		from = sparse_transpose(from, class, 1);
	}
	PROTECT(from);

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (ge) ? 'g' : ((tr) ? 't' : 's');
	cl[2] = (class[2] == 'R') ? 'C' : class[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	dim = GET_SLOT(to, Matrix_DimSym);
	pdim = INTEGER(dim);
	pdim[0] = m;
	pdim[1] = n;

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] != 's' || sy)
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (!ge) {
		ult = (tr && class[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ulf;
		if (ult != 'U') {
			SEXP uplo = PROTECT(mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
		if (di != 'N') {
			SEXP diag = PROTECT(mkString("U"));
			SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}
	}

	/* It remains to set some subset of 'p', 'i', 'j', 'x' ... */

#define BAND_CASES \
	do { \
		switch (class[0]) { \
		case 'l': \
			BAND_SUBCASES(int, LOGICAL, SHOW); \
			break; \
		case 'i': \
			BAND_SUBCASES(int, INTEGER, SHOW); \
			break; \
		case 'd': \
			BAND_SUBCASES(double, REAL, SHOW); \
			break; \
		case 'z': \
			BAND_SUBCASES(Rcomplex, COMPLEX, SHOW); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	if (class[2] != 'T') {
		SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), *pp1 = INTEGER(p1),
			d, j, k, kend, nnz0 = pp0[n], nnz1 = 0;
		pp0++; *(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

		if (class[1] == 's' && !sy) {
			Matrix_memset(pp1, 0, n, sizeof(int));
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((d = j - pi0[k]) >= a && d <= b)
						++pp1[j];
					if (d != 0 && -d >= a && -d <= b)
						++pp1[pi0[k]];
					++k;
				}
			}
			for (j = 0; j < n; ++j) {
				nnz1 += pp1[j];
				pp1[j] = nnz1;
			}
		} else {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((d = j - pi0[k]) >= a && d <= b)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
		}

		if (nnz1 == nnz0 && (class[1] != 's' || sy)) {
			/* No need to allocate in this case: band has all nonzero elements */
			SET_SLOT(to, Matrix_iSym, i0);
			if (class[0] != 'n') {
				SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
				SET_SLOT(to, Matrix_xSym, x0);
				UNPROTECT(1); /* x0 */
			}
			if (class[2] == 'R')
				to = sparse_transpose(to, cl, 1);
			UNPROTECT(5); /* p1, i0, p0, to, from */
			return to;
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, Matrix_iSym, i1);

#undef BAND_SUBCASES
#define BAND_SUBCASES(_CTYPE_, _PTR_, _MASK_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
			if (class[1] == 's' && !sy) { \
				int *pp1_; \
				Matrix_Calloc(pp1_, n, int); \
				Matrix_memcpy(pp1_, pp1 - 1, n, sizeof(int)); \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((d = j - pi0[k]) >= a && d <= b) { \
							pi1[pp1_[j]] = pi0[k]; \
							_MASK_(px1[pp1_[j]] = px0[k]); \
							++pp1_[j]; \
						} \
						if (d != 0 && -d >= a && -d <= b) { \
							pi1[pp1_[pi0[k]]] = j; \
							_MASK_(px1[pp1_[pi0[k]]] = px0[k]); \
							++pp1_[pi0[k]]; \
						} \
						++k; \
					} \
				} \
				Matrix_Free(pp1_, n); \
			} else { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((d = j - pi0[k]) >= a && d <= b) { \
							*(pi1++) = pi0[k]; \
							_MASK_(*(px1++) = px0[k]); \
						} \
						++k; \
					} \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			BAND_SUBCASES(int, LOGICAL, HIDE);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			SET_SLOT(to, Matrix_xSym, x1);
			BAND_CASES;
			UNPROTECT(2); /* x1, x0 */
		}
		if (class[2] == 'R')
			to = sparse_transpose(to, cl, 1);

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), d;
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = 0;

		if (class[1] == 's' && !sy) {
			for (k = 0; k < nnz0; ++k) {
				if ((d = pj0[k] - pi0[k]) >= a && d <= b)
					++nnz1;
				if (d != 0 && -d >= a && -d <= b)
					++nnz1;
			}
		} else {
			for (k = 0; k < nnz0; ++k) {
				if ((d = pj0[k] - pi0[k]) >= a && d <= b)
					++nnz1;
			}
		}

		if (nnz1 == nnz0 && (class[1] != 's' || sy)) {
			/* No need to allocate in this case: band has all nonzero elements */
			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);
			if (class[0] != 'n') {
				SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
				SET_SLOT(to, Matrix_xSym, x0);
				UNPROTECT(1); /* x0 */
			}
			UNPROTECT(4); /* j0, i0, to, from */
			return to;
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);

#undef BAND_SUBCASES
#define BAND_SUBCASES(_CTYPE_, _PTR_, _MASK_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
			if (class[1] == 's' && !sy) { \
				for (k = 0; k < nnz0; ++k) { \
					if ((d = pj0[k] - pi0[k]) >= a && d <= b) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						_MASK_(*(px1++) = px0[k]); \
					} \
					if (d != 0 && -d >= a && -d <= b) { \
						*(pi1++) = pj0[k]; \
						*(pj1++) = pi0[k]; \
						_MASK_(*(px1++) = px0[k]); \
					} \
				} \
			} else { \
				for (k = 0; k < nnz0; ++k) { \
					if ((d = pj0[k] - pi0[k]) >= a && d <= b) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						_MASK_(*(px1++) = px0[k]); \
					} \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			BAND_SUBCASES(int, LOGICAL, HIDE);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			SET_SLOT(to, Matrix_xSym, x1);
			BAND_CASES;
			UNPROTECT(2); /* x1, x0 */
		}

	}

#undef BAND_CASES
#undef BAND_SUBCASES

	UNPROTECT(6);
	return to;
}

/* band(<[CRT]sparseMatrix>, k1, k2), tri[ul](<[CRT]sparseMatrix>, k) */
/* NB: argument validation more or less copied from R_dense_band() */
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1);

	int a, b;
	if (k1 == R_NilValue)
		a = (m > 0) ? 1-m : 0;
	else if ((a = asInteger(k1)) == NA_INTEGER || a < -m || a > n)
		error(_("'%s' must be an integer from %s to %s"),
		      "k1", "-Dim[1]", "Dim[2]");
	if (k2 == R_NilValue)
		b = (n > 0) ? n-1 : 0;
	else if ((b = asInteger(k2)) == NA_INTEGER || b < -m || b > n)
		error(_("'%s' must be an integer from %s to %s"),
		      "k2", "-Dim[1]", "Dim[2]");
	else if (b < a)
		error(_("'%s' must be less than or equal to '%s'"), "k1", "k2");

	return sparse_band(from, valid[ivalid], a, b);
}

SEXP sparse_diag_get(SEXP obj, const char *class, int names)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	UNPROTECT(1); /* dim */

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	}

	SEXP res = PROTECT(allocVector(kind2type(class[0]), r));

#define SDG_CASES \
	do { \
		switch (class[0]) { \
		case 'n': \
		case 'l': \
			SDG_LOOP(int, LOGICAL, SHOW, FIRSTOF, INCREMENT_LOGICAL, 0, 1); \
			break; \
		case 'i': \
			SDG_LOOP(int, INTEGER, SHOW, FIRSTOF, INCREMENT_INTEGER, 0, 1); \
			break; \
		case 'd': \
			SDG_LOOP(double, REAL, SHOW, FIRSTOF, INCREMENT_REAL, 0.0, 1.0); \
			break; \
		case 'z': \
			SDG_LOOP(Rcomplex, COMPLEX, SHOW, FIRSTOF, INCREMENT_COMPLEX, Matrix_zzero, Matrix_zone); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	if (di != 'N') {

		int j;

#undef SDG_LOOP
#define SDG_LOOP(_CTYPE_, _PTR_, _MASK_, _REPLACE_, _INCREMENT_, _ZERO_, _ONE_) \
		do { \
			_CTYPE_ *pres = _PTR_(res); \
			for (j = 0; j < r; ++j) \
				*(pres++) = _ONE_; \
		} while (0)

		SDG_CASES;

	} else if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int j, k, kend, *pp = INTEGER(p), *pi = INTEGER(i);
		pp++;

#undef SDG_LOOP
#define SDG_LOOP(_CTYPE_, _PTR_, _MASK_, _REPLACE_, _INCREMENT_, _ZERO_, _ONE_) \
		do { \
			_MASK_(_CTYPE_ *px = _PTR_(x)); \
			_CTYPE_ *pres = _PTR_(res); \
			if (class[1] == 'g') { \
				for (j = 0, k = 0; j < r; ++j) { \
					pres[j] = _ZERO_; \
					kend = pp[j]; \
					while (k < kend) { \
						if (pi[k] != j) \
							++k; \
						else { \
							pres[j] = _REPLACE_(px[k], 1); \
							k = kend; \
						} \
					} \
				} \
			} else if ((class[2] == 'C') == (ul != 'U')) { \
				for (j = 0, k = 0; j < r; ++j) { \
					kend = pp[j]; \
					pres[j] = (k < kend && pi[k] == j) \
						? _REPLACE_(px[k], 1) : _ZERO_; \
					k = kend; \
				} \
			} else { \
				for (j = 0, k = 0; j < r; ++j) { \
					kend = pp[j]; \
					pres[j] = (k < kend && pi[kend - 1] == j) \
						? _REPLACE_(px[kend - 1], 1) : _ZERO_; \
					k = kend; \
				} \
			} \
		} while (0)

		if (class[0] == 'n') {
			SDG_LOOP(int, LOGICAL, HIDE, SECONDOF, , 0, 1);
		} else {
			SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
			SDG_CASES;
			UNPROTECT(1); /* x */
		}

		UNPROTECT(2); /* i, p */

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, nnz = XLENGTH(i);

#undef SDG_LOOP
#define SDG_LOOP(_CTYPE_, _PTR_, _MASK_, _REPLACE_, _INCREMENT_, _ZERO_, _ONE_) \
		do { \
			_MASK_(_CTYPE_ *px = _PTR_(x)); \
			_CTYPE_ *pres = _PTR_(res); \
			Matrix_memset(pres, 0, r, sizeof(_CTYPE_)); \
			for (k = 0; k < nnz; ++k) { \
				if (*pi == *pj) \
					_INCREMENT_(pres[*pi], (*px)); \
				++pi; ++pj; _MASK_(++px); \
			} \
		} while (0)

		if (class[0] == 'n')
			SDG_LOOP(int, LOGICAL, HIDE, SECONDOF, INCREMENT_PATTERN, 0, 1);
		else {
			SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
			SDG_CASES;
			UNPROTECT(1); /* x */
		}

		UNPROTECT(2); /* j, i */

	}

#undef SDG_CASES
#undef SDG_LOOP

	if (names) {
		/* NB: The logic here must be adjusted once the validity method
		   for 'symmetricMatrix' enforces symmetric 'Dimnames'
		*/
		SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
			rn = VECTOR_ELT(dn, 0),
			cn = VECTOR_ELT(dn, 1);
		if (cn == R_NilValue) {
			if (class[1] == 's' && rn != R_NilValue)
				setAttrib(res, R_NamesSymbol, rn);
		} else {
			if (class[1] == 's')
				setAttrib(res, R_NamesSymbol, cn);
			else if (rn != R_NilValue &&
			         (rn == cn || equal_string_vectors(rn, cn, r)))
				setAttrib(res, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

	UNPROTECT(1); /* res */
	return res;
}

/* diag(<[CRT]sparseMatrix>, names=) */
SEXP R_sparse_diag_get(SEXP obj, SEXP names)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int names_;
	if (TYPEOF(names) != LGLSXP || LENGTH(names) < 1 ||
	    (names_ = LOGICAL(names)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "names", "TRUE", "FALSE");

	return sparse_diag_get(obj, valid[ivalid], names_);
}

SEXP sparse_diag_set(SEXP from, const char *class, SEXP value)
{
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(class));
	int v = LENGTH(value) != 1, full = 0;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	}

#define SDS_CASES \
	do { \
		switch (class[0]) { \
		case 'n': \
		case 'l': \
			SDS_LOOP(int, LOGICAL, SHOW, ISNZ_LOGICAL); \
			break; \
		case 'i': \
			SDS_LOOP(int, INTEGER, SHOW, ISNZ_INTEGER); \
			break; \
		case 'd': \
			SDS_LOOP(double, REAL, SHOW, ISNZ_REAL); \
			break; \
		case 'z': \
			SDS_LOOP(Rcomplex, COMPLEX, SHOW, ISNZ_COMPLEX); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	if (class[2] != 'T') {

		int n_ = (class[2] == 'C') ? n : m;
		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n_ + 1));
		SET_SLOT(to, Matrix_pSym, p1);
		int *pp0 = INTEGER(p0), *pp1 = INTEGER(p1), *pi0 = INTEGER(i0),
			j, k, kend, nd0 = 0, nd1 = 0;
		pp0++; *(pp1++) = 0;

		if (class[1] == 'g') {
			for (j = 0, k = 0; j < r; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] < j)
						++k;
					else {
						if (pi0[k] == j)
							++nd0;
						k = kend;
					}
				}
				pp1[j] = kend - nd0;
			}
			for (j = r; j < n_; ++j)
				pp1[j] = pp0[j] - nd0;
		} else if (di != 'N') {
			for (j = 0; j < n_; ++j)
				pp1[j] = pp0[j];
		} else if ((class[2] == 'C') == (ul == 'U')) {
			for (j = 0, k = 0; j < n_; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[kend - 1] == j)
					++nd0;
				k = kend;
				pp1[j] = kend - nd0;
			}
		} else {
			for (j = 0, k = 0; j < n_; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[k] == j)
					++nd0;
				k = kend;
				pp1[j] = kend - nd0;
			}
		}

#undef SDS_LOOP
#define SDS_LOOP(_CTYPE_, _PTR_, _MASK_, _NZ_) \
		do { \
			_CTYPE_ *pvalue = _PTR_(value); \
			if (v) { \
				for (j = 0; j < r; ++j) { \
					if (_NZ_(pvalue[j])) \
						++nd1; \
					pp1[j] += nd1; \
				} \
				for (j = r; j < n_; ++j) \
					pp1[j] += nd1; \
			} else if (_NZ_(pvalue[0])) { \
				full = 1; \
				for (j = 0; j < r; ++j) \
					pp1[j] += ++nd1; \
				for (j = r; j < n_; ++j) \
					pp1[j] +=   nd1; \
			} \
		} while (0)

		SDS_CASES;

		if (nd1 - nd0 > INT_MAX - pp0[n_ - 1])
			error(_("%s cannot exceed %s"), "p[length(p)]", "2^31-1");

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[n_ - 1]));
		SET_SLOT(to, iSym, i1);
		int *pi1 = INTEGER(i1);

#undef SDS_LOOP
#define SDS_LOOP(_CTYPE_, _PTR_, _MASK_, _NZ_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
			_CTYPE_ *pvalue = _PTR_(value); \
			for (j = 0, k = 0; j < r; ++j) { \
				kend = pp0[j]; \
				while (k < kend && pi0[k] < j) { \
					       *(pi1++) = pi0[k] ; \
					_MASK_(*(px1++) = px0[k]); \
					++k; \
				} \
				if (k < kend && pi0[k] == j) \
					++k; \
				if ((v) ? _NZ_(pvalue[j]) : full) { \
					       *(pi1++) = j                           ; \
					_MASK_(*(px1++) = (v) ? pvalue[j] : pvalue[0]); \
				} \
				while (k < kend) { \
					       *(pi1++) = pi0[k] ; \
					_MASK_(*(px1++) = px0[k]); \
					++k; \
				} \
			} \
			for (j = r; j < n_; ++j) { \
				kend = pp0[j]; \
				while (k < kend) { \
					       *(pi1++) = pi0[k] ; \
					_MASK_(*(px1++) = px0[k]); \
					++k; \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			SDS_LOOP(int, LOGICAL, HIDE, ISNZ_LOGICAL);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), pp1[n_ - 1]));
			SET_SLOT(to, Matrix_xSym, x1);
			SDS_CASES;
			UNPROTECT(2); /* x1, x0 */
		}

		UNPROTECT(4); /* i1, p1, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), j, nd0 = 0, nd1 = 0;
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = nnz0;

		if (di == 'N')
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] == pj0[k])
					++nd0;

#undef SDS_LOOP
#define SDS_LOOP(_CTYPE_, _PTR_, _MASK_, _NZ_) \
		do { \
			_CTYPE_ *pvalue = _PTR_(value); \
			if (v) { \
				for (j = 0; j < r; ++j) \
					if (_NZ_(pvalue[j])) \
						++nd1; \
			} else if (_NZ_(pvalue[0])) { \
				full = 1; \
				nd1 = r; \
			} \
		} while (0)

		SDS_CASES;

		if (nd1 - nd0 > R_XLEN_T_MAX - nnz0)
			error(_("%s cannot exceed %s"), "length(i)", "R_XLEN_T_MAX");
		nnz1 += nd1 - nd0;

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#undef SDS_LOOP
#define SDS_LOOP(_CTYPE_, _PTR_, _MASK_, _NZ_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
			_CTYPE_ *pvalue = _PTR_(value); \
			for (k = 0; k < nnz0; ++k) { \
				if (pi0[k] != pj0[k]) { \
					*(pi1++) = pi0[k]; \
					*(pj1++) = pj0[k]; \
					_MASK_(*(px1++) = px0[k]); \
				} \
			} \
			if (v) { \
				for (j = 0; j < r; ++j) { \
					if (_NZ_(pvalue[j])) { \
						*(pi1++) = *(pj1++) = j; \
						_MASK_((*px1++) = pvalue[j]); \
					} \
				} \
			} else if (full) { \
				for (j = 0; j < r; ++j) { \
					*(pi1++) = *(pj1++) = j; \
					_MASK_((*px1++) = pvalue[0]); \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			SDS_LOOP(int, LOGICAL, HIDE, ISNZ_LOGICAL);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			SET_SLOT(to, Matrix_xSym, x1);
			SDS_CASES;
			UNPROTECT(2); /* x1, x0 */
		}

		UNPROTECT(4); /* j1, i1, j0, i0 */

	}

#undef SDS_CASES
#undef SDS_LOOP

	UNPROTECT(1); /* to */
	return to;
}

/* diag(<[CRT]sparseMatrix>) <- value */
SEXP R_sparse_diag_set(SEXP from, SEXP value)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *class = valid[ivalid];

	SEXPTYPE tx = kind2type(class[0]), tv = TYPEOF(value);

	switch (tv) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
#ifdef MATRIX_ENABLE_ZMATRIX
	case CPLXSXP:
#endif
		break;
	default:
		error(_("replacement diagonal has incompatible type \"%s\""),
		      type2char(tv));
		break;
	}

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	R_xlen_t len = XLENGTH(value);
	if (len != 1 && len != r)
		error(_("replacement diagonal has wrong length"));

	if (tv <= tx) {
		PROTECT(from);
		PROTECT(value = coerceVector(value, tx));
	} else {
		/* defined in ./coerce.c : */
		SEXP sparse_as_kind(SEXP, const char *, char);
		PROTECT(from = sparse_as_kind(from, class, type2kind(tv)));
#ifndef MATRIX_ENABLE_IMATRIX
		if (tv == INTSXP)
			value = coerceVector(value, REALSXP);
#endif
		PROTECT(value);
		class = valid[R_check_class_etc(from, valid)];
	}

	from = sparse_diag_set(from, class, value);
	UNPROTECT(2);
	return from;
}

SEXP sparse_diag_U2N(SEXP from, const char *class)
{
	if (class[1] != 't')
		return from;

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	if (di == 'N')
		return from;

	SEXP val = PROTECT(ScalarLogical(1));
	from = R_sparse_diag_set(from, val);
	UNPROTECT(1); /* val */

	return from;
}

/* diagU2N(<[CRT]sparseMatrix>), parallel to R-level ..diagU2N(),
   though that is more general, working for _all_ Matrix */
SEXP R_sparse_diag_U2N(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_diag_U2N(from, valid[ivalid]);
}

SEXP sparse_diag_N2U(SEXP from, const char *class)
{
	if (class[1] != 't')
		return from;

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	if (di != 'N')
		return from;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */

	if (n == 0)
		PROTECT(from = duplicate(from));
	else {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */
		if (ul == 'U')
			PROTECT(from = sparse_band(from, class,  1, n - 1));
		else
			PROTECT(from = sparse_band(from, class, 1 - n, -1));
	}

	PROTECT(diag = mkString("U"));
	SET_SLOT(from, Matrix_diagSym, diag);
	UNPROTECT(2); /* diag, from */

	return from;
}

/* diagN2U(<[CRT]sparseMatrix>), parallel to R-level ..diagN2U(),
   though that is more general, working for _all_ Matrix */
SEXP R_sparse_diag_N2U(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_diag_N2U(from, valid[ivalid]);
}

SEXP sparse_transpose(SEXP from, const char *class, int lazy)
{
	SEXP to;
	if (class[2] == 'T' || !lazy)
		PROTECT(to = NEW_OBJECT_OF_CLASS(class));
	else {
		char cl[] = "...Matrix";
		cl[0] = class[0];
		cl[1] = class[1];
		cl[2] = (class[2] == 'C') ? 'R' : 'C';
		PROTECT(to = NEW_OBJECT_OF_CLASS(cl));
	}

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n) {
		UNPROTECT(1); /* dim */
		PROTECT(dim = GET_SLOT(to, Matrix_DimSym));
		pdim = INTEGER(dim);
		pdim[0] = n;
		pdim[1] = m;
	} else if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_reversed_DimNames(to, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */
		if (ul == 'U') {
			PROTECT(uplo = mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			char di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		} else {
			SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
			if (LENGTH(factors) > 0)
				SET_SLOT(to, Matrix_factorSym, factors);
			UNPROTECT(1); /* factors */
		}
	}

	/* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

	if (class[2] == 'T') {
		/* No need to allocate in this case: need only reverse 'i' and 'j' */
		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		SET_SLOT(to, Matrix_iSym, j0);
		SET_SLOT(to, Matrix_jSym, i0);
		UNPROTECT(2); /* j0, i0 */
		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			SET_SLOT(to, Matrix_xSym, x0);
			UNPROTECT(1); /* x */
		}
		UNPROTECT(1); /* to */
		return to;
	}

	/* Now dealing only with [CR]sparseMatrix ... */

	SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
		p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(from,        iSym));

	if (lazy) {
		/* No need to allocate in this case: need only reverse 'i' and 'j' */
		SEXP jSym = (class[2] == 'C') ? Matrix_jSym : Matrix_iSym;
		SET_SLOT(to, Matrix_pSym, p0);
		SET_SLOT(to,        jSym, i0);
		UNPROTECT(2); /* i0, p0 */
		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			SET_SLOT(to, Matrix_xSym, x0);
			UNPROTECT(1); /* x */
		}
		UNPROTECT(1); /* to */
		return to;
	}

	int m_ = (class[2] == 'C') ? m : n, n_ = (class[2] == 'C') ? n : m;
	SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) m_ + 1)),
		i1 = PROTECT(allocVector(INTSXP, INTEGER(p0)[n_]));
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to,        iSym, i1);

	/* defined in ./coerce.c : */
	void trans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, int, int);

	if (class[0] == 'n')
		trans(p0, i0, NULL, p1, i1, NULL, m_, n_);
	else {
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
			x1 = PROTECT(allocVector(TYPEOF(x0), INTEGER(p0)[n_]));
		SET_SLOT(to, Matrix_xSym, x1);
		trans(p0, i0, x0, p1, i1, x1, m_, n_);
		UNPROTECT(2); /* x1, x0 */
	}
	UNPROTECT(5); /* i1, p1, i0, p0, to */
	return to;
}

/* t(<[CRT]sparseMatrix>) */
SEXP R_sparse_transpose(SEXP from, SEXP lazy)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	int lazy_;
	if (TYPEOF(lazy) != LGLSXP || LENGTH(lazy) < 1 ||
	    (lazy_ = LOGICAL(lazy)[0]) == NA_LOGICAL)
		error(_("invalid '%s' to %s()"), "lazy", __func__);

	return sparse_transpose(from, valid[ivalid], lazy_);
}

SEXP sparse_force_symmetric(SEXP from, const char *class, char ul)
{
	char ulf = 'U', ult = 'U';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ulf = ult = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */
	}
	if (ul != '\0')
		ult = ul;

	if (class[1] == 's') {
		/* .s[CRT]Matrix */
		if (ulf == ult)
			return from;
		SEXP to = PROTECT(sparse_transpose(from, class, 0));
		if (class[0] == 'z') {
			/* Need _conjugate_ transpose */
			SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
			conjugate(x);
			UNPROTECT(1); /* x */
		}
		UNPROTECT(1) /* to */;
		return to;
	}

	/* Now handling just .[gt][CRT]Matrix ... */

	char cl[] = ".s.Matrix";
	cl[0] = class[0];
	cl[2] = class[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to symmetrize a non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (ult != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	char di = 'N';
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */
	}

	/* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

#define FS_CASES \
	do { \
		switch (class[0]) { \
		case 'l': \
			FS_SUBCASES(int, LOGICAL, SHOW, 1); \
			break; \
		case 'i': \
			FS_SUBCASES(int, INTEGER, SHOW, 1); \
			break; \
		case 'd': \
			FS_SUBCASES(double, REAL, SHOW, 1.0); \
			break; \
		case 'z': \
			FS_SUBCASES(Rcomplex, COMPLEX, SHOW, Matrix_zone); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	if (class[1] == 't' && di == 'N' && ulf == ult) {

		/* No need to allocate in this case: we have the triangle we want */
		if (class[2] != 'T') {
			SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
			SET_SLOT(to, Matrix_pSym, p);
			UNPROTECT(1); /* p */
		}
		if (class[2] != 'R') {
			SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
			SET_SLOT(to, Matrix_iSym, i);
			UNPROTECT(1); /* i */
		}
		if (class[2] != 'C') {
			SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
			SET_SLOT(to, Matrix_jSym, j);
			UNPROTECT(1); /* j */
		}
		if (class[0] != 'n') {
			SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
			SET_SLOT(to, Matrix_xSym, x);
			UNPROTECT(1); /* x */
		}
		UNPROTECT(1); /* to */
		return to;

	} else if (class[2] != 'T') {

		/* Symmetrizing square .[gt][CR]Matrix ... */

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
			i0 = PROTECT(GET_SLOT(from,        iSym));
		int j, k, kend,
			*pp0 = INTEGER(p0),
			*pp1 = INTEGER(p1),
			*pi0 = INTEGER(i0),
			nnz0 = pp0[n],
			nnz1 = 0;
		pp0++; *(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

		/* Counting number of nonzero elements in triangle, by "column" ... */

		if (class[1] == 't') {
			if (di != 'N') {
				/* Have triangular matrix with unit diagonal */
				if (ulf != ult) {
					/* Returning identity matrix */
					for (j = 0; j < n; ++j)
						pp1[j] = ++nnz1;
				} else {
					/* Returning symmetric matrix with unit diagonal */
					for (j = 0; j < n; ++j)
						pp1[j] = ++nnz1 + pp0[j];
					nnz1 += nnz0;
				}
			} else if (ulf == ((class[2] == 'C') ? 'U' : 'L')) {
				/* Have triangular matrix with non-unit "trailing" diagonal
				   and returning diagonal part */
				for (j = 0; j < n; ++j) {
					if (pp0[j - 1] < pp0[j] && pi0[pp0[j] - 1] == j)
						++nnz1;
					pp1[j] = nnz1;
				}
			} else {
				/* Have triangular matrix with non-unit "leading" diagonal
				   and returning diagonal part */
				for (j = 0; j < n; ++j) {
					if (pp0[j - 1] < pp0[j] && pi0[pp0[j - 1]] == j)
						++nnz1;
					pp1[j] = nnz1;
				}
			}
		} else if (ult == ((class[2] == 'C') ? 'U' : 'L')) {
			/* Have general matrix and returning upper triangle */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] <= j)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
		} else {
			/* Have general matrix and returning lower triangle */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] >= j)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
		}

		/* Now allocating and filling out slots ... */

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#undef FS_SUBCASES
#define FS_SUBCASES(_CTYPE_, _PTR_, _MASK_, _ONE_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
			if (class[1] == 't') { \
				if (di != 'N') { \
					/* Have triangular matrix with unit diagonal */ \
					if (ulf != ult) { \
						/* Returning identity matrix */ \
						for (j = 0; j < n; ++j) { \
							*(pi1++) = j; \
							_MASK_(*(px1++) = _ONE_); \
						} \
					} else if (ulf == ((class[2] == 'C') ? 'U' : 'L')) { \
						/* Returning symmetric matrix    */ \
						/* with unit "trailing" diagonal */ \
						for (j = 0, k = 0; j < n; ++j) { \
							kend = pp0[j]; \
							while (k < kend) { \
								*(pi1++) = pi0[k]; \
								_MASK_(*(px1++) = px0[k]); \
								++k; \
							} \
							*(pi1++) = j; \
							_MASK_(*(px1++) = _ONE_); \
						} \
					} else { \
						/* Returning symmetric matrix   */ \
						/* with unit "leading" diagonal */ \
						for (j = 0, k = 0; j < n; ++j) { \
							*(pi1++) = j; \
							_MASK_(*(px1++) = _ONE_); \
							kend = pp0[j]; \
							while (k < kend) { \
								*(pi1++) = pi0[k]; \
								_MASK_(*(px1++) = px0[k]); \
								++k; \
							} \
						} \
					} \
				} else if (ulf == ((class[2] == 'C') ? 'U' : 'L')) { \
					/* Have triangular matrix with non-unit "trailing" */ \
					/* diagonal and returning diagonal part            */ \
					for (j = 0; j < n; ++j) { \
						if (pp0[j - 1] < pp0[j] && pi0[pp0[j] - 1] == j) { \
							*(pi1++) = j; \
							_MASK_(*(px1++) = px0[pp0[j] - 1]); \
						} \
					} \
				} else { \
					/* Have triangular matrix with non-unit "leading" */ \
					/* diagonal and returning diagonal part           */ \
					for (j = 0; j < n; ++j) { \
						if (pp0[j - 1] < pp0[j] && pi0[pp0[j - 1]] == j) { \
							*(pi1++) = j; \
							_MASK_(*(px1++) = px0[pp0[j - 1]]); \
						} \
					} \
				} \
			} else if (ult == ((class[2] == 'C') ? 'U' : 'L')) { \
				/* Have general matrix and returning upper triangle */ \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] <= j) { \
							*(pi1++) = pi0[k]; \
							_MASK_(*(px1++) = px0[k]); \
						} \
						++k; \
					} \
				} \
			} else { \
				/* Have general matrix and returning lower triangle */ \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] >= j) { \
							*(pi1++) = pi0[k]; \
							_MASK_(*(px1++) = px0[k]); \
						} \
						++k; \
					} \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			FS_SUBCASES(int, LOGICAL, HIDE, 1);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			SET_SLOT(to, Matrix_xSym, x1);
			FS_CASES;
			UNPROTECT(2); /* x1, x0 */
		}

	} else {

		/* Symmetrizing square .[gt]TMatrix ... */

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0),
			*pj0 = INTEGER(j0);
		R_xlen_t j, k, nnz0 = XLENGTH(i0), nnz1 = 0;

		/* Counting number of nonzero elements in triangle ... */

		if (class[1] == 't' && di != 'N')
			nnz1 = (ulf == ult) ? n + nnz0 : n;
		else {
			if (ult == 'U') {
				for (k = 0; k < nnz0; ++k)
					if (pi0[k] <= pj0[k])
						++nnz1;
			} else {
				for (k = 0; k < nnz0; ++k)
					if (pi0[k] >= pj0[k])
						++nnz1;
			}
		}

		/* Now allocating and filling out slots ... */

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1),
			*pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);

#undef FS_SUBCASES
#define FS_SUBCASES(_CTYPE_, _PTR_, _MASK_, _ONE_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
			if (class[1] == 't' && di != 'N') { \
				if (ulf == ult) { \
					Matrix_memcpy(pi1, pi0, nnz0, sizeof(int)); \
					Matrix_memcpy(pj1, pj0, nnz0, sizeof(int)); \
					_MASK_(Matrix_memcpy(px1, px0, nnz0, sizeof(_CTYPE_))); \
					pi1 += nnz0; \
					pj1 += nnz0; \
					_MASK_(px1 += nnz0); \
				} \
				for (j = 0; j < n; ++j) { \
					*(pi1++) = *(pj1++) = j; \
					_MASK_(*(px1++) = _ONE_); \
				} \
			} else { \
				if (ult == 'U') { \
					for (k = 0; k < nnz0; ++k) { \
						if (pi0[k] <= pj0[k]) { \
							*(pi1++) = pi0[k]; \
							*(pj1++) = pj0[k]; \
							_MASK_(*(px1++) = px0[k]); \
						} \
					} \
				} else { \
					for (k = 0; k < nnz0; ++k) { \
						if (pi0[k] <= pj0[k]) { \
							*(pi1++) = pi0[k]; \
							*(pj1++) = pj0[k]; \
							_MASK_(*(px1++) = px0[k]); \
						} \
					} \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			FS_SUBCASES(int, LOGICAL, HIDE, 1);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			SET_SLOT(to, Matrix_xSym, x1);
			FS_CASES;
			UNPROTECT(2); /* x1, x0 */
		}

	}

#undef FS_CASES
#undef FS_SUBCASES

	UNPROTECT(5);
	return to;
}

/* forceSymmetric(<[CRT]sparseMatrix>, uplo) */
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char ul = '\0';
	if (uplo != R_NilValue) {
		if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
			error(_("invalid '%s' to %s()"), "uplo", __func__);
	}

	return sparse_force_symmetric(from, valid[ivalid], ul);
}

SEXP sparse_symmpart(SEXP from, const char *class)
{
	if (class[0] != 'z' && class[0] != 'd') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_kind(SEXP, const char *, char);
		from = sparse_as_kind(from, class, 'd');
	}
	if (class[0] != 'z' && class[1] == 's')
		return from;

	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from, &pid);

	char cl[] = ".s.Matrix";
	cl[0] = (class[0] != 'z') ? 'd' : 'z';
	cl[2] = class[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to get symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	} else if (class[2] == 'R') {
		SEXP uplo = PROTECT(mkString("L"));
		ul = 'L';
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

#define ASSIGN_REAL(_X_, _Y_) \
	do { _X_   = _Y_                  ; } while (0)
#define ASSIGN_COMPLEX(_X_, _Y_) \
	do { _X_.r = _Y_.r; _X_.i  = _Y_.i; } while (0)

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int j, k, kend, *pp0 = INTEGER(p0) + 1, *pi0 = INTEGER(i0),
			nnz = pp0[n - 1];

		if (class[1] == 'g') {

			cl[1] = 'g';
			REPROTECT(from = sparse_transpose(from, cl, 0), pid);

			SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
				p0_ = PROTECT(GET_SLOT(from, Matrix_pSym)),
				i0_ = PROTECT(GET_SLOT(from,        iSym)),
				x0_ = PROTECT(GET_SLOT(from, Matrix_xSym));
			int k_, kend_, *pp0_ = INTEGER(p0_) + 1, *pi0_ = INTEGER(i0_),
				*pp1 = INTEGER(p1);
			*(pp1++) = 0;

			for (j = 0, k = 0, k_ = 0; j < n; ++j) {
				kend = pp0[j];
				kend_ = pp0_[j];
				pp1[j] = pp1[j - 1];
				while (k < kend) {
					if (pi0[k] > j)
						k = kend;
					else {
						while (k_ < kend_ && pi0_[k_] < pi0[k]) {
							++pp1[j];
							++k_;
						}
						++pp1[j];
						if (k_ < kend_ && pi0_[k_] == pi0[k])
							++k_;
						++k;
					}
				}
				while (k_ < kend_) {
					if (pi0_[k_] > j)
						k_ = kend_;
					else {
						++pp1[j];
						++k_;
					}
				}
			}

			SEXP i1 = PROTECT(allocVector(INTSXP, pp1[n - 1])),
				x1 = PROTECT(allocVector(kind2type(cl[0]), pp1[n - 1]));
			int *pi1 = INTEGER(i1);

#undef SP_LOOP
#define SP_LOOP(_ASSIGN_, _INCREMENT_) \
			do { \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend = pp0[j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (pi0[k] > j) \
							k = kend; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								*pi1 = pi0_[k_]; \
								_ASSIGN_((*px1), 0.5 * px0_[k_]); \
								++k_; ++pi1; ++px1; \
							} \
							*pi1 = pi0[k]; \
							_ASSIGN_((*px1), 0.5 * px0[k]); \
							if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
								_INCREMENT_((*px1), 0.5 * px0_[k_]); \
								++k_; \
							} \
							++k; ++pi1; ++px1; \
						} \
					} \
					while (k_ < kend_) { \
						if (pi0_[k_] > j) \
							k_ = kend_; \
						else { \
							*pi1 = pi0_[k_]; \
							_ASSIGN_((*px1), 0.5 * px0_[k_]); \
							++k_; ++pi1; ++px1; \
						} \
					} \
				} \
			} while (0)


			if (cl[0] == 'd') {
				double *px0 = REAL(x0), *px1 = REAL(x1),
					*px0_ = REAL(x0_);
				SP_LOOP(ASSIGN_REAL, INCREMENT_REAL);
			} else {
				Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1),
					*px0_ = COMPLEX(x0_);
				SP_LOOP(ASSIGN_COMPLEX, INCREMENT_COMPLEX);
			}

			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to,        iSym, i1);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(6); /* x1, i1, p1, x0_, i0_, p0_ */

		} else if (class[1] == 't') {

			int leading = (class[2] == 'C') == (ul != 'U');

			if (di == 'N') {

				SEXP x1 = PROTECT(allocVector(kind2type(cl[0]), nnz));

#undef SP_LOOP
#define SP_LOOP(_ASSIGN_) \
				do { \
					if (leading) { \
						for (j = 0, k = 0; j < n; ++j) { \
							kend = pp0[j]; \
							if (k < kend) { \
								if (pi0[k] == j) \
									*px1 = *px0; \
								else \
									_ASSIGN_((*px1), 0.5 * (*px0)); \
								++k, ++px0; ++px1; \
								while (k < kend) { \
									_ASSIGN_((*px1), 0.5 * (*px0)); \
									++k; ++px0; ++px1; \
								} \
							} \
						} \
					} else { \
						for (j = 0, k = 0; j < n; ++j) { \
							kend = pp0[j]; \
							if (k < kend) { \
								while (k < kend - 1) { \
									_ASSIGN_((*px1), 0.5 * (*px0)); \
									++k; ++px0; ++px1; \
								} \
								if (pi0[k] == j) \
									*px1 = *px0; \
								else \
									_ASSIGN_((*px1), 0.5 * (*px0)); \
								++k; ++px0; ++px1; \
							} \
						} \
					} \
				} while (0)

				if (cl[0] == 'd') {
					double *px0 = REAL(x0), *px1 = REAL(x1);
					SP_LOOP(ASSIGN_REAL);
				} else {
					Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
					SP_LOOP(ASSIGN_COMPLEX);
				}

				SET_SLOT(to, Matrix_pSym, p0);
				SET_SLOT(to,        iSym, i0);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */

			} else {

				nnz += n;
				SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
					i1 = PROTECT(allocVector(INTSXP, nnz)),
					x1 = PROTECT(allocVector(kind2type(cl[0]), nnz));
				int *pp1 = INTEGER(p1), *pi1 = INTEGER(i1);
				*(pp1++) = 0;

#undef SP_LOOP
#define SP_LOOP(_ASSIGN_, _ONE_) \
				do { \
					if (leading) { \
						for (j = 0, k = 0; j < n; ++j) { \
							kend = pp0[j]; \
							pp1[j] = pp1[j - 1] + kend - k + 1; \
							*pi1 = j; \
							_ASSIGN_((*px1), _ONE_); \
							++pi1, ++px1; \
							while (k < kend) { \
								*pi1 = *pi0; \
								_ASSIGN_((*px1), 0.5 * (*px0)); \
								++k; ++pi0; ++pi1; ++px0; ++px1; \
							} \
						} \
					} else { \
						for (j = 0, k = 0; j < n; ++j) { \
							kend = pp0[j]; \
							pp1[j] = pp1[j - 1] + kend - k + 1; \
							while (k < kend) { \
								*pi1 = *pi0; \
								_ASSIGN_((*px1), 0.5 * (*px0)); \
								++k; ++pi0; ++pi1; ++px0; ++px1; \
							} \
							*pi1 = j; \
							_ASSIGN_((*px1), _ONE_); \
							++pi1; ++px1; \
						} \
					} \
				} while (0)

				if (cl[0] == 'd') {
					double *px0 = REAL(x0), *px1 = REAL(x1);
					SP_LOOP(ASSIGN_REAL, 1.0);
				} else {
					Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
					SP_LOOP(ASSIGN_COMPLEX, Matrix_zone);
				}

				SET_SLOT(to, Matrix_pSym, p1);
				SET_SLOT(to,        iSym, i1);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(3); /* p1, i1, x1 */

			}

		} else {

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

			if (cl[0] == 'd')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(duplicate(x0));
				zeroIm(x1);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}

		}

		UNPROTECT(3); /* x0, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz = XLENGTH(i0);

		if (class[1] == 'g') {

			SEXP i1 = PROTECT(allocVector(INTSXP, nnz)),
				j1 = PROTECT(allocVector(INTSXP, nnz)),
				x1 = PROTECT(allocVector(kind2type(cl[0]), nnz));
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#undef SP_LOOP
#define SP_LOOP(_ASSIGN_) \
			do { \
				for (k = 0; k < nnz; ++k) { \
					if (*pi0 == *pj0) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						*px1 = *px0; \
					} else if (*pi0 < *pj0) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						_ASSIGN_((*px1), 0.5 * (*px0)); \
					} else { \
						*pi1 = *pj0; \
						*pj1 = *pi0; \
						_ASSIGN_((*px1), 0.5 * (*px0)); \
					} \
					++pi0; ++pi1; ++pj0; ++pj1; ++px0; ++px1; \
				} \
			} while (0)

			if (cl[0] == 'd') {
				double *px0 = REAL(x0), *px1 = REAL(x1);
				SP_LOOP(ASSIGN_REAL);
			} else {
				Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
				SP_LOOP(ASSIGN_COMPLEX);
			}

			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(3); /* x1, j1, i1 */

		} else if (class[1] == 't') {

			if (di == 'N') {

				SEXP x1 = PROTECT(allocVector(kind2type(cl[0]), nnz));

#undef SP_LOOP
#define SP_LOOP(_ASSIGN_) \
				do { \
					for (k = 0; k < nnz; ++k) { \
						if (*pi0 == *pj0) \
							*px1 = *px0; \
						else \
							_ASSIGN_((*px1), 0.5 * (*px0)); \
						++px0; ++px1; \
					} \
				} while (0)

				if (cl[0] == 'd') {
					double *px0 = REAL(x0), *px1 = REAL(x1);
					SP_LOOP(ASSIGN_REAL);
				} else {
					Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
					SP_LOOP(ASSIGN_COMPLEX);
				}

				SET_SLOT(to, Matrix_iSym, i0);
				SET_SLOT(to, Matrix_jSym, j0);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */

			} else {

				SEXP i1 = PROTECT(allocVector(INTSXP, nnz + n)),
					j1 = PROTECT(allocVector(INTSXP, nnz + n)),
					x1 = PROTECT(allocVector(kind2type(cl[0]), nnz + n));
				int j, *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#undef SP_LOOP
#define SP_LOOP(_ASSIGN_, _ONE_) \
				do { \
					for (k = 0; k < nnz; ++k) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						_ASSIGN_((*px1), 0.5 * (*px0)); \
						++pi0; ++pi1; ++pj0; ++pj1; ++px0; ++px1; \
					} \
					for (j = 0; j < n; ++j) { \
						*pi1 = *pj1 = j; \
						_ASSIGN_((*px1), _ONE_); \
						++pi1; ++pj1; ++px1; \
					} \
				} while (0)

				if (cl[0] == 'd') {
					double *px0 = REAL(x0), *px1 = REAL(x1);
					SP_LOOP(ASSIGN_REAL, 1.0);
				} else {
					Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
					SP_LOOP(ASSIGN_COMPLEX, Matrix_zone);
				}

				SET_SLOT(to, Matrix_iSym, i1);
				SET_SLOT(to, Matrix_jSym, j1);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(3); /* x1, j1, i1 */

			}

		} else {

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

			if (cl[0] == 'd')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(duplicate(x0));
				zeroIm(x1);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}

		}

		UNPROTECT(3); /* x0, j0, i1 */

	}

	UNPROTECT(2); /* to, from */
	return to;
}

/* symmpart(<[CRT]sparseMatrix>) */
SEXP R_sparse_symmpart(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_symmpart(from, valid[ivalid]);
}

SEXP sparse_skewpart(SEXP from, const char *class)
{
	if (class[0] != 'z' && class[0] != 'd') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_kind(SEXP, const char *, char);
		from = sparse_as_kind(from, class, 'd');
	}

	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from, &pid);

	char cl[] = "...Matrix";
	cl[0] = (class[0] != 'z') ? 'd' : 'z';
	cl[1] = (class[1] != 's') ? 'g' : 's';
	cl[2] = class[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to get skew-symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (class[1] == 's') {

		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		/* Skew-symmetric part of Hermitian matrix is imaginary part */
		if (class[0] != 'z') {
			if (class[2] != 'T') {
				SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1));
				int *pp1 = INTEGER(p1);
				Matrix_memset(pp1, 0, (R_xlen_t) n + 1, sizeof(int));
				SET_SLOT(to, Matrix_pSym, p1);
				UNPROTECT(1); /* p1 */
			}
		} else {
			if (class[2] != 'T') {
				SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym));
				SET_SLOT(to, Matrix_pSym, p0);
				UNPROTECT(1); /* p0 */
			}
			if (class[2] != 'R') {
				SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym));
				SET_SLOT(to, Matrix_iSym, i0);
				UNPROTECT(1); /* i0 */
			}
			if (class[2] != 'C') {
				SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
				SET_SLOT(to, Matrix_jSym, j0);
				UNPROTECT(1); /* j0 */
			}
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(duplicate(x0));
			zeroRe(x1);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(2); /* x1, x0 */
		}

	} else if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1));
		int j, k, kend, k_, kend_, *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			*pp1 = INTEGER(p1);
		pp0++; *(pp1++) = 0;

		REPROTECT(from = sparse_transpose(from, cl, 0), pid);

		SEXP
			p0_ = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0_ = PROTECT(GET_SLOT(from,        iSym)),
			x0_ = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pp0_ = INTEGER(p0_), *pi0_ = INTEGER(i0_), *pp1_;
		pp0_++;
		Matrix_Calloc(pp1_, n, int);

		for (j = 0, k = 0, k_ = 0; j < n; ++j) {
			kend = pp0[j];
			kend_ = pp0_[j];
			while (k < kend) {
				if (pi0[k] >= j)
					k = kend;
				else {
					while (k_ < kend_ && pi0_[k_] < pi0[k]) {
						++pp1_[j];
						++pp1_[pi0_[k_]];
						++k_;
					}
					++pp1_[j];
					++pp1_[pi0[k]];
					if (k_ < kend_ && pi0_[k_] == pi0[k])
						++k_;
					++k;
				}
			}
			while (k_ < kend_) {
				if (pi0_[k_] >= j)
					k_ = kend_;
				else {
					++pp1_[j];
					++pp1_[pi0_[k_]];
					++k_;
				}
			}
		}

		for (j = 0; j < n; ++j) {
			pp1[j] = pp1[j - 1] + pp1_[j];
			pp1_[j] = pp1[j - 1];
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[n - 1])),
			x1 = PROTECT(allocVector(kind2type(cl[0]), pp1[n - 1]));
		int *pi1 = INTEGER(i1);

#undef SP_LOOP
#define SP_LOOP(_ASSIGN_, _INCREMENT_) \
		do { \
			for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
				kend = pp0[j]; \
				kend_ = pp0_[j]; \
				while (k < kend) { \
					if (pi0[k] >= j) \
						k = kend; \
					else { \
						while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
							pi1[pp1_[j]] = pi0_[k_]; \
							_ASSIGN_(px1[pp1_[j]], -0.5 * px0_[k_]); \
							pi1[pp1_[pi0_[k_]]] = j; \
							_ASSIGN_(px1[pp1_[pi0_[k_]]], -px1[pp1_[j]]); \
							++pp1_[j]; \
							++pp1_[pi0_[k_]]; \
							++k_; \
						} \
						pi1[pp1_[j]] = pi0[k]; \
						_ASSIGN_(px1[pp1_[j]], 0.5 * px0[k]); \
						if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
							_INCREMENT_(px1[pp1_[j]], -0.5 * px0_[k_]); \
							++k_; \
						} \
						pi1[pp1_[pi0[k]]] = j; \
						_ASSIGN_(px1[pp1_[pi0[k]]], -px1[pp1_[j]]); \
						++pp1_[j]; \
						++pp1_[pi0[k]]; \
						++k; \
					} \
				} \
				while (k_ < kend_) { \
					if (pi0_[k_] >= j) \
						k_ = kend_; \
					else { \
						pi1[pp1_[j]] = pi0_[k_]; \
						_ASSIGN_(px1[pp1_[j]], -0.5 * px0_[k_]); \
						pi1[pp1_[pi0_[k_]]] = j; \
						_ASSIGN_(px1[pp1_[pi0_[k_]]], -px1[pp1_[j]]); \
						++pp1_[j]; \
						++pp1_[pi0_[k_]]; \
						++k_; \
					} \
				} \
			} \
		} while (0)

		if (cl[0] == 'd') {
			double *px0 = REAL(x0), *px1 = REAL(x1),
				*px0_ = REAL(x0_);
			SP_LOOP(ASSIGN_REAL, INCREMENT_REAL);
		} else {
			Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1),
				*px0_ = COMPLEX(x0_);
			SP_LOOP(ASSIGN_COMPLEX, INCREMENT_COMPLEX);
		}

		Matrix_Free(pp1_, n);
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to,        iSym, i1);
		SET_SLOT(to, Matrix_xSym, x1);
		UNPROTECT(9); /* x1, i1, p1, x0, i0, p0, x0_, i0_, p0_ */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = nnz0;

		for (k = 0; k < nnz0; ++k)
			if (pi0[k] == pj0[k])
				--nnz1;
		nnz1 *= 2;

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1)),
			x1 = PROTECT(allocVector(kind2type(cl[0]), nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#undef SP_LOOP
#define SP_LOOP(_ASSIGN_) \
		do { \
			for (k = 0; k < nnz0; ++k) { \
				if (*pi0 != *pj0) { \
					*pi1 = *pi0; \
					*pj1 = *pj0; \
					_ASSIGN_((*px1),  0.5 * (*px0)); \
					++pi1; ++pj1; ++px1; \
					*pi1 = *pj0; \
					*pj1 = *pi0; \
					_ASSIGN_((*px1), -0.5 * (*px0)); \
					++pi1; ++pj1; ++px1; \
				} \
				++pi0; ++pj0; ++px0; \
			} \
		} while (0)

		if (cl[0] == 'd') {
			double *px0 = REAL(x0), *px1 = REAL(x1);
			SP_LOOP(ASSIGN_REAL);
		} else {
			Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
			SP_LOOP(ASSIGN_COMPLEX);
		}

		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		SET_SLOT(to, Matrix_xSym, x1);
		UNPROTECT(6); /* x1, j1, i1, x0, j0, i0 */

	}

#undef ASSIGN_REAL
#undef ASSIGN_COMPLEX
#undef INCREMENT_REAL
#undef INCREMENT_COMPLEX
#undef SP_LOOP

	UNPROTECT(2); /* to, from */
	return to;
}

/* skewpart(<[CRT]sparseMatrix>) */
SEXP R_sparse_skewpart(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_skewpart(from, valid[ivalid]);
}

SEXP Tsparse_aggregate(SEXP from)
{
	static const char *valid[] = { VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid];

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */

	SEXP to,
		i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
		j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
		i1 = NULL, j1 = NULL;

	/* defined in ./coerce.c : */
	void taggr(SEXP, SEXP, SEXP, SEXP *, SEXP *, SEXP *, int, int);

	if (cl[0] == 'n') {
		taggr(j0, i0, NULL, &j1, &i1, NULL, n, m);
		if (!i1) {
			UNPROTECT(2); /* j0, i0 */
			return from;
		}
		PROTECT(i1);
		PROTECT(j1);
		PROTECT(to = NEW_OBJECT_OF_CLASS(cl));
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		UNPROTECT(5); /* to, j1, i1, j0, i0 */
	} else {
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
			x1 = NULL;
		taggr(j0, i0, x0, &j1, &i1, &x1, n, m);
		if (!i1) {
			UNPROTECT(3); /* x0, j0, i0 */
			return from;
		}
		PROTECT(i1);
		PROTECT(j1);
		PROTECT(x1);
		PROTECT(to = NEW_OBJECT_OF_CLASS(cl));
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		SET_SLOT(to, Matrix_xSym, x1);
		UNPROTECT(7); /* to, x1, j1, i1, x0, j0, i0 */
	}

	PROTECT(to);

	if (m != n || n > 0) {
		PROTECT(dim = GET_SLOT(to, Matrix_DimSym));
		pdim = INTEGER(dim);
		pdim[0] = m;
		pdim[1] = n;
		UNPROTECT(1); /* dim */
	}

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (cl[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (cl[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorSym, factors);
		UNPROTECT(1); /* factors */
	}

	UNPROTECT(1); /* to */
	return to;
}

/* isDiagonal(<[CR]sparseMatrix>) */
#define CR_IS_DIAGONAL(_C_, _I_) \
SEXP _C_ ## sparse_is_diagonal(SEXP obj) \
{ \
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	UNPROTECT(1); /* dim */ \
	if (m != n) \
		return ScalarLogical(0); \
	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)); \
	int *pp = INTEGER(p); \
	if (pp[n] > n) { \
		UNPROTECT(1); /* p */ \
		return ScalarLogical(0); \
	} \
	SEXP i = PROTECT(GET_SLOT(obj, Matrix_ ## _I_ ## Sym)); \
	int d, j, *pi = INTEGER(i); \
	Rboolean res = TRUE; \
	for (j = 0; j < n; ++j) { \
		if ((d = pp[j+1] - pp[j]) > 1 || (d == 1 && *(pi++) != j)) { \
			res = FALSE; \
			break; \
		} \
	} \
	UNPROTECT(2); /* i, p */ \
	return ScalarLogical(res); \
}

/* Csparse_is_diagonal() */
CR_IS_DIAGONAL(C, i)
/* Rsparse_is_diagonal() */
CR_IS_DIAGONAL(R, j)

/* isDiagonal(<TsparseMatrix>) */
SEXP Tsparse_is_diagonal(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */
	if (m != n)
		return ScalarLogical(0);
	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	int *pi = INTEGER(i), *pj = INTEGER(j);
	R_xlen_t k, nnz = XLENGTH(i);
	Rboolean res = TRUE;
	for (k = 0; k < nnz; ++k) {
		if (*(pi++) != *(pj++)) {
			res = FALSE;
			break;
		}
	}
	UNPROTECT(2); /* j, i */
	return ScalarLogical(res);
}

/* isTriangular(<.g[CR]Matrix>, upper) */
#define CR_IS_TRIANGULAR(_C_, _I_, _UPPER_, _LOWER_) \
SEXP _C_ ## sparse_is_triangular(SEXP obj, SEXP upper) \
{ \
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	UNPROTECT(1); /* dim */ \
	if (m != n) \
		return ScalarLogical(0); \
	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)), \
		i = PROTECT(GET_SLOT(obj, Matrix_ ## _I_ ## Sym)); \
	int j, k, kend, *pp = INTEGER(p), *pi = INTEGER(i), \
		need_upper = asLogical(upper); \
	Rboolean res = TRUE; \
	++pp; \
	if (need_upper == NA_LOGICAL) { \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if (_LOWER_) \
					goto opposite; \
				++k; \
			} \
		} \
		UNPROTECT(2); /* i, p */ \
		RETURN_TRUE_OF_KIND("U"); \
	opposite: \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if (_UPPER_) { \
					res = FALSE; \
					goto nokind; \
				} \
				++k; \
			} \
		} \
		UNPROTECT(2); /* i, p */ \
		RETURN_TRUE_OF_KIND("L"); \
	} else if (need_upper != 0) { \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if (_LOWER_) { \
					res = FALSE; \
					goto nokind; \
				} \
				++k; \
			} \
		} \
	} else { \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if (_UPPER_) { \
					res = FALSE; \
					goto nokind; \
				} \
				++k; \
			} \
		} \
	} \
nokind: \
	UNPROTECT(2); /* i, p */ \
	return ScalarLogical(res); \
}

/* Csparse_is_triangular() */
CR_IS_TRIANGULAR(C, i, pi[k] < j, pi[k] > j)
/* Rsparse_is_triangular() */
CR_IS_TRIANGULAR(R, j, pi[k] > j, pi[k] < j)

/* isTriangular(<.gTMatrix>, upper) */
SEXP Tsparse_is_triangular(SEXP obj, SEXP upper)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */
	if (m != n)
		return ScalarLogical(0);
	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	int *pi = INTEGER(i), *pj = INTEGER(j), need_upper = asLogical(upper);
	R_xlen_t k, nnz = XLENGTH(i);
	Rboolean res = TRUE;
	if (need_upper == NA_LOGICAL) {
		for (k = 0; k < nnz; ++k)
			if (pi[k] > pj[k])
				goto opposite;
		UNPROTECT(2); /* j, i */
		RETURN_TRUE_OF_KIND("U");
	opposite:
		for (k = 0; k < nnz; ++k)
			if (pi[k] < pj[k]) {
				res = FALSE;
				goto nokind;
			}
		UNPROTECT(2); /* j, i */
		RETURN_TRUE_OF_KIND("L");
	} else if (need_upper != 0) {
		for (k = 0; k < nnz; ++k)
			if (pi[k] > pj[k]) {
				res = FALSE;
				goto nokind;
			}
	} else {
		for (k = 0; k < nnz; ++k)
			if (pi[k] < pj[k]) {
				res = FALSE;
				goto nokind;
			}
	}
nokind:
	UNPROTECT(2); /* j, i */
	return ScalarLogical(res);
}

#define CR_IS_SYMMETRIC_LOOP(_XCOND_) \
	do { \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if ((i = pi[k]) >= j) { \
					if (i == j) \
						++pp_[j]; \
					k = kend; \
					break; \
				} \
				if (pp_[i] == pp[i] || pi[pp_[i]] != j || (_XCOND_)) { \
					res = FALSE; \
					goto finish; \
				} \
				++pp_[i]; \
				++pp_[j]; \
				++k; \
			} \
		} \
	} while (0)

/* isSymmetric(<.g[CR]Matrix>, tol = 0, checkDN) */
#define CR_IS_SYMMETRIC(_C_, _I_) \
SEXP _C_ ## sparse_is_symmetric(SEXP obj, SEXP checkDN) \
{ \
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), n = pdim[0], s = pdim[1] == n; \
	UNPROTECT(1); /* dim */ \
	if (!s) \
		return ScalarLogical(0); \
	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)); \
	s = asLogical(checkDN) == 0 || DimNames_is_symmetric(dimnames); \
	UNPROTECT(1); /* dimnames */ \
	if (!s) \
		return ScalarLogical(0); \
	SEXP p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)), \
		i0 = PROTECT(GET_SLOT(obj, Matrix_ ## _I_ ## Sym)); \
	int i, j, k, kend, *pp_, *pp = INTEGER(p0), *pi = INTEGER(i0), \
		nprotect = 2; \
	Rboolean res = TRUE; \
	Matrix_Calloc(pp_, n, int); \
	Matrix_memcpy(pp_, pp, n, sizeof(int)); \
	++pp; \
	/* For all X[i,j] in "leading" triangle, */ \
	/* need that X[j,i] exists and X[j,i] == X[i,j] */ \
	if (!HAS_SLOT(obj, Matrix_xSym)) \
		CR_IS_SYMMETRIC_LOOP(0); \
	else { \
		SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)); \
		++nprotect; \
		switch (TYPEOF(x0)) { \
		case LGLSXP: \
		{ \
			int *px = LOGICAL(x0); \
			CR_IS_SYMMETRIC_LOOP( \
				px[pp_[i]] == NA_LOGICAL \
				? (px[k] != NA_LOGICAL) \
				: (px[k] == NA_LOGICAL || px[pp_[i]] != px[k])); \
			break; \
		} \
		case INTSXP: \
		{ \
			int *px = INTEGER(x0); \
			CR_IS_SYMMETRIC_LOOP( \
				px[pp_[i]] == NA_INTEGER \
				? (px[k] != NA_INTEGER) \
				: (px[k] == NA_INTEGER || px[pp_[i]] != px[k])); \
			break; \
		} \
		case REALSXP: \
		{ \
			double *px = REAL(x0); \
			CR_IS_SYMMETRIC_LOOP( \
				ISNAN(px[pp_[i]]) \
				? !ISNAN(px[k]) \
				: (ISNAN(px[k]) || px[pp_[i]] != px[k])); \
			break; \
		} \
		case CPLXSXP: \
		{ \
			Rcomplex *px = COMPLEX(x0); \
			CR_IS_SYMMETRIC_LOOP( \
				ISNAN(px[pp_[i]].r) || ISNAN(px[pp_[i]].i) \
				? !(ISNAN(px[k].r) || ISNAN(px[k].i)) \
				: (ISNAN(px[k].r) || ISNAN(px[k].i) || \
				   px[pp_[i]].r != px[k].r || px[pp_[i]].i != px[k].i)); \
			break; \
		} \
		default: \
			ERROR_INVALID_TYPE(x0, __func__); \
			break; \
		} \
	} \
	/* Need upper, lower triangles to have same number of nonzero elements */ \
	for (j = 0; j < n; ++j) { \
		if (pp_[j] != pp[j]) { \
			res = FALSE; \
			goto finish; \
		} \
	} \
finish: \
	Matrix_Free(pp_, n); \
	UNPROTECT(nprotect); /* x0, i0, p0 */ \
	return ScalarLogical(res); \
}

/* Csparse_is_symmetric() */
/* FIXME: not checking for real diagonal in complex case */
CR_IS_SYMMETRIC(C, i)
/* Rsparse_is_symmetric() */
/* FIXME: not checking for real diagonal in complex case */
CR_IS_SYMMETRIC(R, j)

/* colSums(<CsparseMatrix>), rowSums(<RsparseMatrix>) */
SEXP CRsparse_colSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse)
{
	static const char *valid[] = { VALID_CSPARSE, VALID_RSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);
	const char *cl = valid[ivalid];
	if (cl[1] == 's')
		return CRsparse_rowSums(obj, narm, mean, sparse);

	int doSparse = asLogical(sparse) != 0,
		doNaRm = asLogical(narm) != 0,
		doMean = asLogical(mean) != 0,
		doCount = doNaRm && doMean;

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int margin = (cl[2] == 'C') ? 1 : 0,
		*pdim = INTEGER(dim), m = pdim[!margin], n = pdim[margin];
	UNPROTECT(1); /* dim */

	char di = 'N';
	if (cl[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */
		if (doSparse && di != 'N')
			warning(_("sparseResult=TRUE inefficient for unit triangular 'x'"));
	}

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	++nprotect;
	int *pp = INTEGER(p) + 1, j, k, kend, count = m;

	PROTECT_INDEX pid;
	SEXP res, vl = NULL, vi = NULL, vx = NULL;
	int *pvi = NULL;

	if (!doSparse) {
		PROTECT_WITH_INDEX(
			res = allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, n), &pid);
		++nprotect;
	} else {
		int nnz = n;
		if (di == 'N') {
			nnz = 0;
			for (j = 0; j < n; ++j)
				if (pp[j-1] < pp[j])
					++nnz;
		}

		char cl_[] = ".sparseVector";
		cl_[0] = (((cl[0] == 'n' || cl[0] == 'l') && !doMean)
		          ? 'i' : ((cl[0] != 'z') ? 'd' : 'z'));
		PROTECT(res = NEW_OBJECT_OF_CLASS(cl_));
		PROTECT(vl = ScalarInteger(n));
		PROTECT(vi = allocVector(INTSXP, nnz));
		PROTECT_WITH_INDEX(
			vx = allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, nnz), &pid);
		nprotect += 4;
		pvi = INTEGER(vi);
	}

	if (cl[0] == 'n') {
		double *pres = (doSparse) ? REAL(vx) : REAL(res);
		if (!doSparse) {
			int u = (di == 'N') ? 0 : 1;
			for (j = 0; j < n; ++j) {
				*pres = pp[j] - pp[j-1] + u;
				if (doMean)
					*pres /= count;
				++pres;
			}
		} else if (di == 'N') {
			for (j = 0; j < n; ++j) {
				if (pp[j-1] < pp[j]) {
					*pvi = j + 1;
					*pres = pp[j] - pp[j-1];
					if (doMean)
						*pres /= count;
					++pvi;
					++pres;
				}
			}
		} else {
			for (j = 0; j < n; ++j) {
				*pvi = j + 1;
				*pres = pp[j] - pp[j-1] + 1;
				if (doMean)
					*pres /= count;
				++pvi;
				++pres;
			}
		}
	} else {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));

#define CR_COLSUMS_LOOP \
		do { \
			k = 0; \
			if (!doSparse) { \
				if (di == 'N') { \
					for (j = 0; j < n; ++j) { \
						kend = pp[j]; \
						DO_INIT(ZERO); \
						while (k < kend) { DO_INCR; ++k; } \
						DO_SCALE; \
						++pres; \
					} \
				} else { \
					for (j = 0; j < n; ++j) { \
						kend = pp[j]; \
						DO_INIT(ONE); \
						while (k < kend) { DO_INCR; ++k; } \
						DO_SCALE; \
						++pres; \
					} \
				} \
			} else { \
				if (di == 'N') { \
					for (j = 0; j < n; ++j) { \
						kend = pp[j]; \
						if (k < kend) { \
							*pvi = j + 1; \
							DO_INIT(ZERO); \
							while (k < kend) { DO_INCR; ++k; } \
							DO_SCALE; \
							++pvi; \
							++pres; \
						} \
					} \
				} else { \
					for (j = 0; j < n; ++j) { \
						kend = pp[j]; \
						*pvi = j + 1; \
						DO_INIT(ONE); \
						while (k < kend) { DO_INCR; ++k; } \
						DO_SCALE; \
						++pvi; \
						++pres; \
					} \
				} \
			} \
		} while (0)

#define CR_COLSUMS(_CTYPE1_, _PTR1_, _CTYPE2_, _PTR2_) \
		do { \
			_CTYPE1_ *pres = (doSparse) ? _PTR1_(vx) : _PTR1_(res); \
			_CTYPE2_ *px   = _PTR2_(x); \
			CR_COLSUMS_LOOP; \
		} while (0)

		switch (cl[0]) {
		case 'l':

#define ZERO         0.0
#define ONE          1.0
#define DO_INIT(_U_) \
			do { \
				*pres = _U_; \
				if (doCount) \
					count = m; \
			} while (0)
#define DO_INCR \
			do { \
				if (px[k] != NA_LOGICAL) { \
					if (px[k]) *pres += 1.0; \
				} else if (!doNaRm) \
					*pres = NA_REAL; \
				else if (doMean) \
					--count; \
			} while (0)
#define DO_SCALE     if (doMean) *pres /= count

			CR_COLSUMS(double, REAL, int, LOGICAL);
			break;

#undef DO_INCR

		case 'i':

#define DO_INCR \
			do { \
				if (px[k] != NA_INTEGER) \
					*pres += px[k]; \
				else if (!doNaRm) \
					*pres = NA_REAL; \
				else if (doMean) \
					--count; \
			} while (0)

			CR_COLSUMS(double, REAL, int, INTEGER);
			break;

#undef DO_INCR

		case 'd':

#define DO_INCR \
			do { \
				if (!(doNaRm && ISNAN(px[k]))) \
					*pres += px[k]; \
				else if (doMean) \
					--count; \
			} while (0)

			CR_COLSUMS(double, REAL, double, REAL);
			break;

#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_SCALE

		case 'z':

#define ZERO         Matrix_zzero
#define ONE          Matrix_zone
#define DO_INCR \
			do { \
				if (!(doNaRm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) { \
					(*pres).r += px[k].r; \
					(*pres).i += px[k].i; \
				} else if (doMean) \
					--count; \
			} while (0)
#define DO_SCALE \
			do { \
				if (doMean) { \
					(*pres).r /= count; \
					(*pres).i /= count; \
				} \
			} while (0)

			CR_COLSUMS(Rcomplex, COMPLEX, Rcomplex, COMPLEX);
			break;

#undef ZERO
#undef ONE
#undef DO_INIT
#undef DO_INCR
#undef DO_SCALE

		default:
			break;
		}

#undef CR_COLSUMS
#undef CR_COLSUMS_LOOP

		UNPROTECT(1); /* x */
	}

	if (doSparse) {
		if ((cl[0] == 'n' || cl[0] == 'l') && !doMean)
			REPROTECT(vx = coerceVector(vx, INTSXP), pid);

		SET_SLOT(res, Matrix_lengthSym, vl);
		SET_SLOT(res, Matrix_iSym,      vi);
		SET_SLOT(res, Matrix_xSym,      vx);
	} else {
		if ((cl[0] == 'n' || cl[0] == 'l') && !doMean)
			REPROTECT(res = coerceVector(res, INTSXP), pid);

		SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
			nms = VECTOR_ELT(dimnames, margin);
		if (!isNull(nms))
			setAttrib(res, R_NamesSymbol, nms);
		UNPROTECT(1); /* dimnames */
	}

	UNPROTECT(nprotect);
	return res;
}

/* rowSums(<CsparseMatrix>), colSums(<RsparseMatrix>) */
SEXP CRsparse_rowSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse)
{
	static const char *valid[] = { VALID_CSPARSE, VALID_RSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);
	const char *cl = valid[ivalid];

	int doSparse = asLogical(sparse) != 0,
		doNaRm = asLogical(narm) != 0,
		doMean = asLogical(mean) != 0;

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int margin = (cl[2] == 'C') ? 0 : 1,
		*pdim = INTEGER(dim), m = pdim[margin], n = pdim[!margin];
	UNPROTECT(1); /* dim */

	char ul = 'U', di = 'N';
	if (cl[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		if (cl[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
			if (doSparse && di != 'N')
				warning(_("sparseResult=TRUE inefficient for unit triangular 'x'"));
		}
	}

	SEXP iSym = (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym,
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj,        iSym));
	nprotect += 2;
	int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k, kend, *pcount = NULL;

	SEXP x = NULL;
	if (cl[0] != 'n') {
		PROTECT(x = GET_SLOT(obj, Matrix_xSym));
		++nprotect;
	}

	PROTECT_INDEX pid;
	SEXP res;
	PROTECT_WITH_INDEX(
		res = allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, m), &pid);
	++nprotect;

#define CR_ROWSUMS_LOOP \
	do { \
		k = 0; \
		if (cl[1] != 's') { \
			for (j = 0; j < n; ++j) { \
				kend = pp[j]; \
				while (k < kend) { DO_INCR; ++k; } \
			} \
		} else if (ul == ((cl[2] == 'C') ? 'U' : 'L')) { \
			for (j = 0; j < n; ++j) { \
				kend = pp[j]; \
				if (k < kend) { \
					while (kend - k > 1) { DO_INCR_SYMM; ++k; } \
					if (pi[k] == j) \
						DO_INCR; \
					else \
						DO_INCR_SYMM; \
					++k; \
				} \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				kend = pp[j]; \
				if (k < kend) { \
					if (pi[k] == j) \
						DO_INCR; \
					else \
						DO_INCR_SYMM; \
					++k; \
					while (k < kend) { DO_INCR_SYMM; ++k; } \
				} \
			} \
		} \
	} while (0)

#define CR_ROWSUMS_X(_CTYPE1_, _PTR1_, _CTYPE2_, _PTR2_) \
	do { \
		_CTYPE2_ *px = _PTR2_(x); \
		CR_ROWSUMS_N(_CTYPE1_, _PTR1_); \
	} while (0)

#define CR_ROWSUMS_N(_CTYPE1_, _PTR1_) \
	do { \
		_CTYPE1_ *pres = _PTR1_(res), u = (di == 'N') ? ZERO : ONE; \
		if (doNaRm && doMean && cl[0] != 'n') { \
			Matrix_Calloc(pcount, m, int); \
			for (k = 0; k < m; ++k) { \
				pres[k] = u; \
				pcount[k] = n; \
			} \
		} else { \
			for (k = 0; k < m; ++k) \
				pres[k] = u; \
		} \
		CR_ROWSUMS_LOOP; \
	} while (0)

	switch (cl[0]) {
	case 'n':

#define ZERO         0.0
#define ONE          1.0
#define DO_INCR      pres[pi[k]] += 1.0
#define DO_INCR_SYMM \
		do { \
			pres[pi[k]] += 1.0; \
			pres[j]     += 1.0; \
		} while (0)

		CR_ROWSUMS_N(double, REAL);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'l':

#define DO_INCR \
		do { \
			if (px[k] != NA_LOGICAL) { \
				if (px[k]) \
					pres[pi[k]] += 1.0; \
			} else if (!doNaRm) \
				pres[pi[k]] = NA_REAL; \
			else if (doMean) \
				--pcount[pi[k]]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (px[k] != NA_LOGICAL) { \
				if (px[k]) { \
					pres[pi[k]] += 1.0; \
					pres[j]     += 1.0; \
				} \
			} else if (!doNaRm) { \
				pres[pi[k]] = NA_REAL; \
				pres[j]     = NA_REAL; \
			} else if (doMean) { \
				--pcount[pi[k]]; \
				--pcount[j]; \
			} \
		} while (0)

		CR_ROWSUMS_X(double, REAL, int, LOGICAL);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'i':

#define DO_INCR \
		do { \
			if (px[k] != NA_INTEGER) \
				pres[pi[k]] += px[k]; \
			else if (!doNaRm) \
				pres[pi[k]] = NA_REAL; \
			else if (doMean) \
				--pcount[pi[k]]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (px[k] != NA_INTEGER) { \
				pres[pi[k]] += px[k]; \
				pres[j]     += px[k]; \
			} else if (!doNaRm) { \
				pres[pi[k]] = NA_REAL; \
				pres[j]     = NA_REAL; \
			} else if (doMean) { \
				--pcount[pi[k]]; \
				--pcount[j]; \
			} \
		} while (0)

		CR_ROWSUMS_X(double, REAL, int, INTEGER);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'd':

#define DO_INCR \
		do { \
			if (!(doNaRm && ISNAN(px[k]))) \
				pres[pi[k]] += px[k]; \
			else if (doMean) \
				--pcount[pi[k]]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (!(doNaRm && ISNAN(px[k]))) { \
				pres[pi[k]] += px[k]; \
				pres[j]     += px[k]; \
			} else if (doMean) { \
				--pcount[pi[k]]; \
				--pcount[j]; \
			} \
		} while (0)

		CR_ROWSUMS_X(double, REAL, double, REAL);
		break;

#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_INCR_SYMM

	case 'z':

#define ZERO         Matrix_zzero
#define ONE          Matrix_zone
#define DO_INCR \
		do { \
			if (!(doNaRm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) { \
				pres[pi[k]].r += px[k].r; \
				pres[pi[k]].i += px[k].i; \
			} else if (doMean) \
				--pcount[pi[k]]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (!(doNaRm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) { \
				pres[pi[k]].r += px[k].r; \
				pres[pi[k]].i += px[k].i; \
				pres[j].r     += px[k].r; \
				pres[j].i     += px[k].i; \
			} else if (doMean) { \
				--pcount[pi[k]]; \
				--pcount[j]; \
			} \
		} while (0)

		CR_ROWSUMS_X(Rcomplex, COMPLEX, Rcomplex, COMPLEX);
		break;

#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_INCR_SYMM

	default:
		break;
	}

#undef CR_ROWSUMS
#undef CR_ROWSUMS_LOOP

	if (doMean) {
		if (cl[0] != 'z') {
			double *pres = REAL(res);
			if (doNaRm && cl[0] != 'n') {
				for (k = 0; k < m; ++k)
					pres[k] /= pcount[k];
				Matrix_Free(pcount, m);
			} else {
				for (k = 0; k < m; ++k)
					pres[k] /= n;
			}
		} else {
			Rcomplex *pres = COMPLEX(res);
			if (doNaRm) {
				for (k = 0; k < m; ++k) {
					pres[k].r /= pcount[k];
					pres[k].i /= pcount[k];
				}
				Matrix_Free(pcount, m);
			} else {
				for (k = 0; k < m; ++k) {
					pres[k].r /= n;
					pres[k].i /= n;
				}
			}
		}
	}

	if ((cl[0] == 'n' || cl[0] == 'l') && !doMean)
		REPROTECT(res = coerceVector(res, INTSXP), pid);
	if (doSparse) {
		/* defined in ./sparseVector.c : */
		SEXP v2spV(SEXP);
		REPROTECT(res = v2spV(res), pid);
	} else {
		SEXP dimnames;
		if (cl[1] != 's')
			PROTECT(dimnames = GET_SLOT(obj, Matrix_DimNamesSym));
		else
			PROTECT(dimnames = get_symmetrized_DimNames(obj, -1));
		SEXP nms = VECTOR_ELT(dimnames, margin);
		if (!isNull(nms))
			setAttrib(res, R_NamesSymbol, nms);
		UNPROTECT(1); /* dimnames */
	}

	UNPROTECT(nprotect);
	return res;
}
