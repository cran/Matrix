#include <math.h> /* fabs, hypot */
#include "Mdefines.h"
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
#define DROP0_LOOP1(_CTYPE_, _PTR_, _ISNZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0); \
			for (k = 0; k < nnz0; ++k) { \
				if (_ISNZ_(*px0)) \
					++nnz1; \
				++px0; \
			} \
		} while (0)

		DROP0_CASES(DROP0_LOOP1);
		if (nnz1 == nnz0) {
			UNPROTECT(2); /* p0, x0 */
			return from;
		}
		PROTECT(to = newObject(class));

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
#define DROP0_LOOP2(_CTYPE_, _PTR_, _ISNZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			for (j = 0, k = 0; j < n; ++j) { \
				pp1[j] = pp1[j - 1]; \
				kend = pp0[j]; \
				while (k < kend) { \
					if (_ISNZ_(*px0)) { \
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
		PROTECT(to = newObject(class));

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
#define DROP0_LOOP2(_CTYPE_, _PTR_, _ISNZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			for (k = 0; k < nnz0; ++k) { \
				if (_ISNZ_(*px0)) { \
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
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorsSym, factors);
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
   though that is more general, working for _all_ Matrix
*/
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
   though that is more general, working for _all_ Matrix
*/
SEXP R_sparse_diag_N2U(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_diag_N2U(from, valid[ivalid]);
}

SEXP sparse_band(SEXP from, const char *class, int a, int b)
{
	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */

	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) */
	/* to be triangularMatrix                       */
	if ((m == 0 || n == 0 || (a <= 1 - m && b >= n - 1)) &&
	    (m != n || n > 1 || class[1] == 't'))
		return from;

	int ge = 0, sy = 0, tr = 0;
	ge = m != n || !((tr = a >= 0 || b <= 0 || class[1] == 't') ||
	                 (sy = a == -b && class[1] == 's'));

	char ul0 = 'U', ul1 = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul0 = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		if (class[1] == 't') {
			/* Be fast if band contains entire triangle */
			if ((ul0 == 'U') ? (a <= 0 && b >= n - 1) : (b >= 0 && a <= 1 - m))
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
		ul0 = (ul0 == 'U') ? 'L' : 'U';
		from = sparse_transpose(from, class, 1);
	}
	PROTECT(from);

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (ge) ? 'g' : ((tr) ? 't' : 's');
	cl[2] = (class[2] == 'R') ? 'C' : class[2];
	SEXP to = PROTECT(newObject(cl));

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
		ul1 = (tr && class[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ul0;
		if (ul1 != 'U') {
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
	if (k1 == R_NilValue) // tril()
		a = -m ; // was (m > 0) ? 1 - m : 0;
	else if ((a = asInteger(k1)) == NA_INTEGER || a < -m || a > n)
		error(_("'%s' (%d) must be an integer from %s (%d) to %s (%d)"),
		      "k1", a, "-Dim[1]", -m, "Dim[2]", n);
	if (k2 == R_NilValue) // triu()
		b = n; // was (n > 0) ? n - 1 : 0;
	else if ((b = asInteger(k2)) == NA_INTEGER || b < -m || b > n)
		error(_("'%s' (%d) must be an integer from %s (%d) to %s (%d)"),
		      "k2", b, "-Dim[1]", -m, "Dim[2]", n);
	else if (b < a)
		error(_("'%s' (%d) must be less than or equal to '%s' (%d)"),
		      "k1", a, "k2", b);

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

	SEXP res = PROTECT(allocVector(kindToType(class[0]), r));

#define DG_CASES \
	do { \
		switch (class[0]) { \
		case 'n': \
		case 'l': \
			DG_LOOP(int, LOGICAL, SHOW, FIRSTOF, INCREMENT_LOGICAL, 0, 1); \
			break; \
		case 'i': \
			DG_LOOP(int, INTEGER, SHOW, FIRSTOF, INCREMENT_INTEGER, 0, 1); \
			break; \
		case 'd': \
			DG_LOOP(double, REAL, SHOW, FIRSTOF, INCREMENT_REAL, 0.0, 1.0); \
			break; \
		case 'z': \
			DG_LOOP(Rcomplex, COMPLEX, SHOW, FIRSTOF, INCREMENT_COMPLEX, Matrix_zzero, Matrix_zone); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	if (di != 'N') {

		int j;

#undef DG_LOOP
#define DG_LOOP(_CTYPE_, _PTR_, _MASK_, _REPLACE_, _INCREMENT_, _ZERO_, _ONE_) \
		do { \
			_CTYPE_ *pres = _PTR_(res); \
			for (j = 0; j < r; ++j) \
				*(pres++) = _ONE_; \
		} while (0)

		DG_CASES;

	} else if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(obj,        iSym));
		int j, k, kend, *pp0 = INTEGER(p0), *pi0 = INTEGER(i0);
		pp0++;

#undef DG_LOOP
#define DG_LOOP(_CTYPE_, _PTR_, _MASK_, _REPLACE_, _INCREMENT_, _ZERO_, _ONE_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0)); \
			_CTYPE_ *pres = _PTR_(res); \
			if (class[1] == 'g') { \
				for (j = 0, k = 0; j < r; ++j) { \
					pres[j] = _ZERO_; \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] != j) \
							++k; \
						else { \
							pres[j] = _REPLACE_(px0[k], 1); \
							k = kend; \
						} \
					} \
				} \
			} else if ((class[2] == 'C') == (ul != 'U')) { \
				for (j = 0, k = 0; j < r; ++j) { \
					kend = pp0[j]; \
					pres[j] = (k < kend && pi0[k] == j) \
						? _REPLACE_(px0[k], 1) : _ZERO_; \
					k = kend; \
				} \
			} else { \
				for (j = 0, k = 0; j < r; ++j) { \
					kend = pp0[j]; \
					pres[j] = (k < kend && pi0[kend - 1] == j) \
						? _REPLACE_(px0[kend - 1], 1) : _ZERO_; \
					k = kend; \
				} \
			} \
		} while (0)

		if (class[0] == 'n') {
			DG_LOOP(int, LOGICAL, HIDE, SECONDOF, , 0, 1);
		} else {
			SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym));
			DG_CASES;
			UNPROTECT(1); /* x0 */
		}

		UNPROTECT(2); /* i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0);

#undef DG_LOOP
#define DG_LOOP(_CTYPE_, _PTR_, _MASK_, _REPLACE_, _INCREMENT_, _ZERO_, _ONE_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0)); \
			_CTYPE_ *pres = _PTR_(res); \
			Matrix_memset(pres, 0, r, sizeof(_CTYPE_)); \
			for (k = 0; k < nnz0; ++k) { \
				if (*pi0 == *pj0) \
					_INCREMENT_(pres[*pi0], (*px0)); \
				++pi0; ++pj0; _MASK_(++px0); \
			} \
		} while (0)

		if (class[0] == 'n')
			DG_LOOP(int, LOGICAL, HIDE, SECONDOF, INCREMENT_PATTERN, 0, 1);
		else {
			SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym));
			DG_CASES;
			UNPROTECT(1); /* x0 */
		}

		UNPROTECT(2); /* j0, i0 */

	}

#undef DG_CASES
#undef DG_LOOP

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
			         (rn == cn || equal_character_vectors(rn, cn, r)))
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
	SEXP to = PROTECT(newObject(class));
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

#define DS_CASES \
	do { \
		switch (class[0]) { \
		case 'n': \
		case 'l': \
			DS_LOOP(int, LOGICAL, SHOW, ISNZ_LOGICAL); \
			break; \
		case 'i': \
			DS_LOOP(int, INTEGER, SHOW, ISNZ_INTEGER); \
			break; \
		case 'd': \
			DS_LOOP(double, REAL, SHOW, ISNZ_REAL); \
			break; \
		case 'z': \
			DS_LOOP(Rcomplex, COMPLEX, SHOW, ISNZ_COMPLEX); \
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

#undef DS_LOOP
#define DS_LOOP(_CTYPE_, _PTR_, _MASK_, _ISNZ_) \
		do { \
			_CTYPE_ *pvalue = _PTR_(value); \
			if (v) { \
				for (j = 0; j < r; ++j) { \
					if (_ISNZ_(pvalue[j])) \
						++nd1; \
					pp1[j] += nd1; \
				} \
				for (j = r; j < n_; ++j) \
					pp1[j] += nd1; \
			} else if (_ISNZ_(pvalue[0])) { \
				full = 1; \
				for (j = 0; j < r; ++j) \
					pp1[j] += ++nd1; \
				for (j = r; j < n_; ++j) \
					pp1[j] +=   nd1; \
			} \
		} while (0)

		DS_CASES;

		if (nd1 - nd0 > INT_MAX - pp0[n_ - 1])
			error(_("%s cannot exceed %s"), "p[length(p)]", "2^31-1");

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[n_ - 1]));
		SET_SLOT(to, iSym, i1);
		int *pi1 = INTEGER(i1);

#undef DS_LOOP
#define DS_LOOP(_CTYPE_, _PTR_, _MASK_, _ISNZ_) \
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
				if ((v) ? _ISNZ_(pvalue[j]) : full) { \
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
			DS_LOOP(int, LOGICAL, HIDE, ISNZ_LOGICAL);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), pp1[n_ - 1]));
			SET_SLOT(to, Matrix_xSym, x1);
			DS_CASES;
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

#undef DS_LOOP
#define DS_LOOP(_CTYPE_, _PTR_, _MASK_, _ISNZ_) \
		do { \
			_CTYPE_ *pvalue = _PTR_(value); \
			if (v) { \
				for (j = 0; j < r; ++j) \
					if (_ISNZ_(pvalue[j])) \
						++nd1; \
			} else if (_ISNZ_(pvalue[0])) { \
				full = 1; \
				nd1 = r; \
			} \
		} while (0)

		DS_CASES;

		if (nd1 - nd0 > R_XLEN_T_MAX - nnz0)
			error(_("%s cannot exceed %s"), "length(i)", "R_XLEN_T_MAX");
		nnz1 += nd1 - nd0;

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#undef DS_LOOP
#define DS_LOOP(_CTYPE_, _PTR_, _MASK_, _ISNZ_) \
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
					if (_ISNZ_(pvalue[j])) { \
						*(pi1++) = *(pj1++) = j; \
						_MASK_(*(px1++) = pvalue[j]); \
					} \
				} \
			} else if (full) { \
				for (j = 0; j < r; ++j) { \
					*(pi1++) = *(pj1++) = j; \
					_MASK_(*(px1++) = pvalue[0]); \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			DS_LOOP(int, LOGICAL, HIDE, ISNZ_LOGICAL);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			SET_SLOT(to, Matrix_xSym, x1);
			DS_CASES;
			UNPROTECT(2); /* x1, x0 */
		}

		UNPROTECT(4); /* j1, i1, j0, i0 */

	}

#undef DS_CASES
#undef DS_LOOP

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

	SEXPTYPE tx = kindToType(class[0]), tv = TYPEOF(value);

	switch (tv) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
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
#ifndef MATRIX_ENABLE_IMATRIX
		if (tv == INTSXP) {
		PROTECT(from = sparse_as_kind(from, class, 'd'));
		PROTECT(value = coerceVector(value, REALSXP));
		} else {
#endif
		PROTECT(from = sparse_as_kind(from, class, typeToKind(tv)));
		PROTECT(value);
#ifndef MATRIX_ENABLE_IMATRIX
		}
#endif
		class = valid[R_check_class_etc(from, valid)];
	}

	from = sparse_diag_set(from, class, value);
	UNPROTECT(2);
	return from;
}

SEXP sparse_transpose(SEXP from, const char *class, int lazy)
{
	SEXP to;
	if (class[2] == 'T' || !lazy)
		PROTECT(to = newObject(class));
	else {
		char cl[] = "...Matrix";
		cl[0] = class[0];
		cl[1] = class[1];
		cl[2] = (class[2] == 'C') ? 'R' : 'C';
		PROTECT(to = newObject(cl));
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
			SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
			if (LENGTH(factors) > 0)
				SET_SLOT(to, Matrix_factorsSym, factors);
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
		error(_("invalid '%s' to '%s'"), "lazy", __func__);

	return sparse_transpose(from, valid[ivalid], lazy_);
}

SEXP sparse_force_symmetric(SEXP from, const char *class, char ul)
{
	char ul0 = 'U', ul1 = 'U';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul0 = ul1 = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */
	}
	if (ul != '\0')
		ul1 = ul;

	if (class[1] == 's') {
		/* .s[CRT]Matrix */
		if (ul0 == ul1)
			return from;
		SEXP to = PROTECT(sparse_transpose(from, class, 0));
		if (class[0] == 'z') {
			/* Need _conjugate_ transpose */
			SEXP x1 = PROTECT(GET_SLOT(to, Matrix_xSym));
			conjugate(x1);
			UNPROTECT(1); /* x1 */
		}
		UNPROTECT(1) /* to */;
		return to;
	}

	/* Now handling just .[gt][CRT]Matrix ... */

	char cl[] = ".s.Matrix";
	cl[0] = class[0];
	cl[2] = class[2];
	SEXP to = PROTECT(newObject(cl));

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

	if (ul1 != 'U') {
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

	if (class[1] == 't' && di == 'N' && ul0 == ul1) {

		/* No need to allocate in this case: we have the triangle we want */
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
		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			SET_SLOT(to, Matrix_xSym, x0);
			UNPROTECT(1); /* x0 */
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
				if (ul0 != ul1) {
					/* Returning identity matrix */
					for (j = 0; j < n; ++j)
						pp1[j] = ++nnz1;
				} else {
					/* Returning symmetric matrix with unit diagonal */
					for (j = 0; j < n; ++j)
						pp1[j] = ++nnz1 + pp0[j];
					nnz1 += nnz0;
				}
			} else if (ul0 == ((class[2] == 'C') ? 'U' : 'L')) {
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
		} else if (ul1 == ((class[2] == 'C') ? 'U' : 'L')) {
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
					if (ul0 != ul1) { \
						/* Returning identity matrix */ \
						for (j = 0; j < n; ++j) { \
							*(pi1++) = j; \
							_MASK_(*(px1++) = _ONE_); \
						} \
					} else if (ul0 == ((class[2] == 'C') ? 'U' : 'L')) { \
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
				} else if (ul0 == ((class[2] == 'C') ? 'U' : 'L')) { \
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
			} else if (ul1 == ((class[2] == 'C') ? 'U' : 'L')) { \
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
			nnz1 = (ul0 == ul1) ? n + nnz0 : n;
		else {
			if (ul1 == 'U') {
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
				if (ul0 == ul1) { \
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
				if (ul1 == 'U') { \
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
			error(_("invalid '%s' to '%s'"), "uplo", __func__);
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
	SEXP to = PROTECT(newObject(cl));

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
				x1 = PROTECT(allocVector(kindToType(cl[0]), pp1[n - 1]));
			int *pi1 = INTEGER(i1);

#undef SP_LOOP
#define SP_LOOP(_CTYPE_, _PTR_, _ASSIGN_, _INCREMENT_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1), \
					*px0_ = _PTR_(x0_); \
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

			if (cl[0] == 'd')
				SP_LOOP(double, REAL, ASSIGN_REAL, INCREMENT_REAL);
			else
				SP_LOOP(Rcomplex, COMPLEX, ASSIGN_COMPLEX, INCREMENT_COMPLEX);

			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to,        iSym, i1);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(6); /* x1, i1, p1, x0_, i0_, p0_ */

		} else if (class[1] == 't') {

			int leading = (class[2] == 'C') == (ul != 'U');

			if (di == 'N') {

				SEXP x1 = PROTECT(allocVector(kindToType(cl[0]), nnz));

#undef SP_LOOP
#define SP_LOOP(_CTYPE_, _PTR_, _ASSIGN_) \
				do { \
					_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
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

				if (cl[0] == 'd')
					SP_LOOP(double, REAL, ASSIGN_REAL);
				else
					SP_LOOP(Rcomplex, COMPLEX, ASSIGN_COMPLEX);

				SET_SLOT(to, Matrix_pSym, p0);
				SET_SLOT(to,        iSym, i0);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */

			} else {

				nnz += n;
				SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
					i1 = PROTECT(allocVector(INTSXP, nnz)),
					x1 = PROTECT(allocVector(kindToType(cl[0]), nnz));
				int *pp1 = INTEGER(p1), *pi1 = INTEGER(i1);
				*(pp1++) = 0;

#undef SP_LOOP
#define SP_LOOP(_CTYPE_, _PTR_, _ASSIGN_, _ONE_) \
				do { \
					_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
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

				if (cl[0] == 'd')
					SP_LOOP(double, REAL, ASSIGN_REAL, 1.0);
				else
					SP_LOOP(Rcomplex, COMPLEX, ASSIGN_COMPLEX, Matrix_zone);

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
				/* Symmetric part of Hermitian matrix is real part */
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
				x1 = PROTECT(allocVector(kindToType(cl[0]), nnz));
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#undef SP_LOOP
#define SP_LOOP(_CTYPE_, _PTR_, _ASSIGN_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
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

			if (cl[0] == 'd')
				SP_LOOP(double, REAL, ASSIGN_REAL);
			else
				SP_LOOP(Rcomplex, COMPLEX, ASSIGN_COMPLEX);

			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(3); /* x1, j1, i1 */

		} else if (class[1] == 't') {

			if (di == 'N') {

				SEXP x1 = PROTECT(allocVector(kindToType(cl[0]), nnz));

#undef SP_LOOP
#define SP_LOOP(_CTYPE_, _PTR_, _ASSIGN_) \
				do { \
					_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
					for (k = 0; k < nnz; ++k) { \
						if (*pi0 == *pj0) \
							*px1 = *px0; \
						else \
							_ASSIGN_((*px1), 0.5 * (*px0)); \
						++px0; ++px1; \
					} \
				} while (0)

				if (cl[0] == 'd')
					SP_LOOP(double, REAL, ASSIGN_REAL);
				else
					SP_LOOP(Rcomplex, COMPLEX, ASSIGN_COMPLEX);

				SET_SLOT(to, Matrix_iSym, i0);
				SET_SLOT(to, Matrix_jSym, j0);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */

			} else {

				SEXP i1 = PROTECT(allocVector(INTSXP, nnz + n)),
					j1 = PROTECT(allocVector(INTSXP, nnz + n)),
					x1 = PROTECT(allocVector(kindToType(cl[0]), nnz + n));
				int j, *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#undef SP_LOOP
#define SP_LOOP(_CTYPE_, _PTR_, _ASSIGN_, _ONE_) \
				do { \
					_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
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

				if (cl[0] == 'd')
					SP_LOOP(double, REAL, ASSIGN_REAL, 1.0);
				else
					SP_LOOP(Rcomplex, COMPLEX, ASSIGN_COMPLEX, Matrix_zone);

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
				/* Symmetric part of Hermitian matrix is real part */
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
	SEXP to = PROTECT(newObject(cl));

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
			x1 = PROTECT(allocVector(kindToType(cl[0]), pp1[n - 1]));
		int *pi1 = INTEGER(i1);

#undef SP_LOOP
#define SP_LOOP(_CTYPE_, _PTR_, _ASSIGN_, _INCREMENT_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1), \
				*px0_ = _PTR_(x0_); \
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

		if (cl[0] == 'd')
			SP_LOOP(double, REAL, ASSIGN_REAL, INCREMENT_REAL);
		else
			SP_LOOP(Rcomplex, COMPLEX, ASSIGN_COMPLEX, INCREMENT_COMPLEX);

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
			x1 = PROTECT(allocVector(kindToType(cl[0]), nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#undef SP_LOOP
#define SP_LOOP(_CTYPE_, _PTR_, _ASSIGN_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
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

		if (cl[0] == 'd')
			SP_LOOP(double, REAL, ASSIGN_REAL);
		else
			SP_LOOP(Rcomplex, COMPLEX, ASSIGN_COMPLEX);

		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		SET_SLOT(to, Matrix_xSym, x1);
		UNPROTECT(6); /* x1, j1, i1, x0, j0, i0 */

	}

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

int sparse_is_symmetric(SEXP obj, const char *class, int checkDN)
{
	if (class[1] == 's')
		return 1;

	if (checkDN) {
		SEXP dimnames = GET_SLOT(obj, Matrix_DimNamesSym);
		if (!DimNames_is_symmetric(dimnames))
			return 0;
	}

	if (class[1] == 't')
		return sparse_is_diagonal(obj, class);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return 1;

	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		obj = sparse_as_Csparse(obj, class);
	}
	PROTECT(obj);

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(obj,        iSym));
	int i, j, k, kend, *pp_, *pp0 = INTEGER(p0) + 1, *pi0 = INTEGER(i0);
	Matrix_Calloc(pp_, n, int);
	Matrix_memcpy(pp_, pp0 - 1, n, sizeof(int));

	int ans = 0;

#define IS_LOOP(_CTYPE_, _PTR_, _MASK_, _NOTREAL_, _NOTCONJ_) \
	do { \
		_MASK_(_CTYPE_ *px0 = _PTR_(x0)); \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp0[j]; \
			while (k < kend) { \
				i = pi0[k]; \
				if (i >= j) { \
					if (i == j) { \
						if (_NOTREAL_(px0[k])) \
							goto finish; \
						++pp_[j]; \
					} \
					k = kend; \
				} else { \
					if (pp_[i] == pp0[i] || pi0[pp_[i]] != j || \
						_NOTCONJ_(px0[k], px0[pp_[i]])) \
						goto finish; \
					++pp_[i]; \
					++pp_[j]; \
					++k; \
				} \
			} \
		} \
	} while (0)

#undef NOTCONJ_PATTERN
#define NOTCONJ_PATTERN(_X_, _Y_) 0

	/* For all X[i,j], i >= j, we require:
	     o  that X[j,i] exists
	     o  that X[j,i] == Conj(X[i,j])
	     o  that Im(X[j,i]) == 0 if i == j
	*/
	if (class[0] == 'n')
		IS_LOOP(int, LOGICAL, HIDE, NOTREAL_PATTERN, NOTCONJ_PATTERN);
	else {
		SEXP x0 = GET_SLOT(obj, Matrix_xSym);
		switch (class[0]) {
		case 'l':
			IS_LOOP(int, LOGICAL, SHOW, NOTREAL_LOGICAL, NOTCONJ_LOGICAL);
			break;
		case 'i':
			IS_LOOP(int, INTEGER, SHOW, NOTREAL_INTEGER, NOTCONJ_INTEGER);
			break;
		case 'd':
			IS_LOOP(double, REAL, SHOW, NOTREAL_REAL, NOTCONJ_REAL);
			break;
		case 'z':
			IS_LOOP(Rcomplex, COMPLEX, SHOW, NOTREAL_COMPLEX, NOTCONJ_COMPLEX);
			break;
		default:
			break;
		}
	}

	/* We further require that the upper and lower triangles
	   have the same number of entries ...
	*/
	for (j = 0; j < n; ++j)
		if (pp_[j] != pp0[j])
			goto finish;

	ans = 1;

finish:
	Matrix_Free(pp_, n);
	UNPROTECT(3); /* i0, p0, obj */

#undef IS_LOOP

	return ans;
}

/* isSymmetric(<[CRT]sparseMatrix>, checkDN, tol = 0)
   NB: requires symmetric nonzero pattern
   TODO: support 'tol', 'scale' arguments and bypass all.equal ??
*/
SEXP R_sparse_is_symmetric(SEXP obj, SEXP checkDN)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int checkDN_;
	if (TYPEOF(checkDN) != LGLSXP || LENGTH(checkDN) < 1 ||
	    (checkDN_ = LOGICAL(checkDN)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "checkDN", "TRUE", "FALSE");

	return ScalarLogical(sparse_is_symmetric(obj, valid[ivalid], checkDN_));
}

int sparse_is_triangular(SEXP obj, const char *class, int upper)
{
	if (class[1] == 't') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (upper == NA_LOGICAL || (upper != 0) == (ul == 'U'))
			return (ul == 'U') ? 1 : -1;
		else if (sparse_is_diagonal(obj, class))
			return (ul == 'U') ? -1 : 1;
		else
			return 0;
	}

	if (class[1] == 's') {
		if (!sparse_is_diagonal(obj, class))
			return 0;
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (upper == NA_LOGICAL)
			return (ul == 'U') ? 1 : -1;
		else
			return (upper != 0) ? 1 : -1;
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return (upper != 0) ? 1 : -1;

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(obj,        iSym));
		UNPROTECT(2); /* i0, p0 */
		int j, k, kend, *pp0 = INTEGER(p0) + 1, *pi0 = INTEGER(i0);

		if (upper == NA_LOGICAL) {
			/* Examine  last entry in each "column" */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[kend - 1] > j)
					break;
				k = kend;
			}
			if (j == n)
				return (class[2] == 'C') ? 1 : -1;
			/* Examine first entry in each "column" */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[k] < j)
					break;
				k = kend;
			}
			if (j == n)
				return (class[2] == 'C') ? -1 : 1;
			return 0;
		} else if ((class[2] == 'C') == (upper != 0)) {
			/* Examine  last entry in each "column" */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[kend - 1] > j)
					return 0;
				k = kend;
			}
			return (class[2] == 'C') ? 1 : -1;
		} else {
			/* Examine first entry in each "column" */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[k] < j)
					return 0;
				k = kend;
			}
			return (class[2] == 'C') ? -1 : 1;
		}

	} else {

		SEXP i0 = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(obj, Matrix_jSym));
		UNPROTECT(2); /* i0, j0 */
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0);

		if (upper == NA_LOGICAL) {
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] > pj0[k])
					break;
			if (k == nnz0)
				return  1;
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] < pj0[k])
					break;
			if (k == nnz0)
				return -1;
			return 0;
		} else if (upper != 0) {
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] > pj0[k])
					return 0;
			return  1;
		} else {
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] < pj0[k])
					return 0;
			return -1;
		}

	}
}

/* isTriangular(<[CRT]sparseMatrix>, upper)
   NB: requires triangular nonzero pattern
*/
SEXP R_sparse_is_triangular(SEXP obj, SEXP upper)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	if (TYPEOF(upper) != LGLSXP || LENGTH(upper) < 1)
		error(_("'%s' must be %s or %s or %s"), "upper", "TRUE", "FALSE", "NA");
	int upper_ = LOGICAL(upper)[0];

	int ans_ = sparse_is_triangular(obj, valid[ivalid], upper_);
	SEXP ans = allocVector(LGLSXP, 1);
	LOGICAL(ans)[0] = ans_ != 0;
	if (upper_ == NA_LOGICAL && ans_ != 0) {
		PROTECT(ans);
		static
		SEXP kindSym = NULL;
		SEXP kindVal = PROTECT(mkString((ans_ > 0) ? "U" : "L"));
		if (!kindSym) kindSym = install("kind");
		setAttrib(ans, kindSym, kindVal);
		UNPROTECT(2);
	}
	return ans;
}

int sparse_is_diagonal(SEXP obj, const char *class)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return 1;

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(obj,        iSym));
		UNPROTECT(2); /* i0, p0 */
		int j, k, kend, *pp0 = INTEGER(p0) + 1, *pi0 = INTEGER(i0);

		for (j = 0, k = 0; j < n; ++j) {
			kend = pp0[j];
			if (kend - k > 1 || (kend - k == 1 && pi0[k] != j))
				return 0;
			k = kend;
		}
		return 1;

	} else {

		SEXP i0 = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(obj, Matrix_jSym));
		UNPROTECT(2); /* i0, j0 */
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0);

		for (k = 0; k < nnz0; ++k)
			if (*(pi0++) != *(pj0++))
				return 0;
		return 1;

	}
}

/* isDiagonal(<[CRT]sparseMatrix>)
   NB: requires diagonal nonzero pattern
*/
SEXP R_sparse_is_diagonal(SEXP obj)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	return ScalarLogical(sparse_is_diagonal(obj, valid[ivalid]));
}

#define   MAP(_I_) work[_I_]
#define NOMAP(_I_)      _I_

#define CAST_PATTERN(_X_)   1
#define CAST_LOGICAL(_X_) (_X_ != 0)
#define CAST_INTEGER(_X_)  _X_
#define CAST_REAL(_X_)     _X_
#define CAST_COMPLEX(_X_)  _X_

#define SUM_CASES(_MAP_) \
do { \
	if (class[0] == 'n') { \
		if (mean) \
		SUM_LOOP(int, LOGICAL, double, REAL, HIDE, \
		         0.0, 1.0, NA_REAL, ISNA_PATTERN, \
		         _MAP_, CAST_PATTERN, INCREMENT_REAL, SCALE2_REAL); \
		else \
		SUM_LOOP(int, LOGICAL, int, INTEGER, HIDE, \
		         0, 1, NA_INTEGER, ISNA_PATTERN, \
		         _MAP_, CAST_PATTERN, INCREMENT_INTEGER, SCALE2_REAL); \
	} else { \
		SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)); \
		switch (class[0]) { \
		case 'l': \
			if (mean) \
			SUM_LOOP(int, LOGICAL, double, REAL, SHOW, \
			         0.0, 1.0, NA_REAL, ISNA_LOGICAL, \
			         _MAP_, CAST_LOGICAL, INCREMENT_REAL, SCALE2_REAL); \
			else \
			SUM_LOOP(int, LOGICAL, int, INTEGER, SHOW, \
			         0, 1, NA_INTEGER, ISNA_LOGICAL, \
			         _MAP_, CAST_LOGICAL, INCREMENT_INTEGER, SCALE2_REAL); \
			break; \
		case 'i': \
			SUM_LOOP(int, INTEGER, double, REAL, SHOW, \
			         0.0, 1.0, NA_REAL, ISNA_INTEGER, \
			         _MAP_, CAST_INTEGER, INCREMENT_REAL, SCALE2_REAL); \
			break; \
		case 'd': \
			SUM_LOOP(double, REAL, double, REAL, SHOW, \
			         0.0, 1.0, NA_REAL, ISNA_REAL, \
			         _MAP_, CAST_REAL, INCREMENT_REAL, SCALE2_REAL); \
			break; \
		case 'z': \
			SUM_LOOP(Rcomplex, COMPLEX, Rcomplex, COMPLEX, SHOW, \
			         Matrix_zzero, Matrix_zone, Matrix_zna, ISNA_COMPLEX, \
			         _MAP_, CAST_COMPLEX, INCREMENT_COMPLEX, SCALE2_COMPLEX); \
			break; \
		default: \
			break; \
		} \
		UNPROTECT(1); /* x0 */ \
	} \
} while (0)

#define SUM_TYPEOF(c) (c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

static
void Csparse_colsum(SEXP obj, const char *class,
                    int m, int n, char di, int narm, int mean,
                    SEXP res)
{
	int narm_ = narm && mean && class[0] != 'n';

	SEXP p0 = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp0 = INTEGER(p0) + 1, j, k, kend, nnz1 = n, count = -1;

	if (isS4(res)) {

		if (di == 'N') {
			nnz1 = 0;
			for (j = 0; j < n; ++j)
				if (pp0[j - 1] < pp0[j])
					++nnz1;
		}

		SEXP
		j1 = PROTECT(allocVector(INTSXP, nnz1)),
		x1 = PROTECT(allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(res, Matrix_iSym, j1);
		SET_SLOT(res, Matrix_xSym, x1);

		int *pj1 = INTEGER(j1);
		if (di != 'N')
			for (j = 0; j < n; ++j)
				*(pj1++) = j + 1;
		else
			for (j = 0; j < n; ++j)
				if (pp0[j - 1] < pp0[j])
					*(pj1++) = j + 1;

#define SUM_LOOP(_CTYPE0_, _PTR0_, _CTYPE1_, _PTR1_, _MASK_, \
		         _ZERO_, _ONE_, _NA_, _ISNA_, \
		         _MAP_, _CAST_, _INCREMENT_, _SCALE2_) \
		do { \
			_MASK_(_CTYPE0_ *px0 = _PTR0_(x0)); \
			       _CTYPE1_ *px1 = _PTR1_(x1) , tmp; \
			for (j = 0, k = 0; j < n; ++j) { \
				kend = pp0[j]; \
				if (k < kend || nnz1 == n) { \
					*px1 = (di != 'N') ? _ONE_ : _ZERO_; \
					if (mean) \
						count = m; \
					while (k < kend) { \
						if (_ISNA_(*px0)) { \
							if (!narm) \
								*px1 = _NA_; \
							else if (narm_) \
								--count; \
						} else { \
							tmp = _CAST_(*px0); \
							_INCREMENT_((*px1), tmp); \
						} \
						_MASK_(++px0); \
						++k; \
					} \
					if (mean) \
						_SCALE2_((*px1), count); \
					++px1; \
				} \
			} \
		} while (0)

		SUM_CASES(MAP);
		UNPROTECT(2); /* x1, j1 */

	} else {

		SEXP x1 = res;
		SUM_CASES(NOMAP);

	}

#undef SUM_LOOP

	UNPROTECT(1); /* p0 */
	return;
}

static
void Csparse_rowsum(SEXP obj, const char *class,
                    int m, int n, char di, int narm, int mean,
                    SEXP res, SEXP iSym)
{
	int narm_ = narm && mean && class[0] != 'n';

	SEXP p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(obj, iSym));
	int *pp0 = INTEGER(p0) + 1, *pi0 = INTEGER(i0), i, j, k, kend,
		nnz0 = pp0[n - 1], nnz1 = m;

	if (isS4(res)) {

		int *work;
		Matrix_Calloc(work, m, int);
		if (di != 'N')
			for (i = 0; i < m; ++i)
				work[i] = i;
		else {
			if (class[1] != 's') {
				for (k = 0; k < nnz0; ++k)
					++work[pi0[k]];
			} else {
				for (j = 0, k = 0; j < n; ++j) {
					kend = pp0[j];
					while (k < kend) {
						++work[pi0[k]];
						if (pi0[k] != j)
						++work[     j];
						++k;
					}
				}
			}
			nnz1 = 0;
			for (i = 0; i < m; ++i)
				work[i] = (work[i]) ? nnz1++ : -1;
		}

		SEXP
		i1 = PROTECT(allocVector(INTSXP, nnz1)),
		x1 = PROTECT(allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(res, Matrix_iSym, i1);
		SET_SLOT(res, Matrix_xSym, x1);
		int *pi1 = INTEGER(i1);
		if (narm_)
			for (i = 0; i < nnz1; ++i)
				pi1[i] = n;

#define SUM_LOOP(_CTYPE0_, _PTR0_, _CTYPE1_, _PTR1_, _MASK_, \
		         _ZERO_, _ONE_, _NA_, _ISNA_, \
		         _MAP_, _CAST_, _INCREMENT_, _SCALE2_) \
		do { \
			_MASK_(_CTYPE0_ *px0 = _PTR0_(x0)); \
			       _CTYPE1_ *px1 = _PTR1_(x1) ; \
			       _CTYPE1_  tmp = (di != 'N') ? _ONE_ : _ZERO_; \
			for (i = 0; i < nnz1; ++i) \
				px1[i] = tmp; \
			if (class[1] != 's') { \
				for (k = 0; k < nnz0; ++k) { \
					if (_ISNA_(px0[k])) { \
						if (!narm) \
							px1[_MAP_(pi0[k])] = _NA_; \
						else if (narm_) \
							--pi1[_MAP_(pi0[k])]; \
					} else { \
						tmp = _CAST_(px0[k]); \
						_INCREMENT_(px1[_MAP_(pi0[k])], tmp); \
					} \
				} \
			} else { \
				int off; \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						off = pi0[k] != j; \
						if (_ISNA_(px0[k])) { \
							if (!narm) { \
								px1[_MAP_(pi0[k])] = _NA_; \
								if (off) \
								px1[_MAP_(     j)] = _NA_; \
							} else if (narm_) { \
								--pi1[_MAP_(pi0[k])]; \
								if (off) \
								--pi1[_MAP_(     j)]; \
							} \
						} else { \
							tmp = _CAST_(px0[k]); \
							_INCREMENT_(px1[_MAP_(pi0[k])], tmp); \
							if (off) \
							_INCREMENT_(px1[_MAP_(     j)], tmp); \
						} \
						++k; \
					} \
				} \
			} \
			if (mean) { \
				if (narm_) \
					for (i = 0; i < nnz1; ++i) \
						_SCALE2_(px1[i], pi1[i]); \
				else \
					for (i = 0; i < nnz1; ++i) \
						_SCALE2_(px1[i], n); \
			} \
		} while (0)

		SUM_CASES(MAP);
		for (i = 0; i < m; ++i)
			if (work[i] >= 0)
				*(pi1++) = i + 1;
		Matrix_Free(work, m);
		UNPROTECT(2); /* x1, i1 */

	} else {

		SEXP x1 = res;
		int *pi1 = NULL;
		if (narm_) {
			Matrix_Calloc(pi1, m, int);
			for (i = 0; i < m; ++i)
				pi1[i] = n;
		}
		SUM_CASES(NOMAP);
		if (narm_)
			Matrix_Free(pi1, m);

	}

#undef SUM_LOOP

	UNPROTECT(2); /* i0, p0 */
	return;
}

static
void Tsparse_colsum(SEXP obj, const char *class,
                    int m, int n, char di, int narm, int mean,
                    SEXP res, SEXP iSym, SEXP jSym)
{
	int narm_ = narm && mean && class[0] != 'n';
	if (narm_)
		obj = Tsparse_aggregate(obj);
	PROTECT(obj);

	SEXP i0 = PROTECT(GET_SLOT(obj, iSym)),
		j0 = PROTECT(GET_SLOT(obj, jSym));
	int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), j, nnz1 = n;
	R_xlen_t k, nnz0 = XLENGTH(i0);

	if (isS4(res)) {

		int *work;
		Matrix_Calloc(work, n, int);
		if (di != 'N')
			for (j = 0; j < n; ++j)
				work[j] = j;
		else {
			if (class[1] != 's') {
				for (k = 0; k < nnz0; ++k)
					++work[pj0[k]];
			} else {
				for (k = 0; k < nnz0; ++k) {
					++work[pj0[k]];
					if (pi0[k] != pj0[k])
					++work[pi0[k]];
				}
			}
			nnz1 = 0;
			for (j = 0; j < n; ++j)
				work[j] = (work[j]) ? nnz1++ : -1;
		}

		SEXP
		j1 = PROTECT(allocVector(INTSXP, nnz1)),
		x1 = PROTECT(allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(res, Matrix_iSym, j1);
		SET_SLOT(res, Matrix_xSym, x1);
		int *pj1 = INTEGER(j1);
		if (narm_)
			for (j = 0; j < nnz1; ++j)
				pj1[j] = m;

#define SUM_LOOP(_CTYPE0_, _PTR0_, _CTYPE1_, _PTR1_, _MASK_, \
		         _ZERO_, _ONE_, _NA_, _ISNA_, \
		         _MAP_, _CAST_, _INCREMENT_, _SCALE2_) \
		do { \
			_MASK_(_CTYPE0_ *px0 = _PTR0_(x0)); \
			       _CTYPE1_ *px1 = _PTR1_(x1) ; \
			       _CTYPE1_  tmp = (di != 'N') ? _ONE_ : _ZERO_; \
			for (j = 0; j < nnz1; ++j) \
				px1[j] = tmp; \
			if (class[1] != 's') { \
				for (k = 0; k < nnz0; ++k) { \
					if (_ISNA_(px0[k])) { \
						if (!narm) \
							px1[_MAP_(pj0[k])] = _NA_; \
						else if (narm_) \
							--pj1[_MAP_(pj0[k])]; \
					} else { \
						tmp = _CAST_(px0[k]); \
						_INCREMENT_(px1[_MAP_(pj0[k])], tmp); \
					} \
				} \
			} else { \
				int off; \
				for (k = 0; k < nnz0; ++k) { \
					off = pi0[k] != pj0[k]; \
					if (_ISNA_(px0[k])) { \
						if (!narm) { \
							px1[_MAP_(pj0[k])] = _NA_; \
							if (off) \
							px1[_MAP_(pi0[k])] = _NA_; \
						} else if (narm_) { \
							--pj1[_MAP_(pj0[k])]; \
							if (off) \
							--pj1[_MAP_(pi0[k])]; \
						} \
					} else { \
						tmp = _CAST_(px0[k]); \
						_INCREMENT_(px1[_MAP_(pj0[k])], tmp); \
						if (off) \
						_INCREMENT_(px1[_MAP_(pi0[k])], tmp); \
					} \
				} \
			} \
			if (mean) { \
				if (narm_) \
					for (j = 0; j < nnz1; ++j) \
						_SCALE2_(px1[j], pj1[j]); \
				else \
					for (j = 0; j < nnz1; ++j) \
						_SCALE2_(px1[j], m); \
			} \
		} while (0)

		SUM_CASES(MAP);
		for (j = 0; j < n; ++j)
			if (work[j] >= 0)
				*(pj1++) = j + 1;
		Matrix_Free(work, n);
		UNPROTECT(2); /* x1, j1 */

	} else {

		SEXP x1 = res;
		int *pj1 = NULL;
		if (narm_)
			Matrix_Calloc(pj1, n, int);
		SUM_CASES(NOMAP);
		if (narm_)
			Matrix_Free(pj1, n);

	}

#undef SUM_LOOP

	UNPROTECT(3); /* j0, i0, obj */
	return;
}

SEXP sparse_marginsum(SEXP obj, const char *class, int margin,
                      int narm, int mean, int sparse)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1],
		r = (margin == 0) ? m : n;

	SEXP res;
	SEXPTYPE type = SUM_TYPEOF(class[0]);
	if (sparse) {
		char cl[] = ".sparseVector";
		cl[0] = (type == CPLXSXP) ? 'z' : ((type == REALSXP) ? 'd' : 'i');
		PROTECT(res = newObject(cl));

		SEXP length = PROTECT(ScalarInteger(r));
		SET_SLOT(res, Matrix_lengthSym, length);
		UNPROTECT(1); /* length */
	} else {
		PROTECT(res = allocVector(type, r));

		SEXP dimnames = (class[1] != 's')
			? GET_SLOT(obj, Matrix_DimNamesSym)
			: get_symmetrized_DimNames(obj, -1),
			marnames = VECTOR_ELT(dimnames, margin);
		if (marnames != R_NilValue) {
			PROTECT(marnames);
			setAttrib(res, R_NamesSymbol, marnames);
			UNPROTECT(1); /* marnames */
		}
	}

	char di = 'N';
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(obj, Matrix_diagSym);
		di = *CHAR(STRING_ELT(diag, 0));
	}

	if (margin == 0) {
		if (class[2] == 'C') {
			Csparse_rowsum(obj, class, m, n, di, narm, mean, res,
			               Matrix_iSym);
		} else if (class[2] == 'R') {
			if (class[1] != 's')
			Csparse_colsum(obj, class, n, m, di, narm, mean, res);
			else
			Csparse_rowsum(obj, class, n, m, di, narm, mean, res,
			               Matrix_jSym);
		} else {
			Tsparse_colsum(obj, class, n, m, di, narm, mean, res,
			               Matrix_jSym, Matrix_iSym);
		}
	} else {
		if (class[2] == 'C') {
			if (class[1] != 's')
			Csparse_colsum(obj, class, m, n, di, narm, mean, res);
			else
			Csparse_rowsum(obj, class, m, n, di, narm, mean, res,
			               Matrix_iSym);
		} else if (class[2] == 'R') {
			Csparse_rowsum(obj, class, n, m, di, narm, mean, res,
			               Matrix_jSym);
		} else {
			Tsparse_colsum(obj, class, m, n, di, narm, mean, res,
			               Matrix_iSym, Matrix_jSym);
		}
	}

	UNPROTECT(1); /* res */
	return res;
}

/* (row|col)(Sums|Means)(<[CRT]sparseMatrix>) */
SEXP R_sparse_marginsum(SEXP obj, SEXP margin,
                        SEXP narm, SEXP mean, SEXP sparse)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int margin_;
	if (TYPEOF(margin) != INTSXP || LENGTH(margin) < 1 ||
	    ((margin_ = INTEGER(margin)[0]) != 0 && margin_ != 1))
		error(_("'%s' must be %d or %d"), "margin", 0, 1);

	int narm_;
	if (TYPEOF(narm) != LGLSXP || LENGTH(narm) < 1 ||
	    (narm_ = LOGICAL(narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	int mean_;
	if (TYPEOF(mean) != LGLSXP || LENGTH(mean) < 1 ||
	    (mean_ = LOGICAL(mean)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "mean", "TRUE", "FALSE");

	int sparse_;
	if (TYPEOF(sparse) != LGLSXP || LENGTH(sparse) < 1 ||
	    (sparse_ = LOGICAL(sparse)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "sparse", "TRUE", "FALSE");

	return sparse_marginsum(obj, valid[ivalid],
	                        margin_, narm_, mean_, sparse_);
}

#undef SUM_CASES
#undef SUM_TYPEOF

#define TRY_INCREMENT(_LABEL_) \
	do { \
		if ((s >= 0) \
			? ( t <= MATRIX_INT_FAST64_MAX - s) \
			: (-t <= s - MATRIX_INT_FAST64_MIN)) { \
			s += t; \
			t = 0; \
			count = 0; \
		} else { \
			over = 1; \
			goto _LABEL_; \
		} \
	} while (0)

#define LONGDOUBLE_AS_DOUBLE(v) \
	(v > DBL_MAX) ? R_PosInf : ((v < -DBL_MAX) ? R_NegInf : (double) v);

SEXP sparse_sum(SEXP obj, const char *class, int narm)
{
	if (class[2] == 'T')
		obj = Tsparse_aggregate(obj);
	PROTECT(obj);

	SEXP res;

	if (!narm && (class[0] == 'l' || class[0] == 'i')) {
		SEXP x = GET_SLOT(obj, Matrix_xSym);
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
		R_xlen_t nx = XLENGTH(x);
		while (nx--) {
			if (*px == NA_INTEGER) {
				res = allocVector(INTSXP, 1);
				INTEGER(res)[0] = NA_INTEGER;
				UNPROTECT(1); /* obj */
				return res;
			}
			++px;
		}
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char di = 'N';
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(obj, Matrix_diagSym);
		di = *CHAR(STRING_ELT(diag, 0));
	}

	int symmetric = class[1] == 's';

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j_, k = 0, kend,
			n_ = (class[2] == 'C') ? n : m;

		if (class[0] == 'n') {
			Matrix_int_fast64_t nnz = pp[n_ - 1];
			if (di != 'N')
				nnz += n;
			if (symmetric) {
				SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
				char ul = *CHAR(STRING_ELT(uplo, 0));

				nnz *= 2;
				for (j_ = 0; j_ < n_; ++j_) {
					kend = pp[j_];
					if (k < kend && pi[(ul == 'U') ? kend - 1 : k] == j_)
						--nnz;
					k = kend;
				}
			}
			if (nnz <= INT_MAX) {
				res = allocVector(INTSXP, 1);
				INTEGER(res)[0] = (int) nnz;
			} else {
				res = allocVector(REALSXP, 1);
				REAL(res)[0] = (double) nnz;
			}
			UNPROTECT(3); /* i, p, obj */
			return res;
		}

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		UNPROTECT(2); /* i, p */

		if (class[0] == 'z') {
			Rcomplex *px = COMPLEX(x);
			long double zr = (di == 'N') ? 0.0L : n, zi = 0.0L;
			for (j_ = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				while (k < kend) {
					if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
						zr += (symmetric && pi[k] != j_)
							? 2.0L * px[k].r : px[k].r;
						zi += (symmetric && pi[k] != j_)
							? 2.0L * px[k].i : px[k].i;
					}
					++k;
				}
			}
			res = allocVector(CPLXSXP, 1);
			COMPLEX(res)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
			COMPLEX(res)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
		} else if (class[0] == 'd') {
			double *px = REAL(x);
			long double zr = (di == 'N') ? 0.0L : n;
			for (j_ = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				while (k < kend) {
					if (!(narm && ISNAN(px[k])))
						zr += (symmetric && pi[k] != j_)
							? 2.0L * px[k] : px[k];
					++k;
				}
			}
			res = allocVector(REALSXP, 1);
			REAL(res)[0] = LONGDOUBLE_AS_DOUBLE(zr);
		} else {
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			Matrix_int_fast64_t s = (di == 'N') ? 0LL : n, t = 0LL;
			unsigned int count = 0;
			int over = 0;
			for (j_ = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				while (k < kend) {
					if (!narm || px[k] != NA_INTEGER) {
						int d = (symmetric && pi[k] != j_) ? 2 : 1;
						if (count > UINT_MAX - d)
							TRY_INCREMENT(ifoverC);
						t += (d == 2) ? 2LL * px[k] : px[k];
						count += d;
					}
					++k;
				}
			}
			TRY_INCREMENT(ifoverC);
		ifoverC:
			if (over) {
				long double zr = (long double) s + (long double) t;
				for (; j_ < n_; ++j_) {
					kend = pp[j_];
					while (k < kend) {
						if (!narm || px[k] != NA_INTEGER)
							zr += (symmetric && pi[k] != j_)
								? 2.0L * px[k] : px[k];
						++k;
					}
				}
				res = allocVector(REALSXP, 1);
				REAL(res)[0] = LONGDOUBLE_AS_DOUBLE(zr);
			} else if (s > INT_MIN && s <= INT_MAX) {
				res = allocVector(INTSXP, 1);
				INTEGER(res)[0] = (int) s;
			} else {
				res = allocVector(REALSXP, 1);
				REAL(res)[0] = (double) s;
			}
		}

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);

		if (class[0] == 'n') {
			Matrix_int_fast64_t nnz = (Matrix_int_fast64_t) kend;
			if (di != 'N')
				nnz += n;
			if (symmetric) {
				nnz *= 2;
				for (k = 0; k < kend; ++k)
					if (pi[k] == pj[k])
						--nnz;
			}
			if (nnz <= INT_MAX) {
				res = allocVector(INTSXP, 1);
				INTEGER(res)[0] = (int) nnz;
			} else {
				res = allocVector(REALSXP, 1);
				REAL(res)[0] = (double) nnz;
			}
			UNPROTECT(3); /* j, i, obj */
			return res;
		}

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		UNPROTECT(2); /* j, i */

		if (class[0] == 'z') {
			Rcomplex *px = COMPLEX(x);
			long double zr = (di == 'N') ? 0.0L : n, zi = 0.0L;
			for (k = 0; k < kend; ++k)
				if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
					zr += (symmetric && pi[k] != pj[k])
						? 2.0L * px[k].r : px[k].r;
					zi += (symmetric && pi[k] != pj[k])
						? 2.0L * px[k].i : px[k].i;
				}
			res = allocVector(CPLXSXP, 1);
			COMPLEX(res)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
			COMPLEX(res)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
		} else if (class[0] == 'd') {
			double *px = REAL(x);
			long double zr = (di == 'N') ? 0.0L : n;
			for (k = 0; k < kend; ++k)
				if (!(narm && ISNAN(px[k])))
					zr += (symmetric && pi[k] != pj[k])
						? 2.0L * px[k] : px[k];
			res = allocVector(REALSXP, 1);
			REAL(res)[0] = LONGDOUBLE_AS_DOUBLE(zr);
		} else {
			int *px = (class[0] == 'i') ? INTEGER(x) : LOGICAL(x);
			Matrix_int_fast64_t s = (di == 'N') ? 0LL : n, t = 0LL;
			unsigned int count = 0;
			int over = 0;
			for (k = 0; k < kend; ++k) {
				if (!narm || px[k] != NA_INTEGER) {
					int d = (symmetric && pi[k] != pj[k]) ? 2 : 1;
					if (count > UINT_MAX - d)
						TRY_INCREMENT(ifoverT);
					t += (d == 2) ? 2LL * px[k] : px[k];
					count += d;
				}
			}
			TRY_INCREMENT(ifoverT);
		ifoverT:
			if (over) {
				long double zr = (long double) s + (long double) t;
				for (; k < kend; ++k)
					if (!(narm && px[k] == NA_INTEGER))
						zr += (symmetric && pi[k] != pj[k])
							? 2.0L * px[k] : px[k];
				res = allocVector(REALSXP, 1);
				REAL(res)[0] = LONGDOUBLE_AS_DOUBLE(zr);
			} else if (s > INT_MIN && s <= INT_MAX) {
				res = allocVector(INTSXP, 1);
				INTEGER(res)[0] = (int) s;
			} else {
				res = allocVector(REALSXP, 1);
				REAL(res)[0] = (double) s;
			}
		}

	}

	UNPROTECT(1); /* obj */
	return res;
}

/* sum(<[CRT]sparseMatrix>) */
SEXP R_sparse_sum(SEXP obj, SEXP narm)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int narm_;
	if (TYPEOF(narm) != LGLSXP || LENGTH(narm) < 1 ||
	    (narm_ = LOGICAL(narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	return sparse_sum(obj, valid[ivalid], narm_);
}

SEXP sparse_prod(SEXP obj, const char *class, int narm)
{
	if (class[2] == 'T')
		obj = Tsparse_aggregate(obj);
	PROTECT(obj);

	SEXP res = PROTECT(allocVector((class[0] == 'z') ? CPLXSXP : REALSXP, 1));

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (class[1] == 't') {
			SEXP diag = GET_SLOT(obj, Matrix_diagSym);
			di = *CHAR(STRING_ELT(diag, 0));
		}
	}

	int symmetric = (class[1] != 's')
		? 0 : (((class[2] == 'C') == (ul == 'U')) ? 1 : -1);
	long double zr = 1.0L, zi = 0.0L;

	Matrix_int_fast64_t mn = (Matrix_int_fast64_t) m * n,
		nnz, nnzmax = (symmetric) ? (mn + n) / 2 : mn;

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p) + 1, *pi = INTEGER(i), i_, j_, k = 0, kend,
			m_ = (class[2] == 'C') ? m : n, n_ = (class[2] == 'C') ? n : m,
			seen0 = 0;

		nnz = pp[n_ - 1];
		if (di != 'N')
			nnz += n;
		if (class[0] == 'n') {
			REAL(res)[0] = (nnz == nnzmax) ? 1.0 : 0.0;
			UNPROTECT(4); /* i, p, res, obj */
			return res;
		}

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		UNPROTECT(2); /* i, p */

		if (class[0] == 'z') {
			Rcomplex *px = COMPLEX(x);
			long double zr0, zi0;
			for (j_ = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				if (seen0 || kend - k == m_) {
					while (k < kend) {
						if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
							zr0 = zr; zi0 = zi;
							zr = zr0 * px[k].r - zi0 * px[k].i;
							zi = zr0 * px[k].i + zi0 * px[k].r;
							if (symmetric && pi[k] != j_) {
							zr0 = zr; zi0 = zi;
							zr = zr0 * px[k].r - zi0 * px[k].i;
							zi = zr0 * px[k].i + zi0 * px[k].r;
							}
						}
						++k;
					}
				} else {
					int i0 = (symmetric >= 0) ? 0 : j_,
						i1 = (symmetric <= 0) ? m_ : j_ + 1;
					for (i_ = i0; i_ < i1; ++i_) {
						if (seen0 || (k < kend && pi[k] == i_)) {
						if (k >= kend)
							break;
						if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
							zr0 = zr; zi0 = zi;
							zr = zr0 * px[k].r - zi0 * px[k].i;
							zi = zr0 * px[k].i + zi0 * px[k].r;
							if (symmetric && pi[k] != j_) {
							zr0 = zr; zi0 = zi;
							zr = zr0 * px[k].r - zi0 * px[k].i;
							zi = zr0 * px[k].i + zi0 * px[k].r;
							}
						}
						++k;
						} else if (di == 'N' || i_ != j_) {
						zr *= 0.0L;
						zi *= 0.0L;
						seen0 = 1;
						}
					}
				}
			}
		} else if (class[0] == 'd') {
			double *px = REAL(x);
			for (j_ = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				if (seen0 || kend - k == m_) {
					while (k < kend) {
						if (!(narm && ISNAN(px[k])))
							zr *= (symmetric && pi[k] != j_)
								? (long double) px[k] * px[k] : px[k];
						++k;
					}
				} else {
					int i0 = (symmetric >= 0) ? 0 : j_,
						i1 = (symmetric <= 0) ? m_ : j_ + 1;
					for (i_ = i0; i_ < i1; ++i_) {
						if (seen0 || (k < kend && pi[k] == i_)) {
						if (k >= kend)
							break;
						if (!(narm && ISNAN(px[k])))
							zr *= (symmetric && pi[k] != j_)
								? (long double) px[k] * px[k] : px[k];
						++k;
						} else if (di == 'N' || i_ != j_) {
						zr *= 0.0L;
						seen0 = 1;
						}
					}
				}
			}
		} else {
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			for (j_ = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				if (seen0 || kend - k == m_) {
					while (k < kend) {
						if (px[k] != NA_INTEGER)
							zr *= (symmetric && pi[k] != j_)
								? (long double) px[k] * px[k] : px[k];
						else if (!narm)
							zr *= NA_REAL;
						++k;
					}
				} else {
					int i0 = (symmetric >= 0) ? 0 : j_,
						i1 = (symmetric <= 0) ? m_ : j_ + 1;
					for (i_ = i0; i_ < i1; ++i_) {
						if (seen0 || (k < kend && pi[k] == i_)) {
						if (k >= kend)
							break;
						if (px[k] != NA_INTEGER)
							zr *= (symmetric && pi[k] != j_)
								? (long double) px[k] * px[k] : px[k];
						else if (!narm)
							zr *= NA_REAL;
						++k;
						} else if (di == 'N' || i_ != j_) {
						zr *= 0.0L;
						seen0 = 1;
						}
					}
				}
			}
		}

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);

		nnz = (Matrix_int_fast64_t) kend;
		if (di != 'N')
			nnz += n;
		if (class[0] == 'n') {
			REAL(res)[0] = (nnz == nnzmax) ? 1.0 : 0.0;
			UNPROTECT(4); /* j, i, res, obj */
			return res;
		}
		if (nnz < nnzmax)
			zr = 0.0;

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		UNPROTECT(2); /* j, i */

		if (class[0] == 'z') {
			Rcomplex *px = COMPLEX(x);
			long double zr0, zi0;
			for (k = 0; k < kend; ++k)
				if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
					zr0 = zr; zi0 = zi;
					zr = zr0 * px[k].r - zi0 * px[k].i;
					zi = zr0 * px[k].i + zi0 * px[k].r;
					if (symmetric && pi[k] != pj[k]) {
					zr0 = zr; zi0 = zi;
					zr = zr0 * px[k].r - zi0 * px[k].i;
					zi = zr0 * px[k].i + zi0 * px[k].r;
					}
				}
		} else if (class[0] == 'd') {
			double *px = REAL(x);
			for (k = 0; k < kend; ++k)
				if (!(narm && ISNAN(px[k])))
					zr *= (symmetric && pi[k] != pj[k])
						? (long double) px[k] * px[k] : px[k];
		} else {
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			for (k = 0; k < kend; ++k)
				if (px[k] != NA_INTEGER)
					zr *= (symmetric && pi[k] != pj[k])
						? (long double) px[k] * px[k] : px[k];
				else if (!narm)
					zr *= NA_REAL;
		}

	}

	if (class[0] == 'z') {
		COMPLEX(res)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
		COMPLEX(res)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
	} else
		   REAL(res)[0]   = LONGDOUBLE_AS_DOUBLE(zr);
	UNPROTECT(2); /* res, obj */
	return res;
}

/* prod(<[CRT]sparseMatrix>) */
SEXP R_sparse_prod(SEXP obj, SEXP narm)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int narm_;
	if (TYPEOF(narm) != LGLSXP || LENGTH(narm) < 1 ||
	    (narm_ = LOGICAL(narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	return sparse_prod(obj, valid[ivalid], narm_);
}

#undef TRY_INCREMENT
#undef LONGDOUBLE_AS_DOUBLE

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
		PROTECT(to = newObject(cl));
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
		PROTECT(to = newObject(cl));
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
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorsSym, factors);
		UNPROTECT(1); /* factors */
	}

	UNPROTECT(1); /* to */
	return to;
}
