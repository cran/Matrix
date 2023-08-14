#include "dense.h"

SEXP dense_band(SEXP from, const char *class, int a, int b, int new)
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

#define BAND_CASES(_DO_) \
	do { \
		switch (class[0]) { \
		case 'n': \
		case 'l': \
			_DO_(int, LOGICAL, i); \
			break; \
		case 'i': \
			_DO_(int, INTEGER, i); \
			break; \
		case 'd': \
			_DO_(double, REAL, d); \
			break; \
		case 'z': \
			_DO_(Rcomplex, COMPLEX, z); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define IF_UNPACKED_ELSE(_IF_, _ELSE_) \
	do { \
		if (class[2] != 'p') \
			BAND_CASES(_IF_); \
		else \
			BAND_CASES(_ELSE_); \
	} while (0)

#define UNPACKED_MAKE_BANDED(_CTYPE_, _PTR_, _PREFIX_) \
	_PREFIX_ ## dense_unpacked_make_banded(_PTR_(x1), m, n, a, b, di)

#define PACKED_MAKE_BANDED(_CTYPE_, _PTR_, _PREFIX_) \
	_PREFIX_ ## dense_packed_make_banded(_PTR_(x1), n, a, b, ult, di)

#define UNPACKED_COPY_DIAGONAL(_CTYPE_, _PTR_, _PREFIX_) \
	do { \
		Matrix_memset(_PTR_(x1), 0, len, sizeof(_CTYPE_)); \
		if (a <= 0 && b >= 0) \
			_PREFIX_ ## dense_unpacked_copy_diagonal( \
				_PTR_(x1), _PTR_(x0), n, len, 'U', di); \
	} while (0)

#define PACKED_COPY_DIAGONAL(_CTYPE_, _PTR_, _PREFIX_)\
	do { \
		Matrix_memset(_PTR_(x1), 0, len, sizeof(_CTYPE_)); \
		if (a <= 0 && b >= 0) \
			_PREFIX_ ## dense_packed_copy_diagonal( \
				_PTR_(x1), _PTR_(x0), n, len, ult, ulf, di); \
	} while (0)

	char ulf = 'U', ult = 'U', di = 'N';
	if (class[1] != 'g') {
		if (ge) {
			/* defined in ./coerce.c : */
			SEXP dense_as_general(SEXP, const char *, int);
			SEXP to = PROTECT(dense_as_general(from, class, 1)),
				x1 = PROTECT(GET_SLOT(to, Matrix_xSym));
			BAND_CASES(UNPACKED_MAKE_BANDED);
			UNPROTECT(2);
			return to;
		}

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

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (ge) ? 'g' :                      ((sy) ? 's' : 't')       ;
	cl[2] = (ge) ? 'e' : ((class[2] != 'p') ? ((sy) ? 'y' : 'r') : 'p');
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

	if (ge) {
		SEXP x1 = PROTECT(GET_SLOT(from, Matrix_xSym));
		if (new) {
			x1 = duplicate(x1);
			UNPROTECT(1);
			PROTECT(x1);
		}
		SET_SLOT(to, Matrix_xSym, x1);
		BAND_CASES(UNPACKED_MAKE_BANDED);
		UNPROTECT(2); /* x, to */
		return to;
	}

	/* Returning .(sy|sp|tr|tp)Matrix ... */

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

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1;

	if (tr && class[1] == 't') {
		/* Result is either a diagonal matrix or a zero matrix */
		R_xlen_t len = XLENGTH(x0);
		PROTECT(x1 = allocVector(TYPEOF(x0), len));
		IF_UNPACKED_ELSE(UNPACKED_COPY_DIAGONAL, PACKED_COPY_DIAGONAL);
	} else {
		if (!new)
			PROTECT(x1 = x0);
		else if (sy || (tr && (class[1] == 'g' || ulf == ult || n <= 1)))
			PROTECT(x1 = duplicate(x0));
		else if (class[2] != 'p')
			/* band is "opposite" the stored triangle: */
			PROTECT(x1 = unpacked_force(x0, n, ulf, '\0'));
		else
			/* band is "opposite" the stored triangle: */
			PROTECT(x1 = packed_transpose(x0, n, ulf));
		IF_UNPACKED_ELSE(UNPACKED_MAKE_BANDED, PACKED_MAKE_BANDED);
	}
	SET_SLOT(to, Matrix_xSym, x1);

#undef BAND_CASES
#undef IF_UNPACKED_ELSE
#undef UNPACKED_MAKE_BANDED
#undef UNPACKED_COPY_DIAGONAL
#undef PACKED_MAKE_BANDED
#undef PACKED_COPY_DIAGONAL

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* band(<denseMatrix>, k1, k2), tri[ul](<denseMatrix>, k) */
/* band(     <matrix>, k1, k2), tri[ul](     <matrix>, k) */
/* NB: argument validation more or less copied by R_sparse_band() */
SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2)
{
	static const char *valid[] = {
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
	int ivalid = R_check_class_etc(from, valid), isS4 = ivalid >= 0;
	if (!isS4) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, int, int);
		from = matrix_as_dense(from, ".ge", '\0', '\0', 0, 1);
	}
	PROTECT(from);
	if (!isS4)
		ivalid = R_check_class_etc(from, valid);

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
		error(_("'%s' must be less than or equal to '%s'"),
		      "k1", "k2");

	from = dense_band(from, valid[ivalid], a, b, isS4);
	UNPROTECT(1);
	return from;
}

/* colSums(<denseMatrix>) */
SEXP R_dense_colSums(SEXP obj, SEXP narm, SEXP mean)
{
	static const char *valid[] = {
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);
	const char *cl = valid[ivalid];
	if (cl[1] == 's')
		return R_dense_rowSums(obj, narm, mean);

	int doNaRm = asLogical(narm) != 0,
		doMean = asLogical(mean) != 0,
		doCount = doNaRm && doMean;

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */

	char ul = 'U', di = 'N';
	if (cl[1] == 't') {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */
	}

	SEXP res = PROTECT(allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, n)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int i, j, count = m;

#define DENSE_COLSUMS_LOOP \
	do { \
		if (cl[1] == 'g') { /* general */ \
			for (j = 0; j < n; ++j, ++pres) { \
				DO_INIT(ZERO); \
				for (i = 0; i < m; ++i, ++px) \
					DO_INCR; \
				DO_SCALE; \
			} \
		} else if (cl[2] != 'p') { \
			if (ul == 'U') { /* unpacked upper triangular */ \
				if (di == 'N') { \
					for (j = 0; j < n; ++j, ++pres) { \
						DO_INIT(ZERO); \
						for (i = 0; i <= j; ++i, ++px) \
							DO_INCR; \
						DO_SCALE; \
						px += n-j-1; \
					} \
				} else { \
					for (j = 0; j < n; ++j, ++pres) { \
						DO_INIT(ONE); \
						for (i = 0; i < j; ++i, ++px) \
							DO_INCR; \
						DO_SCALE; \
						px += n-j; \
					} \
				} \
			} else { /* unpacked lower triangular */ \
				if (di == 'N') { \
					for (j = 0; j < n; ++j, ++pres) { \
						px += j; \
						DO_INIT(ZERO); \
						for (i = j; i < n; ++i, ++px) \
							DO_INCR; \
						DO_SCALE; \
					} \
				} else { \
					for (i = j = 0; j < n; ++j, i = j, ++pres) { \
						px += j+1; \
						DO_INIT(ONE); \
						for (i = j+1; i < n; ++i, ++px) \
							DO_INCR; \
						DO_SCALE; \
					} \
				} \
			} \
		} else { \
			if (ul == 'U') { /* packed upper triangular */ \
				if (di == 'N') { \
					for (j = 0; j < n; ++j, ++pres) { \
						DO_INIT(ZERO); \
						for (i = 0; i <= j; ++i, ++px) \
							DO_INCR; \
						DO_SCALE; \
					} \
				} else { \
					for (j = 0; j < n; ++j, ++pres) { \
						DO_INIT(ONE); \
						for (i = 0; i < j; ++i, ++px) \
							DO_INCR; \
						DO_SCALE; \
						++px; \
					} \
				} \
			} else { /* packed lower triangular */ \
				if (di == 'N') { \
					for (j = 0; j < n; ++j, ++pres) { \
						DO_INIT(ZERO); \
						for (i = j; i < n; ++i, ++px) \
							DO_INCR; \
						DO_SCALE; \
					} \
				} else { \
					for (i = j = 0; j < n; ++j, i = j, ++pres) { \
						++px; \
						DO_INIT(ONE); \
						for (i = j+1; i < n; ++i, ++px) \
							DO_INCR; \
						DO_SCALE; \
					} \
				} \
			} \
		} \
	} while (0)

#define DENSE_COLSUMS(_CTYPE1_, _PTR1_, _CTYPE2_, _PTR2_) \
	do { \
		_CTYPE1_ *pres = _PTR1_(res); \
		_CTYPE2_ *px = _PTR2_(x); \
		DENSE_COLSUMS_LOOP; \
	} while (0)

	switch (cl[0]) {
	case 'n':

#define ZERO         0.0
#define ONE          1.0
#define DO_INIT(_U_) *pres = _U_
#define DO_INCR      if (*px) *pres += 1.0
#define DO_SCALE     if (doMean) *pres /= count

		DENSE_COLSUMS(double, REAL, int, LOGICAL);
		break;

#undef DO_INIT
#undef DO_INCR

	case 'l':

#define DO_INIT(_U_) \
		do { \
			*pres = _U_; \
			if (doCount) \
				count = m; \
		} while (0)
#define DO_INCR \
		do { \
			if (*px != NA_LOGICAL) { \
				if (*px) *pres += 1.0; \
			} else if (!doNaRm) \
				*pres = NA_REAL; \
			else if (doMean) \
				--count; \
		} while (0)

	DENSE_COLSUMS(double, REAL, int, LOGICAL);
	break;

#undef DO_INCR

	case 'i':

#define DO_INCR \
		do { \
			if (*px != NA_INTEGER) \
				*pres += *px; \
			else if (!doNaRm) \
				*pres = NA_REAL; \
			else if (doMean) \
				--count; \
		} while (0)

		DENSE_COLSUMS(double, REAL, int, INTEGER);
		break;

#undef DO_INCR

	case 'd':

#define DO_INCR \
		do { \
			if (!(doNaRm && ISNAN(*px))) \
				*pres += *px; \
			else if (doMean) \
				--count; \
		} while (0)

		DENSE_COLSUMS(double, REAL, double, REAL);
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
			if (!(doNaRm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
				(*pres).r += (*px).r; \
				(*pres).i += (*px).i; \
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

		DENSE_COLSUMS(Rcomplex, COMPLEX, Rcomplex, COMPLEX);
		break;

#undef ZERO
#undef ONE
#undef DO_INIT
#undef DO_INCR
#undef DO_SCALE

	default:
		break;
	}

#undef DENSE_COLSUMS
#undef DENSE_COLSUMS_LOOP

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		nms = VECTOR_ELT(dimnames, 1);
	if (!isNull(nms))
		setAttrib(res, R_NamesSymbol, nms);

	UNPROTECT(3); /* dimnames, x, res */
	return res;
}

/* rowSums(<denseMatrix>) */
SEXP R_dense_rowSums(SEXP obj, SEXP narm, SEXP mean)
{
	static const char *valid[] = {
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);
	const char *cl = valid[ivalid];

	int doNaRm = asLogical(narm) != 0,
		doMean = asLogical(mean) != 0;

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
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
		}
	}

	SEXP res = PROTECT(allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, m)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int i, j, *pcount = NULL;

#define DENSE_ROWSUMS_LOOP \
	do { \
		if (cl[1] == 'g') { /* general */ \
			for (j = 0; j < n; ++j) \
				for (i = 0; i < m; ++i, ++px) \
					DO_INCR; \
		} else if (cl[1] == 't') { \
			if (cl[2] != 'p') { \
				if (ul == 'U') { /* unpacked upper triangular */ \
					if (di == 'N') { \
						for (j = 0; j < n; ++j) { \
							for (i = 0; i <= j; ++i, ++px) \
								DO_INCR; \
							px += n-j-1; \
						} \
					} else { \
						for (j = 0; j < n; ++j) { \
							for (i = 0; i < j; ++i, ++px) \
								DO_INCR; \
							px += n-j; \
						} \
					} \
				} else { /* unpacked lower triangular */ \
					if (di == 'N') { \
						for (j = 0; j < n; ++j) { \
							px += j; \
							for (i = j; i < n; ++i, ++px) \
								DO_INCR; \
						} \
					} else { \
						for (i = j = 0; j < n; ++j, i = j) { \
							px += j+1; \
							for (i = j+1; i < n; ++i, ++px) \
								DO_INCR; \
						} \
					} \
				} \
			} else { \
				if (ul == 'U') { /* packed upper triangular */ \
					if (di == 'N') { \
						for (j = 0; j < n; ++j) \
							for (i = 0; i <= j; ++i, ++px) \
								DO_INCR; \
					} else { \
						for (j = 0; j < n; ++j) { \
							for (i = 0; i < j; ++i, ++px) \
								DO_INCR; \
							++px; \
						} \
					} \
				} else { /* packed lower triangular */ \
					if (di == 'N') { \
						for (j = 0; j < n; ++j) \
							for (i = j; i < n; ++i, ++px) \
								DO_INCR; \
					} else { \
						for (i = j = 0; j < n; ++j, i = j) { \
							++px; \
							for (i = j+1; i < n; ++i, ++px) \
								DO_INCR; \
						} \
					} \
				} \
			} \
		} else { \
			if (cl[2] != 'p') { \
				if (ul == 'U') { /* unpacked upper symmetric */ \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i, ++px) \
							DO_INCR_SYMM; \
						DO_INCR; \
						px += n-j; \
					} \
				} else { /* unpacked lower symmetric */ \
					for (i = j = 0; j < n; ++j, i = j) { \
						px += j; \
						DO_INCR; \
						++px; \
						for (i = j+1; i < n; ++i, ++px) \
							DO_INCR_SYMM; \
					} \
				} \
			} else { \
				if (ul == 'U') { /* packed upper symmetric */ \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i, ++px) \
							DO_INCR_SYMM; \
						DO_INCR; \
						++px; \
					} \
				} else { /* packed lower symmetric */ \
					for (i = j = 0; j < n; ++j, i = j) { \
						DO_INCR; \
						++px; \
						for (i = j+1; i < n; ++i, ++px) \
							DO_INCR_SYMM; \
					} \
				} \
			} \
		} \
	} while (0)

#define DENSE_ROWSUMS(_CTYPE1_, _PTR1_, _CTYPE2_, _PTR2_) \
	do { \
		_CTYPE1_ *pres = _PTR1_(res), u = (di == 'N') ? ZERO : ONE; \
		_CTYPE2_ *px   = _PTR2_(x); \
		if (doNaRm && doMean && cl[0] != 'n') { \
			Matrix_Calloc(pcount, m, int); \
			for (i = 0; i < m; ++i) { \
				pres[i] = u; \
				pcount[i] = n; \
			} \
		} else { \
			for (i = 0; i < m; ++i) \
				pres[i] = u; \
		} \
		DENSE_ROWSUMS_LOOP; \
	} while (0)

	switch (cl[0]) {
	case 'n':

#define ZERO         0.0
#define ONE          1.0
#define DO_INCR      if (*px) pres[i] += 1.0
#define DO_INCR_SYMM \
		do { \
			if (*px) { \
				pres[i] += 1.0; \
				pres[j] += 1.0; \
			} \
		} while (0)

		DENSE_ROWSUMS(double, REAL, int, LOGICAL);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'l':

#define DO_INCR \
		do { \
			if (*px != NA_LOGICAL) { \
				if (*px) \
					pres[i] += 1.0; \
			} else if (!doNaRm) \
				pres[i] = NA_REAL; \
			else if (doMean) \
				--pcount[i]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (*px != NA_LOGICAL) { \
				if (*px) { \
					pres[i] += 1.0; \
					pres[j] += 1.0; \
				} \
			} else if (!doNaRm) { \
				pres[i] = NA_REAL; \
				pres[j] = NA_REAL; \
			} else if (doMean) { \
				--pcount[i]; \
				--pcount[j]; \
			} \
		} while (0)

		DENSE_ROWSUMS(double, REAL, int, LOGICAL);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'i':

#define DO_INCR \
		do { \
			if (*px != NA_INTEGER) \
				pres[i] += *px; \
			else if (!doNaRm) \
				pres[i] = NA_REAL; \
			else if (doMean) \
				--pcount[i]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (*px != NA_INTEGER) { \
				pres[i] += *px; \
				pres[j] += *px; \
			} else if (!doNaRm) { \
				pres[i] = NA_REAL; \
				pres[j] = NA_REAL; \
			} else if (doMean) { \
				--pcount[i]; \
				--pcount[j]; \
			} \
		} while (0)

		DENSE_ROWSUMS(double, REAL, int, INTEGER);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'd':

#define DO_INCR \
		do { \
			if (!(doNaRm && ISNAN(*px))) \
				pres[i] += *px; \
			else if (doMean) \
				--pcount[i]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (!(doNaRm && ISNAN(*px))) { \
				pres[i] += *px; \
				pres[j] += *px; \
			} else if (doMean) { \
				--pcount[i]; \
				--pcount[j]; \
			} \
		} while (0)

		DENSE_ROWSUMS(double, REAL, double, REAL);
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
			if (!(doNaRm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
				pres[i].r += (*px).r; \
				pres[i].i += (*px).i; \
			} else if (doMean) \
				--pcount[i]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (!(doNaRm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
				pres[i].r += (*px).r; \
				pres[i].i += (*px).i; \
				pres[j].r += (*px).r; \
				pres[j].i += (*px).i; \
			} else if (doMean) { \
				--pcount[i]; \
				--pcount[j]; \
			} \
		} while (0)

		DENSE_ROWSUMS(Rcomplex, COMPLEX, Rcomplex, COMPLEX);
		break;

#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_INCR_SYMM

	default:
		break;
	}

#undef DENSE_ROWSUMS
#undef DENSE_ROWSUMS_LOOP

	if (doMean) {
		if (cl[0] != 'z') {
			double *pres = REAL(res);
			if (doNaRm && cl[0] != 'n') {
				for (i = 0; i < m; ++i)
					pres[i] /= pcount[i];
				Matrix_Free(pcount, m);
			} else {
				for (i = 0; i < m; ++i)
					pres[i] /= n;
			}
		} else {
			Rcomplex *pres = COMPLEX(res);
			if (doNaRm) {
				for (i = 0; i < m; ++i) {
					pres[i].r /= pcount[i];
					pres[i].i /= pcount[i];
				}
				Matrix_Free(pcount, m);
			} else {
				for (i = 0; i < m; ++i) {
					pres[i].r /= n;
					pres[i].i /= n;
				}
			}
		}
	}

	SEXP dimnames;
	if (cl[1] != 's')
		PROTECT(dimnames = GET_SLOT(obj, Matrix_DimNamesSym));
	else
		PROTECT(dimnames = get_symmetrized_DimNames(obj, -1));
	SEXP nms = VECTOR_ELT(dimnames, 0);
	if (!isNull(nms))
		setAttrib(res, R_NamesSymbol, nms);

	UNPROTECT(3); /* dimnames, x, res */
	return res;
}

/* MJ: unused */
#if 0

/**
 * Perform a left cyclic shift of columns j to k in the upper triangular
 * matrix x, then restore it to upper triangular form with Givens rotations.
 * The algorithm is based on the Fortran routine DCHEX from Linpack.
 *
 * The lower triangle of x is not modified.
 *
 * @param x Matrix stored in column-major order
 * @param ldx leading dimension of x
 * @param j column number (0-based) that will be shifted to position k
 * @param k last column number (0-based) to be shifted
 * @param cosines cosines of the Givens rotations
 * @param sines sines of the Givens rotations
 *
 * @return 0 for success
 */
static
int left_cyclic(double *x, int ldx, int j, int k,
                double *cosines, double *sines)
{
	if (j < 0)
		error(_("incorrect left cyclic shift, j (%d) < 0"),
		      j);
	if (j >= k)
		error(_("incorrect left cyclic shift, j (%d) >= k (%d)"),
		      j, k);
	if (ldx < k)
		error(_("incorrect left cyclic shift, k (%d) > ldx (%d)"),
		      k, ldx);

	double *lastcol = (double *) R_alloc((size_t) k + 1, sizeof(double));
	int i;
	/* keep a copy of column j */
	for (i = 0; i <= j; i++)
		lastcol[i] = x[i + j*ldx];
	/* For safety, zero the rest */
	for (i = j+1; i <= k; i++)
		lastcol[i] = 0.;
	for (int jj = j+1, ind = 0; jj <= k; jj++, ind++) {
		/* columns to be shifted */
		int diagind = jj*(ldx+1); //  ind == (jj-j) - 1
		double tmp = x[diagind], cc, ss;
		/* Calculate the Givens rotation. */
		/* This modified the super-diagonal element */
		F77_CALL(drotg)(x+diagind-1, &tmp, cosines+ind, sines+ind);
		cc = cosines[ind];
		ss = sines[ind];
		/* Copy column jj+1 to column jj. */
		for (i = 0; i < jj; i++)
			x[i + (jj-1)*ldx] = x[i+jj*ldx];
		/* Apply rotation to columns up to k */
		for (i = jj; i < k; i++) {
			tmp = cc*x[(jj-1)+i*ldx] + ss*x[jj+i*ldx];
			x[jj+i*ldx] = cc*x[jj+i*ldx] - ss*x[(jj-1)+i*ldx];
			x[(jj-1)+i*ldx] = tmp;
		}
		/* Apply rotation to lastcol */
		lastcol[jj] = -ss*lastcol[jj-1]; lastcol[jj-1] *= cc;
	}
	/* Copy lastcol to column k */
	for(i = 0; i <= k; i++)
		x[i+k*ldx] = lastcol[i];
	return 0;
}

static
SEXP getGivens(double *x, int ldx, int jmin, int rank)
{
	int shiftlen = (rank - jmin) - 1;
	SEXP ans = PROTECT(allocVector(VECSXP, 4)), nms, cosines, sines;
	SET_VECTOR_ELT(ans, 0, ScalarInteger(jmin));
	SET_VECTOR_ELT(ans, 1, ScalarInteger(rank));
	SET_VECTOR_ELT(ans, 2, cosines = allocVector(REALSXP, shiftlen));
	SET_VECTOR_ELT(ans, 3,   sines = allocVector(REALSXP, shiftlen));
	setAttrib(ans, R_NamesSymbol, nms = allocVector(STRSXP, 4));
	SET_STRING_ELT(nms, 0, mkChar("jmin"));
	SET_STRING_ELT(nms, 1, mkChar("rank"));
	SET_STRING_ELT(nms, 2, mkChar("cosines"));
	SET_STRING_ELT(nms, 3, mkChar("sines"));
	if (left_cyclic(x, ldx, jmin, rank - 1, REAL(cosines), REAL(sines)))
		error(_("unknown error in getGivens"));
	UNPROTECT(1);
	return ans;
}

static
SEXP checkGivens(SEXP X, SEXP jmin, SEXP rank)
{
	if (!(isReal(X) && isMatrix(X)))
		error(_("X must be a numeric (double precision) matrix"));
	SEXP ans = PROTECT(allocVector(VECSXP, 2)),
		Xcp = PROTECT(duplicate(X));
	int Xdims = INTEGER(getAttrib(X, R_DimSymbol));
	SET_VECTOR_ELT(ans, 0, Xcp);
	SET_VECTOR_ELT(ans, 1, getGivens(REAL(Xcp), Xdims[0],
	                                 asInteger(jmin), asInteger(rank)));
	UNPROTECT(2);
	return ans;
}

SEXP lsq_dense_chol(SEXP X, SEXP y)
{
	if (!(isReal(X) && isMatrix(X)))
		error(_("X must be a numeric (double precision) matrix"));
	if (!(isReal(y) && isMatrix(y)))
		error(_("y must be a numeric (double precision) matrix"));
	int *Xdims = INTEGER(getAttrib(X, R_DimSymbol)),
		*ydims = INTEGER(getAttrib(y, R_DimSymbol));
	if (Xdims[0] != ydim[0])
		error(_("number of rows in y (%d) does not match "
		        "number of rows in X (%d)"),
		      ydims[0], Xdims[0]);
	int n = Xdims[0], p = Xdims[1], k = ydims[1];
	if (p < 1 || k < 1)
		return allocMatrix(REALSXP, p, k);
	SEXP ans = PROTECT(allocMatrix(REALSXP, p, k));
	double d_one = 1.0, d_zero = 0.0,
		*xpx = (double *) R_alloc((size_t) p * p, sizeof(double));
	int info;
	F77_CALL(dgemm)("T", "N", &p, &k, &n, &d_one, REAL(X), &n, REAL(y),
	                &n, &d_zero, REAL(ans), &p FCONE FCONE);
	F77_CALL(dsyrk)("U", "T", &p, &n, &d_one, REAL(X), &n, &d_zero,
	                xpx, &p FCONE FCONE);
	F77_CALL(dposv)("U", &p, &k, xpx, &p, REAL(ans), &p, &info FCONE);
	if (info)
		error(_("LAPACK dposv returned error code %d"), info);
	UNPROTECT(1);
	return ans;
}

SEXP lsq_dense_qr(SEXP X, SEXP y)
{
	if (!(isReal(X) && isMatrix(X)))
		error(_("X must be a numeric (double precision) matrix"));
	if (!(isReal(y) && isMatrix(y)))
		error(_("y must be a numeric (double precision) matrix"));
	int *Xdims = INTEGER(getAttrib(X, R_DimSymbol)),
		*ydims = INTEGER(getAttrib(y, R_DimSymbol));
	if (Xdims[0] != ydim[0])
		error(_("number of rows in y (%d) does not match "
		        "number of rows in X (%d)"),
		      ydims[0], Xdims[0]);
	int n = Xdims[0], p = Xdims[1], k = ydims[1];
	if (p < 1 || k < 1)
		return allocMatrix(REALSXP, p, k);
	SEXP ans = PROTECT(duplicate(y));
	double *xvals = (double *) R_alloc((size_t) n * p, sizeof(double)),
		*work, tmp;
	int lwork = -1, info;
	Memcpy(xvals, REAL(X), (size_t) n * p);
	F77_CALL(dgels)("N", &n, &p, &k, xvals, &n, REAL(ans), &n,
	                &tmp, &lwork, &info FCONE);
	if (info)
		error(_("LAPACK dgels returned error code %d"), info);
	lwork = (int) tmp;
	work = (double *) R_alloc((size_t) lwork, sizeof(double));
	F77_CALL(dgels)("N", &n, &p, &k, xvals, &n, REAL(ans), &n,
	                work, &lwork, &info FCONE);
	if (info)
		error(_("LAPACK dgels returned error code %d"), info);
	UNPROTECT(1);
	return ans;
}

/* Rank-Correcting/Adapting LAPACK  QR Decomposition
 * From Doug Bates' initial import; __unused__
 *
 * Provides a qr() with 'rcond' and rank reduction while(rcond < tol),
 * possibly via Givens rotations but WITHOUT PIVOTING
 *
 * .Call(Matrix:::lapack_qr, A, 1e-17)
 *  --> ~/R/MM/Pkg-ex/Matrix/qr-rank-deficient.R
 *
 * TODO: export as Matrix::qrNoPiv() or qr1()  or similar
 */
SEXP lapack_qr(SEXP Xin, SEXP tl)
{
	if (!(isReal(Xin) && isMatrix(Xin)))
		error(_("X must be a real (numeric) matrix"));
	double tol = asReal(tl);
	if (tol < 0.0)
		error(_("tol, given as %g, must be >= 0"), tol);
	if (tol > 1.0)
		error(_("tol, given as %g, must be <= 1"), tol);
	SEXP ans = PROTECT(allocVector(VECSXP, 5)), X, qraux, pivot;
	int *Xdims = INTEGER(getAttrib(Xin, R_DimSymbol)),
		n = Xdims[0],
		p = Xdims[1],
		/* size of triangular part of decomposition : */
		trsz = (n < p) ? n : p,
		i;
	SET_VECTOR_ELT(ans, 0, X = duplicate(Xin));
	SET_VECTOR_ELT(ans, 2, qraux = allocVector(REALSXP, trsz));
	SET_VECTOR_ELT(ans, 3, pivot = allocVector(INTSXP, p));
	for (i = 0; i < p; i++)
		INTEGER(pivot)[i] = i + 1;
	SEXP nms, Givens = PROTECT(allocVector(VECSXP, trsz - 1));
	setAttrib(ans, R_NamesSymbol, nms = allocVector(STRSXP, 5));
	SET_STRING_ELT(nms, 0, mkChar("qr"));
	SET_STRING_ELT(nms, 1, mkChar("rank"));
	SET_STRING_ELT(nms, 2, mkChar("qraux"));
	SET_STRING_ELT(nms, 3, mkChar("pivot"));
	SET_STRING_ELT(nms, 4, mkChar("Givens"));
	int rank = trsz, nGivens = 0;
	double rcond = 0.0;
	if (n > 0 && p > 0) {
		double *xpt = REAL(X), *work, tmp;
		int *iwork, lwork, info;
		lwork = -1;
		F77_CALL(dgeqrf)(&n, &p, xpt, &n, REAL(qraux), &tmp, &lwork,
		                 &info);
		if (info)
			error(_("LAPACK dgeqrf returned error code %d"), info);
		lwork = (int) tmp;
		work = (double *) R_alloc(((size_t) lwork < (size_t) 3 * trsz)
		                          ? (size_t) 3 * trsz : (size_t) lwork,
		                          sizeof(double));
		F77_CALL(dgeqrf)(&n, &p, xpt, &n, REAL(qraux), work, &lwork,
		                 &info);
		if (info)
			error(_("LAPACK dgeqrf returned error code %d"), info);
		iwork = (int *) R_alloc((size_t) trsz, sizeof(int));
		F77_CALL(dtrcon)("1", "U", "N", &rank, xpt, &n, &rcond,
		                 work, iwork, &info FCONE FCONE FCONE);
		if (info)
			error(_("LAPACK dtrcon returned error code %d"), info);
		while (rcond < tol) { /* check diagonal elements */
			double el, minabs = (xpt[0] < 0.0) ? -xpt[0]: xpt[0];
			int jmin = 0;
			for (i = 1; i < rank; i++) {
				el = xpt[i*n]; // had  i*(n+1)  which looks wrong to MM
				if (el < 0.0)
					el = -el;
				if (el < minabs) {
					jmin = i;
					minabs = el;
				}
			}
			if (jmin < (rank - 1)) {
				SET_VECTOR_ELT(Givens, nGivens,
				               getGivens(xpt, n, jmin, rank));
				nGivens++;
			} // otherwise jmin == (rank - 1), so just "drop that column"
			rank--;
			// new  rcond := ... for reduced rank
			F77_CALL(dtrcon)("1", "U", "N", &rank, xpt, &n, &rcond,
			                 work, iwork, &info FCONE FCONE FCONE);
			if (info)
				error(_("LAPACK dtrcon returned error code %d"), info);
		}
	}
	SEXP Gcpy;
	SET_VECTOR_ELT(ans, 1, ScalarInteger(rank));
	SET_VECTOR_ELT(ans, 4, Gcpy = allocVector(VECSXP, nGivens));
	for (i = 0; i < nGivens; i++)
		SET_VECTOR_ELT(Gcpy, i, VECTOR_ELT(Givens, i));
	setAttrib(ans, install("useLAPACK"), ScalarLogical(1));
	setAttrib(ans, install("rcond"), ScalarReal(rcond));
	UNPROTECT(2);
	return ans;
}

#endif /* MJ */
