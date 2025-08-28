#include "Mdefines.h"
#include "idz.h"
#include "dense.h"

SEXP dense_band(SEXP from, const char *class, int a, int b)
{
	/* defined in ./coerce.c : */
	SEXP dense_as_general(SEXP, const char *, int);

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */

	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) */
	/* to be triangularMatrix                       */
	/* Don't return 'from' if from@x has attributes */
	/* (notably 'class' or 'dim'), which can happen */
	/* if from = matrix_as_dense(..., new=0).       */
	if ((m == 0 || n == 0 || (a <= 1 - m && b >= n - 1)) &&
	    (m != n || n > 1 || class[1] == 't') &&
	    !ANY_ATTRIB(GET_SLOT(from, Matrix_xSym)))
		return from;

	int ge = 0, sy = 0, tr = 0;
	ge = m != n || !((tr = a >= 0 || b <= 0 || class[1] == 't') ||
	                 (sy = a == -b && class[1] == 's'));

#define BAND_CASES(_DO_) \
	do { \
		switch (class[0]) { \
		case 'n': \
		case 'l': \
			_DO_(i, int, LOGICAL); \
			break; \
		case 'i': \
			_DO_(i, int, INTEGER); \
			break; \
		case 'd': \
			_DO_(d, double, REAL); \
			break; \
		case 'z': \
			_DO_(z, Rcomplex, COMPLEX); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define BAND2(_PREFIX_, _CTYPE_, _PTR_) \
	_PREFIX_ ## band2(_PTR_(x1), m, n, a, b, di)

#define BAND1(_PREFIX_, _CTYPE_, _PTR_) \
	_PREFIX_ ## band1(_PTR_(x1), n, a, b, ul1, di)

#define DCPY2(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		Matrix_memset(px1, 0, XLENGTH(x1), sizeof(_CTYPE_)); \
		if (a <= 0 && b >= 0) \
			_PREFIX_ ## dcpy2(px1, px0, n, XLENGTH(x1),      'U', di); \
	} while (0)

#define DCPY1(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		Matrix_memset(px1, 0, XLENGTH(x1), sizeof(_CTYPE_)); \
		if (a <= 0 && b >= 0) \
			_PREFIX_ ## dcpy1(px1, px0, n, XLENGTH(x1), ul1, ul0, di); \
	} while (0)

	char ul0 = 'U', ul1 = 'U', di = 'N';
	if (class[1] != 'g') {
		if (ge) {
			PROTECT(from = dense_as_general(from, class, 1));
			SEXP x1 = PROTECT(GET_SLOT(from, Matrix_xSym));
			BAND_CASES(BAND2);
			UNPROTECT(2); /* x1, from */
			return from;
		}

		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul0 = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		if (class[1] == 't') {
			/* Be fast if band contains entire triangle */
			if ((ul0 == 'U')
			    ? (a <= 0 && b >= n - 1) : (b >= 0 && a <= 1 - m))
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

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1;

	if (ge) {
		PROTECT(x1 = duplicate(x0));
		if (ANY_ATTRIB(x1))
			CLEAR_ATTRIB(x1);
		SET_SLOT(to, Matrix_xSym, x1);
		BAND_CASES(BAND2);
		UNPROTECT(3); /* x1, x0, to */
		return to;
	}

	/* Returning .(sy|sp|tr|tp)Matrix ... */

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

	if (tr && class[1] == 't') {
		if ((ul0 == 'U') ? (b <= 0) : (a >= 0)) {
			/* Result is either a diagonal matrix or a zero matrix : */
			PROTECT(x1 = allocVector(TYPEOF(x0), XLENGTH(x0)));
			if (class[2] != 'p')
				BAND_CASES(DCPY2);
			else
				BAND_CASES(DCPY1);
		} else {
			PROTECT(x1 = duplicate(x0));
			if (class[2] != 'p')
				BAND_CASES(BAND2);
			else
				BAND_CASES(BAND1);
		}
	} else {
		if (sy || (tr && (class[1] == 'g' || ul0 == ul1 || n <= 1))) {
			PROTECT(x1 = duplicate(x0));
			if (ANY_ATTRIB(x1))
				CLEAR_ATTRIB(x1);
		} else {
			/* Band is "opposite" the stored triangle : */
			PROTECT(from = dense_transpose(from, class));
			x1 = GET_SLOT(from, Matrix_xSym);
			UNPROTECT(1);
			PROTECT(x1);
		}
		if (class[2] != 'p')
			BAND_CASES(BAND2);
		else
			BAND_CASES(BAND1);
	}
	SET_SLOT(to, Matrix_xSym, x1);

#undef BAND_CASES
#undef BAND2
#undef BAND1
#undef DCPY2
#undef DCPY1

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* band(<denseMatrix>, k1, k2), tri[ul](<denseMatrix>, k) */
/* NB: argument validation more or less copied by R_sparse_band() */
SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2)
{
	if (!isS4(from)) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, int, int);
		from = matrix_as_dense(from, ".ge", '\0', '\0', 0, 0);
	}
	PROTECT(from);
	static const char *valid[] = { VALID_DENSE, "" };
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

	from = dense_band(from, valid[ivalid], a, b);
	UNPROTECT(1);
	return from;
}

SEXP dense_diag_get(SEXP obj, const char *class, int names)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n, j;
	UNPROTECT(1); /* dim */

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		if (class[2] == 'p') {
			SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
			ul = *CHAR(STRING_ELT(uplo, 0));
			UNPROTECT(1); /* uplo */
		}
		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	}

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		res = PROTECT(allocVector(TYPEOF(x), r));

#define DG_LOOP(_CTYPE_, _PTR_, _ONE_) \
	do { \
		_CTYPE_ *pres = _PTR_(res), *px = _PTR_(x); \
		if (di == 'U') \
			for (j = 0; j < r; ++j) \
				*(pres++) = _ONE_; \
		else if (class[2] != 'p') { \
			R_xlen_t m1a = (R_xlen_t) m + 1; \
			for (j = 0; j < r; ++j, px += m1a) \
				*(pres++) = *px; \
		} \
		else if (ul == 'U') \
			for (j = 0; j < n; px += (++j) + 1) \
				*(pres++) = *px; \
		else \
			for (j = 0; j < n; px += n - (j++)) \
				*(pres++) = *px; \
	} while (0)

	switch (class[0]) {
	case 'n':
	case 'l':
		DG_LOOP(int, LOGICAL, 1);
		break;
	case 'i':
		DG_LOOP(int, INTEGER, 1);
		break;
	case 'd':
		DG_LOOP(double, REAL, 1.0);
		break;
	case 'z':
		DG_LOOP(Rcomplex, COMPLEX, Matrix_zone);
		break;
	default:
		break;
	}

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

#undef DG_LOOP

	UNPROTECT(2); /* x, res */
	return res;
}

SEXP R_dense_diag_get(SEXP obj, SEXP names)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int names_;
	if (TYPEOF(names) != LGLSXP || LENGTH(names) < 1 ||
	    (names_ = LOGICAL(names)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "names", "TRUE", "FALSE");

	return dense_diag_get(obj, valid[ivalid], names_);
}

SEXP dense_diag_set(SEXP from, const char *class, SEXP value, int new)
{
	SEXP to = PROTECT(newObject(class));
	int v = LENGTH(value) != 1;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n, j;
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (new) {
		x = duplicate(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
	}
	SET_SLOT(to, Matrix_xSym, x);

#define DS_LOOP(_CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pvalue = _PTR_(value); \
		if (class[2] != 'p') { \
			R_xlen_t m1a = (R_xlen_t) m + 1; \
			if (v) \
				for (j = 0; j < r; ++j, px += m1a) \
					*px = *(pvalue++); \
			else \
				for (j = 0; j < r; ++j, px += m1a) \
					*px = *pvalue; \
		} else if (ul == 'U') { \
			if (v) \
				for (j = 0; j < n; px += (++j) + 1) \
					*px = *(pvalue++); \
			else \
				for (j = 0; j < n; px += (++j) + 1) \
					*px = *pvalue; \
		} else { \
			if (v) \
				for (j = 0; j < n; px += n - (j++)) \
					*px = *(pvalue++); \
			else \
				for (j = 0; j < n; px += n - (j++)) \
					*px = *pvalue; \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
	case 'l':
		DS_LOOP(int, LOGICAL);
		break;
	case 'i':
		DS_LOOP(int, INTEGER);
		break;
	case 'd':
		DS_LOOP(double, REAL);
		break;
	case 'z':
		DS_LOOP(Rcomplex, COMPLEX);
		break;
	default:
		break;
	}

#undef DS_LOOP

	UNPROTECT(2); /* x, to */
	return to;
}

SEXP R_dense_diag_set(SEXP from, SEXP value)
{
	static const char *valid[] = { VALID_DENSE, "" };
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

	int new = 1;
	if (tv <= tx) {
		PROTECT(from);
		PROTECT(value = coerceVector(value, tx));
	} else {
		/* defined in ./coerce.c : */
		SEXP dense_as_kind(SEXP, const char *, char, int);
#ifndef MATRIX_ENABLE_IMATRIX
		if (tv == INTSXP) {
		PROTECT(from = dense_as_kind(from, class, 'd', 0));
		PROTECT(value = coerceVector(value, REALSXP));
		} else {
#endif
		PROTECT(from = dense_as_kind(from, class, typeToKind(tv), 0));
		PROTECT(value);
#ifndef MATRIX_ENABLE_IMATRIX
		}
#endif
		class = valid[R_check_class_etc(from, valid)];
		new = 0;
	}

	from = dense_diag_set(from, class, value, new);
	UNPROTECT(2);
	return from;
}

SEXP dense_transpose(SEXP from, const char *class)
{
	SEXP to = PROTECT(newObject(class));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], i, j;
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
	if (class[1] == 's' || class[1] == 'p' || class[1] == 'o')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_reversed_DimNames(to, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
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

			if (class[1] == 'o' && n > 0) {
				SEXP sd = PROTECT(GET_SLOT(from, Matrix_sdSym));
				SET_SLOT(to, Matrix_sdSym, sd);
				UNPROTECT(1); /* sd */
			}
		}
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), XLENGTH(x0)));
	SET_SLOT(to, Matrix_xSym, x1);

#define TRANS_LOOP(_CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (class[2] != 'p') { \
			R_xlen_t mn1s = XLENGTH(x0) - 1; \
			for (j = 0; j < m; ++j, px0 -= mn1s) \
				for (i = 0; i < n; ++i, px0 += m) \
					*(px1++) = *px0; \
		} else if (ul == 'U') { \
			for (j = 0; j < n; ++j) \
				for (i = j; i < n; ++i) \
					*(px1++) = *(px0 + PACKED_AR21_UP(j, i)); \
		} else { \
			R_xlen_t n2 = (R_xlen_t) n * 2; \
			for (j = 0; j < n; ++j) \
				for (i = 0; i <= j; ++i) \
					*(px1++) = *(px0 + PACKED_AR21_LO(j, i, n2)); \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
	case 'l':
		TRANS_LOOP(int, LOGICAL);
		break;
	case 'i':
		TRANS_LOOP(int, INTEGER);
		break;
	case 'c':
	case 'd':
		TRANS_LOOP(double, REAL);
		break;
	case 'z':
		TRANS_LOOP(Rcomplex, COMPLEX);
		break;
	default:
		break;
	}

#undef TRANS_LOOP

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

SEXP R_dense_transpose(SEXP from)
{
	static const char *valid[] = {
		"dpoMatrix", "dppMatrix", "corMatrix", "copMatrix",
		VALID_DENSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return dense_transpose(from, valid[ivalid]);
}

SEXP dense_force_symmetric(SEXP from, const char *class, char ul)
{
	char ul0 = 'U', ul1 = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul0 = ul1 = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */
		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	}

	if (ul != '\0')
		ul1 = ul;

	if (class[1] == 's') {
		/* .s[yp]Matrix */
		if (ul0 == ul1)
			return from;
		SEXP to = PROTECT(dense_transpose(from, class));
		if (class[0] == 'z') {
			/* Need _conjugate_ transpose */
			SEXP x1 = PROTECT(GET_SLOT(to, Matrix_xSym));
			conjugate(x1);
			UNPROTECT(1); /* x1 */
		}
		UNPROTECT(1) /* to */;
		return to;
	}

	/* Now handling just .(ge|tr|tp)Matrix ... */

	char cl[] = ".s.Matrix";
	cl[0] = class[0];
	cl[2] = (class[2] != 'p') ? 'y' : 'p';
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

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));

	if (class[1] == 'g' || ul0 == ul1)
		SET_SLOT(to, Matrix_xSym, x0);
	else {
		SEXP x1 = PROTECT(allocVector(TYPEOF(x0), XLENGTH(x0)));
		SET_SLOT(to, Matrix_xSym, x1);

		R_xlen_t len = XLENGTH(x1);

#define DCPY(_PREFIX_, _CTYPE_, _PTR_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			Matrix_memset(px1, 0, len, sizeof(_CTYPE_)); \
			if (class[2] != 'p') \
				_PREFIX_ ## dcpy2(px1, px0, n, len,     '\0', di); \
			else \
				_PREFIX_ ## dcpy1(px1, px0, n, len, ul1, ul0, di); \
		} while (0)

		switch (class[0]) {
		case 'n':
		case 'l':
			DCPY(i, int, LOGICAL);
			break;
		case 'i':
			DCPY(i, int, INTEGER);
			break;
		case 'd':
			DCPY(d, double, REAL);
			break;
		case 'z':
			DCPY(z, Rcomplex, COMPLEX);
			break;
		default:
			break;
		}

#undef DCPY

		UNPROTECT(1); /* x1 */
	}

	UNPROTECT(2); /* x0, to */
	return to;
}

SEXP R_dense_force_symmetric(SEXP from, SEXP uplo)
{
	static const char *valid[] = { VALID_DENSE, "" };
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

	return dense_force_symmetric(from, valid[ivalid], ul);
}

SEXP dense_symmpart(SEXP from, const char *class)
{
	if (class[0] != 'z' && class[0] != 'd') {
		/* defined in ./coerce.c : */
		SEXP dense_as_kind(SEXP, const char *, char, int);
		from = dense_as_kind(from, class, 'd', 0);
	}
	if (class[0] != 'z' && class[1] == 's')
		return from;
	PROTECT(from);

	char cl[] = ".s.Matrix";
	cl[0] = (class[0] != 'z') ? 'd' : 'z';
	cl[2] = (class[2] != 'p') ? 'y' : 'p';
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
	}

	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (class[0] == 'z' || class[0] == 'd') {
		x = duplicate(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
	}
	SET_SLOT(to, Matrix_xSym, x);

	if (class[1] == 's') {
		/* Symmetric part of Hermitian matrix is real part */
		zeroIm(x);
		UNPROTECT(3); /* x, to, from */
		return to;
	}

	int i, j;

#define SP_LOOP(_CTYPE_, _PTR_, _INCREMENT_, _SCALE1_, _ONE_) \
	do { \
		_CTYPE_ *px = _PTR_(x); \
		if (class[1] == 'g') { \
			_CTYPE_ *py = px; \
			for (j = 0; j < n; ++j) { \
				for (i = j + 1; i < n; ++i) { \
					px += n; \
					py += 1; \
					_INCREMENT_((*px), (*py)); \
					_SCALE1_((*px), 0.5); \
				} \
				px = (py += j + 2); \
			} \
		} else if (class[2] != 'p') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) { \
						_SCALE1_((*px), 0.5); \
						px += 1; \
					} \
					px += n - j; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					px += j + 1; \
					for (i = j + 1; i < n; ++i) { \
						_SCALE1_((*px), 0.5); \
						px += 1; \
					} \
				} \
			} \
			if (di != 'N') { \
				R_xlen_t n1a = (R_xlen_t) n + 1; \
				px = _PTR_(x); \
				for (j = 0; j < n; ++j, px += n1a) \
					*px = _ONE_; \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) { \
						_SCALE1_((*px), 0.5); \
						px += 1; \
					} \
					px += 1; \
				} \
				if (di != 'N') { \
					px = _PTR_(x); \
					for (j = 0; j < n; px += (++j) + 1) \
						*px = _ONE_; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					px += 1; \
					for (i = j + 1; i < n; ++i) { \
						_SCALE1_((*px), 0.5); \
						px += 1; \
					} \
				} \
				if (di != 'N') { \
					px = _PTR_(x); \
					for (j = 0; j < n; px += n - (j++)) \
						*px = _ONE_; \
				} \
			} \
		} \
	} while (0)

	if (cl[0] == 'd')
		SP_LOOP(double, REAL, INCREMENT_REAL, SCALE1_REAL, 1.0);
	else
		SP_LOOP(Rcomplex, COMPLEX, INCREMENT_COMPLEX, SCALE1_COMPLEX, Matrix_zone);

#undef SP_LOOP

	UNPROTECT(3); /* x, to, from */
	return to;
}

SEXP R_dense_symmpart(SEXP from)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return dense_symmpart(from, valid[ivalid]);
}

SEXP dense_skewpart(SEXP from, const char *class)
{
	if (class[0] != 'z' && class[0] != 'd') {
		/* defined in ./coerce.c : */
		SEXP dense_as_kind(SEXP, const char *, char, int);
		from = dense_as_kind(from, class, 'd', 0);
	}
	PROTECT(from);

	char cl[] = "...Matrix";
	cl[0] = (class[0] != 'z') ? 'd' : 'z';
	cl[1] = (class[1] != 's') ? 'g' : 's';
	cl[2] = (class[1] != 's') ? 'e' :
		((class[0] != 'z') ? 'C' : ((class[2] != 'p') ? 'y' : 'p'));
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

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (class[1] == 's' && ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	if (class[1] == 's' && class[0] != 'z') {
		/* Skew-symmetric part of Hermitian matrix is imaginary part */
		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1));
		int *pp = INTEGER(p);
		Matrix_memset(pp, 0, (R_xlen_t) n + 1, sizeof(int));
		SET_SLOT(to, Matrix_pSym, p);
		UNPROTECT(3); /* p, to, from */
		return to;
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1 = x0;

	if (class[1] == 's') {
		/* Skew-symmetric part of Hermitian matrix is imaginary part */
		x1 = duplicate(x1);
		UNPROTECT(1); /* x1 */
		PROTECT(x1);
		SET_SLOT(to, Matrix_xSym, x1);
		zeroRe(x1);
		UNPROTECT(3); /* x1, to, from */
		return to;
	}

	if (class[2] == 'p' || class[0] == 'z' || class[0] == 'd') {
		if ((Matrix_int_fast64_t) n * n > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");
		x1 = allocVector(TYPEOF(x0), (R_xlen_t) n * n);
	}
	PROTECT(x1);
	SET_SLOT(to, Matrix_xSym, x1);

	int i, j;
	R_xlen_t upos = 0, lpos = 0;

#define SP_LOOP(_CTYPE_, _PTR_, _INCREMENT_, _ASSIGN_, _ZERO_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				lpos = j; \
				for (i = 0; i < j; ++i) { \
					_ASSIGN_(px1[upos], 0.5 * px0[upos]); \
					_INCREMENT_(px1[upos], -0.5 * px0[lpos]); \
					_ASSIGN_(px1[lpos], -px1[upos]); \
					upos += 1; \
					lpos += n; \
				} \
				px1[upos] = _ZERO_; \
				upos += n - j; \
			} \
		} else if (class[2] != 'p') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					lpos = j; \
					for (i = 0; i < j; ++i) { \
						_ASSIGN_(px1[upos], 0.5 * px0[upos]); \
						_ASSIGN_(px1[lpos], -px1[upos]); \
						upos += 1; \
						lpos += n; \
					} \
					px1[upos] = _ZERO_; \
					upos += n - j; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					upos = lpos; \
					px1[lpos] = _ZERO_; \
					for (i = j + 1; i < n; ++i) { \
						upos += n; \
						lpos += 1; \
						_ASSIGN_(px1[lpos], 0.5 * px0[lpos]); \
						_ASSIGN_(px1[upos], -px1[lpos]); \
					} \
					lpos += j + 2; \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j, ++px0) { \
					lpos = j; \
					for (i = 0; i < j; ++i, ++px0) { \
						_ASSIGN_(px1[upos], 0.5 * (*px0)); \
						_ASSIGN_(px1[lpos], -px1[upos]); \
						upos += 1; \
						lpos += n; \
					} \
					px1[upos] = _ZERO_; \
					upos += n - j; \
				} \
			} else { \
				for (j = 0; j < n; ++j, ++px0) { \
					upos = lpos; \
					px1[lpos] = _ZERO_; \
					for (i = j + 1; i < n; ++i, ++px0) { \
						upos += n; \
						lpos += 1; \
						_ASSIGN_(px1[lpos], 0.5 * (*px0)); \
						_ASSIGN_(px1[upos], -px1[lpos]); \
					} \
					lpos += j + 2; \
				} \
			} \
		} \
	} while (0)

	if (cl[0] == 'd')
		SP_LOOP(double, REAL, INCREMENT_REAL, ASSIGN_REAL, 0.0);
	else
		SP_LOOP(Rcomplex, COMPLEX, INCREMENT_COMPLEX, ASSIGN_COMPLEX, Matrix_zzero);

#undef SP_LOOP

	UNPROTECT(4); /* x1, x0, to, from */
	return to;
}

SEXP R_dense_skewpart(SEXP from)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return dense_skewpart(from, valid[ivalid]);
}

int dense_is_symmetric(SEXP obj, const char *class, int checkDN)
{
	if (class[1] == 's')
		return 1;

	if (checkDN) {
		SEXP dimnames = GET_SLOT(obj, Matrix_DimNamesSym);
		if (!DimNames_is_symmetric(dimnames))
			return 0;
	}

	if (class[1] == 't')
		return dense_is_diagonal(obj, class);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return 1;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j;

#define IS_LOOP(_CTYPE_, _PTR_, _NOTREAL_, _NOTCONJ_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *py = px; \
		for (j = 0; j < n; px = (py += (++j) + 1)) { \
			if (_NOTREAL_((*px))) \
				return 0; \
			for (i = j + 1; i < n; ++i) { \
				px += n; \
				py += 1; \
				if (_NOTCONJ_((*px), (*py))) \
					return 0; \
			} \
		} \
		return 1; \
	} while (0)

	switch (class[0]) {
	case 'n':
		IS_LOOP(int, LOGICAL, NOTREAL_PATTERN, NOTCONJ_PATTERN);
		break;
	case 'l':
		IS_LOOP(int, LOGICAL, NOTREAL_LOGICAL, NOTCONJ_LOGICAL);
		break;
	case 'i':
		IS_LOOP(int, INTEGER, NOTREAL_INTEGER, NOTCONJ_INTEGER);
		break;
	case 'd':
		IS_LOOP(double, REAL, NOTREAL_REAL, NOTCONJ_REAL);
		break;
	case 'z':
		IS_LOOP(Rcomplex, COMPLEX, NOTREAL_COMPLEX, NOTCONJ_COMPLEX);
		break;
	default:
		break;
	}

#undef IS_LOOP

	return 0;
}

SEXP R_dense_is_symmetric(SEXP obj, SEXP checkDN)
{
	if (!isS4(obj)) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, int, int);
		obj = matrix_as_dense(obj, ".ge", '\0', '\0', 0, 0);
	}
	PROTECT(obj);
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int checkDN_;
	if (TYPEOF(checkDN) != LGLSXP || LENGTH(checkDN) < 1 ||
	    (checkDN_ = LOGICAL(checkDN)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "checkDN", "TRUE", "FALSE");

	SEXP ans = ScalarLogical(dense_is_symmetric(obj, valid[ivalid], checkDN_));
	UNPROTECT(1);
	return ans;
}

int dense_is_triangular(SEXP obj, const char *class, int upper)
{
	if (class[1] == 't') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (upper == NA_LOGICAL || (upper != 0) == (ul == 'U'))
			return (ul == 'U') ? 1 : -1;
		else if (dense_is_diagonal(obj, class))
			return (ul == 'U') ? -1 : 1;
		else
			return 0;
	}

	if (class[1] == 's') {
		if (!dense_is_diagonal(obj, class))
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

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j;

#define IT_LOOP(_CTYPE_, _PTR_, _ISNZ_) \
	do { \
		_CTYPE_ *px; \
		if (upper == NA_LOGICAL) { \
			px = _PTR_(x); \
			for (j = 0; j < n; px += (++j)) { \
				px += 1; \
				for (i = j + 1; i < n; ++i, px += 1) { \
					if (_ISNZ_(*px)) { \
						j = n; \
						break; \
					} \
				} \
			} \
			if (j == n) \
				return  1; \
			px = _PTR_(x); \
			for (j = 0; j < n; px += n - (++j)) { \
				for (i = 0; i < j; ++i, px += 1) { \
					if (_ISNZ_(*px)) { \
						j = n; \
						break; \
					} \
				} \
				px += 1; \
			} \
			if (j == n) \
				return -1; \
			return 0; \
		} else if (upper != 0) { \
			px = _PTR_(x); \
			for (j = 0; j < n; px += (++j)) { \
				px += 1; \
				for (i = j + 1; i < n; ++i, px += 1) \
					if (_ISNZ_(*px)) \
						return 0; \
			} \
			return  1; \
		} else { \
			px = _PTR_(x); \
			for (j = 0; j < n; px += n - (++j)) { \
				for (i = 0; i < j; ++i, px += 1) \
					if (_ISNZ_(*px)) \
						return 0; \
				px += 1; \
			} \
			return -1; \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
		IT_LOOP(int, LOGICAL, ISNZ_PATTERN);
		break;
	case 'l':
		IT_LOOP(int, LOGICAL, ISNZ_LOGICAL);
		break;
	case 'i':
		IT_LOOP(int, INTEGER, ISNZ_INTEGER);
		break;
	case 'd':
		IT_LOOP(double, REAL, ISNZ_REAL);
		break;
	case 'z':
		IT_LOOP(Rcomplex, COMPLEX, ISNZ_COMPLEX);
		break;
	default:
		break;
	}

#undef IT_LOOP

	return 0;
}

SEXP R_dense_is_triangular(SEXP obj, SEXP upper)
{
	if (!isS4(obj)) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, int, int);
		obj = matrix_as_dense(obj, ".ge", '\0', '\0', 0, 0);
	}
	PROTECT(obj);
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	if (TYPEOF(upper) != LGLSXP || LENGTH(upper) < 1)
		error(_("'%s' must be %s or %s or %s"), "upper", "TRUE", "FALSE", "NA");
	int upper_ = LOGICAL(upper)[0];

	int ans_ = dense_is_triangular(obj, valid[ivalid], upper_);
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
	UNPROTECT(1);
	return ans;
}

int dense_is_diagonal(SEXP obj, const char *class)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return 1;

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = *CHAR(STRING_ELT(uplo, 0));
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j;

#define ID_LOOP(_CTYPE_, _PTR_, _ISNZ_) \
	do { \
		_CTYPE_ *px = _PTR_(x); \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (_ISNZ_(*px)) \
						return 0; \
					px += 1; \
				} \
				px += 1; \
				for (i = j + 1; i < n; ++i) { \
					if (_ISNZ_(*px)) \
						return 0; \
					px += 1; \
				} \
			} \
		} else if (class[2] != 'p') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) { \
						if (_ISNZ_(*px)) \
							return 0; \
						px += 1; \
					} \
					px += n - j; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					px += j + 1; \
					for (i = j + 1; i < n; ++i) { \
						if (_ISNZ_(*px)) \
							return 0; \
						px += 1; \
					} \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) { \
						if (_ISNZ_(*px)) \
							return 0; \
						px += 1; \
					} \
					px += 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					px += 1; \
					for (i = j + 1; i < n; ++i) { \
						if (_ISNZ_(*px)) \
							return 0; \
						px += 1; \
					} \
				} \
			} \
		} \
		return 1; \
	} while (0)

	switch (class[0]) {
	case 'n':
		ID_LOOP(int, LOGICAL, ISNZ_PATTERN);
		break;
	case 'l':
		ID_LOOP(int, LOGICAL, ISNZ_LOGICAL);
		break;
	case 'i':
		ID_LOOP(int, INTEGER, ISNZ_INTEGER);
		break;
	case 'd':
		ID_LOOP(double, REAL, ISNZ_REAL);
		break;
	case 'z':
		ID_LOOP(Rcomplex, COMPLEX, ISNZ_COMPLEX);
		break;
	default:
		break;
	}

#undef ID_LOOP

	return 0;
}

SEXP R_dense_is_diagonal(SEXP obj)
{
	if (!isS4(obj)) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, int, int);
		obj = matrix_as_dense(obj, ".ge", '\0', '\0', 0, 0);
	}
	PROTECT(obj);
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	SEXP ans = ScalarLogical(dense_is_diagonal(obj, valid[ivalid]));
	UNPROTECT(1);
	return ans;
}

#define CAST_PATTERN(_X_) (_X_ != 0)
#define CAST_LOGICAL(_X_) (_X_ != 0)
#define CAST_INTEGER(_X_)  _X_
#define CAST_REAL(_X_)     _X_
#define CAST_COMPLEX(_X_)  _X_

#define SUM_CASES \
do { \
	switch (class[0]) { \
	case 'n': \
		if (mean) \
		SUM_LOOP(int, LOGICAL, double, REAL, \
		         0.0, 1.0, NA_REAL, ISNA_PATTERN, \
		         CAST_PATTERN, INCREMENT_REAL, SCALE2_REAL); \
		else \
		SUM_LOOP(int, LOGICAL, int, INTEGER, \
		         0, 1, NA_INTEGER, ISNA_PATTERN, \
		         CAST_PATTERN, INCREMENT_INTEGER, SCALE2_REAL); \
		break; \
	case 'l': \
		if (mean) \
		SUM_LOOP(int, LOGICAL, double, REAL, \
		         0.0, 1.0, NA_REAL, ISNA_LOGICAL, \
		         CAST_LOGICAL, INCREMENT_REAL, SCALE2_REAL); \
		else \
		SUM_LOOP(int, LOGICAL, int, INTEGER, \
		         0, 1, NA_INTEGER, ISNA_LOGICAL, \
		         CAST_LOGICAL, INCREMENT_INTEGER, SCALE2_REAL); \
		break; \
	case 'i': \
		SUM_LOOP(int, INTEGER, double, REAL, \
		         0.0, 1.0, NA_REAL, ISNA_INTEGER, \
		         CAST_INTEGER, INCREMENT_REAL, SCALE2_REAL); \
		break; \
	case 'd': \
		SUM_LOOP(double, REAL, double, REAL, \
		         0.0, 1.0, NA_REAL, ISNA_REAL, \
		         CAST_REAL, INCREMENT_REAL, SCALE2_REAL); \
		break; \
	case 'z': \
		SUM_LOOP(Rcomplex, COMPLEX, Rcomplex, COMPLEX, \
		         Matrix_zzero, Matrix_zone, Matrix_zna, ISNA_COMPLEX, \
		         CAST_COMPLEX, INCREMENT_COMPLEX, SCALE2_COMPLEX); \
		break; \
	default: \
		break; \
	} \
} while (0)

#define SUM_TYPEOF(c) (c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

static
void dense_colsum(SEXP x, const char *class,
                  int m, int n, char ul, char di, int narm, int mean,
                  SEXP res)
{
	int i, j, count = -1, narm_ = narm && mean && class[0] != 'n',
		unpacked = class[2] != 'p';

#define SUM_LOOP(_CTYPE0_, _PTR0_, _CTYPE1_, _PTR1_, \
	             _ZERO_, _ONE_, _NA_, _ISNA_, \
	             _CAST_, _INCREMENT_, _SCALE2_) \
	do { \
		_CTYPE0_ *px0 = _PTR0_(  x); \
		_CTYPE1_ *px1 = _PTR1_(res), tmp; \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				*px1 = _ZERO_; \
				SUM_KERNEL(for (i = 0; i < m; ++i), _NA_, _ISNA_, \
				           _CAST_, _INCREMENT_, _SCALE2_); \
				px1 += 1; \
			} \
		} else if (di == 'N') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					*px1 = _ZERO_; \
					SUM_KERNEL(for (i = 0; i <= j; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_, _SCALE2_); \
					if (unpacked) \
						px0 += n - j - 1; \
					px1 += 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (unpacked) \
						px0 += j; \
					*px1 = _ZERO_; \
					SUM_KERNEL(for (i = j; i < n; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_, _SCALE2_); \
					px1 += 1; \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					*px1 = _ONE_; \
					SUM_KERNEL(for (i = 0; i < j; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_, _SCALE2_); \
					++px0; \
					if (unpacked) \
						px0 += n - j - 1; \
					px1 += 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (unpacked) \
						px0 += j; \
					++px0; \
					*px1 = _ONE_; \
					SUM_KERNEL(for (i = j + 1; i < n; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_, _SCALE2_); \
					px1 += 1; \
				} \
			} \
		} \
	} while (0)

#define SUM_KERNEL(_FOR_, _NA_, _ISNA_, _CAST_, _INCREMENT_, _SCALE2_) \
	do { \
		if (mean) \
			count = m; \
		_FOR_ { \
			if (_ISNA_(*px0)) { \
				if (!narm) \
					*px1 = _NA_; \
				else if (narm_) \
					--count; \
			} else { \
				tmp = _CAST_(*px0); \
				_INCREMENT_((*px1), tmp); \
			} \
			++px0; \
		} \
		if (mean) \
			_SCALE2_((*px1), count); \
	} while (0)

	SUM_CASES;

#undef SUM_LOOP
#undef SUM_KERNEL

	return;
}

static
void dense_rowsum(SEXP x, const char *class,
                  int m, int n, char ul, char di, int narm, int mean,
                  SEXP res)
{
	int i, j, *count = NULL, narm_ = narm && mean && class[0] != 'n',
		unpacked = class[2] != 'p', symmetric = class[1] == 's';
	if (narm_) {
		Matrix_Calloc(count, m, int);
		for (i = 0; i < m; ++i)
			count[i] = n;
	}

#define SUM_LOOP(_CTYPE0_, _PTR0_, _CTYPE1_, _PTR1_, \
		         _ZERO_, _ONE_, _NA_, _ISNA_, \
		         _CAST_, _INCREMENT_, _SCALE2_) \
	do { \
		_CTYPE0_ *px0 = _PTR0_(  x); \
		_CTYPE1_ *px1 = _PTR1_(res), tmp = (di == 'N') ? _ZERO_ : _ONE_; \
		for (i = 0; i < m; ++i) \
			px1[i] = tmp; \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(for (i = 0; i < m; ++i), _NA_, _ISNA_, \
				           _CAST_, _INCREMENT_); \
		} else if (class[1] == 's' || di == 'N') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					SUM_KERNEL(for (i = 0; i <= j; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_); \
					if (unpacked) \
						px0 += n - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (unpacked) \
						px0 += j; \
					SUM_KERNEL(for (i = j; i < n; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_); \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					SUM_KERNEL(for (i = 0; i < j; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_); \
					++px0; \
					if (unpacked) \
						px0 += n - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (unpacked) \
						px0 += j; \
					++px0; \
					SUM_KERNEL(for (i = j + 1; i < n; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_); \
				} \
			} \
		} \
		if (mean) { \
			if (narm_) \
				for (i = 0; i < m; ++i) \
					_SCALE2_(px1[i], count[i]); \
			else \
				for (i = 0; i < m; ++i) \
					_SCALE2_(px1[i], n); \
		} \
	} while (0)

#define SUM_KERNEL(_FOR_, _NA_, _ISNA_, _CAST_, _INCREMENT_) \
	do { \
		_FOR_ { \
			int again = symmetric && i != j; \
			if (_ISNA_(*px0)) { \
				if (!narm) { \
					px1[i] = _NA_; \
					if (again) \
					px1[j] = _NA_; \
				} else if (narm_) { \
					--count[i]; \
					if (again) \
					--count[j]; \
				} \
			} else { \
				tmp = _CAST_(*px0); \
				_INCREMENT_(px1[i], tmp); \
				if (again) \
				_INCREMENT_(px1[j], tmp); \
			} \
			++px0; \
		} \
	} while (0)

	SUM_CASES;

#undef SUM_LOOP
#undef SUM_KERNEL

	if (narm_)
		Matrix_Free(count, m);
	return;
}

SEXP dense_marginsum(SEXP obj, const char *class, int margin,
                     int narm, int mean)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1],
		r = (margin == 0) ? m : n;

	SEXP res = PROTECT(allocVector(SUM_TYPEOF(class[0]), r)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));

	SEXP dimnames = (class[1] != 's')
		? GET_SLOT(obj, Matrix_DimNamesSym)
		: get_symmetrized_DimNames(obj, -1),
		marnames = VECTOR_ELT(dimnames, margin);
	if (marnames != R_NilValue) {
		PROTECT(marnames);
		setAttrib(res, R_NamesSymbol, marnames);
		UNPROTECT(1); /* marnames */
	}

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (class[1] == 't') {
			SEXP diag = GET_SLOT(obj, Matrix_diagSym);
			di = *CHAR(STRING_ELT(diag, 0));
		}
	}

	if (margin == 0 || class[1] == 's')
		dense_rowsum(x, class, m, n, ul, di, narm, mean, res);
	else
		dense_colsum(x, class, m, n, ul, di, narm, mean, res);

	UNPROTECT(2); /* x, res */
	return res;
}

/* (row|col)(Sums|Means)(<denseMatrix>) */
SEXP R_dense_marginsum(SEXP obj, SEXP margin,
                       SEXP narm, SEXP mean)
{
	static const char *valid[] = { VALID_DENSE, "" };
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

	return dense_marginsum(obj, valid[ivalid], margin_, narm_, mean_);
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

SEXP dense_sum(SEXP obj, const char *class, int narm)
{
	SEXP res;

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

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, unpacked = class[2] != 'p', symmetric = class[1] == 's';

#define SUM_LOOP \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's' || di == 'N') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					SUM_KERNEL(for (i = 0; i <= j; ++i)); \
					if (unpacked) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (unpacked) \
						px += j; \
					SUM_KERNEL(for (i = j; i < m; ++i)); \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					SUM_KERNEL(for (i = 0; i < j; ++i)); \
					++px; \
					if (unpacked) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (unpacked) \
						px += j; \
					++px; \
					SUM_KERNEL(for (i = j + 1; i < m; ++i)); \
				} \
			} \
		} \
	} while (0)

	if (class[0] == 'n') {
		int *px = LOGICAL(x);
		Matrix_int_fast64_t s = (di == 'N') ? 0LL : n;

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (*px != 0) \
					s += (symmetric && i != j) ? 2 : 1; \
				++px; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

		if (s <= INT_MAX) {
			res = allocVector(INTSXP, 1);
			INTEGER(res)[0] = (int) s;
		} else {
			res = allocVector(REALSXP, 1);
			REAL(res)[0] = (double) s;
		}
		return res;
	}

	if (!narm && (class[0] == 'l' || class[0] == 'i')) {
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (*px == NA_INTEGER) { \
					res = allocVector(INTSXP, 1); \
					INTEGER(res)[0] = NA_INTEGER; \
					return res; \
				} \
				++px; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

	}

	if (class[0] == 'z') {
		Rcomplex *px = COMPLEX(x);
		long double zr = (di == 'N') ? 0.0L : n, zi = 0.0L;

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
					zr += (symmetric && i != j) \
						? 2.0L * (*px).r : (*px).r; \
					zi += (symmetric && i != j) \
						? 2.0L * (*px).i : (*px).i; \
				} \
				++px; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

		res = allocVector(CPLXSXP, 1);
		COMPLEX(res)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
		COMPLEX(res)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
	} else if (class[0] == 'd') {
		double *px = REAL(x);
		long double zr = (di == 'N') ? 0.0L : n;

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && ISNAN(*px))) \
					zr += (symmetric && i != j) \
						? 2.0L * *px : *px; \
				++px; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

		res = allocVector(REALSXP, 1);
		REAL(res)[0] = LONGDOUBLE_AS_DOUBLE(zr);
	} else {
		int *px = (class[0] == 'i') ? INTEGER(x) : LOGICAL(x);
		Matrix_int_fast64_t s = (di == 'N') ? 0LL : n, t = 0LL;
		unsigned int count = 0;
		int over = 0;

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && *px == NA_INTEGER)) { \
					int d = (symmetric && i != j) ? 2 : 1; \
					if (count > UINT_MAX - d) \
						TRY_INCREMENT(ifover); \
					t += (d == 2) ? 2LL * *px : *px; \
					count += d; \
				} \
				++px; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

		TRY_INCREMENT(ifover);
	ifover:
		if (over) {
			long double zr = (di == 'N') ? 0.0L : n; /* FIXME: wasteful */
			px = (class[0] == 'i') ? INTEGER(x) : LOGICAL(x);

#define SUM_KERNEL(_FOR_) \
			do { \
				_FOR_ { \
					if (!(narm && *px == NA_INTEGER)) \
						zr += (symmetric && i != j) \
							? 2.0L * *px : *px; \
					++px; \
				} \
			} while (0)

			SUM_LOOP;

#undef SUM_KERNEL

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

#undef SUM_LOOP

	return res;
}

/* sum(<denseMatrix>) */
SEXP R_dense_sum(SEXP obj, SEXP narm)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int narm_;
	if (TYPEOF(narm) != LGLSXP || LENGTH(narm) < 1 ||
	    (narm_ = LOGICAL(narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	return dense_sum(obj, valid[ivalid], narm_);
}

SEXP dense_prod(SEXP obj, const char *class, int narm)
{
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

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, unpacked = class[2] != 'p', symmetric = class[1] == 's';
	long double zr = 1.0L, zi = 0.0L;

#define PROD_LOOP \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				PROD_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					PROD_KERNEL(for (i = 0; i <= j; ++i)); \
					if (unpacked) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (unpacked) \
						px += j; \
					PROD_KERNEL(for (i = j; i < m; ++i)); \
				} \
			} \
		} else if (di == 'N') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					if (j == 1) { zr *= 0.0L; zi *= 0.0L; } \
					PROD_KERNEL(for (i = 0; i <= j; ++i)); \
					if (unpacked) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (j == 1) { zr *= 0.0L; zi *= 0.0L; } \
					if (unpacked) \
						px += j; \
					PROD_KERNEL(for (i = j; i < m; ++i)); \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					if (j == 1) { zr *= 0.0L; zi *= 0.0L; } \
					PROD_KERNEL(for (i = 0; i < j; ++i)); \
					++px; \
					if (unpacked) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (j == 1) { zr *= 0.0L; zi *= 0.0L; } \
					if (unpacked) \
						px += j; \
					++px; \
					PROD_KERNEL(for (i = j + 1; i < m; ++i)); \
				} \
			} \
		} \
	} while (0)

	if (class[0] == 'n') {
		int *px = LOGICAL(x);
		if (class[1] == 't')
			REAL(res)[0] = (n > 1 || (n == 1 && *px == 0)) ? 0.0 : 1.0;
		else {

#define PROD_KERNEL(_FOR_) \
			do { \
				_FOR_ { \
					if (*px == 0) { \
						REAL(res)[0] = 0.0; \
						UNPROTECT(1); /* res */ \
						return res; \
					} \
					++px; \
				} \
			} while (0)

			PROD_LOOP;

#undef PROD_KERNEL

			REAL(res)[0] = 1.0;
		}
		UNPROTECT(1); /* res */
		return res;
	}

	if (class[0] == 'z') {
		Rcomplex *px = COMPLEX(x);
		long double zr0, zi0;

#define PROD_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
					zr0 = zr; zi0 = zi; \
					zr = zr0 * (*px).r - zi0 * (*px).i; \
					zi = zr0 * (*px).i + zi0 * (*px).r; \
					if (symmetric && i != j) { \
					zr0 = zr; zi0 = zi; \
					zr = zr0 * (*px).r - zi0 * (*px).i; \
					zi = zr0 * (*px).i + zi0 * (*px).r; \
					} \
				} \
				++px; \
			} \
		} while (0)

		PROD_LOOP;

#undef PROD_KERNEL

	} else if (class[0] == 'd') {
		double *px = REAL(x);

#define PROD_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && ISNAN(*px))) \
					zr *= (symmetric && i != j) \
						? (long double) *px * *px : *px; \
				++px; \
			} \
		} while (0)

		PROD_LOOP;

#undef PROD_KERNEL

	} else {
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);

#define PROD_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (*px != NA_INTEGER) \
					zr *= (symmetric && i != j) \
						? (long double) *px * *px : *px; \
				else if (!narm) \
					zr *= NA_REAL; \
				++px; \
			} \
		} while (0)

		PROD_LOOP;

#undef PROD_KERNEL

	}

#undef PROD_LOOP

	if (class[0] == 'z') {
		COMPLEX(res)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
		COMPLEX(res)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
	} else
		   REAL(res)[0]   = LONGDOUBLE_AS_DOUBLE(zr);
	UNPROTECT(1); /* res */
	return res;
}

/* prod(<denseMatrix>) */
SEXP R_dense_prod(SEXP obj, SEXP narm)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int narm_;
	if (TYPEOF(narm) != LGLSXP || LENGTH(narm) < 1 ||
	    (narm_ = LOGICAL(narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	return dense_prod(obj, valid[ivalid], narm_);
}

#undef TRY_INCREMENT
#undef LONGDOUBLE_AS_DOUBLE

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
