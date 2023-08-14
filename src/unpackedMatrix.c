#include "unpackedMatrix.h"

#define PACK(_PREFIX_, _CTYPE_, _ONE_) \
void _PREFIX_ ## dense_pack(_CTYPE_ *dest, const _CTYPE_ *src, int n, \
                            char uplo, char diag) \
{ \
	int i, j; \
	R_xlen_t dpos = 0, spos = 0; \
	if (uplo == 'U') { \
		for (j = 0; j < n; spos += n-(++j)) \
			for (i = 0; i <= j; ++i) \
				dest[dpos++] = src[spos++]; \
		if (diag != 'N') { \
			dpos = 0; \
			for (j = 0; j < n; dpos += (++j)+1) \
				dest[dpos] = _ONE_; \
		} \
	} else { \
		for (j = 0; j < n; spos += (++j)) \
			for (i = j; i < n; ++i) \
				dest[dpos++] = src[spos++]; \
		if (diag != 'N') { \
			dpos = 0; \
			for (j = 0; j < n; dpos += n-(j++)) \
				dest[dpos] = _ONE_; \
		} \
	} \
	return; \
}

/**
 * @brief Pack a square `unpackedMatrix`.
 *
 * Copies the upper or lower triangular part of `src` to `dest`,
 * where it is stored contiguously ("packed"). Optionally resets
 * the diagonal elements to 1.
 *
 * @param dest,src Pointers to the first elements of length-`(n*(n+1))/2`
 *     and length-`n*n` (resp.) arrays, usually the "data" of the `x`
 *     slot of an `n`-by-`n` `packedMatrix` and `unpackedMatrix` (resp.).
 * @param n Size of matrix being packed.
 * @param uplo,diag `char` specifying whether the "nontrivial part"
 *     is upper (`'U'`) or lower (`'L'`) and whether to "force" a
 *     unit diagonal (`'U'`) or not (`'N'`).
 */
/* ddense_pack() */
PACK(d, double, 1.0)
/* idense_pack() */
PACK(i, int, 1)
/* zdense_pack() */
PACK(z, Rcomplex, Matrix_zone)

#undef PACK

#define MAKE_TRIANGULAR(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## dense_unpacked_make_triangular(_CTYPE_ *x, int m, int n, \
                                                char uplo, char diag) \
{ \
	int i, j, r = (m < n) ? m : n; \
	R_xlen_t pos = 0; \
	if (uplo == 'U') { \
		for (j = 0; j < r; pos += (++j)+1) \
			for (i = j+1; i < m; ++i) \
				x[++pos] = _ZERO_; \
	} else { \
		for (j = 0; j < r; pos += m-(j++)) \
			for (i = 0; i < j; ++i) \
				x[pos++] = _ZERO_; \
		for (j = r; j < n; ++j) \
			for (i = 0; i < m; ++i) \
				x[pos++] = _ZERO_; \
	} \
	if (diag != 'N') { \
		pos = 0; \
		R_xlen_t m1a = (R_xlen_t) m + 1; \
		for (j = 0; j < r; ++j, pos += m1a) \
			x[pos] = _ONE_; \
	} \
	return; \
}

/**
 * @brief Make triangular an `unpackedMatrix`.
 *
 * "Triangularizes" the elements of an `m`-by-`n` `unpackedMatrix`,
 * which need not be square (though all `triangularMatrix` _are_).
 *
 * @param x A pointer to the first element of a length-`m*n` array,
 *     usually the "data" of the `x` slot of an `m`-by-`n` `unpackedMatrix`.
 * @param m,n Dimensions of matrix being triangularized.
 * @param uplo,diag `char` constants specifying whether the matrix
 *     should be upper (`'U'`) or lower (`'L'`) triangular and
 *     whether it should have a unit diagonal (`'U'`) or not (`'N'`).
 */
/* ddense_unpacked_make_triangular() */
MAKE_TRIANGULAR(d, double, 0.0, 1.0)
/* idense_unpacked_make_triangular() */
MAKE_TRIANGULAR(i, int, 0, 1)
/* zdense_unpacked_make_triangular() */
MAKE_TRIANGULAR(z, Rcomplex, Matrix_zzero, Matrix_zone)

#undef MAKE_TRIANGULAR

#define MAKE_SYMMETRIC_SET_EQ(_X_, _DEST_, _SRC_) \
	_X_[_DEST_] = _X_[_SRC_]

#define MAKE_SYMMETRIC_SET_CJ(_X_, _DEST_, _SRC_) \
	do { \
		_X_[_DEST_].r =  _X_[_SRC_].r; \
		_X_[_DEST_].i = -_X_[_SRC_].i; \
	} while (0)

#define MAKE_SYMMETRIC(_PREFIX_, _CTYPE_, _SET_) \
void _PREFIX_ ## dense_unpacked_make_symmetric(_CTYPE_ *x, int n, \
                                               char uplo) \
{ \
	int i, j, n1s = n - 1; \
	R_xlen_t upos = n, lpos = 1; \
	if (uplo == 'U') { \
		for (j = 0; j < n; upos = (lpos += (++j)+1) + n1s) \
			for (i = j+1; i < n; ++i, upos += n, ++lpos) \
				_SET_(x, lpos, upos); \
	} else { \
		for (j = 0; j < n; upos = (lpos += (++j)+1) + n1s) \
			for (i = j+1; i < n; ++i, upos += n, ++lpos) \
				_SET_(x, upos, lpos); \
	} \
	return; \
}

/**
 * @brief Make symmetric a square `unpackedMatrix`.
 *
 * "Symmetrizes" the elements of an `n`-by-`n` `unpackedMatrix`.
 *
 * @param x A pointer to the first element of a length-`n*n` array,
 *     usually the "data" of the `x` slot of an `n`-by-`n` `unpackedMatrix`.
 * @param n Size of matrix being symmetrized.
 * @param uplo A `char` specifying whether to copy the upper triangle
 *     to the lower triangle (`'U'`) or to do the reverse (`'L'`).
 */
/* ddense_unpacked_make_symmetric() */
MAKE_SYMMETRIC(d, double, MAKE_SYMMETRIC_SET_EQ)
/* idense_unpacked_make_symmetric() */
MAKE_SYMMETRIC(i, int, MAKE_SYMMETRIC_SET_EQ)
/* zdense_unpacked_make_symmetric() */
MAKE_SYMMETRIC(z, Rcomplex, MAKE_SYMMETRIC_SET_CJ)

#undef MAKE_SYMMETRIC
#undef MAKE_SYMMETRIC_SET_CJ
#undef MAKE_SYMMETRIC_SET_EQ

#define MAKE_BANDED(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## dense_unpacked_make_banded(_CTYPE_ *x, int m, int n, \
                                            int a, int b, char diag) \
{ \
	if (m == 0 || n == 0) \
		return; \
	if (a > b || a >= n || b <= -m) { \
		Matrix_memset(x, 0, (R_xlen_t) m * n, sizeof(_CTYPE_)); \
		return; \
	} \
	if (a <= -m) a = 1-m; \
	if (b >=  n) b = n-1; \
 \
	int i, j, i0, i1, \
		j0 = (a < 0) ? 0 : a, \
		j1 = (b < n-m) ? m+b : n; \
 \
	if (j0 > 0) { \
		R_xlen_t dx; \
		Matrix_memset(x, 0, dx = (R_xlen_t) m * j0, sizeof(_CTYPE_)); \
		x += dx; \
	} \
	for (j = j0; j < j1; ++j, x += m) { \
		i0 = j - b; \
		i1 = j - a + 1; \
		for (i = 0; i < i0; ++i) \
			*(x + i) = _ZERO_; \
		for (i = i1; i < m; ++i) \
			*(x + i) = _ZERO_; \
	} \
	if (j1 < n) \
		Matrix_memset(x, 0, (R_xlen_t) m * (n - j1), sizeof(_CTYPE_)); \
	if (diag != 'N' && a <= 0 && b >= 0) { \
		x -= m * (R_xlen_t) j; \
		R_xlen_t m1a = (R_xlen_t) m + 1; \
		for (j = 0; j < n; ++j, x += m1a) \
			*x = _ONE_; \
	} \
	return; \
}
MAKE_BANDED(d, double, 0.0, 1.0)
MAKE_BANDED(i, int, 0, 1)
MAKE_BANDED(z, Rcomplex, Matrix_zzero, Matrix_zone)

#undef MAKE_BANDED

#define COPY_DIAGONAL(_PREFIX_, _CTYPE_, _ONE_) \
void _PREFIX_ ## dense_unpacked_copy_diagonal(_CTYPE_ *dest, \
                                              const _CTYPE_ *src, \
                                              int n, R_xlen_t len, \
                                              char uplo, char diag) \
{ \
	int j; \
	R_xlen_t n1a = (R_xlen_t) n + 1; \
	if (diag != 'N') { \
		for (j = 0; j < n; ++j, dest += n1a) \
			*dest = _ONE_; \
	} else if (len == n) { \
		/* copying from diagonalMatrix */ \
		for (j = 0; j < n; ++j, dest += n1a, ++src) \
			*dest = *src; \
	} else if (len == (n * n1a) / 2) { \
		/* copying from packedMatrix */ \
		if (uplo == 'U') { \
			for (j = 0; j < n; dest += n1a, src += (++j)+1) \
				*dest = *src; \
		} else { \
			for (j = 0; j < n; dest += n1a, src += n-(j++)) \
				*dest = *src; \
		} \
	} else if (len == (R_xlen_t) n * n) { \
		/* copying from square unpackedMatrix */ \
		for (j = 0; j < n; ++j, dest += n1a, src += n1a) \
			*dest = *src; \
	} else { \
		error(_("incompatible '%s' and '%s' in %s()"), "n", "len", __func__); \
	} \
	return; \
}

/**
 * Copy a length-`n` diagonal to a length-`n*n` array.
 *
 * @param dest A pointer to the first element of a length-`n*n` array,
 *     usually the "data" of the `x` slot of an `n`-by-`n` `unpackedMatrix`.
 * @param src A pointer to the first element of a length-`n`,
 *     length-`(n*(n+1))/2`, or length-`n*n` array, usually the "data"
 *     of the `x` slot of an `n`-by-`n` `diagonalMatrix`, `packedMatrix`,
 *     or `unpackedMatrix`, respectively.
 * @param n Size of matrix being copied from and to.
 * @param len Length of `src` array.
 * @param uplo,diag `char` constants specifying whether `src` stores an
 *     upper (`'U'`) or lower (`'L'`) triangle when `len == (n*(n+1))/2` and
 *     whether the matrix should have a unit diagonal (`'U'`) or not (`'N'`).
 */
/* ddense_unpacked_copy_diagonal() */
COPY_DIAGONAL(d, double, 1.0)
/* idense_unpacked_copy_diagonal() */
COPY_DIAGONAL(i, int, 1)
/* zdense_unpacked_copy_diagonal() */
COPY_DIAGONAL(z, Rcomplex, Matrix_zone)

#undef COPY_DIAGONAL

#define IS_TRIANGULAR(_PREFIX_, _CTYPE_, \
                      _L_IS_NOT_ZERO_, _U_IS_NOT_ZERO_) \
static Rboolean _PREFIX_ ## dense_unpacked_is_triangular(const _CTYPE_ *x, \
                                                         int n, \
                                                         char uplo) \
{ \
	int i, j; \
	if (uplo == 'U') { \
		for (j = 0; j < n; x += (++j)+1) \
			for (i = j+1; i < n; ++i) \
				if (_L_IS_NOT_ZERO_) \
					return FALSE; \
	} else { \
		for (j = 0; j < n; x += n-(j++)) \
			for (i = 0; i < j; ++i) \
				if (_U_IS_NOT_ZERO_) \
					return FALSE; \
	} \
	return TRUE; \
}

/* ddense_unpacked_is_triangular() */
IS_TRIANGULAR(d, double,
              ISNAN(*(++x)) || *x != 0.0,
              ISNAN(*x) || *(x++) != 0.0)
/* idense_unpacked_is_triangular() */
IS_TRIANGULAR(i, int,
              *(++x) != 0,
              *(x++) != 0)
/* zdense_unpacked_is_triangular() */
IS_TRIANGULAR(z, Rcomplex,
              ISNAN((*(++x)).r) || (*x).r != 0.0 ||
              ISNAN((*x).i) || (*x).i != 0.0,
              ISNAN((*x).r) || (*x).r != 0.0 ||
              ISNAN((*x).i) || (*(x++)).i != 0.0)

#undef IS_TRIANGULAR

#define IS_SYMMETRIC(_PREFIX_, _CTYPE_, \
                     _U_IS_NA_, _L_IS_NOT_NA_, _L_IS_NOT_EQUAL_) \
static Rboolean _PREFIX_ ## dense_unpacked_is_symmetric(const _CTYPE_ *x, \
                                                        int n) \
{ \
	int i, j; \
	R_xlen_t upos = 0, lpos = 0; \
	for (j = 0; j < n; upos = (lpos += (++j)+1)) { \
		for (i = j+1; i < n; ++i) { \
			upos += n; ++lpos; \
			if (_U_IS_NA_) { \
				if (_L_IS_NOT_NA_) \
					return FALSE; \
			} else { \
				if (_L_IS_NOT_EQUAL_) \
					return FALSE; \
			} \
		} \
	} \
	return TRUE; \
}

/* ddense_unpacked_is_symmetric() */
IS_SYMMETRIC(d, double,
             ISNAN(x[upos]),
             !ISNAN(x[lpos]),
             ISNAN(x[lpos]) || x[lpos] != x[upos])
/* ldense_unpacked_is_symmetric() */
IS_SYMMETRIC(l, int,
             x[upos] == NA_LOGICAL,
             x[lpos] != NA_LOGICAL,
             (x[upos] == 0) ? (x[lpos] != 0) : (x[lpos] == 0))
/* ndense_unpacked_is_symmetric() ... macro still works, if we cheat */
IS_SYMMETRIC(n, int,
             x[upos] == 0,
             x[lpos] != 0,
             x[lpos] == 0)
/* idense_unpacked_is_symmetric() */
IS_SYMMETRIC(i, int,
             x[upos] == NA_INTEGER,
             x[lpos] != NA_INTEGER,
             x[lpos] != x[upos])
/* zdense_unpacked_is_symmetric() */
IS_SYMMETRIC(z, Rcomplex,
             ISNAN(x[upos].r) || ISNAN(x[upos].i) ,
             !(ISNAN(x[lpos].r) || ISNAN(x[lpos].i)),
             ISNAN(x[lpos].r) || ISNAN(x[lpos].i) ||
             x[upos].r != x[lpos].r || x[upos].i != -x[lpos].i)

#undef IS_SYMMETRIC

#define IS_DIAGONAL(_PREFIX_, _CTYPE_, _OD_IS_NOT_ZERO_) \
static Rboolean _PREFIX_ ## dense_unpacked_is_diagonal(const _CTYPE_ *x, \
                                                       int n) \
{ \
	int i, j; \
	for (j = 0; j < n; ++j) { \
		for (i = 0; i < j; ++i) \
			if (_OD_IS_NOT_ZERO_) \
				return FALSE; \
		++x; /* skip over diagonal */ \
		for (i = j+1; i < n; ++i) \
			if (_OD_IS_NOT_ZERO_) \
				return FALSE; \
	} \
	return TRUE; \
}

/* ddense_unpacked_is_diagonal() */
IS_DIAGONAL(d, double, ISNAN(*x) || *(x++) != 0.0)
/* idense_unpacked_is_diagonal() */
IS_DIAGONAL(i, int, *(x++) != 0)
/* zdense_unpacked_is_diagonal() */
IS_DIAGONAL(z, Rcomplex,
            ISNAN((*x).r) || (*(x)).r != 0.0 ||
            ISNAN((*x).i) || (*(x++)).i != 0.0)

#undef IS_DIAGONAL

SEXP unpacked_force(SEXP x, int n, char uplo, char diag)
{
	SEXPTYPE tx = TYPEOF(x);
	if (tx < LGLSXP || tx > CPLXSXP)
		ERROR_INVALID_TYPE(x, __func__);
	R_xlen_t nx = XLENGTH(x);
	SEXP y = PROTECT(allocVector(tx, nx));

	if (diag == '\0') {

#define FORCE_SYMMETRIC(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *py = _PTR_(y); \
		Matrix_memcpy(py, px, nx, sizeof(_CTYPE_)); \
		_PREFIX_ ## dense_unpacked_make_symmetric(py, n, uplo); \
	} while (0)

	switch (tx) {
	case LGLSXP:
		FORCE_SYMMETRIC(i, int, LOGICAL);
		break;
	case INTSXP:
		FORCE_SYMMETRIC(i, int, INTEGER);
		break;
	case REALSXP:
		FORCE_SYMMETRIC(d, double, REAL);
		break;
	case CPLXSXP:
		FORCE_SYMMETRIC(z, Rcomplex, COMPLEX);
		break;
	default:
		break;
	}

#undef FORCE_SYMMETRIC

	} else {

#define FORCE_TRIANGULAR(_PREFIX_, _CTYPE_, _PTR_, _ONE_) \
		do { \
			_CTYPE_ *px = _PTR_(x), *py = _PTR_(y); \
			Matrix_memcpy(py, px, nx, sizeof(_CTYPE_)); \
			_PREFIX_ ## dense_unpacked_make_triangular( \
				py, n, n, uplo, diag); \
			if (diag != 'N') { \
				R_xlen_t n1a = (R_xlen_t) n + 1; \
				for (int j = 0; j < n; ++j, py += n1a) \
					*py = _ONE_; \
			} \
		} while (0)

		switch (tx) {
		case LGLSXP:
			FORCE_TRIANGULAR(i, int, LOGICAL, 1);
			break;
		case INTSXP:
			FORCE_TRIANGULAR(i, int, INTEGER, 1);
			break;
		case REALSXP:
			FORCE_TRIANGULAR(d, double, REAL, 1.0);
			break;
		case CPLXSXP:
			FORCE_TRIANGULAR(z, Rcomplex, COMPLEX, Matrix_zone);
			break;
		default:
			break;
		}

#undef FORCE_TRIANGULAR

	}

	UNPROTECT(1);
	return y;
}

/* forceSymmetric(x, uplo), returning .syMatrix */
SEXP unpackedMatrix_force_symmetric(SEXP from, SEXP uplo_to)
{
	static const char *valid[] = {
	/* 0 */ "dgeMatrix", "lgeMatrix", "ngeMatrix",
	/* 3 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *clf = valid[ivalid];

	char ulf = 'U', ult = 'U';
	if (clf[1] != 'g') {
		/* .(tr|sy)Matrix */
		SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ulf = ult = *CHAR(STRING_ELT(uplo_from, 0));
		UNPROTECT(1); /* uplo_from */
	}

	if (!isNull(uplo_to) &&
	    (TYPEOF(uplo_to) != STRSXP || LENGTH(uplo_to) < 1 ||
	    (uplo_to = STRING_ELT(uplo_to, 0)) == NA_STRING ||
	    ((ult = *CHAR(uplo_to)) != 'U' && ult != 'L')))
		error(_("invalid '%s' to %s()"), "uplo", __func__);

	if (clf[1] == 's') {
		/* .syMatrix */
		if (ulf == ult)
			return from;
		SEXP to = PROTECT(unpackedMatrix_transpose(from));
		if (clf[0] == 'z') {
			/* Need _conjugate_ transpose */
			SEXP x_to = PROTECT(GET_SLOT(to, Matrix_xSym));
			conjugate(x_to);
			UNPROTECT(1); /* x_to */
		}
		UNPROTECT(1); /* to */
		return to;
	}

	/* Now handling just .(ge|tr)Matrix ... */

	char clt[] = ".syMatrix";
	clt[0] = clf[0];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
		x_from = PROTECT(GET_SLOT(from, Matrix_xSym));

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
		PROTECT(uplo_to = mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo_to);
		UNPROTECT(1); /* uplo_to */
	}

	if (clf[1] == 'g' || ulf == ult) {
		/* .geMatrix or .trMatrix with correct uplo */
		SET_SLOT(to, Matrix_xSym, x_from);
	} else {
		/* .trMatrix with incorrect uplo */
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */

		SEXPTYPE tx = TYPEOF(x_from);
		R_xlen_t nx = XLENGTH(x_from);
		SEXP x_to = PROTECT(allocVector(tx, nx));

#define COPY_DIAGONAL(_PREFIX_, _CTYPE_, _PTR_) \
		do { \
			Matrix_memset(_PTR_(x_to), 0, nx, sizeof(_CTYPE_)); \
			_PREFIX_ ## dense_unpacked_copy_diagonal( \
				_PTR_(x_to), _PTR_(x_from), n, nx, 'U' /* unused */, di); \
		} while (0)

		switch (tx) {
		case LGLSXP: /* [nl]..Matrix */
			COPY_DIAGONAL(i, int, LOGICAL);
			break;
		case INTSXP: /* i..Matrix */
			COPY_DIAGONAL(i, int, INTEGER);
			break;
		case REALSXP: /* d..Matrix */
			COPY_DIAGONAL(d, double, REAL);
			break;
		case CPLXSXP: /* z..Matrix */
			COPY_DIAGONAL(z, Rcomplex, COMPLEX);
			break;
		default:
			ERROR_INVALID_TYPE(x_from, __func__);
			break;
		}

#undef COPY_DIAGONAL

		SET_SLOT(to, Matrix_xSym, x_to);
		UNPROTECT(1); /* x_to */
	}

	UNPROTECT(2); /* x_from, to */
	return to;
}

#define UPM_IS_TR(_RES_, _X_, _N_, _UPLO_) \
do { \
	switch (TYPEOF(_X_)) { \
	case LGLSXP: \
		_RES_ = idense_unpacked_is_triangular(LOGICAL(_X_), _N_, _UPLO_); \
		break; \
	case INTSXP: \
		_RES_ = idense_unpacked_is_triangular(INTEGER(_X_), _N_, _UPLO_); \
		break; \
	case REALSXP: \
		_RES_ = ddense_unpacked_is_triangular(REAL(_X_), _N_, _UPLO_); \
		break; \
	case CPLXSXP: \
		_RES_ = zdense_unpacked_is_triangular(COMPLEX(_X_), _N_, _UPLO_); \
		break; \
	default: \
		ERROR_INVALID_TYPE(_X_, __func__); \
		_RES_ = FALSE; \
		break; \
	} \
} while (0)

#define UPM_IS_SY(_RES_, _X_, _N_, _LDENSE_) \
do { \
	switch (TYPEOF(_X_)) { \
	case LGLSXP: \
		_RES_ = (_LDENSE_ \
		         ? ldense_unpacked_is_symmetric(LOGICAL(_X_), _N_) \
		         : ndense_unpacked_is_symmetric(LOGICAL(_X_), _N_)); \
		break; \
	case INTSXP: \
		_RES_ = idense_unpacked_is_symmetric(INTEGER(_X_), _N_); \
		break; \
	case REALSXP: \
		_RES_ = ddense_unpacked_is_symmetric(REAL(_X_), _N_); \
		break; \
	case CPLXSXP: \
		_RES_ = zdense_unpacked_is_symmetric(COMPLEX(_X_), _N_); \
		break; \
	default: \
		ERROR_INVALID_TYPE(_X_, __func__); \
		_RES_ = FALSE; \
		break; \
	} \
} while (0)

#define UPM_IS_DI(_RES_, _X_, _N_) \
do { \
	switch (TYPEOF(_X_)) { \
	case LGLSXP: \
		_RES_ = idense_unpacked_is_diagonal(LOGICAL(_X_), _N_); \
		break; \
	case INTSXP: \
		_RES_ = idense_unpacked_is_diagonal(INTEGER(_X_), _N_); \
		break; \
	case REALSXP: \
		_RES_ = ddense_unpacked_is_diagonal(REAL(_X_), _N_); \
		break; \
	case CPLXSXP: \
		_RES_ = zdense_unpacked_is_diagonal(COMPLEX(_X_), _N_); \
		break; \
	default: \
		ERROR_INVALID_TYPE(_X_, __func__); \
		_RES_ = FALSE; \
		break; \
	} \
} while (0)

#define RETURN_GE_IS_TR(_X_, _N_, _UPPER_, _NPROT_) \
do { \
	Rboolean res = FALSE; \
	if (_UPPER_ == NA_LOGICAL) { \
		UPM_IS_TR(res, _X_, _N_, 'U'); \
		if (res) { \
			UNPROTECT(_NPROT_); \
			RETURN_TRUE_OF_KIND("U"); \
		} \
		UPM_IS_TR(res, _X_, _N_, 'L'); \
		if (res) { \
			UNPROTECT(_NPROT_); \
			RETURN_TRUE_OF_KIND("L"); \
		} \
	} else { \
		UPM_IS_TR(res, _X_, _N_, (_UPPER_ != 0) ? 'U' : 'L'); \
	} \
	UNPROTECT(_NPROT_); \
	return ScalarLogical(res); \
} while (0)

/* isTriangular(x, upper) */
SEXP unpackedMatrix_is_triangular(SEXP obj, SEXP upper)
{
	static const char *valid[] = {
	/* 0 */ "dgeMatrix", "lgeMatrix", "ngeMatrix",
	/* 3 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	int need_upper = asLogical(upper);

	if (ivalid < 3) {
		/* .geMatrix: need to do a complete triangularity check */
		SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
		int *pdim = INTEGER(dim), n = pdim[0], s = pdim[1] == n;
		UNPROTECT(1); /* dim */
		if (!s)
			return ScalarLogical(0);
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		RETURN_GE_IS_TR(x, n, need_upper, /* unprotect this many: */ 1);
	} else {
		/* .(tr|sy)Matrix */
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

#define IF_DIAGONAL \
		Rboolean res = FALSE; \
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)), \
			dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)); \
		int n = INTEGER(dim)[0]; \
		UPM_IS_TR(res, x, n, (ul == 'U') ? 'L' : 'U'); \
		UNPROTECT(2); /* dim, x */ \
		if (res)

		if (ivalid < 6) {
			/* .trMatrix: fast if 'upper', 'uplo' agree; else need diagonal */
			if (need_upper == NA_LOGICAL)
				RETURN_TRUE_OF_KIND((ul == 'U') ? "U" : "L");
			else if ((need_upper) ? ul == 'U' : ul != 'U')
				return ScalarLogical(1);
			else {
				IF_DIAGONAL {
					return ScalarLogical(1);
				}
			}
		} else {
			/* .syMatrix: triangular iff diagonal (upper _and_ lower tri.) */
			IF_DIAGONAL {
				if (need_upper == NA_LOGICAL)
					RETURN_TRUE_OF_KIND("U");
				else
					return ScalarLogical(1);
			}
		}

#undef IF_DIAGONAL

		return ScalarLogical(0);
	}
}

/* isTriangular(x, upper) */
SEXP matrix_is_triangular(SEXP obj, SEXP upper)
{
	SEXP dim = PROTECT(getAttrib(obj, R_DimSymbol));
	int *pdim = INTEGER(dim), n = pdim[0], s = pdim[1] == n;
	UNPROTECT(1); /* dim */
	if (!s)
		return ScalarLogical(0);
	int need_upper = asLogical(upper);
	RETURN_GE_IS_TR(obj, n, need_upper, /* unprotect this many: */ 0);
}

#undef RETURN_GE_IS_TR

/* isSymmetric(x, tol = 0, checkDN) */
/* FIXME: not checking for real diagonal in complex case */
SEXP unpackedMatrix_is_symmetric(SEXP obj, SEXP checkDN)
{
	static const char *valid[] = {
	/* 0 */ "dgeMatrix", "lgeMatrix", "ngeMatrix",
	/* 3 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0) {
		ERROR_INVALID_CLASS(obj, __func__);
		return R_NilValue;
	} else if (ivalid < 6) {
		/* .(ge|tr)Matrix */
		SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
		int *pdim = INTEGER(dim), n = pdim[0], s = pdim[1] == n;
		UNPROTECT(1); /* dim */
		if (s && asLogical(checkDN) != 0) {
			SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
			s = DimNames_is_symmetric(dimnames);
			UNPROTECT(1); /* dimnames */
		}
		if (!s)
			return ScalarLogical(0);
		Rboolean res = FALSE;
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		if (ivalid < 3) {
			/* .geMatrix: need to do a complete symmetry check */
			UPM_IS_SY(res, x, n, ivalid == 1);
		} else {
			/* .trMatrix: symmetric iff diagonal (upper _and_ lower tri.) */
			SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
			char ul = (*CHAR(STRING_ELT(uplo, 0)) == 'U') ? 'L' : 'U';
			UNPROTECT(1); /* uplo */
			UPM_IS_TR(res, x, n, ul);
		}
		UNPROTECT(1); /* x */
		return ScalarLogical(res);
	} else {
		/* .syMatrix: symmetric by definition */
		return ScalarLogical(1);
	}
}

/* isSymmetric(x, tol = 0) */
SEXP matrix_is_symmetric(SEXP obj, SEXP checkDN)
{
	SEXP dim = PROTECT(getAttrib(obj, R_DimSymbol));
	int *pdim = INTEGER(dim), n = pdim[0], s = pdim[1] == n;
	UNPROTECT(1); /* dim */
	if (s && asLogical(checkDN) != 0) {
		SEXP dimnames = PROTECT(getAttrib(obj, R_DimNamesSymbol));
		s = isNull(dimnames) || DimNames_is_symmetric(dimnames);
		UNPROTECT(1); /* dimnames */
	}
	if (!s)
		return ScalarLogical(0);
	Rboolean res = FALSE;
	UPM_IS_SY(res, obj, n, 1);
	return ScalarLogical(res);
}

/* isDiagonal(x) */
SEXP unpackedMatrix_is_diagonal(SEXP obj)
{
	static const char *valid[] = {
	/* 0 */ "dgeMatrix", "lgeMatrix", "ngeMatrix",
	/* 3 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0], s = pdim[1] == n;
	UNPROTECT(1); /* dim */
	if (!s)
		return ScalarLogical(0);
	Rboolean res = FALSE;
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	if (ivalid < 3) {
		/* .geMatrix: need to do a complete diagonality check */
		UPM_IS_DI(res, x, n);
	} else {
		/* .(tr|sy)Matrix: diagonal iff stored triangle is zero off diagonal */
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = (*CHAR(STRING_ELT(uplo, 0)) == 'U') ? 'L' : 'U';
		UNPROTECT(1); /* uplo */
		UPM_IS_TR(res, x, n, ul);
	}
	UNPROTECT(1); /* x */
	return ScalarLogical(res);
}

/* isDiagonal(x) */
SEXP matrix_is_diagonal(SEXP obj)
{
	SEXP dim = PROTECT(getAttrib(obj, R_DimSymbol));
	int *pdim = INTEGER(dim), n = pdim[0], s = pdim[1] == n;
	UNPROTECT(1); /* dim */
	if (!s)
		return ScalarLogical(0);
	Rboolean res = FALSE;
	UPM_IS_DI(res, obj, n);
	return ScalarLogical(res);
}

#undef UPM_IS_DI
#undef UPM_IS_TR
#undef UPM_IS_SY

/* t(x) */
/* MJ: Technically no need to do full transpose of .(tr|sy)Matrix ...  */
/*     but then identical(.@x, t(t(.))@x) can be FALSE ...             */
SEXP unpackedMatrix_transpose(SEXP from)
{
	static const char *valid[] = {
	/* 0 */ "dgeMatrix", "lgeMatrix", "ngeMatrix",
	/* 3 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 6 */ "corMatrix", "dpoMatrix",
	/* 8 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid];

	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

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
	if (ivalid < 6)
		set_reversed_DimNames(to, dimnames);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (ivalid >= 3) {
		SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ulf = *CHAR(STRING_ELT(uplo_from, 0));
		UNPROTECT(1); /* uplo_from */
		if (ulf == 'U') {
			SEXP uplo_to = PROTECT(mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo_to);
			UNPROTECT(1); /* uplo_to */
		}
		if (ivalid < 6) {
			/* .trMatrix */
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			char di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		} else {
			/* .syMatrix */
			SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
			if (LENGTH(factors) > 0)
				SET_SLOT(to, Matrix_factorSym, factors);
			UNPROTECT(1); /* factors */

			if (ivalid == 6) {
				/* corMatrix */
				SEXP sd = PROTECT(GET_SLOT(from, Matrix_sdSym));
				if (LENGTH(sd) > 0)
					SET_SLOT(to, Matrix_sdSym, sd);
				UNPROTECT(1); /* sd */
			}
		}
	}

	SEXPTYPE tx;
	R_xlen_t nx;
	SEXP x_from = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x_to = PROTECT(allocVector(tx = TYPEOF(x_from), nx = XLENGTH(x_from)));

#define UPM_T(_CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x_from), *px1 = _PTR_(x_to); \
		int i, j; \
		R_xlen_t nx1s = nx - 1; \
		for (j = 0; j < m; ++j, px0 -= nx1s) \
			for (i = 0; i < n; ++i, px0 += m) \
				*(px1++) = *px0; \
	} while (0)

	switch (tx) {
	case LGLSXP: /* [nl]..Matrix */
		UPM_T(int, LOGICAL);
		break;
	case INTSXP: /* i..Matrix */
		UPM_T(int, INTEGER);
		break;
	case REALSXP: /* d..Matrix */
		UPM_T(double, REAL);
		break;
	case CPLXSXP: /* z..Matrix */
		UPM_T(Rcomplex, COMPLEX);
		break;
	default:
		ERROR_INVALID_TYPE(x_from, __func__);
		break;
	}

#undef UPM_T

	SET_SLOT(to, Matrix_xSym, x_to);

	UNPROTECT(3); /* x_to, x_from, to */
	return to;
}

/* diag(x, names) */
SEXP unpackedMatrix_diag_get(SEXP obj, SEXP nms)
{
	int do_nms = asLogical(nms);
	if (do_nms == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "names", "TRUE", "FALSE");

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim),
		m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	UNPROTECT(1); /* dim */

	char ul = '\0', di = '\0';
	if (HAS_SLOT(obj, Matrix_uploSym)) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		if (HAS_SLOT(obj, Matrix_diagSym)) {
			SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	}

	SEXPTYPE tx;
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	res = PROTECT(allocVector(tx = TYPEOF(x), r));

#define UPM_D_G(_CTYPE_, _PTR_, _ONE_) \
	do { \
		_CTYPE_ *pres = _PTR_(res); \
		int j; \
		if (di == 'U') { \
			for (j = 0; j < r; ++j) \
				*(pres++) = _ONE_; \
		} else { \
			_CTYPE_ *px = _PTR_(x); \
			R_xlen_t m1a = (R_xlen_t) m + 1; \
			for (j = 0; j < r; ++j, px += m1a) \
				*(pres++) = *px; \
		} \
	} while (0)

	switch (tx) {
	case LGLSXP: /* [nl]..Matrix */
		UPM_D_G(int, LOGICAL, 1);
		break;
	case INTSXP: /* i..Matrix */
		UPM_D_G(int, INTEGER, 1);
		break;
	case REALSXP: /* d..Matrix */
		UPM_D_G(double, REAL, 1.0);
		break;
	case CPLXSXP: /* z..Matrix */
		UPM_D_G(Rcomplex, COMPLEX, Matrix_zone);
		break;
	default:
		ERROR_INVALID_TYPE(x, __func__);
		break;
	}

#undef UPM_D_G

	if (do_nms) {
		/* NB: The logic here must be adjusted once the validity method
		   for 'symmetricMatrix' enforces symmetric 'Dimnames' */
		SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
			rn = VECTOR_ELT(dn, 0),
			cn = VECTOR_ELT(dn, 1);
		if (isNull(cn)) {
			if (ul != '\0' && di == '\0' && !isNull(rn))
				setAttrib(res, R_NamesSymbol, rn);
		} else {
			if (ul != '\0' && di == '\0')
				setAttrib(res, R_NamesSymbol, cn);
			else if (!isNull(rn) &&
			         (rn == cn || equal_string_vectors(rn, cn, r)))
				setAttrib(res, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

	UNPROTECT(2); /* res, x */
	return res;
}

/* diag(x) <- value */
SEXP unpackedMatrix_diag_set(SEXP obj, SEXP val)
{
	static const char *valid[] = {
	/* 0 */ "dgeMatrix", "lgeMatrix", "ngeMatrix",
	/* 3 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim),
		m = pdim[0], n = pdim[1], r = (m < n) ? m : n;

	SEXPTYPE tv = TYPEOF(val);
	if (tv < LGLSXP || tv > REALSXP)
		/* Upper bound can become CPLXSXP once we have proper zMatrix */
		error(_("replacement diagonal has incompatible type \"%s\""),
		      type2char(tv));

	R_xlen_t nv = XLENGTH(val);
	if (nv != 1 && nv != r)
		error(_("replacement diagonal has wrong length"));

	SEXP x;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pid);
	SEXPTYPE tx = TYPEOF(x);

	/* Allocate and coerce as necessary */
	SEXP res;
	if (tv <= tx) {
		PROTECT(val = coerceVector(val, tv = tx));
		PROTECT(res = NEW_OBJECT_OF_CLASS(valid[ivalid]));
		REPROTECT(x = duplicate(x), pid);
	} else { /* tv > tx */
		/* dMatrix is only possibility until we have proper [iz]Matrix */
		PROTECT(val = coerceVector(val, tv = REALSXP));
		char cl[] = "d..Matrix";
		cl[1] = valid[ivalid][1];
		cl[2] = valid[ivalid][2];
		PROTECT(res = NEW_OBJECT_OF_CLASS(cl));
		REPROTECT(x = coerceVector(x, tx = tv), pid);
	}

	if (m != n || n > 0)
		SET_SLOT(res, Matrix_DimSym, dim);

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (valid[ivalid][1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(res, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

#define UPM_D_S(_CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pval = _PTR_(val); \
		R_xlen_t m1a = (R_xlen_t) m + 1; \
		int j; \
		if (nv == 1) \
			for (j = 0; j < r; ++j, px += m1a) \
				*px = *pval; \
		else \
			for (j = 0; j < r; ++j, px += m1a) \
				*px = *(pval++); \
	} while (0)

	switch (tx) {
	case LGLSXP:
		UPM_D_S(int, LOGICAL);
		break;
	case INTSXP:
		UPM_D_S(int, INTEGER);
		break;
	case REALSXP:
		UPM_D_S(double, REAL);
		break;
	case CPLXSXP:
		UPM_D_S(Rcomplex, COMPLEX);
		break;
	default:
		ERROR_INVALID_TYPE(x, __func__);
		break;
	}

#undef UPM_D_S

	SET_SLOT(res, Matrix_xSym, x);

	UNPROTECT(4); /* x, res, val, dim */
	return res;
}

/* symmpart(x) */
SEXP unpackedMatrix_symmpart(SEXP from)
{
	static const char *valid[] = {
	/* 0 */ "dgeMatrix", "lgeMatrix", "ngeMatrix",
	/* 3 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	const char *clf = valid[ivalid];
	if (clf[0] == 'd' && clf[1] == 's')
		return from;

	char clt[] = ".syMatrix";
	clt[0] = (clf[0] != 'z') ? 'd' : 'z';
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to get symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (clf[1] != 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP x;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);
	REPROTECT(x = (clf[0] == clt[0]) ? duplicate(x) : coerceVector(x, REALSXP),
	          pid);
	if (clf[0] == 'n')
		na2one(x);

	if (clf[1] == 'g') {

		int i, j;
		R_xlen_t upos = 0, lpos = 0;

#define UPM_SYMMPART_GE(_CTYPE_, _PTR_, _ASSIGN_OFFDIAG_) \
		do { \
			_CTYPE_ *px = _PTR_(x); \
			for (j = 0; j < n; ++j) { \
				for (i = j+1; i < n; ++i) { \
					upos += n; ++lpos; \
					_ASSIGN_OFFDIAG_(upos, lpos); \
				} \
				upos = (lpos += j+2); \
			} \
		} while (0)

#define ASSIGN_OFFDIAG_DGE(_UPOS_, _LPOS_) \
		do { \
			px[_UPOS_] += px[_LPOS_]; \
			px[_UPOS_] *= 0.5; \
		} while (0)

#define ASSIGN_OFFDIAG_ZGE(_UPOS_, _LPOS_) \
		do { \
			px[_UPOS_].r += px[_LPOS_].r; \
			px[_UPOS_].i += px[_LPOS_].i; \
			px[_UPOS_].r *= 0.5; \
			px[_UPOS_].i *= 0.5; \
		} while (0)

			if (clf[0] != 'z')
				UPM_SYMMPART_GE(double, REAL, ASSIGN_OFFDIAG_DGE);
			else
				UPM_SYMMPART_GE(Rcomplex, COMPLEX, ASSIGN_OFFDIAG_ZGE);

	} else {

		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		if (clf[1] != 's') {

			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			char di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */

			int i, j;

#define UPM_SYMMPART_TR(_CTYPE_, _PTR_, _ASSIGN_OFFDIAG_, _ASSIGN_ONDIAG_) \
			do { \
				_CTYPE_ *px = _PTR_(x); \
				if (ul == 'U') { \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i, ++px) \
							_ASSIGN_OFFDIAG_; \
						px += n-j; \
					} \
				} else { \
					for (j = 0; j < n; ++j) { \
						px += j+1; \
						for (i = j+1; i < n; ++i, ++px) \
							_ASSIGN_OFFDIAG_; \
					} \
				} \
				if (di != 'N') { \
					R_xlen_t n1a = (R_xlen_t) n + 1; \
					px = _PTR_(x); \
					for (j = 0; j < n; ++j, px += n1a) \
						_ASSIGN_ONDIAG_; \
				} \
			} while (0)

			if (clt[0] != 'z')
				UPM_SYMMPART_TR(double, REAL,
				                *px *= 0.5,
				                *px  = 1.0);
			else
				UPM_SYMMPART_TR(Rcomplex, COMPLEX,
				                do { (*px).r *= 0.5; (*px).i *= 0.5; } while (0),
				                do { (*px).r  = 1.0; (*px).i  = 0.0; } while (0));

		} else { /* clf[1] == 's' */

			if (clt[0] == 'z')
				/* Symmetric part of Hermitian matrix is real part */
				zeroIm(x);

		}

	}

	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(2); /* x, to */
	return to;
}

/* symmpart(x) */
SEXP matrix_symmpart(SEXP from)
{
	SEXP dim = PROTECT(getAttrib(from, R_DimSymbol));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to get symmetric part of non-square matrix"));

	SEXP to, x;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(x = from, &pid);

	int i, j;
	R_xlen_t upos = 0, lpos = 0, nn = (R_xlen_t) n * n;

	switch (TYPEOF(x)) {
	case LGLSXP:
	case INTSXP:
		REPROTECT(x = coerceVector(x, REALSXP), pid);
	case REALSXP:
		PROTECT(to = NEW_OBJECT_OF_CLASS("dsyMatrix"));
		if (!MAYBE_REFERENCED(x))
			SET_ATTRIB(x, R_NilValue);
		else {
			REPROTECT(x = allocVector(REALSXP, nn), pid);
			Matrix_memcpy(REAL(x), REAL(from), nn, sizeof(double));
		}
		UPM_SYMMPART_GE(double, REAL, ASSIGN_OFFDIAG_DGE);
		break;
#ifdef MATRIX_ENABLE_ZMATRIX
	case CPLXSXP:
		PROTECT(to = NEW_OBJECT_OF_CLASS("zsyMatrix"));
		if (!MAYBE_REFERENCED(x))
			SET_ATTRIB(x, R_NilValue);
		else {
			REPROTECT(x = allocVector(CPLXSXP, nn), pid);
			Matrix_memcpy(COMPLEX(x), COMPLEX(from), nn, sizeof(Rcomplex));
		}
		UPM_SYMMPART_GE(Rcomplex, COMPLEX, ASSIGN_OFFDIAG_ZGE);
		break;
#endif
	default:
		ERROR_INVALID_TYPE(x, __func__);
		break;
	}

	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	SET_SLOT(to, Matrix_xSym, x);

	SEXP dimnames = PROTECT(getAttrib(from, R_DimNamesSymbol));
	if (!isNull(dimnames))
		set_symmetrized_DimNames(to, dimnames, -1);

	UNPROTECT(4); /* dimnames, to, x, dim */
	return to;
}

#undef ASSIGN_OFFDIAG_DGE
#undef ASSIGN_OFFDIAG_ZGE
#undef UPM_SYMMPART_GE
#undef UPM_SYMMPART_TR

/* skewpart(x) */
SEXP unpackedMatrix_skewpart(SEXP from)
{
	static const char *valid[] = {
	/* 0 */ "dgeMatrix", "lgeMatrix", "ngeMatrix",
	/* 3 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *clf = valid[ivalid];

	char clt[] = "...Matrix";
	clt[0] = (clf[0] != 'z') ? 'd' : 'z';
	clt[1] = (clf[1] != 's') ? 'g' : 's';
	clt[2] = (clf[1] != 's') ? 'e' : ((clf[0] != 'z') ? 'C' : 'y');
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to get skew-symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1);

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (clf[1] != 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1);

	char ul = 'U';
	if (clf[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (clf[1] == 's' && ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	SEXP x;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);

	if (clf[1] != 's') {

		SEXP y;
		int i, j;
		R_xlen_t upos = 0, lpos = 0;

#define UPM_SKEWPART(_CTYPE_, _PTR_, _X_, _Y_, \
		             _ASSIGN_OFFDIAG_, _ASSIGN_ONDIAG_) \
		do { \
			_CTYPE_ *px = _PTR_(_X_), *py = _PTR_(_Y_); \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					lpos = j; \
					for (i = 0; i < j; ++i) { \
						_ASSIGN_OFFDIAG_(upos, lpos); \
						++upos; lpos += n; \
					} \
					_ASSIGN_ONDIAG_(upos); \
					upos += n-j; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					upos = lpos; \
					_ASSIGN_ONDIAG_(lpos); \
					for (i = j+1; i < n; ++i) { \
						upos += n; ++lpos; \
						_ASSIGN_OFFDIAG_(lpos, upos); \
					} \
					lpos += j+2; \
				} \
			} \
		} while (0)

#define ASSIGN_ONDIAG_DGE(_DPOS_) \
		py[_DPOS_] = 0.0

#define ASSIGN_OFFDIAG_DGE(_UPOS_, _LPOS_) \
		do { \
			py[_UPOS_] = 0.5 * (px[_UPOS_] - px[_LPOS_]); \
			py[_LPOS_] = -py[_UPOS_]; \
		} while (0)

#define ASSIGN_OFFDIAG_DTR(_UPOS_, _LPOS_) \
		do { \
			py[_UPOS_] = 0.5 * px[_UPOS_]; \
			py[_LPOS_] = -py[_UPOS_]; \
		} while (0)

#define ASSIGN_ONDIAG_ZGE(_DPOS_) \
		py[_DPOS_].r = py[_DPOS_].i = 0.0

#define ASSIGN_OFFDIAG_ZGE(_UPOS_, _LPOS_) \
		do { \
			py[_UPOS_].r = 0.5 * (px[_UPOS_].r - px[_LPOS_].r); \
			py[_UPOS_].i = 0.5 * (px[_UPOS_].i - px[_LPOS_].i); \
			py[_LPOS_].r = -py[upos].r; \
			py[_LPOS_].i = -py[upos].i; \
		} while (0)

#define ASSIGN_OFFDIAG_ZTR(_UPOS_, _LPOS_) \
		do { \
			py[_UPOS_].r = 0.5 * px[_UPOS_].r; \
			py[_UPOS_].i = 0.5 * px[_UPOS_].i; \
			py[_LPOS_].r = -py[upos].r; \
			py[_LPOS_].i = -py[upos].i; \
		} while (0)

		if (clf[0] != 'z') {
			if (clf[0] == 'd')
				PROTECT(y = allocVector(REALSXP, (R_xlen_t) n * n));
			else
				PROTECT(x = y = coerceVector(x, REALSXP));
			if (clf[1] == 'g')
				UPM_SKEWPART(double, REAL, x, y,
				             ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
			else
				UPM_SKEWPART(double, REAL, x, y,
				             ASSIGN_OFFDIAG_DTR, ASSIGN_ONDIAG_DGE);
			} else { /* clf[0] == 'z' */
			PROTECT(y = allocVector(CPLXSXP, (R_xlen_t) n * n));
			if (clf[1] == 'g')
				UPM_SKEWPART(Rcomplex, COMPLEX, x, y,
				             ASSIGN_OFFDIAG_ZGE, ASSIGN_ONDIAG_ZGE);
			else
				UPM_SKEWPART(Rcomplex, COMPLEX, x, y,
				             ASSIGN_OFFDIAG_ZTR, ASSIGN_ONDIAG_ZGE);
		}

		SET_SLOT(to, Matrix_xSym, y);
		UNPROTECT(1); /* y */

	} else { /* clf[1] == 's' */

		if (clf[0] != 'z') {
			/* Skew-symmetric part of symmetric matrix is zero matrix */
			R_xlen_t n1a = (R_xlen_t) n + 1;
			SEXP p = PROTECT(allocVector(INTSXP, n1a));
			int *pp = INTEGER(p);
			Matrix_memset(pp, 0, n1a, sizeof(int));
			SET_SLOT(to, Matrix_pSym, p);
			UNPROTECT(1); /* p */
		} else {
			/* Skew-symmetric part of Hermitian matrix is imaginary part */
			REPROTECT(x = duplicate(x), pid);
			zeroRe(x);
			SET_SLOT(to, Matrix_xSym, x);
		}

	}

	UNPROTECT(2); /* x, to */
	return to;
}

/* skewpart(x) */
SEXP matrix_skewpart(SEXP from)
{
	SEXP dim = PROTECT(getAttrib(from, R_DimSymbol));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to get skew-symmetric part of non-square matrix"));

	SEXP to, x;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(x = from, &pid);

	char ul = 'U';
	int i, j;
	R_xlen_t upos = 0, lpos = 0;

	switch (TYPEOF(x)) {
	case LGLSXP:
	case INTSXP:
		REPROTECT(x = coerceVector(x, REALSXP), pid);
	case REALSXP:
		PROTECT(to = NEW_OBJECT_OF_CLASS("dgeMatrix"));
		if (!MAYBE_REFERENCED(x)) {
			SET_ATTRIB(x, R_NilValue);
			UPM_SKEWPART(double, REAL, x, x,
			             ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
		} else {
			REPROTECT(x = allocVector(REALSXP, (R_xlen_t) n * n), pid);
			UPM_SKEWPART(double, REAL, from, x,
			             ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
		}
		break;
#ifdef MATRIX_ENABLE_ZMATRIX
	case CPLXSXP:
		PROTECT(to = NEW_OBJECT_OF_CLASS("zgeMatrix"));
		if (!MAYBE_REFERENCED(from)) {
			SET_ATTRIB(x, R_NilValue);
			UPM_SKEWPART(Rcomplex, COMPLEX, x, x,
			             ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
		} else {
			REPROTECT(x = allocVector(CPLXSXP, (R_xlen_t) n * n), pid);
			UPM_SKEWPART(Rcomplex, COMPLEX, from, x,
			             ASSIGN_OFFDIAG_ZGE, ASSIGN_ONDIAG_ZGE);
		}
		break;
#endif
	default:
		ERROR_INVALID_TYPE(x, __func__);
		break;
	}

	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	SET_SLOT(to, Matrix_xSym, x);

	SEXP dimnames = PROTECT(getAttrib(from, R_DimNamesSymbol));
	if (!isNull(dimnames))
		set_symmetrized_DimNames(to, dimnames, -1);

	UNPROTECT(4); /* dimnames, to, x, dim */
	return to;
}

#undef ASSIGN_ONDIAG_DGE
#undef ASSIGN_ONDIAG_ZGE
#undef ASSIGN_OFFDIAG_DGE
#undef ASSIGN_OFFDIAG_DTR
#undef ASSIGN_OFFDIAG_ZGE
#undef ASSIGN_OFFDIAG_ZTR
#undef UPM_SKEWPART
