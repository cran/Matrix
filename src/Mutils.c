#include "Mutils.h"

/**
 * A safe version of `NEW_OBJECT(MAKE_CLASS(what))`, protecting the
 * intermediate R object.  The caller must protect the return value
 * of this function.
 *
 * @param A string specifying the name of a defined S4 class.
 */
SEXP NEW_OBJECT_OF_CLASS(const char *what)
{
	SEXP class = PROTECT(MAKE_CLASS(what)), obj = NEW_OBJECT(class);
	UNPROTECT(1);
	return obj;
}

/* memset() but passing length and size rather than their product
   which can overflow size_t ... hence _safer_ than Memzero()
*/
void *Matrix_memset(void *dest, int ch, R_xlen_t length, size_t size)
{
	if (dest && length > 0 && size > 0) {

		char *dest_ = (char *) dest;
		size_t N = SIZE_MAX / size;

#if (SIZE_MAX < R_XLEN_T_MAX)
		R_xlen_t S_M = (R_xlen_t) SIZE_MAX;
		if (length <= S_M) {
#endif

			/* 'length' is representable as size_t : */

			size_t n = (size_t) length;
			if (n <= N)
				memset(dest_, ch, n * size);
			else {
				size_t d = N * size;
				while (n > N) {
					memset(dest_, ch, d);
					dest_ += d;
					n -= d;
				}
				memset(dest_, ch, n * size);
			}

#if (SIZE_MAX < R_XLEN_T_MAX)
		} else {

			/* 'length' would overflow size_t : */

			size_t n, d = N * size;
			while (length > S_M) {
				n = SIZE_MAX;
				while (n > N) {
					memset(dest_, ch, d);
					dest_ += d;
					n -= d;
				}
				memset(dest_, ch, n * size);
				length -= S_M;
			}
			n = (size_t) length;
			while (n > N) {
				memset(dest_, ch, d);
				dest_ += d;
				n -= d;
			}
			memset(dest_, ch, n * size);

		}
#endif

	}

	return dest;
}

/* memcpy() but passing length and size rather than their product
   which can overflow size_t ... hence _safer_ than Memcpy()
*/
void *Matrix_memcpy(void *dest, const void *src, R_xlen_t length, size_t size)
{
	if (dest && src && length > 0 && size > 0) {

		char *dest_ = (char *) dest;
		const char *src_ = (const char *) src;

		size_t N = SIZE_MAX / size;

#if (SIZE_MAX < R_XLEN_T_MAX)
		R_xlen_t S_M = (R_xlen_t) SIZE_MAX;
		if (length <= S_M) {
#endif

			/* 'length' is representable as size_t : */

			size_t n = (size_t) length;
			if (n <= N)
				memcpy(dest_, src_, n * size);
			else {
				size_t d = N * size;
				while (n > N) {
					memcpy(dest_, src_, d);
					dest_ += d;
					src_ += d;
					n -= d;
				}
				memcpy(dest_, src_, n * size);
			}

#if (SIZE_MAX < R_XLEN_T_MAX)
		} else {

			/* 'length' would overflow size_t : */

			size_t n, d = N * size;
			while (length > S_M) {
				n = SIZE_MAX;
				while (n > N) {
					memcpy(dest_, src_, d);
					dest_ += d;
					src_ += d;
					n -= d;
				}
				memcpy(dest_, src_, n * size);
				length -= S_M;
			}
			n = (size_t) length;
			while (n > N) {
				memcpy(dest_, src_, d);
				dest_ += d;
				n -= d;
			}
			memcpy(dest_, src_, n * size);

		}
#endif

	}

	return dest;
}


/* More for 'Dimnames' ============================================== */

Rboolean DimNames_is_trivial(SEXP dn)
{
	if (!(isNull(VECTOR_ELT(dn, 0)) &&
	    isNull(VECTOR_ELT(dn, 1))))
		return FALSE;
	Rboolean res = TRUE;
	SEXP ndn = PROTECT(getAttrib(dn, R_NamesSymbol));
	if (!isNull(ndn))
		res = FALSE;
	UNPROTECT(1);
	return res;
}

Rboolean DimNames_is_symmetric(SEXP dn)
{
	/* NB: Assuming here that we have the 'Dimnames' slot
	   of a _valid_ Matrix, so that the elements are either
	   NULL or character vectors

	   Keep synchronized with symmetricMatrix_validate() above,
	   (which must do slightly more)! */

	SEXP rn, cn;
	int n;
	if (!isNull(rn = VECTOR_ELT(dn, 0)) &&
	    !isNull(cn = VECTOR_ELT(dn, 1)) &&
	    rn != cn &&
	    ((n = LENGTH(rn)) != LENGTH(cn) || !equal_string_vectors(rn, cn, n)))
		return FALSE;
	Rboolean res = TRUE;
	SEXP ndn = PROTECT(getAttrib(dn, R_NamesSymbol));
	const char *ndn0, *ndn1;
	if (!isNull(ndn) &&
	    *(ndn0 = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
	    *(ndn1 = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
	    strcmp(ndn0, ndn1) != 0)
		res = FALSE;
	UNPROTECT(1);
	return res;
}

SEXP R_DimNames_is_symmetric(SEXP dn)
{
	return ScalarLogical(DimNames_is_symmetric(dn));
}

/**
 * @brief Produce symmetric `Dimnames` from possibly asymmetric ones.
 *
 * Roughly `dest[1:2] <- rep(src[j], 2)`, where `j` is either 1 or 2
 * depending on `J`.  If `J` is 0 or 1, then `j = J+1`.  If `J` is -1,
 * then `j = 1` if and only if `src[[2]]` is `NULL` and `src[[1]]`
 * is not.  For speed, it is assumed that `dest` is newly allocated,
 * i.e., that it is `list(NULL, NULL)`.
 *
 * @param dest,src Lists of length 2, typically the `Dimnames` slots
 *     of two square `Matrix` of equal size.
 * @param J An integer, one of -1, 0, and 1.
 */
void symmDN(SEXP dest, SEXP src, int J /* -1|0|1 */)
{
	SEXP s;
	if (J < 0) {
		if (!isNull(s = VECTOR_ELT(src, J = 1)) ||
			!isNull(s = VECTOR_ELT(src, J = 0))) {
			SET_VECTOR_ELT(dest, 0, s);
			SET_VECTOR_ELT(dest, 1, s);
		} else {
			J = 1;
		}
	} else {
		if (!isNull(s = VECTOR_ELT(src, J))) {
			SET_VECTOR_ELT(dest, 0, s);
			SET_VECTOR_ELT(dest, 1, s);
		}
	}
	/* names(dimnames(.)) */
	PROTECT(s = getAttrib(src, R_NamesSymbol));
	if (!isNull(s)) {
		SEXP destnms = PROTECT(allocVector(STRSXP, 2));
		if (*CHAR(s = STRING_ELT(s, J)) != '\0') {
			SET_STRING_ELT(destnms, 0, s);
			SET_STRING_ELT(destnms, 1, s);
		}
		setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return;
}

/**
 * @brief Reverse (or "transpose") `Dimnames`.
 *
 * Roughly `dest[1:2] <- src[2:1]`.  For speed, it is assumed that
 * `dest` is newly allocated, i.e., that it is `list(NULL, NULL)`.
 *
 * @param dest,src Lists of length 2, typically the `Dimnames` slots
 *     of two square `Matrix` of equal size.
 */
void revDN(SEXP dest, SEXP src) {
	SEXP s;
	if (!isNull(s = VECTOR_ELT(src, 0)))
		SET_VECTOR_ELT(dest, 1, s);
	if (!isNull(s = VECTOR_ELT(src, 1)))
		SET_VECTOR_ELT(dest, 0, s);
	PROTECT(s = getAttrib(src, R_NamesSymbol));
	if (!isNull(s)) {
		SEXP srcnms = s, destnms = PROTECT(allocVector(STRSXP, 2));
		if (*CHAR(s = STRING_ELT(srcnms, 0)) != '\0')
			SET_STRING_ELT(destnms, 1, s);
		if (*CHAR(s = STRING_ELT(srcnms, 1)) != '\0')
			SET_STRING_ELT(destnms, 0, s);
		setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return;
}

SEXP R_symmDN(SEXP dn)
{
	/* Be fast (do nothing!) when dimnames = list(NULL, NULL) */
	if (DimNames_is_trivial(dn))
		return dn;
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	symmDN(newdn, dn, -1);
	UNPROTECT(1);
	return newdn;
}

SEXP R_revDN(SEXP dn)
{
	/* Be fast (do nothing!) when dimnames = list(NULL, NULL) */
	if (DimNames_is_trivial(dn))
		return dn;
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	revDN(newdn, dn);
	UNPROTECT(1);
	return newdn;
}

SEXP get_symmetrized_DimNames(SEXP obj, int J) {
	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	if (DimNames_is_trivial(dn)) {
		UNPROTECT(1);
		return dn;
	}
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	symmDN(newdn, dn, J);
	UNPROTECT(2);
	return newdn;
}

SEXP get_reversed_DimNames(SEXP obj) {
	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	if (DimNames_is_trivial(dn)) {
		UNPROTECT(1);
		return dn;
	}
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	revDN(newdn, dn);
	UNPROTECT(2);
	return newdn;
}

void set_symmetrized_DimNames(SEXP obj, SEXP dn, int J) {
	if (!DimNames_is_trivial(dn)) {
		SEXP newdn = PROTECT(allocVector(VECSXP, 2));
		symmDN(newdn, dn, J);
		SET_SLOT(obj, Matrix_DimNamesSym, newdn);
		UNPROTECT(1);
	}
	return;
}

void set_reversed_DimNames(SEXP obj, SEXP dn) {
	if (!DimNames_is_trivial(dn)) {
		SEXP newdn = PROTECT(allocVector(VECSXP, 2));
		revDN(newdn, dn);
		SET_SLOT(obj, Matrix_DimNamesSym, newdn);
		UNPROTECT(1);
	}
	return;
}

void set_DimNames(SEXP obj, SEXP dn)
{
	if (!DimNames_is_trivial(dn)) {
		SEXP s, newdn = PROTECT(allocVector(VECSXP, 2));
		if (!isNull(s = VECTOR_ELT(dn, 0)))
			SET_VECTOR_ELT(newdn, 0, s);
		if (!isNull(s = VECTOR_ELT(dn, 1)))
			SET_VECTOR_ELT(newdn, 1, s);
		PROTECT(s = getAttrib(dn, R_NamesSymbol));
		if (!isNull(s))
			setAttrib(newdn, R_NamesSymbol, s);
		SET_SLOT(obj, Matrix_DimNamesSym, newdn);
		UNPROTECT(2);
	}
	return;
}


/* For 'factors' ==================================================== */

SEXP get_factor(SEXP obj, const char *nm)
{
	SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorSym)), val = R_NilValue;
	if (LENGTH(factors) > 0) {
		SEXP valid = PROTECT(getAttrib(factors, R_NamesSymbol));
		int i = strmatch2(nm, valid);
		if (i >= 0)
			val = VECTOR_ELT(factors, i);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return val;
}

void set_factor(SEXP obj, const char *nm, SEXP val)
{
	PROTECT(val);
	SEXP factors;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(factors = GET_SLOT(obj, Matrix_factorSym), &pid);
	if (LENGTH(factors) > 0) {
		SEXP valid = PROTECT(getAttrib(factors, R_NamesSymbol));
		int i = strmatch2(nm, valid);
		UNPROTECT(1);
		if (i >= 0) {
			SET_VECTOR_ELT(factors, i, val);
			UNPROTECT(2);
			return;
		}
	}
	REPROTECT(factors = append_to_named_list(factors, nm, val), pid);
	SET_SLOT(obj, Matrix_factorSym, factors);
	UNPROTECT(2);
	return;
}

/**
 * @brief Subassign by name to the `factors` slot of a `compMatrix`.
 *
 * Like `obj\@factors[[nm]] <- val`, but modifying `obj` (rather than a copy)
 * _even if_ `obj` is referenced elsewhere, supporting "automagic" caching of
 * factorizations by R functions taking `compMatrix` as an argument.
 * _Use with care!_
 *
 * @param obj A `compMatrix`.
 * @param nm A length-1 `STRSXP` giving a factor name.
 * @param val A `SEXP`, usually a `MatrixFactorization`.
 * @param warn A length-1 `LGLSXP`. Warn if `obj` has no `factors` slot
 *     (in which case `obj` is untouched)?
 *
 * @return `val`.
 */
SEXP R_set_factor(SEXP obj, SEXP nm, SEXP val, SEXP warn)
{
	if (TYPEOF(nm) != STRSXP || LENGTH(nm) < 1 ||
		(nm = STRING_ELT(nm, 0)) == NA_STRING)
	error(_("invalid factor name"));
	else if (HAS_SLOT(obj, Matrix_factorSym))
		set_factor(obj, CHAR(nm), val);
	else if (asLogical(warn) != 0)
		warning(_("attempt to set factor on %s without '%s' slot"),
		        "Matrix", "factors");
	return val;
}

/**
 * @brief Empty the 'factors' slot of a 'compMatrix'.
 *
 * Like `obj\@factors <- list()`, but modifying `obj` (rather than a copy)
 * _even if_ `obj` is referenced elsewhere, supporting "automagic" clearing
 * of the `factors` slot by R functions taking `compMatrix` as an argument.
 * _Use with care!_
 *
 * @param obj A `compMatrix`.
 * @param warn A length-1 LGLSXP. Warn if `obj` has no `factors` slot
 *     (in which case `obj` is untouched)?
 *
 * @return `TRUE` if `obj` has a nonempty `factors` slot, `FALSE` otherwise.
 */
SEXP R_empty_factors(SEXP obj, SEXP warn)
{
	/* If there is a nonempty 'factors' slot, then replace it with list() */
	if (HAS_SLOT(obj, Matrix_factorSym)) {
		SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorSym));
		if (LENGTH(factors) > 0) {
			PROTECT(factors = allocVector(VECSXP, 0));
			SET_SLOT(obj, Matrix_factorSym, factors);
			UNPROTECT(2);
			return ScalarLogical(1); /* slot was reset */
		}
		UNPROTECT(1);
	} else if (asLogical(warn) != 0)
		warning(_("attempt to discard factors from %s without '%s' slot"),
		        "Matrix", "factors");
	return ScalarLogical(0); /* no-op */
}


/* For permutations ================================================= */

/* Poor man's C translation of LAPACK dswap */
static void dswap(int n, double *x, int incx, double *y, int incy)
{
	double tmp;
	while (n--) {
		tmp = *x;
		*x = *y;
		*y = tmp;
		x += incx;
		y += incy;
	}
	return;
}

/* Poor man's C translation of LAPACK dsyswapr */
static void dsyswapr(char uplo, int n, double *x, int k0, int k1)
{
	double tmp, *x0 = x + (R_xlen_t) k0 * n, *x1 = x + (R_xlen_t) k1 * n;
	if (uplo == 'U') {
		dswap(k0, x0, 1, x1, 1);
		tmp = x0[k0];
		x0[k0] = x1[k1];
		x1[k1] = tmp;
		dswap(k1 - k0 - 1, x0 + k0 + n, n, x1 + k0 + 1, 1);
		dswap(n - k1 - 1, x1 + k0 + n, n, x1 + k1 + n, n);
	} else {
		dswap(k0, x + k0, n, x + k1, n);
		tmp = x0[k0];
		x0[k0] = x1[k1];
		x1[k1] = tmp;
		dswap(k1 - k0 - 1, x0 + k0 + 1, 1, x0 + k1 + n, n);
		dswap(n - k1 - 1, x0 + k1 + 1, 1, x1 + k1 + 1, 1);
	}
	return;
}

void rowPerm(double *x, int m, int n, int *p, int off, int invert)
{
	int i, k0, k1;
	for (i = 0; i < m; ++i)
		p[i] = -(p[i] - off + 1);
	if (!invert) {
		for (i = 0; i < m; ++i) {
			if (p[i] > 0)
				continue;
			k0 = i;
			p[k0] = -p[k0];
			k1 = p[k0] - 1;
			while (p[k1] < 0) {
				dswap(n, x + k0, m, x + k1, m);
				k0 = k1;
				p[k0] = -p[k0];
				k1 = p[k0] - 1;
			}
		}
	} else {
		for (i = 0; i < m; ++i) {
			if (p[i] > 0)
				continue;
			k0 = i;
			p[k0] = -p[k0];
			k1 = p[k0] - 1;
			while (k1 != k0) {
				dswap(n, x + k0, m, x + k1, m);
				p[k1] = -p[k1];
				k1 = p[k1] - 1;
			}
		}
	}
	for (i = 0; i < m; ++i)
		p[i] = p[i] + off - 1;
	return;
}

void symPerm(double *x, int n, char uplo, int *p, int off, int invert)
{
	int i, k0, k1;
	for (i = 0; i < n; ++i)
		p[i] = -(p[i] - off + 1);
	if (!invert) {
		for (i = 0; i < n; ++i) {
			if (p[i] > 0)
				continue;
			k0 = i;
			p[k0] = -p[k0];
			k1 = p[k0] - 1;
			while (p[k1] < 0) {
				if (k0 < k1)
					dsyswapr(uplo, n, x, k0, k1);
				else
					dsyswapr(uplo, n, x, k1, k0);
				k0 = k1;
				p[k0] = -p[k0];
				k1 = p[k0] - 1;
			}
		}
	} else {
		for (i = 0; i < n; ++i) {
			if (p[i] > 0)
				continue;
			k0 = i;
			p[k0] = -p[k0];
			k1 = p[k0] - 1;
			while (k1 != k0) {
				if (k0 < k1)
					dsyswapr(uplo, n, x, k0, k1);
				else
					dsyswapr(uplo, n, x, k1, k0);
				p[k1] = -p[k1];
				k1 = p[k1] - 1;
			}
		}
	}
	for (i = 0; i < n; ++i)
		p[i] = p[i] + off - 1;
	return;
}

int isPerm(const int *p, int n, int off)
{
	int res = 1;
	if (n <= 0)
		return res;
	int i, j;
	char *work;
	Matrix_Calloc(work, n, char);
	for (i = 0; i < n; ++i) {
		if (p[i] == NA_INTEGER || (j = p[i] - off) < 0 || j >= n || work[j]) {
			res = 0;
			break;
		}
		work[j] = 1;
	}
	Matrix_Free(work, n);
	return res;
}

int signPerm(const int *p, int n, int off)
{
	if (!isPerm(p, n, off))
		error(_("attempt to get sign of non-permutation"));
	int sign = 1;
	if (n <= 0)
		return sign;
	int i, pos = 0;
	char *work;
	Matrix_Calloc(work, n, char);
	while (pos < n) {
		work[pos] = 1;
		i = p[pos] - off;
		while (!work[i]) { /* transposition */
			sign = -sign;
			work[i] = 1;
			i = p[i] - off;
		}
		while (pos < n && work[pos])
			++pos;
	}
	Matrix_Free(work, n);
	return sign;
}

void invertPerm(const int *p, int *ip, int n, int off, int ioff)
{
	if (!isPerm(p, n, off))
		error(_("attempt to invert non-permutation"));
	int j;
	for (j = 0; j < n; ++j)
		ip[p[j] - off] = j + ioff;
	return;
}

void asPerm(const int *p, int *ip, int m, int n, int off, int ioff)
{
	int i, j, tmp;
	for (i = 0; i < n; ++i)
		ip[i] = i + ioff;
	for (i = 0; i < m; ++i) {
		j = p[i] - off;
		if (j < 0 || j >= n)
			error(_("invalid transposition vector"));
		if (j != i) {
			tmp = ip[j];
			ip[j] = ip[i];
			ip[i] = tmp;
		}
	}
	return;
}

SEXP R_isPerm(SEXP p, SEXP off)
{
	if (TYPEOF(p) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "p", "integer");
	if (TYPEOF(off) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "off", "integer");
	if (XLENGTH(off) != 1)
		error(_("'%s' does not have length %d"), "off", 1);
	int off_ = INTEGER(off)[0];
	if (off_ == NA_INTEGER)
		error(_("'%s' is NA"), "off");
	R_xlen_t n_ = XLENGTH(p);
	if (n_ > INT_MAX)
		return ScalarLogical(0);
	return ScalarLogical(isPerm(INTEGER(p), (int) n_, off_));
}

SEXP R_signPerm(SEXP p, SEXP off)
{
	if (TYPEOF(p) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "p", "integer");
	if (TYPEOF(off) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "off", "integer");
	if (XLENGTH(off) != 1)
		error(_("'%s' does not have length %d"), "off", 1);
	int off_ = INTEGER(off)[0];
	if (off_ == NA_INTEGER)
		error(_("'%s' is NA"), "off");
	R_xlen_t n_ = XLENGTH(p);
	if (n_ > INT_MAX)
		error(_("attempt to get sign of non-permutation"));
	return ScalarInteger(signPerm(INTEGER(p), (int) n_, off_));
}

SEXP R_invertPerm(SEXP p, SEXP off, SEXP ioff)
{
	if (TYPEOF(p) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "p", "integer");
	if (TYPEOF(off) != INTSXP || TYPEOF(ioff) != INTSXP)
		error(_("'%s' or '%s' is not of type \"%s\""), "off", "ioff", "integer");
	if (XLENGTH(off) != 1 || XLENGTH(ioff) != 1)
		error(_("'%s' or '%s' does not have length %d"), "off", "ioff", 1);
	int off_ = INTEGER(off)[0], ioff_ = INTEGER(ioff)[0];
	if (off_ == NA_INTEGER || ioff_ == NA_INTEGER)
		error(_("'%s' or '%s' is NA"), "off", "ioff");
	R_xlen_t n_ = XLENGTH(p);
	if (n_ > INT_MAX)
		error(_("attempt to invert non-permutation"));
	SEXP ip = PROTECT(allocVector(INTSXP, n_));
	invertPerm(INTEGER(p), INTEGER(ip), (int) n_, off_, ioff_);
	UNPROTECT(1);
	return ip;
}

SEXP R_asPerm(SEXP p, SEXP off, SEXP ioff, SEXP n)
{
	if (TYPEOF(p) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "p", "integer");
	R_xlen_t m_ = XLENGTH(p);
	if (m_ > INT_MAX)
		error(_("'%s' has length exceeding %s"), "p", "2^31-1");
	if (TYPEOF(off) != INTSXP || TYPEOF(ioff) != INTSXP)
		error(_("'%s' or '%s' is not of type \"%s\""), "off", "ioff", "integer");
	if (XLENGTH(off) != 1 || XLENGTH(ioff) != 1)
		error(_("'%s' or '%s' does not have length %d"), "off", "ioff", 1);
	int off_ = INTEGER(off)[0], ioff_ = INTEGER(ioff)[0];
	if (off_ == NA_INTEGER || ioff_ == NA_INTEGER)
		error(_("'%s' or '%s' is NA"), "off", "ioff");
	if (TYPEOF(n) != INTSXP)
		error(_("'%s' is not of type \"%s\""), "n", "integer");
	if (XLENGTH(n) != 1)
		error(_("'%s' does not have length %d"), "n", 1);
	int n_ = INTEGER(n)[0];
	if (n_ == NA_INTEGER || n_ < m_)
		error(_("'%s' is NA or less than %s"), "n", "length(p)");
	SEXP ip = PROTECT(allocVector(INTSXP, n_));
	asPerm(INTEGER(p), INTEGER(ip), (int) m_, n_, off_, ioff_);
	UNPROTECT(1);
	return ip;
}


/* For inheritance ================================================== */

char type2kind(SEXPTYPE type)
{
	switch (type) {
	case LGLSXP:
		return 'l';
	case INTSXP:
#ifdef MATRIX_ENABLE_IMATRIX
		return 'i';
#endif
	case REALSXP:
		return 'd';
#ifdef MATRIX_ENABLE_ZMATRIX
	case CPLXSXP:
		return 'z';
#endif
	default:
		error(_("unexpected type \"%s\" in %s()"), type2char(type), __func__);
		return '\0';
	}
}

SEXPTYPE kind2type(char kind)
{
	switch (kind) {
	case 'n':
	case 'l':
		return LGLSXP;
#ifdef MATRIX_ENABLE_IMATRIX
	case 'i':
		return INTSXP;
#endif
	case 'd':
		return REALSXP;
#ifdef MATRIX_ENABLE_ZMATRIX
	case 'z':
		return CPLXSXP;
#endif
	default:
		error(_("unexpected kind \"%c\" in %s()"), kind, __func__);
		return NILSXP;
	}
}

size_t kind2size(char kind)
{
	switch (kind) {
	case 'n':
	case 'l':
#ifdef MATRIX_ENABLE_IMATRIX
	case 'i':
#endif
		return sizeof(int);
	case 'd':
		return sizeof(double);
#ifdef MATRIX_ENABLE_ZMATRIX
	case 'z':
		return sizeof(Rcomplex);
#endif
	default:
		error(_("unexpected kind \"%c\" in %s()"), kind, __func__);
		return 0;
	}
}

const char *Matrix_nonvirtual(SEXP obj, int strict)
{
	if (!IS_S4_OBJECT(obj))
		return "";
	static const char *valid[] = { VALID_NONVIRTUAL, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		return "";
	if (!strict)
		ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
	return valid[ivalid];
}

SEXP R_Matrix_nonvirtual(SEXP obj, SEXP strict)
{
	return mkString(Matrix_nonvirtual(obj, asLogical(strict)));
}

#define RETURN_AS_STRSXP(_C_) \
do { \
	char c = _C_; \
	if (!c) \
		return mkString(""); \
	else { \
		char s[] = { c, '\0' }; \
		return mkString(s); \
	} \
} while (0)

char Matrix_kind(SEXP obj, int i2d)
{
	if (IS_S4_OBJECT(obj)) {
		static const char *valid[] = { VALID_NONVIRTUAL, "" };
		int ivalid = R_check_class_etc(obj, valid);
		if (ivalid < 0)
			return '\0';
		ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
		const char *cl = valid[ivalid];
		return (cl[2] == 'd') ? 'n' : cl[0];
	} else {
		switch (TYPEOF(obj)) {
		case LGLSXP:
			return 'l';
		case INTSXP:
			return (i2d) ? 'd' : 'i';
		case REALSXP:
			return 'd';
		case CPLXSXP:
			return 'z';
		default:
			return '\0';
		}
	}
}

SEXP R_Matrix_kind(SEXP obj, SEXP i2d)
{
	RETURN_AS_STRSXP(Matrix_kind(obj, asLogical(i2d)));
}

char Matrix_shape(SEXP obj)
{
	if (!IS_S4_OBJECT(obj))
		return '\0';
	static const char *valid[] = { VALID_NONVIRTUAL, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		return '\0';
	ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
	const char *cl = valid[ivalid];
	return (cl[2] == 'd' || cl[3] != 'M') ? 'g' : cl[1];
}

SEXP R_Matrix_shape(SEXP obj)
{
	RETURN_AS_STRSXP(Matrix_shape(obj));
}

char Matrix_repr(SEXP obj)
{
	if (!IS_S4_OBJECT(obj))
		return '\0';
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		return '\0';
	ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
	const char *cl = valid[ivalid];
	switch (cl[2]) {
	case 'e':
	case 'r':
	case 'y':
		return 'u'; /* unpackedMatrix */
	case 'p':
		return 'p'; /* packedMatrix */
	case 'C':
	case 'R':
	case 'T':
		return cl[2]; /* [CRT]sparseMatrix */
	case 'i':
		return 'd'; /* diagonalMatrix */
	case 'd':
		return 'i'; /* indMatrix */
	default:
		return '\0';
	}
}

SEXP R_Matrix_repr(SEXP obj)
{
	RETURN_AS_STRSXP(Matrix_repr(obj));
}

#undef RETURN_AS_STRSXP


/* For indexing ===================================================== */

SEXP R_index_triangle(SEXP n, SEXP packed, SEXP upper, SEXP diag)
{
	SEXP r;
	int i, j, n_ = asInteger(n), packed_ = asLogical(packed),
		upper_ = asLogical(upper), diag_ = asLogical(diag);
	Matrix_int_fast64_t
		nn = (Matrix_int_fast64_t) n_ * n_,
		nx = (packed_) ? n_ + (nn - n_) / 2 : nn,
		nr = (diag_) ? n_ + (nn - n_) / 2 : (nn - n_) / 2;
	if (nx > 0x1.0p+53)
		error(_("indices would exceed %s"), "2^53");
	if (nr > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	if (nx > INT_MAX) {

		PROTECT(r = allocVector(REALSXP, (R_xlen_t) nr));
		double k = 1.0, nr_ = (double) nr, *pr = REAL(r);

#define DO_INDEX \
		do { \
			if (packed_) { \
				if (diag_) { \
					while (k <= nr_) \
						*(pr++) = k++; \
				} else if (upper_) { \
					for (j = 0; j < n_; ++j) { \
						for (i = 0; i < j; ++i) \
							*(pr++) = k++; \
						k++; \
					} \
				} else { \
					for (j = 0; j < n_; ++j) { \
						k++; \
						for (i = j+1; i < n_; ++i) \
							*(pr++) = k++; \
					} \
				} \
			} else if (diag_) { \
				if (upper_) { \
					for (j = 0; j < n_; ++j) { \
						for (i = 0; i <= j; ++i) \
							*(pr++) = k++; \
						k += n_-j-1; \
					} \
				} else { \
					for (j = 0; j < n_; ++j) { \
						k += j; \
						for (i = j; i < n_; ++i) \
							*(pr++) = k++; \
					} \
				} \
			} else { \
				if (upper_) { \
					for (j = 0; j < n_; ++j) { \
						for (i = 0; i < j; ++i) \
							*(pr++) = k++; \
						k += n_-j; \
					} \
				} else { \
					for (j = 0; j < n_; ++j) { \
						k += j+1; \
						for (i = j+1; i < n_; ++i) \
							*(pr++) = k++; \
					} \
				} \
			} \
		} while (0)

		DO_INDEX;

	} else {

		PROTECT(r = allocVector(INTSXP, (R_xlen_t) nr));
		int k = 1, nr_ = (int) nr, *pr = INTEGER(r);

		DO_INDEX;

#undef DO_INDEX

	}

	UNPROTECT(1);
	return r;
}

SEXP R_index_diagonal(SEXP n, SEXP packed, SEXP upper)
{
	SEXP r;
	int j, n_ = asInteger(n), packed_ = asLogical(packed),
	upper_ = asLogical(upper);
	Matrix_int_fast64_t
		nn = (Matrix_int_fast64_t) n_ * n_,
		nx = (packed_) ? n_ + (nn - n_) / 2 : nn;
	if (nx > 0x1.0p+53)
		error(_("indices would exceed %s"), "2^53");
	if (nx > INT_MAX) {

		PROTECT(r = allocVector(REALSXP, n_));
		double k = 1.0, *pr = REAL(r);

#define DO_INDEX \
		do { \
			if (!packed_) { \
				for (j = 0; j < n_; ++j) { \
					*(pr++) = k++; \
					k += n_; \
				} \
			} else if (upper_) { \
				for (j = 0; j < n_; ++j) { \
					*(pr++) = k; \
					k += j+2; \
				} \
			} else { \
				for (j = 0; j < n_; ++j) { \
					*(pr++) = k; \
					k += n_-j; \
				} \
			} \
		} while (0)

		DO_INDEX;

	} else {

		PROTECT(r = allocVector(INTSXP, n_));
		int k = 1, *pr = INTEGER(r);
		DO_INDEX;

#undef DO_INDEX

	}

	UNPROTECT(1);
	return r;
}


/* "Miscellaneous" ================================================== */

SEXP R_nnz(SEXP x, SEXP countNA, SEXP nnzmax)
{
	int do_countNA = asLogical(countNA);
	R_xlen_t n = XLENGTH(x), nnz = 0;
	double n_ = asReal(nnzmax);
	if (!ISNAN(n_) && n_ >= 0.0 && n_ < (double) n)
		n = (R_xlen_t) n_;

#define DO_NNZ(_CTYPE_, _PTR_, _NA_, _NZ_, _STRICTLY_NZ_) \
	do { \
		_CTYPE_ *px = _PTR_(x); \
		if (do_countNA == NA_LOGICAL) { \
			while (n-- > 0) { \
				if (_NA_(*px)) \
					return ScalarInteger(NA_INTEGER); \
				if (_NZ_(*px)) \
					++nnz; \
				++px; \
			} \
		} else if (do_countNA != 0) { \
			while (n-- > 0) { \
				if (_NZ_(*px)) \
					++nnz; \
				++px; \
			} \
		} else { \
			while (n-- > 0) { \
				if (_STRICTLY_NZ_(*px)) \
					++nnz; \
				++px; \
			} \
		} \
	} while (0)

	switch (TYPEOF(x)) {
	case LGLSXP:
		DO_NNZ(int, LOGICAL,
		       ISNA_LOGICAL, ISNZ_LOGICAL, STRICTLY_ISNZ_LOGICAL);
		break;
	case INTSXP:
		DO_NNZ(int, INTEGER,
		       ISNA_INTEGER, ISNZ_INTEGER, STRICTLY_ISNZ_INTEGER);
	break;
	case REALSXP:
		DO_NNZ(double, REAL,
		       ISNA_REAL, ISNZ_REAL, STRICTLY_ISNZ_REAL);
	break;
	case CPLXSXP:
		DO_NNZ(Rcomplex, COMPLEX,
		       ISNA_COMPLEX, ISNZ_COMPLEX, STRICTLY_ISNZ_COMPLEX);
	break;
	default:
		ERROR_INVALID_TYPE(x, __func__);
	}

#undef DO_NNZ

	return (nnz <= INT_MAX)
		? ScalarInteger((int) nnz) : ScalarReal((double) nnz);
}

void conjugate(SEXP x)
{
	Rcomplex *px = COMPLEX(x);
	R_xlen_t nx = XLENGTH(x);
	while (nx--) {
		(*px).i = -(*px).i;
		++px;
	}
	return;
}

void zeroRe(SEXP x)
{
	Rcomplex *px = COMPLEX(x);
	R_xlen_t nx = XLENGTH(x);
	while (nx--) {
		if (!ISNAN((*px).r))
			(*px).r = 0.0;
		++px;
	}
	return;
}

void zeroIm(SEXP x)
{
	Rcomplex *px = COMPLEX(x);
	R_xlen_t nx = XLENGTH(x);
	while (nx--) {
		if (!ISNAN((*px).i))
			(*px).i = 0.0;
		++px;
	}
	return;
}

void na2one(SEXP x)
{
	R_xlen_t i, n = XLENGTH(x);
	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *px = LOGICAL(x);
		for (i = 0; i < n; ++i, ++px)
			if (*px == NA_LOGICAL)
				*px = 1;
		break;
	}
	case INTSXP:
	{
		int *px = INTEGER(x);
		for (i = 0; i < n; ++i, ++px)
			if (*px == NA_INTEGER)
				*px = 1;
		break;
	}
	case REALSXP:
	{
		double *px = REAL(x);
		for (i = 0; i < n; ++i, ++px)
			if (ISNAN(*px))
				*px = 1.0;
		break;
	}
	case CPLXSXP:
	{
		Rcomplex *px = COMPLEX(x);
		for (i = 0; i < n; ++i, ++px)
			if (ISNAN((*px).r) || ISNAN((*px).i))
				*px = Matrix_zone;
		break;
	}
	default:
		ERROR_INVALID_TYPE(x, __func__);
		break;
	}
	return;
}

/* That both 's1' and 's2' are STRSXP of length at least 'n' must be
   checked by the caller ... see, e.g., symmetricMatrix_validate() above
*/
Rboolean equal_string_vectors(SEXP s1, SEXP s2, int n)
{
	/* Only check the first 'n' elements, even if 's1' or 's2' is longer ...

	   Note that 'R_compute_identical()' in src/main/identical.c
	   is careful to distinguish between NA_STRING and "NA" in STRSXP,
	   but we need not be here ...

	   MJ: Why not?
	*/

	for (int i = 0; i < n; ++i)
		if (strcmp(CHAR(STRING_ELT(s1, i)), CHAR(STRING_ELT(s2, i))) != 0)
			return FALSE;
	return TRUE;
}

SEXP append_to_named_list(SEXP x, const char *nm, SEXP val)
{
	PROTECT(val);
	R_xlen_t n = XLENGTH(x);
	SEXP y = PROTECT(allocVector(VECSXP, n + 1)),
	ny = PROTECT(allocVector(STRSXP, n + 1)),
	nval = PROTECT(mkChar(nm));
	if (n > 0) {
		SEXP nx = PROTECT(getAttrib(x, R_NamesSymbol));
		R_xlen_t i;
		for (i = 0; i < n; ++i) {
			SET_VECTOR_ELT( y, i, VECTOR_ELT( x, i));
			SET_STRING_ELT(ny, i, STRING_ELT(nx, i));
		}
		UNPROTECT(1);
	}
	SET_VECTOR_ELT( y, n,  val);
	SET_STRING_ELT(ny, n, nval);
	setAttrib(y, R_NamesSymbol, ny);
	UNPROTECT(4);
	return y;
}


/* ================================================================== */
/* ================================================================== */

SEXP Matrix_expand_pointers(SEXP pP)
{
	int n = length(pP) - 1;
	int *p = INTEGER(pP);
	SEXP ans = PROTECT(allocVector(INTSXP, p[n]));

	expand_cmprPt(n, p, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param ij: 2-column integer matrix
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd(SEXP ij, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
	SEXP ans;
	int *ij_di = NULL, n, nprot=1;
	Rboolean check_bounds = asLogical(chk_bnds), one_ind = asLogical(orig_1);

	if (TYPEOF(di) != INTSXP) {
		di = PROTECT(coerceVector(di, INTSXP));
		nprot++;
	}
	if (TYPEOF(ij) != INTSXP) {
		ij = PROTECT(coerceVector(ij, INTSXP));
		nprot++;
	}
	if (!isMatrix(ij) ||
	    (ij_di = INTEGER(getAttrib(ij, R_DimSymbol)))[1] != 2)
		error(_("Argument ij must be 2-column integer matrix"));
	n = ij_di[0];
	int *Di = INTEGER(di), *IJ = INTEGER(ij),
		*j_ = IJ+n;/* pointer offset! */

	if ((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
		ans = PROTECT(allocVector(REALSXP, n));
		double *ii = REAL(ans), nr = (double) Di[0];

#define do_ii_FILL(_i_, _j_) \
		int i; \
		if (check_bounds) { \
			for (i = 0; i < n; i++) { \
				if (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER) \
					ii[i] = NA_INTEGER; \
				else { \
					register int i_i, j_i; \
					if (one_ind) { \
						i_i = _i_[i]-1; \
						j_i = _j_[i]-1; \
					} else { \
						i_i = _i_[i]; \
						j_i = _j_[i]; \
					} \
					if (i_i < 0 || i_i >= Di[0]) \
						error(_("subscript 'i' out of bounds in M[ij]")); \
					if (j_i < 0 || j_i >= Di[1]) \
						error(_("subscript 'j' out of bounds in M[ij]")); \
					ii[i] = i_i + j_i * nr; \
				} \
			} \
		} else { \
			for (i = 0; i < n; i++) \
				ii[i] = (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER) \
					? NA_INTEGER \
					: ((one_ind) \
					   ? ((_i_[i]-1) + (_j_[i]-1) * nr) \
					   :   _i_[i] + _j_[i] * nr); \
		}

		do_ii_FILL(IJ, j_);
	} else {
	ans = PROTECT(allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];

	do_ii_FILL(IJ, j_);
	}
	UNPROTECT(nprot);
	return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param i: integer vector
 * @param j: integer vector of same length as 'i'
 * @param orig_1: logical: if TRUE, "1-origin" otherwise "0-origin"
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
	SEXP ans;
	int n = LENGTH(i), nprot = 1;
	Rboolean check_bounds = asLogical(chk_bnds), one_ind = asLogical(orig_1);

	if (TYPEOF(di)!= INTSXP) {
		di = PROTECT(coerceVector(di,INTSXP));
		nprot++;
	}
	if (TYPEOF(i) != INTSXP) {
		i = PROTECT(coerceVector(i, INTSXP));
		nprot++;
	}
	if (TYPEOF(j) != INTSXP) {
		j = PROTECT(coerceVector(j, INTSXP));
		nprot++;
	}
	if (LENGTH(j) != n)
		error(_("i and j must be integer vectors of the same length"));

	int *Di = INTEGER(di), *i_ = INTEGER(i), *j_ = INTEGER(j);

	if ((Di[0] * (double) Di[1]) >= 1 + (double) INT_MAX) { /* need double */
		ans = PROTECT(allocVector(REALSXP, n));
		double *ii = REAL(ans), nr = (double) Di[0];

		do_ii_FILL(i_, j_);
	} else {
		ans = PROTECT(allocVector(INTSXP, n));
		int *ii = INTEGER(ans), nr = Di[0];

		do_ii_FILL(i_, j_);
	}
	UNPROTECT(nprot);
	return ans;
}
#undef do_ii_FILL

// Almost "Cut n Paste" from ...R../src/main/array.c  do_matrix() :
// used in ../R/Matrix.R as
//
// .External(Mmatrix,
//		 data, nrow, ncol, byrow, dimnames,
//		 missing(nrow), missing(ncol))
SEXP Mmatrix(SEXP args)
{
	SEXP vals, ans, snr, snc, dimnames;
	int nr = 1, nc = 1, byrow, miss_nr, miss_nc;
	R_xlen_t lendat;

	args = CDR(args); /* skip 'name' */
	vals = CAR(args); args = CDR(args);
	/* Supposedly as.vector() gave a vector type, but we check */
	switch (TYPEOF(vals)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
	case STRSXP:
	case RAWSXP:
	case EXPRSXP:
	case VECSXP:
		break;
	default:
		error(_("'data' must be of a vector type"));
	}
	lendat = XLENGTH(vals);
	snr = CAR(args); args = CDR(args);
	snc = CAR(args); args = CDR(args);
	byrow = asLogical(CAR(args)); args = CDR(args);
	if (byrow == NA_INTEGER)
		error(_("invalid '%s' argument"), "byrow");
	dimnames = CAR(args);
	args = CDR(args);
	miss_nr = asLogical(CAR(args)); args = CDR(args);
	miss_nc = asLogical(CAR(args));

	if (!miss_nr) {
		if (!isNumeric(snr)) error(_("non-numeric matrix extent"));
		nr = asInteger(snr);
		if (nr == NA_INTEGER)
			error(_("invalid 'nrow' value (too large or NA)"));
		if (nr < 0)
			error(_("invalid 'nrow' value (< 0)"));
	}
	if (!miss_nc) {
		if (!isNumeric(snc)) error(_("non-numeric matrix extent"));
		nc = asInteger(snc);
		if (nc == NA_INTEGER)
			error(_("invalid 'ncol' value (too large or NA)"));
		if (nc < 0)
			error(_("invalid 'ncol' value (< 0)"));
	}
	if (miss_nr && miss_nc) {
		if (lendat > INT_MAX) error("data is too long");
		nr = (int) lendat;
	} else if (miss_nr) {
		if (lendat > (double) nc * INT_MAX) error("data is too long");
		nr = (int) ceil((double) lendat / (double) nc);
	} else if (miss_nc) {
		if (lendat > (double) nr * INT_MAX) error("data is too long");
		nc = (int) ceil((double) lendat / (double) nr);
	}

	if (lendat > 0) {
		R_xlen_t nrc = (R_xlen_t) nr * nc;
		if (lendat > 1 && nrc % lendat != 0) {
			if (((lendat > nr) && (lendat / nr) * nr != lendat) ||
				((lendat < nr) && (nr / lendat) * lendat != nr))
				warning(_("data length [%d] is not a sub-multiple "
				          "or multiple of the number of rows [%d]"),
				        lendat, nr);
			else if (((lendat > nc) && (lendat / nc) * nc != lendat) ||
					 ((lendat < nc) && (nc / lendat) * lendat != nc))
				warning(_("data length [%d] is not a sub-multiple "
				          "or multiple of the number of columns [%d]"),
						lendat, nc);
		} else if ((lendat > 1) && (nrc == 0))
			warning(_("data length exceeds size of matrix"));
	}

#ifndef LONG_VECTOR_SUPPORT
	if ((double) nr * (double) nc > INT_MAX)
		error(_("too many elements specified"));
#endif

	PROTECT(ans = allocMatrix(TYPEOF(vals), nr, nc));
	if (lendat) {
		if (isVector(vals))
			copyMatrix(ans, vals, byrow);
		else
			copyListMatrix(ans, vals, byrow);
	} else if (isVector(vals)) { /* fill with NAs */
		R_xlen_t N = (R_xlen_t) nr * nc, i;
		switch (TYPEOF(vals)) {
		case STRSXP:
			for (i = 0; i < N; i++)
				SET_STRING_ELT(ans, i, NA_STRING);
			break;
		case LGLSXP:
			for (i = 0; i < N; i++)
				LOGICAL(ans)[i] = NA_LOGICAL;
			break;
		case INTSXP:
			for (i = 0; i < N; i++)
				INTEGER(ans)[i] = NA_INTEGER;
			break;
		case REALSXP:
			for (i = 0; i < N; i++)
				REAL(ans)[i] = NA_REAL;
			break;
		case CPLXSXP:
		{
			/* Initialization must work whether Rcomplex is typedef-ed
			   to a struct { R < 4.3.0 } or to a union { R >= 4.3.0 }
			*/
			Rcomplex zna = { .r = NA_REAL, .i = 0.0 };
			for (i = 0; i < N; i++)
				COMPLEX(ans)[i] = zna;
			break;
		}
		case RAWSXP:
			// FIXME:  N may overflow size_t !!
			memset(RAW(ans), 0, N);
			break;
		default:
			/* don't fill with anything */
			;
		}
	}
	if (!isNull(dimnames)&& length(dimnames) > 0)
		ans = dimnamesgets(ans, dimnames);
	UNPROTECT(1);
	return ans;
}

/**
 * From the two 'x' slots of two dense matrices a and b,
 * compute the 'x' slot of rbind(a, b)
 *
 * Currently, an auxiliary only for setMethod rbind2(<denseMatrix>, <denseMatrix>)
 * in ../R/bind2.R
 *
 * @param a
 * @param b
 *
 * @return
 */
SEXP R_rbind2_vector(SEXP a, SEXP b)
{
	int *d_a = INTEGER(GET_SLOT(a, Matrix_DimSym)),
		*d_b = INTEGER(GET_SLOT(b, Matrix_DimSym)),
		n1 = d_a[0], m = d_a[1],
		n2 = d_b[0];
	if (d_b[1] != m)
		error(_("the number of columns differ in R_rbind2_vector: %d != %d"),
		      m, d_b[1]);
	SEXP
		a_x = GET_SLOT(a, Matrix_xSym),
		b_x = GET_SLOT(b, Matrix_xSym);
	int nprot = 1;
	// Care: can have "ddenseMatrix" "ldenseMatrix" or "ndenseMatrix"
	if (TYPEOF(a_x) != TYPEOF(b_x)) { // choose the "common type"
		// Now know: either LGLSXP or REALSXP. FIXME for iMatrix, zMatrix,..
		if (TYPEOF(a_x) != REALSXP) {
			a_x = PROTECT(duplicate(coerceVector(a_x, REALSXP)));
			nprot++;
		} else if (TYPEOF(b_x) != REALSXP) {
			b_x = PROTECT(duplicate(coerceVector(b_x, REALSXP)));
			nprot++;
		}
	}

	SEXP ans = PROTECT(allocVector(TYPEOF(a_x), m * (n1 + n2)));
	int ii = 0;
	switch (TYPEOF(a_x)) {
	case LGLSXP:
	{
		int
			*r = LOGICAL(ans),
			*ax= LOGICAL(a_x),
			*bx= LOGICAL(b_x);

#define COPY_a_AND_b_j \
		for (int j = 0; j < m; j++) { \
			Memcpy(r+ii, ax+ j*n1, n1); ii += n1; \
			Memcpy(r+ii, bx+ j*n2, n2); ii += n2; \
		}; break

		COPY_a_AND_b_j;
	}
	case REALSXP: {
		double
			*r = REAL(ans),
			*ax= REAL(a_x),
			*bx= REAL(b_x);

		COPY_a_AND_b_j;
	}
	} // switch
	UNPROTECT(nprot);
	return ans;
}

#define TRUE_  ScalarLogical(1)
#define FALSE_ ScalarLogical(0)

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// all0	 <- function(x) !any(is.na(x)) && all(!x) ## ~= allFalse
// allFalse <- function(x) !any(x) && !any(is.na(x)) ## ~= all0
SEXP R_all0(SEXP x) {
	if (!isVectorAtomic(x)) {
		if (length(x) == 0) return TRUE_;
		// Typically S4.  TODO: Call the R code above, instead!
		error(_("Argument must be numeric-like atomic vector"));
	}
	R_xlen_t i, n = XLENGTH(x);
	if (n == 0) return TRUE_;

	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *xx = LOGICAL(x);
		for (i = 0; i < n; i++)
			if (xx[i] == NA_LOGICAL || xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	case INTSXP:
	{
		int *xx = INTEGER(x);
		for (i = 0; i < n; i++)
			if (xx[i] == NA_INTEGER || xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	case REALSXP:
	{
		double *xx = REAL(x);
		for (i = 0; i < n; i++)
			if (ISNAN(xx[i]) || xx[i] != 0.) return FALSE_;
		return TRUE_;
	}
	case RAWSXP:
	{
		unsigned char *xx = RAW(x);
		for (i = 0; i < n; i++)
			if (xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	}
	error(_("Argument must be numeric-like atomic vector"));
	return R_NilValue; // -Wall
}

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// any0 <- function(x) isTRUE(any(x == 0)) ## ~= anyFalse
// anyFalse <- function(x) isTRUE(any(!x)) ## ~= any0
SEXP R_any0(SEXP x) {
	if (!isVectorAtomic(x)) {
		if (length(x) == 0) return FALSE_;
		// Typically S4.  TODO: Call the R code above, instead!
		error(_("Argument must be numeric-like atomic vector"));
	}
	R_xlen_t i, n = XLENGTH(x);
	if (n == 0) return FALSE_;

	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *xx = LOGICAL(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	case INTSXP:
	{
		int *xx = INTEGER(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	case REALSXP:
	{
		double *xx = REAL(x);
		for (i = 0; i < n; i++) if (xx[i] == 0.) return TRUE_;
		return FALSE_;
	}
	case RAWSXP:
	{
		unsigned char *xx = RAW(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	}
	error(_("Argument must be numeric-like atomic vector"));
	return R_NilValue; // -Wall
}

#undef TRUE_
#undef FALSE_
