#include "validity.h"

/* Slot validity methods ===============================================
   Called by various class validity methods (see below).
*/

/**
 * Test that `dim` is a length-2, non-negative integer vector.
 *
 * @param dim A `SEXP`,
 *     typically the `Dim` slot of a (to be validated) `Matrix`.
 *
 * @return A string containing an error message, empty if `dim`
 *     is valid.
 */
char *Dim_validate(SEXP dim)
{
	/* TODO? coerce from REALSXP to INTSXP?
	   // if (TYPEOF(dim) != INTSXP && TYPEOF(dim) != REALSXP)
	   //     return _("'Dim' slot is not numeric");
	   though above is not enough as we must prohibit Dim[i] > INT_MAX

	   FIXME? Prohibit is.object(dim) or maybe just inherits(dim, "factor")
	   and return a different error message in that case?
	*/
	if (TYPEOF(dim) != INTSXP)
		return _("'Dim' slot is not of type \"integer\"");
	if (XLENGTH(dim) != 2)
		return _("'Dim' slot does not have length 2");
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m == NA_INTEGER || n == NA_INTEGER)
		return _("'Dim' slot contains NA");
	if (m < 0 || n < 0)
		return dngettext(Matrix_Domain,
		                 "'Dim' slot contains negative value",
		                 "'Dim' slot contains negative values",
		                 (m < 0 && n < 0) ? 2 : 1);
	return "";
}

SEXP R_Dim_validate(SEXP dim)
{
	char *msg = Dim_validate(dim);
	return (msg[0] == '\0') ? ScalarLogical(1) : mkString(msg);
}

/**
 * Test that `dimnames` is a valid length-2 list.
 *
 * @param dimnames A `SEXP`,
 *     typically the `Dimnames` slot of a (to be validated) `Matrix`.
 * @param pdim Pointer to a length-2, non-negative `int` array,
 *     typically from the `Dim` slot of a (to be validated) `Matrix`.
 *     Array validity _must_ be checked by the caller.
 *
 * @return A string containing an error message, empty if `dimnames`
 *     is valid.
 */
char *DimNames_validate(SEXP dimnames, int *pdim)
{
	if (TYPEOF(dimnames) != VECSXP)
		return _("'Dimnames' slot is not a list");
	if (XLENGTH(dimnames) != 2)
		return _("'Dimnames' slot does not have length 2");

	/* Behave as do_matrix() from src/main/array.c:
	   Dimnames[[i]] must be NULL or _coercible to_ character
	   of length Dim[i] or 0 ... see R_Dimnames_fixup() below */
	for (int i = 0; i < 2; ++i) {
		SEXP s = VECTOR_ELT(dimnames, i);
		if (!isNull(s)) {
			if (!isVector(s)) {
				char *buf;
				SNPRINTF(buf, _("Dimnames[[%d]] is not NULL or a vector"),
				         i+1);
				return buf;
			}
			R_xlen_t ns = XLENGTH(s);
			if (ns != pdim[i] && ns != 0) {
				char *buf;
				SNPRINTF(buf, _("length of Dimnames[[%d]] (%lld) is not equal to Dim[%d] (%d)"),
				         i+1, (long long) ns, i+1, pdim[i]);
				return buf;
			}
		}
	}
	return "";
}

SEXP R_DimNames_validate(SEXP dimnames, SEXP dim)
{
	char *msg = DimNames_validate(dimnames, INTEGER(dim));
	return (msg[0] == '\0') ? ScalarLogical(1) : mkString(msg);
}

/**
 * @brief Sanitize user-supplied `[dD]imnames`.
 *
 * Replaces length-0 vectors with `NULL` and non-character vectors
 * with the result of coercing to character. Intended to emulate the
 * behaviour of `do_matrix()` from `src/main/array.c`.
 *
 * @param dn A list of length 2 passing `DimNames_validate()`.
 *
 * @return A modified copy of `dn`, or `dn` if no modification is
 *    necessary.
 */
SEXP R_DimNames_fixup(SEXP dn)
{
	SEXP s;
	int i;
	Rboolean do_fixup = FALSE;
	for (i = 0; i < 2; ++i) {
		if (!isNull(s = VECTOR_ELT(dn, i)) &&
			(LENGTH(s) == 0 || TYPEOF(s) != STRSXP)) {
			do_fixup = TRUE;
			break;
		}
	}
	if (do_fixup) {
		PROTECT(dn = duplicate(dn));
		for (i = 0; i < 2; ++i) {
			if (isNull(s = VECTOR_ELT(dn, i)))
				continue;
			if (LENGTH(s) == 0)
				SET_VECTOR_ELT(dn, i, R_NilValue);
			else if (TYPEOF(s) != STRSXP) {
				if (inherits(s, "factor"))
					PROTECT(s = asCharacterFactor(s));
				else {
					PROTECT(s = coerceVector(s, STRSXP));
					SET_ATTRIB(s, R_NilValue);
					SET_OBJECT(s, 0);
				}
				SET_VECTOR_ELT(dn, i, s);
				UNPROTECT(1); /* s */
			}
		}
		UNPROTECT(1); /* dn */
	}
	return dn;
}


/* Class validity methods ==============================================
   NB: These assume that validity methods for superclasses
   have already been called via validObject() ...
*/

SEXP Matrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	char *msg = Dim_validate(dim);
	if (msg[0] == '\0') {
		SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
		msg = DimNames_validate(dimnames, INTEGER(dim));
		UNPROTECT(1); /* dimnames */
	}
	UNPROTECT(1); /* dim */
	return (msg[0] == '\0') ? ScalarLogical(1) : mkString(msg);
}

SEXP MatrixFactorization_validate(SEXP obj)
{
	return Matrix_validate(obj);
}

#define TYPEMATRIX_VALIDATE(_PREFIX_, _SEXPTYPE_, _T2C_SEXPTYPE_) \
SEXP _PREFIX_ ## Matrix_validate(SEXP obj) \
{ \
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)); \
	if (TYPEOF(x) != _SEXPTYPE_) \
		UPRET(1, "'x' slot is not of type \"" #_T2C_SEXPTYPE_ "\""); \
	UNPROTECT(1); /* x */ \
	return ScalarLogical(1); \
}
/* dMatrix_validate() */
TYPEMATRIX_VALIDATE(     d, REALSXP,  double)
/* lMatrix_validate() */
TYPEMATRIX_VALIDATE(     l,  LGLSXP, logical)
/* ndenseMatrix_validate() */
/* NB: "nsparseMatrix" has no 'x' slot, only "ndenseMatrix" ... */
TYPEMATRIX_VALIDATE(ndense,  LGLSXP, logical)
/* iMatrix_validate() */
TYPEMATRIX_VALIDATE(     i,  INTSXP, integer)
/* zMatrix_validate() */
TYPEMATRIX_VALIDATE(     z, CPLXSXP, complex)
#undef TYPEMATRIX_VALIDATE

SEXP compMatrix_validate(SEXP obj)
{
	SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorSym));
	if (TYPEOF(factors) != VECSXP)
		UPRET(1, "'factors' slot is not a list");
	if (XLENGTH(factors) > 0) {
		SEXP nms = PROTECT(getAttrib(factors, R_NamesSymbol));
		if (isNull(nms))
			UPRET(2, "'factors' slot has no 'names' attribute");
		UNPROTECT(1); /* nms */
	}
	UNPROTECT(1); /* factors */
	return ScalarLogical(1);
}

SEXP symmetricMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		UPRET(1, "Dim[1] != Dim[2] (matrix is not square)");
	UNPROTECT(1); /* dim */

#ifdef ENFORCE_SYMMETRIC_DIMNAMES
	/* This check can be expensive when both rownames and colnames have
	   nonzero length, and even more so when coercions to character are
	   required ... Users can avoid the expense by setting at least one
	   of rownames and colnames to NULL or by ensuring that they are the
	   same object, as testing for pointer equality is fast ... */

# define ANY_TO_STRING(x) \
	(TYPEOF(x) == STRSXP \
	 ? x \
	 : (inherits(x, "factor") \
	    ? asCharacterFactor(x) \
	    : coerceVector(x, STRSXP)))

	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	ndn = PROTECT(getAttrib(dn, R_NamesSymbol));
	const char *ndn0, *ndn1;
	if (!isNull(ndn) &&
		*(ndn0 = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
		*(ndn1 = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
		strcmp(ndn0, ndn1) != 0)
		UPRET(2, "Dimnames[1] differs from Dimnames[2]");
	if (n > 0) {
		/* NB: It is already known that the length of 'dn[[i]]' is 0 or 'n' */
		SEXP rn, cn;
		if (!isNull(rn = VECTOR_ELT(dn, 0)) &&
			!isNull(cn = VECTOR_ELT(dn, 1)) &&
			LENGTH(rn) == n && LENGTH(cn) == n && rn != cn) {
			PROTECT(rn = ANY_TO_STRING(rn));
			PROTECT(cn = ANY_TO_STRING(cn));
			if (!equal_string_vectors(rn, cn, n))
				UPRET(4, "Dimnames[1] differs from Dimnames[2]");
			UNPROTECT(2); /* cn, rn */
		}
	}
	UNPROTECT(2); /* ndn, dn */

# undef ANY_TO_STRING
#endif

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	if (TYPEOF(uplo) != STRSXP)
		UPRET(1, "'uplo' slot is not of type \"character\"");
	if (XLENGTH(uplo) != 1)
		UPRET(1, "'uplo' slot does not have length 1");
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
		UPRET(1, "'uplo' slot is not \"U\" or \"L\"");
	UNPROTECT(1); /* uplo */

	return ScalarLogical(1);
}

SEXP triangularMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		UPRET(1, "Dim[1] != Dim[2] (matrix is not square)");
	UNPROTECT(1); /* dim */

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	if (TYPEOF(uplo) != STRSXP)
		UPRET(1, "'uplo' slot is not of type \"character\"");
	if (XLENGTH(uplo) != 1)
		UPRET(1, "'uplo' slot does not have length 1");
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
		UPRET(1, "'uplo' slot is not \"U\" or \"L\"");
	UNPROTECT(1); /* uplo */

	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	if (TYPEOF(diag) != STRSXP)
		UPRET(1, "'diag' slot is not of type \"character\"");
	if (XLENGTH(diag) != 1)
		UPRET(1, "'diag' slot does not have length 1");
	const char *di = CHAR(STRING_ELT(diag, 0));
	if (di[0] == '\0' || di[1] != '\0' || (di[0] != 'N' && di[0] != 'U'))
		UPRET(1, "'diag' slot is not \"N\" or \"U\"");
	UNPROTECT(1); /* diag */

	return ScalarLogical(1);
}

SEXP diagonalMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		UPRET(1, "Dim[1] != Dim[2] (matrix is not square)");
	UNPROTECT(1); /* dim */

	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	if (TYPEOF(diag) != STRSXP)
		UPRET(1, "'diag' slot is not of type \"character\"");
	if (XLENGTH(diag) != 1)
		UPRET(1, "'diag' slot does not have length 1");
	const char *di = CHAR(STRING_ELT(diag, 0));
	if (di[0] == '\0' || di[1] != '\0' || (di[0] != 'N' && di[0] != 'U'))
		UPRET(1, "'diag' slot is not \"N\" or \"U\"");
	int nonunit = di[0] == 'N';
	UNPROTECT(1); /* diag */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	if (nonunit) {
		if (XLENGTH(x) != n)
			UPRET(1, "'diag' slot is \"N\" but 'x' slot "
				  "does not have length n=Dim[1]");
	} else {
		if (XLENGTH(x) != 0)
			UPRET(1, "'diag' slot is \"U\" (identity matrix) but 'x' slot "
				  "does not have length 0");
	}
	UNPROTECT(1); /* x */

	return ScalarLogical(1);
}

SEXP indMatrix_validate(SEXP obj)
{
	SEXP margin = PROTECT(GET_SLOT(obj, Matrix_marginSym));
	if (XLENGTH(margin) != 1)
		UPRET(1, "'margin' slot does not have length 1");
	int mg = INTEGER(margin)[0] - 1;
	if (mg != 0 && mg != 1)
		UPRET(1, "'margin' slot is not 1 or 2");
	UNPROTECT(1); /* margin */

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[mg], n = pdim[!mg];
	if (m > 0 && n == 0)
		UPRET(1, (mg == 0) ? "m-by-0 indMatrix invalid for positive 'm' when margin=1" : "0-by-n indMatrix invalid for positive 'n' when margin=2");
	UNPROTECT(1); /* dim */

	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	if (TYPEOF(perm) != INTSXP)
		UPRET(1, "'perm' slot is not of type \"integer\"");
	if (XLENGTH(perm) != m)
		UPRET(1, "'perm' slot does not have length Dim[margin]");
	int *pperm = INTEGER(perm);
	while (m--) {
		if (*pperm == NA_INTEGER)
			UPRET(1, "'perm' slot contains NA");
		if (*pperm < 1 || *pperm > n)
			UPRET(1, "'perm' slot has elements not in {1,...,Dim[1+margin%%2]}");
		++pperm;
	}
	UNPROTECT(1); /* perm */

	return ScalarLogical(1);
}

SEXP pMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		UPRET(1, "Dim[1] != Dim[2] (matrix is not square)");
	UNPROTECT(1); /* dim */

	if (n > 1) {
		SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
		int i, *pperm = INTEGER(perm);
		char *work;
		Matrix_Calloc(work, n, char);
		--work;
		for (i = 0; i < n; ++i) {
			if (work[*pperm])
				break;
			work[*(pperm++)] = 1;
		}
		++work;
		Matrix_Free(work, n);
		UNPROTECT(1); /* perm */
		if (i < n)
			return mkString(_("'perm' slot contains duplicates"));
	}

	return ScalarLogical(1);
}

SEXP CsparseMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	if (TYPEOF(p) != INTSXP)
		UPRET(1, "'p' slot is not of type \"integer\"");
	if (XLENGTH(p) - 1 != n)
		UPRET(1, "'p' slot does not have length Dim[2]+1");
	int *pp = INTEGER(p);
	if (pp[0] != 0)
		UPRET(1, "first element of 'p' slot is not 0");
	int j;
	for (j = 1; j <= n; ++j) {
		if (pp[j] == NA_INTEGER)
			UPRET(1, "'p' slot contains NA");
		if (pp[j] < pp[j-1])
			UPRET(1, "'p' slot is not nondecreasing");
		if (pp[j] - pp[j-1] > m)
			UPRET(1, "first differences of 'p' slot exceed Dim[1]");
	}

	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	if (TYPEOF(i) != INTSXP)
		UPRET(2, "'i' slot is not of type \"integer\"");
	if (XLENGTH(i) < pp[n])
		UPRET(2, "'i' slot has length less than p[length(p)]");
	int *pi = INTEGER(i), k = 0, kend, ik, i0;
	for (j = 1; j <= n; ++j) {
		kend = pp[j];
		i0 = -1;
		while (k < kend) {
			ik = pi[k];
			if (ik == NA_INTEGER)
				UPRET(2, "'i' slot contains NA");
			if (ik < 0 || ik >= m)
				UPRET(2, "'i' slot has elements not in {0,...,Dim[1]-1}");
			if (ik <= i0)
				UPRET(2, "'i' slot is not increasing within columns");
			i0 = ik;
			++k;
		}
	}
	UNPROTECT(2); /* i, p */

	return ScalarLogical(1);
}

SEXP RsparseMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	if (TYPEOF(p) != INTSXP)
		UPRET(1, "'p' slot is not of type \"integer\"");
	if (XLENGTH(p) - 1 != m)
		UPRET(1, "'p' slot does not have length Dim[1]+1");
	int *pp = INTEGER(p);
	if (pp[0] != 0)
		UPRET(1, "first element of 'p' slot is not 0");
	int i;
	for (i = 1; i <= m; ++i) {
		if (pp[i] == NA_INTEGER)
			UPRET(1, "'p' slot contains NA");
		if (pp[i] < pp[i-1])
			UPRET(1, "'p' slot is not nondecreasing");
		if (pp[i] - pp[i-1] > n)
			UPRET(1, "first differences of 'p' slot exceed Dim[2]");
	}

	SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	if (TYPEOF(j) != INTSXP)
		UPRET(2, "'j' slot is not of type \"integer\"");
	if (XLENGTH(j) < pp[m])
		UPRET(2, "'j' slot has length less than p[length(p)]");
	int *pj = INTEGER(j), k = 0, kend, jk, j0;
	for (i = 1; i <= m; ++i) {
		kend = pp[i];
		j0 = -1;
		while (k < kend) {
			jk = pj[k];
			if (jk == NA_INTEGER)
				UPRET(2, "'j' slot contains NA");
			if (jk < 0 || jk >= n)
				UPRET(2, "'j' slot has elements not in {0,...,Dim[2]-1}");
			if (jk <= j0)
				UPRET(2, "'j' slot is not increasing within rows");
			j0 = jk;
			++k;
		}
	}
	UNPROTECT(2); /* j, p */

	return ScalarLogical(1);
}

SEXP TsparseMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */

	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
	j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	if (TYPEOF(i) != INTSXP)
		UPRET(2, "'i' slot is not of type \"integer\"");
	if (TYPEOF(j) != INTSXP)
	UPRET(2, "'j' slot is not of type \"integer\"");
	R_xlen_t nnz;
	if ((nnz = XLENGTH(i)) != XLENGTH(j))
	UPRET(2, "'i' and 'j' slots do not have equal length");
	if (nnz > 0) {
		if (m == 0 || n == 0)
			UPRET(2, "'i' slot has nonzero length but prod(Dim) is 0");
		int *pi = INTEGER(i), *pj = INTEGER(j);
		while (nnz--) {
			if (*pi == NA_LOGICAL)
		UPRET(2, "'i' slot contains NA");
			if (*pj == NA_LOGICAL)
				UPRET(2, "'j' slot contains NA");
			if (*pi < 0 || *pi >= m)
				UPRET(2, "'i' slot has elements not in {0,...,Dim[1]-1}");
			if (*pj < 0 || *pj >= n)
				UPRET(2, "'j' slot has elements not in {0,...,Dim[2]-1}");
			++pi;
			++pj;
		}
	}
	UNPROTECT(2); /* j, i */

	return ScalarLogical(1);
}

SEXP sCMatrix_validate(SEXP obj)
{
	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
	if (pp[n] > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
		int *pi = INTEGER(i), j, k = 0, kend;
		if (ul == 'U') {
			for (j = 0; j < n; ++j) {
				kend = *(++pp);
				while (k < kend) {
					if (pi[k] > j)
						UPRET(2, "uplo=\"U\" but there are entries below the diagonal");
					++k;
				}
			}
		} else {
			for (j = 0; j < n; ++j) {
				kend = *(++pp);
				while (k < kend) {
					if (pi[k] < j)
						UPRET(2, "uplo=\"L\" but there are entries above the diagonal");
					++k;
				}
			}
		}
		UNPROTECT(1); /* i */
	}
	UNPROTECT(1); /* p */

	return ScalarLogical(1);
}

SEXP tCMatrix_validate(SEXP obj)
{
	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	if (di == 'N')
		return sCMatrix_validate(obj);

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
	if (pp[n] > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
		int *pi = INTEGER(i), j, k = 0, kend;
		if (ul == 'U') {
			for (j = 0; j < n; ++j) {
				kend = *(++pp);
				while (k < kend) {
					if (pi[k] >= j) {
						int s = pi[k] == j;
						UPRET(2, (s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal");
					}
					++k;
				}
			}
		} else {
			for (j = 0; j < n; ++j) {
				kend = *(++pp);
				while (k < kend) {
					if (pi[k] <= j) {
						int s = pi[k] == j;
						UPRET(2, (s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal");
					}
					++k;
				}
			}
		}
		UNPROTECT(1); /* i */
	}
	UNPROTECT(1); /* p */

	return ScalarLogical(1);
}

SEXP sRMatrix_validate(SEXP obj)
{
	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
	if (pp[m] > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pj = INTEGER(j), i, k = 0, kend;
		if (ul == 'U') {
			for (i = 0; i < m; ++i) {
				kend = *(++pp);
				while (k < kend) {
					if (pj[k] < i)
						UPRET(2, "uplo=\"U\" but there are entries below the diagonal");
					++k;
				}
			}
		} else {
			for (i = 0; i < m; ++i) {
				kend = *(++pp);
				while (k < kend) {
					if (pj[k] > i)
						UPRET(2, "uplo=\"L\" but there are entries above the diagonal");
					++k;
				}
			}
		}
		UNPROTECT(1); /* j */
	}
	UNPROTECT(1); /* p */

	return ScalarLogical(1);
}

SEXP tRMatrix_validate(SEXP obj)
{
	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	if (di == 'N')
		return sRMatrix_validate(obj);

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
	if (pp[m] > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pj = INTEGER(j), i, k = 0, kend;
		if (ul == 'U') {
			for (i = 0; i < m; ++i) {
				kend = *(++pp);
				while (k < kend) {
					if (pj[k] <= i) {
						int s = pj[k] == i;
						UPRET(2, (s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal");
					}
					++k;
				}
			}
		} else {
			for (i = 0; i < m; ++i) {
				kend = *(++pp);
				while (k < kend) {
					if (pj[k] >= i) {
						int s = pj[k] == i;
						UPRET(2, (s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal");
					}
					++k;
				}
			}
		}
		UNPROTECT(1); /* j */
	}
	UNPROTECT(1); /* p */

	return ScalarLogical(1);
}

SEXP sTMatrix_validate(SEXP obj)
{
	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	R_xlen_t nnz = XLENGTH(i);
	if (nnz > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		if (ul == 'U') {
			while (nnz--)
				if (*(pi++) > *(pj++))
					UPRET(2, "uplo=\"U\" but there are entries below the diagonal");
		} else {
			while (nnz--)
				if (*(pi++) < *(pj++))
					UPRET(2, "uplo=\"L\" but there are entries above the diagonal");
		}
		UNPROTECT(1); /* j */
	}
	UNPROTECT(1); /* i */

	return ScalarLogical(1);
}

SEXP tTMatrix_validate(SEXP obj)
{
	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	if (di == 'N')
		return sTMatrix_validate(obj);

	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	R_xlen_t nnz = XLENGTH(i);
	if (nnz > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		if (ul == 'U') {
			while (nnz--) {
				if (*pi >= *pj) {
					int s = *pi == *pj;
					UPRET(2, (s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal");
				}
				++pi;
				++pj;
			}
		} else {
			while (nnz--) {
				if (*pi <= *pj) {
					int s = *pi == *pj;
					UPRET(2, (s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal");
				}
				++pi;
				++pj;
			}
		}
		UNPROTECT(1); /* j */
	}
	UNPROTECT(1); /* i */

	return ScalarLogical(1);
}

SEXP xgCMatrix_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	if (XLENGTH(x) != XLENGTH(i))
		UPRET(2, "'i' and 'x' slots do not have equal length");
	UNPROTECT(2); /* i, x */
	return ScalarLogical(1);
}

SEXP xsCMatrix_validate(SEXP obj)
{
	SEXP val;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(val = xgCMatrix_validate(obj), &pid);
	if (TYPEOF(val) != STRSXP)
		REPROTECT(val = sCMatrix_validate(obj), pid);
	UNPROTECT(1); /* val */
	return val;
}

SEXP xtCMatrix_validate(SEXP obj)
{
	SEXP val;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(val = xgCMatrix_validate(obj), &pid);
	if (TYPEOF(val) != STRSXP)
		REPROTECT(val = tCMatrix_validate(obj), pid);
	UNPROTECT(1); /* val */
	return val;
}

SEXP xgRMatrix_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	if (XLENGTH(x) != XLENGTH(j))
		UPRET(2, "'j' and 'x' slots do not have equal length");
	UNPROTECT(2); /* j, x */
	return ScalarLogical(1);
}

SEXP xsRMatrix_validate(SEXP obj)
{
	SEXP val;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(val = xgRMatrix_validate(obj), &pid);
	if (TYPEOF(val) != STRSXP)
		REPROTECT(val = sRMatrix_validate(obj), pid);
	UNPROTECT(1); /* val */
	return val;
}

SEXP xtRMatrix_validate(SEXP obj)
{
	SEXP val;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(val = xgRMatrix_validate(obj), &pid);
	if (TYPEOF(val) != STRSXP)
		REPROTECT(val = tRMatrix_validate(obj), pid);
	UNPROTECT(1); /* val */
	return val;
}

SEXP xgTMatrix_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	if (XLENGTH(x) != XLENGTH(i))
		UPRET(2, "'i' and 'x' slots do not have equal length");
	UNPROTECT(2); /* i, x */
	return ScalarLogical(1);
}

SEXP xsTMatrix_validate(SEXP obj)
{
	SEXP val;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(val = xgTMatrix_validate(obj), &pid);
	if (TYPEOF(val) != STRSXP)
		REPROTECT(val = sTMatrix_validate(obj), pid);
	UNPROTECT(1); /* val */
	return val;
}

SEXP xtTMatrix_validate(SEXP obj)
{
	SEXP val;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(val = xgTMatrix_validate(obj), &pid);
	if (TYPEOF(val) != STRSXP)
		REPROTECT(val = tTMatrix_validate(obj), pid);
	UNPROTECT(1); /* val */
	return val;
}

SEXP unpackedMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int *pdim = INTEGER(dim);
	if (XLENGTH(x) != (double) pdim[0] * pdim[1])
		UPRET(2, "'x' slot does not have length prod(Dim)");
	UNPROTECT(2); /* x, dim */
	return ScalarLogical(1);
}

SEXP packedMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int *pdim = INTEGER(dim);
	if (XLENGTH(x) != 0.5 * pdim[0] * (pdim[0] + 1.0))
		UPRET(2, "'x' slot does not have length n*(n+1)/2, n=Dim[1]");
	UNPROTECT(2); /* x, dim */
	return ScalarLogical(1);
}

SEXP dpoMatrix_validate(SEXP obj)
{
	/* NB: Non-finite entries are "valid" because we consider
	   crossprod(x) and tcrossprod(x) to be positive semidefinite
	   even if 'x' contains non-finite entries (for speed) ...
	*/

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int j, n = INTEGER(dim)[0];
	R_xlen_t np1 = (R_xlen_t) n + 1;
	UNPROTECT(1); /* dim */

	/* Non-negative diagonal elements are necessary but _not_ sufficient */
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double *px = REAL(x);
	for (j = 0; j < n; ++j, px += np1)
		if (!ISNAN(*px) && *px < 0.0)
			UPRET(1, "matrix has negative diagonal elements");
	UNPROTECT(1); /* x */

	return ScalarLogical(1);
}

SEXP dppMatrix_validate(SEXP obj)
{
	/* NB: Non-finite entries are "valid" because we consider
	   crossprod(x) and tcrossprod(x) to be positive semidefinite
	   even if 'x' contains non-finite entries (for speed) ...
	*/

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int j, n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */

	/* Non-negative diagonal elements are necessary but _not_ sufficient */
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double *px = REAL(x);
	if (ul == 'U') {
		for (j = 0; j < n; px += (++j)+1)
			if (!ISNAN(*px) && *px < 0.0)
				UPRET(1, "matrix has negative diagonal elements");
	} else {
		for (j = 0; j < n; px += n-(j++))
			if (!ISNAN(*px) && *px < 0.0)
				UPRET(1, "matrix has negative diagonal elements");
	}
	UNPROTECT(1); /* x */

	return ScalarLogical(1);
}

SEXP corMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int j, n = INTEGER(dim)[0];
	R_xlen_t np1 = (R_xlen_t) n + 1;
	UNPROTECT(1); /* dim */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double *px = REAL(x);
	for (j = 0; j < n; ++j, px += np1)
		if (ISNAN(*px) || *px != 1.0)
			UPRET(1, "matrix has nonunit diagonal elements");
	UNPROTECT(1); /* x */

	SEXP sd = PROTECT(GET_SLOT(obj, Matrix_sdSym));
	if (TYPEOF(sd) != REALSXP)
		UPRET(1, "'sd' slot is not of type \"double\"");
	if (XLENGTH(sd) != n)
		UPRET(1, "'sd' slot does not have length n=Dim[1]");
	double *psd = REAL(sd);
	for (j = 0; j < n; ++j)
		if (!ISNAN(psd[j]) && psd[j] < 0.0)
			UPRET(1, "'sd' slot has negative elements");
	UNPROTECT(1); /* sd */

	return ScalarLogical(1);
}

SEXP pcorMatrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int j, n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double *px = REAL(x);
	if (ul == 'U') {
		for (j = 0; j < n; px += (++j)+1)
			if (ISNAN(*px) || *px != 1.0)
				UPRET(1, "matrix has nonunit diagonal elements");
	} else {
		for (j = 0; j < n; px += n-(j++))
			if (ISNAN(*px) || *px != 1.0)
				UPRET(1, "matrix has nonunit diagonal elements");
	}
	UNPROTECT(1); /* x */

	SEXP sd = PROTECT(GET_SLOT(obj, Matrix_sdSym));
	if (TYPEOF(sd) != REALSXP)
		UPRET(1, "'sd' slot is not of type \"double\"");
	if (XLENGTH(sd) != n)
		UPRET(1, "'sd' slot does not have length n=Dim[1]");
	double *psd = REAL(sd);
	for (j = 0; j < n; ++j)
		if (!ISNAN(psd[j]) && psd[j] < 0.0)
			UPRET(1, "'sd' slot has negative elements");
	UNPROTECT(1); /* sd */

	return ScalarLogical(1);
}

SEXP denseLU_validate(SEXP obj)
{
	/* In R, we start by checking that 'obj' would be a valid dgeMatrix */

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	UNPROTECT(1); /* dim */

	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	if (TYPEOF(perm) != INTSXP)
		UPRET(1, "'perm' slot is not of type \"integer\"");
	if (XLENGTH(perm) != r)
		UPRET(1, "'perm' slot does not have length min(Dim)");
	int *pperm = INTEGER(perm);
	while (r--) {
		if (*pperm == NA_INTEGER)
			UPRET(1, "'perm' slot contains NA");
		if (*pperm < 1 || *pperm > m)
			UPRET(1, "'perm' slot has elements not in {1,...,Dim[1]}");
		++pperm;
	}
	UNPROTECT(1); /* perm */

	return ScalarLogical(1);
}

SEXP sparseLU_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		UPRET(1, "Dim[1] != Dim[2] (matrix is not square)");
	UNPROTECT(1); /* dim */

	SEXP L = PROTECT(GET_SLOT(obj, Matrix_LSym)), uplo, diag;
	PROTECT(dim = GET_SLOT(L, Matrix_DimSym));
	pdim = INTEGER(dim);
	if (pdim[0] != n || pdim[1] != n)
		UPRET(2, "dimensions of 'L' slot are not identical to 'Dim'");

	PROTECT(uplo = GET_SLOT(L, Matrix_uploSym));
	if (*CHAR(STRING_ELT(uplo, 0)) == 'U')
	UPRET(3, "'L' slot is upper (not lower) triangular");

	PROTECT(diag = GET_SLOT(L, Matrix_diagSym));
	if (*CHAR(STRING_ELT(diag, 0)) == 'N') {
		SEXP p = PROTECT(GET_SLOT(L, Matrix_pSym)),
			i = PROTECT(GET_SLOT(L, Matrix_iSym)),
			x = PROTECT(GET_SLOT(L, Matrix_xSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k = 0, kend;
		double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			kend = *(++pp);
			if (kend == k || pi[k] != j || px[k] != 1.0)
				UPRET(7, "'L' slot has nonunit diagonal elements");
			k = kend;
		}
		UNPROTECT(3); /* x, i, p */
	}
	UNPROTECT(4); /* diag, uplo, dim, L */

	SEXP U = PROTECT(GET_SLOT(obj, Matrix_USym));
	PROTECT(dim = GET_SLOT(U, Matrix_DimSym));
	pdim = INTEGER(dim);
	if (pdim[0] != n || pdim[1] != n)
		UPRET(2, "dimensions of 'U' slot are not identical to 'Dim'");

	PROTECT(uplo = GET_SLOT(U, Matrix_uploSym));
	if (*CHAR(STRING_ELT(uplo, 0)) != 'U')
		UPRET(3, "'U' slot is lower (not upper) triangular");
	UNPROTECT(3); /* uplo, dim, U */

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		q = PROTECT(GET_SLOT(obj, Matrix_qSym));
	if (TYPEOF(p) != INTSXP)
		UPRET(2, "'p' slot is not of type \"integer\"");
	if (TYPEOF(q) != INTSXP)
		UPRET(2, "'q' slot is not of type \"integer\"");
	if (XLENGTH(p) != n)
		UPRET(2, "'p' slot does not have length Dim[1]");
	if (XLENGTH(q) != n && XLENGTH(q) != 0)
		UPRET(2, "'q' slot does not have length Dim[1] or length 0");
	int i, *pp = INTEGER(p);
	char *work;
	Matrix_Calloc(work, n, char);
	for (i = 0; i < n; ++i) {
		if (*pp == NA_INTEGER)
			FRUPRET(work, n, 2, "'p' slot contains NA");
		if (*pp < 0 || *pp >= n)
			FRUPRET(work, n, 2, "'p' slot has elements not in {0,...,Dim[1]-1}");
		if (work[*pp])
			FRUPRET(work, n, 2, "'p' slot contains duplicates");
		work[*(pp++)] = 1;
	}
	if (LENGTH(q) == n) {
		int *pq = INTEGER(q);
		for (i = 0; i < n; ++i) {
			if (*pq == NA_INTEGER)
				FRUPRET(work, n, 2, "'q' slot contains NA");
			if (*pq < 0 || *pq >= n)
				FRUPRET(work, n, 2, "'q' slot has elements not in {0,...,Dim[2]-1}");
			if (!work[*pq])
				FRUPRET(work, n, 2, "'q' slot contains duplicates");
			work[*(pq++)] = 0;
		}
	}
	Matrix_Free(work, n);
	UNPROTECT(2); /* q, p */

	return ScalarLogical(1);
}

SEXP sparseQR_validate(SEXP obj)
{
	/* MJ: assuming for simplicity that 'V' and 'R' slots are formally valid */

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m < n)
		UPRET(1, "matrix has more columns than rows");
	UNPROTECT(1); /* dim */

	SEXP beta = PROTECT(GET_SLOT(obj, Matrix_betaSym));
	if (TYPEOF(beta) != REALSXP)
		UPRET(1, "'beta' slot is not of type \"double\"");
	if (XLENGTH(beta) != n)
		UPRET(1, "'beta' slot does not have length Dim[2]");
	UNPROTECT(1); /* beta */

	SEXP V = PROTECT(GET_SLOT(obj, Matrix_VSym));
	PROTECT(dim = GET_SLOT(V, Matrix_DimSym));
	pdim = INTEGER(dim);
	int m2 = pdim[0];
	if (m2 < m)
		UPRET(2, "'V' slot has fewer than Dim[1] rows");
	if (m2 > m + n)
		UPRET(2, "'V' slot has more than Dim[1]+Dim[2] rows");
	if (pdim[1] != n)
		UPRET(2, "'V' slot does not have Dim[2] columns");

	SEXP V_p = PROTECT(GET_SLOT(V, Matrix_pSym)),
	V_i = PROTECT(GET_SLOT(V, Matrix_iSym));
	int *V_pp = INTEGER(V_p), *V_pi = INTEGER(V_i), j, k, kend;
	for (j = 0, k = 0; j < n; ++j) {
		kend = *(++V_pp);
		if (k < kend) {
			if (V_pi[k] < j)
				UPRET(4, "'V' slot must be lower trapezoidal but has entries above the diagonal");
		}
		k = kend;
	}
	UNPROTECT(4); /* V_i, V_p, dim, V */

	SEXP R = PROTECT(GET_SLOT(obj, Matrix_RSym));
	PROTECT(dim = GET_SLOT(R, Matrix_DimSym));
	pdim = INTEGER(dim);
	if (pdim[0] != m2)
		UPRET(2, "'R' slot does not have nrow(V) rows");
	if (pdim[1] != n)
		UPRET(2, "'R' slot does not have Dim[2] columns");

	SEXP R_p = PROTECT(GET_SLOT(R, Matrix_pSym)),
	R_i = PROTECT(GET_SLOT(R, Matrix_iSym)),
	R_x = PROTECT(GET_SLOT(R, Matrix_xSym));
	int *R_pp = INTEGER(R_p), *R_pi = INTEGER(R_i);
	double *R_px = REAL(R_x);
	for (j = 0, k = 0; j < n; ++j) {
		kend = *(++R_pp);
		if (k < kend) {
			if (R_pi[kend - 1] > j)
				UPRET(5, "'R' slot must be upper trapezoidal but has entries below the diagonal");
			if (R_pi[kend - 1] == j &&
				!ISNAN(R_px[kend - 1]) && R_px[kend - 1] < 0.0)
				UPRET(5, "'R' slot has negative diagonal elements");
		}
		k = kend;
	}
	UNPROTECT(5); /* R_x, R_i, R_p, dim, R */

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		q = PROTECT(GET_SLOT(obj, Matrix_qSym));
	if (TYPEOF(p) != INTSXP)
		UPRET(2, "'p' slot is not of type \"integer\"");
	if (TYPEOF(q) != INTSXP)
		UPRET(2, "'q' slot is not of type \"integer\"");
	if (XLENGTH(p) != m2)
		UPRET(2, "'p' slot does not have length nrow(V)");
	if (XLENGTH(q) != n && XLENGTH(q) != 0)
		UPRET(2, "'q' slot does not have length Dim[2] or length 0");
	int i, *pp = INTEGER(p);
	char *work;
	Matrix_Calloc(work, m2, char); /* n <= m <= m2 */
	for (i = 0; i < m2; ++i) {
		if (*pp == NA_INTEGER)
			FRUPRET(work, m2, 2, "'p' slot contains NA");
		if (*pp < 0 || *pp >= m2)
			FRUPRET(work, m2, 2, "'p' slot has elements not in {0,...,nrow(V)-1}");
		if (work[*pp])
			FRUPRET(work, m2, 2, "'p' slot contains duplicates");
		work[*(pp++)] = 1;
	}
	if (LENGTH(q) == n) {
		int *pq = INTEGER(q);
		for (i = 0; i < n; ++i) {
			if (*pq == NA_INTEGER)
				FRUPRET(work, m2, 2, "'q' slot contains NA");
			if (*pq < 0 || *pq >= n)
				FRUPRET(work, m2, 2, "'q' slot has elements not in {0,...,Dim[2]-1}");
			if (!work[*pq])
				FRUPRET(work, m2, 2, "'q' slot contains duplicates");
			work[*(pq++)] = 0;
		}
	}
	Matrix_Free(work, m2);
	UNPROTECT(2); /* q, p */

	return ScalarLogical(1);
}

SEXP BunchKaufman_validate(SEXP obj)
{
	/* In R, we start by checking that 'obj' would be a valid dtrMatrix */

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */

	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	if (TYPEOF(perm) != INTSXP)
		UPRET(1, "'perm' slot is not of type \"integer\"");
	if (XLENGTH(perm) != n)
		UPRET(1, "'perm' slot does not have length n=Dim[1]");
	int n_ = n, *pperm = INTEGER(perm);
	while (n_) {
		if (*pperm == NA_INTEGER)
			UPRET(1, "'perm' slot contains NA");
		if (*pperm < -n || *pperm == 0 || *pperm > n)
			UPRET(1, "'perm' slot has elements not in {-n,...,n}\\{0}, n=Dim[1]");
		if (*pperm > 0) {
			pperm += 1;
			n_ -= 1;
		} else if (n_ > 1 && *(pperm + 1) == *pperm) {
			pperm += 2;
			n_ -= 2;
		} else
			UPRET(1, "'perm' slot has an unpaired negative element");
	}
	UNPROTECT(1); /* perm */

	return ScalarLogical(1);
}

SEXP pBunchKaufman_validate(SEXP obj)
{
	/* In R, we start by checking that 'obj' would be a valid dtpMatrix */

	return BunchKaufman_validate(obj);
}

SEXP Cholesky_validate(SEXP obj)
{
	/* In R, we start by checking that 'obj' would be a valid dtrMatrix */

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int j, n = INTEGER(dim)[0];
	R_xlen_t np1 = (R_xlen_t) n + 1;
	UNPROTECT(1); /* dim */

	/* Non-negative diagonal elements are necessary _and_ sufficient */
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double *px = REAL(x);
	for (j = 0; j < n; ++j, px += np1)
		if (!ISNAN(*px) && *px < 0.0)
			UPRET(1, "Cholesky factor has negative diagonal elements");
	UNPROTECT(1); /* x */

#if 0
	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	if (TYPEOF(perm) != INTSXP)
		UPRET(1, "'perm' slot is not of type \"integer\"");
	if (XLENGTH(perm) != n && XLENGTH(perm) != 0)
		UPRET(1, "'perm' slot does not have length Dim[1] or length 0");
	if (LENGTH(perm) == n) {
		int *pperm = INTEGER(perm);
		char *work;
		Matrix_Calloc(work, n, char);
		for (i = 0; i < n; ++i) {
			if (*pperm == NA_INTEGER)
				FRUPRET(work, n, 1, "'perm' slot contains NA");
			if (*pperm < 0 || *pperm >= n)
				FRUPRET(work, n, 1, "'perm' slot has elements not in {0,...,Dim[1]-1}");
			if (work[*pperm])
				FRUPRET(work, n, 1, "'perm' slot contains duplicates");
			work[*(pperm++)] = 1;
		}
		Matrix_Free(work, n);
	}
	UNPROTECT(1); /* perm */
#endif

	return ScalarLogical(1);
}

SEXP pCholesky_validate(SEXP obj)
{
	/* In R, we start by checking that 'obj' would be a valid dtpMatrix */

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int j, n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */

	/* Non-negative diagonal elements are necessary _and_ sufficient */
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double *px = REAL(x);
	if (ul == 'U') {
		for (j = 0; j < n; px += (++j)+1)
			if (!ISNAN(*px) && *px < 0.0)
				UPRET(1, "Cholesky factor has negative diagonal elements");
	} else {
		for (j = 0; j < n; px += n-(j++))
			if (!ISNAN(*px) && *px < 0.0)
				UPRET(1, "Cholesky factor has negative diagonal elements");
	}
	UNPROTECT(1); /* x */

#if 0
	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	if (TYPEOF(perm) != INTSXP)
		UPRET(1, "'perm' slot is not of type \"integer\"");
	if (XLENGTH(perm) != n && XLENGTH(perm) != 0)
		UPRET(1, "'perm' slot does not have length Dim[1] or length 0");
	if (LENGTH(perm) == n) {
		int *pperm = INTEGER(perm);
		char *work;
		Matrix_Calloc(work, n, char);
		for (i = 0; i < n; ++i) {
			if (*pperm == NA_INTEGER)
				FRUPRET(work, n, 1, "'perm' slot contains NA");
			if (*pperm < 0 || *pperm >= n)
				FRUPRET(work, n, 1, "'perm' slot has elements not in {0,...,Dim[1]-1}");
			if (work[*pperm])
				FRUPRET(work, n, 1, "'perm' slot contains duplicates");
			work[*(pperm++)] = 1;
		}
		Matrix_Free(work, n);
	}
	UNPROTECT(1); /* perm */
#endif

	return ScalarLogical(1);
}

SEXP CHMfactor_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		UPRET(1, "Dim[1] != Dim[2] (matrix is not square)");
	UNPROTECT(1); /* dim */

	SEXP type = PROTECT(GET_SLOT(obj, install("type")));
	if (TYPEOF(type) != INTSXP)
		UPRET(1, "'type' slot is not of type \"integer\"");
	if (XLENGTH(type) != 6)
		UPRET(1, "'type' slot does not have length 6");
	int order = INTEGER(type)[0];
	if (order < 0 || order > 4)
		UPRET(1, "type[1] (cholmod_factor.ordering) is not in 0:4");
	UNPROTECT(1); /* type */

	SEXP colcount = PROTECT(GET_SLOT(obj, install("colcount")));
	if (TYPEOF(colcount) != INTSXP)
		UPRET(1, "'colcount' slot is not of type \"integer\"");
	if (XLENGTH(colcount) != n)
		UPRET(1, "'colcount' slot does not have length Dim[2]");
	int j, *pcolcount = INTEGER(colcount);
	for (j = 0; j < n; ++j) {
		if (pcolcount[j] == NA_INTEGER)
			UPRET(1, "'colcount' slot contains NA");
		if (pcolcount[j] < 0 || pcolcount[j] > n - j)
			UPRET(1, "colcount[j] is not in {0,...,Dim[2]-j+1)}");
	}
	UNPROTECT(1); /* colcount */

	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	if (TYPEOF(perm) != INTSXP)
		UPRET(1, "'perm' slot is not of type \"integer\"");
	if (order == 0) {
		if (XLENGTH(perm) != 0)
			UPRET(1, "'perm' slot does not have length 0");
	} else {
		if (XLENGTH(perm) != n)
			UPRET(1, "'perm' slot does not have length Dim[1]");
		int *pperm = INTEGER(perm);
		char *work;
		Matrix_Calloc(work, n, char);
		for (j = 0; j < n; ++j) {
			if (*pperm == NA_INTEGER)
				FRUPRET(work, n, 1, "'perm' slot contains NA");
			if (*pperm < 0 || *pperm >= n)
				FRUPRET(work, n, 1, "'perm' slot has elements not in {0,...,Dim[1]-1}");
			if (work[*pperm])
				FRUPRET(work, n, 1, "'perm' slot contains duplicates");
			work[*(pperm++)] = 1;
		}
		Matrix_Free(work, n);
	}
	UNPROTECT(1); /* perm */

	return ScalarLogical(1);
}

SEXP CHMsimpl_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n == INT_MAX)
		UPRET(1, "Dim[1]+1 is not representable as \"integer\"");
	UNPROTECT(1); /* dim */

	SEXP type = PROTECT(GET_SLOT(obj, install("type")));
	int *ptype = INTEGER(type), mono = ptype[3];
	if (ptype[1] != 0 && ptype[1] != 1)
		UPRET(1, "type[2] (cholmod_factor.is_ll) is not 0 or 1");
	if (ptype[2] != 0)
		UPRET(1, "type[3] (cholmod_factor.is_super) is not 0");
	if (ptype[3] != 0 && ptype[3] != 1)
		UPRET(1, "type[4] (cholmod_factor.is_monotonic) is not 0 or 1");
	UNPROTECT(1); /* type */

	SEXP nxt = PROTECT(GET_SLOT(obj, install("nxt"))),
	prv = PROTECT(GET_SLOT(obj, install("prv")));
	if (TYPEOF(nxt) != INTSXP)
		UPRET(2, "'nxt' slot is not of type \"integer\"");
	if (TYPEOF(prv) != INTSXP)
		UPRET(2, "'prv' slot is not of type \"integer\"");
	if (XLENGTH(nxt) - 2 != n)
		UPRET(2, "'nxt' slot does not have length Dim[2]+2");
	if (XLENGTH(prv) - 2 != n)
		UPRET(2, "'prv' slot does not have length Dim[2]+2");
	int *pnxt = INTEGER(nxt), *pprv = INTEGER(prv),
	j1 = pnxt[n + 1], j2 = pprv[n], count = n + 1;
	while (count--) {
		if (j1 < 0 || j1 > n)
			UPRET(2, "nxt[-(n+1)] has elements not in {0,...,n}, n=Dim[2]");
		if (j2 < 0 || j2 > n + 1 || j2 == n)
			UPRET(2, "prv[-(n+2)] has elements not in {0,...,n+1}\\{n}, n=Dim[2]");
		if ((count >  1) && mono && (pnxt[j1] != j1 + 1 || pprv[j2] != j2 - 1))
			UPRET(2, "type[4] is 1 but columns are not stored in increasing order");
		if ((count >= 1) ? j1 == n : j1 != n)
			UPRET(2, "traversal of 'nxt' slot does not complete in exactly length(nxt) steps");
		if ((count >= 1) ? j2 == n + 1 : j2 != n + 1)
			UPRET(2, "traversal of 'prv' slot does not complete in exactly length(prv) steps");
		j1 = pnxt[j1];
		j2 = pprv[j2];
	}
	if (j1 != -1)
		UPRET(2, "nxt[Dim[2]+1] is not -1");
	if (j2 != -1)
		UPRET(2, "prv[Dim[2]+2] is not -1");

	SEXP nz = PROTECT(GET_SLOT(obj, install("nz")));
	if (TYPEOF(nz) != INTSXP)
		UPRET(3, "'nz' slot is not of type \"integer\"");
	if (XLENGTH(nz) != n)
		UPRET(3, "'nz' slot does not have length Dim[2]");
	int j, *pnz = INTEGER(nz);
	for (j = 0; j < n; ++j) {
		if (pnz[j] == NA_INTEGER)
			UPRET(3, "'nz' slot contains NA");
		if (pnz[j] < 1 || pnz[j] > n - j)
			UPRET(3, "nz[j] is not in {1,...,Dim[2]-j+1)}");
	}

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	if (TYPEOF(p) != INTSXP)
		UPRET(4, "'p' slot is not of type \"integer\"");
	if (XLENGTH(p) - 1 != n)
		UPRET(4, "'p' slot does not have length Dim[2]+1");
	j1 = pnxt[n + 1];
	int *pp = INTEGER(p);
	if (pp[j1] != 0)
		UPRET(4, "column 'j' is stored first but p[j] is not 0");
	for (j = 0; j < n; ++j) {
		j2 = pnxt[j1];
		if (pp[j2] == NA_INTEGER)
			UPRET(4, "'p' slot contains NA");
		if (pp[j2] < pp[j1])
			UPRET(4, "'p' slot is not increasing when traversed in stored column order");
		if (pp[j2] - pp[j1] < pnz[j1])
			UPRET(4, "'i' slot allocates fewer than nz[j] elements for column 'j'");
		if (pp[j2] - pp[j1] > n - j1)
			UPRET(4, "'i' slot allocates more than Dim[2]-j+1 elements for column 'j'");
		j1 = j2;
	}

	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	if (TYPEOF(i) != INTSXP)
		UPRET(5, "'i' slot is not of type \"integer\"");
	if (XLENGTH(i) != pp[n])
		UPRET(5, "'i' slot does not have length p[length(p)]");
	int *pi = INTEGER(i), *pi_, k;
	j1 = pnxt[n + 1];
	for (j = 0; j < n; ++j) {
		pi_ = pi + pp[j1];
		if (pi_[0] != j1)
			UPRET(5, "first entry in column 'j' does not have row index 'j'");
		for (k = 1; k < pnz[j1]; ++k) {
			if (pi_[k] == NA_INTEGER)
				UPRET(5, "'i' slot contains NA");
			if (pi_[k] < 0 || pi_[k] >= n)
				UPRET(5, "'i' slot has elements not in {0,...,Dim[1]-1}");
			if (pi_[k] <= pi_[k - 1])
				UPRET(5, "'i' slot is not increasing within columns");
		}
		j1 = pnxt[j1];
	}
	UNPROTECT(5); /* i, p, nz, prv, nxt */

	return ScalarLogical(1);
}

SEXP CHMsuper_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */

	SEXP type = PROTECT(GET_SLOT(obj, install("type")));
	int *ptype = INTEGER(type);
	if (ptype[1] != 1)
		UPRET(1, "type[2] (cholmod_factor.is_ll) is not 1");
	if (ptype[2] != 1)
		UPRET(1, "type[3] (cholmod_factor.is_super) is not 1");
	if (ptype[3] != 1)
		UPRET(1, "type[4] (cholmod_factor.is_monotonic) is not 1");
	if (ptype[4] < 0)
		UPRET(1, "type[5] (cholmod_factor.maxcsize) is negative");
	if (ptype[5] < 0)
		UPRET(1, "type[6] (cholmod_factor.maxesize) is negative");
	if (n > 0 && ptype[5] >= n)
		UPRET(1, "type[6] (cholmod_factor.maxesize) is not less than Dim[1]");
	UNPROTECT(1); /* type */

	/* FIXME: maxcsize and maxesize are well-defined properties of the
	   factorization, so we should also test that the values are
	   _correct_ ... see ./CHOLMOD/Supernodal/cholmod_super_symbolic.c
	*/

	SEXP super = PROTECT(GET_SLOT(obj, install("super")));
	if (TYPEOF(super) != INTSXP)
		UPRET(1, "'super' slot is not of type \"integer\"");
	R_xlen_t nsuper1a = XLENGTH(super);
	if (nsuper1a - 1 < ((n > 0) ? 1 : 0))
		UPRET(1, "'super' slot has length less than 2");
	if (nsuper1a - 1 > n)
		UPRET(1, "'super' slot has length greater than Dim[2]+1");
	int k, nsuper = (int) (nsuper1a - 1), *psuper = INTEGER(super);
	if (psuper[0] != 0)
		UPRET(1, "first element of 'super' slot is not 0");
	if (psuper[nsuper] != n)
		UPRET(1, "last element of 'super' slot is not Dim[2]");
	for (k = 1; k <= nsuper; ++k) {
		if (psuper[k] == NA_INTEGER)
			UPRET(1, "'super' slot contains NA");
		if (psuper[k] <= psuper[k-1])
			UPRET(1, "'super' slot is not increasing");
	}

	SEXP pi = PROTECT(GET_SLOT(obj, install("pi"))),
		px = PROTECT(GET_SLOT(obj, install("px")));
	if (TYPEOF(pi) != INTSXP)
		UPRET(3, "'pi' slot is not of type \"integer\"");
	if (TYPEOF(px) != INTSXP)
		UPRET(3, "'px' slot is not of type \"integer\"");
	if (XLENGTH(pi) != nsuper1a)
		UPRET(3, "'pi' and 'super' slots do not have equal length");
	if (XLENGTH(px) != nsuper1a)
		UPRET(3, "'px' and 'super' slots do not have equal length");
	int *ppi = INTEGER(pi), *ppx = INTEGER(px), nr, nc;
	if (ppi[0] != 0)
		UPRET(3, "first element of 'pi' slot is not 0");
	if (ppx[0] != 0)
		UPRET(3, "first element of 'px' slot is not 0");
	for (k = 1; k <= nsuper; ++k) {
		if (ppi[k] == NA_INTEGER)
			UPRET(3, "'pi' slot contains NA");
		if (ppx[k] == NA_INTEGER)
			UPRET(3, "'px' slot contains NA");
		if (ppi[k] <= ppi[k-1])
			UPRET(3, "'pi' slot is not increasing");
		if (ppx[k] <= ppx[k-1])
			UPRET(3, "'px' slot is not increasing");
		nr = ppi[k] - ppi[k-1];
		nc = psuper[k] - psuper[k-1];
		if (nr < nc)
			UPRET(3, "first differences of 'pi' slot are less than those of 'super' slot");
		if ((double) nr * nc > INT_MAX)
			UPRET(3, "supernode lengths exceed 2^31-1");
		if (ppx[k] - ppx[k-1] != nr * nc)
			UPRET(3, "first differences of 'px' slot are not equal to supernode lengths");
	}

	SEXP s = PROTECT(GET_SLOT(obj, install("s")));
	if (TYPEOF(s) != INTSXP)
		UPRET(4, "'s' slot is not of type \"integer\"");
	if (XLENGTH(s) != ppi[nsuper])
		UPRET(4, "'s' slot does not have length pi[length(pi)]");
	int i, j, *ps = INTEGER(s);
	for (k = 1; k <= nsuper; ++k) {
		nr = ppi[k] - ppi[k-1];
		nc = psuper[k] - (j = psuper[k-1]);
		for (i = 0; i < nr; ++i) {
			if (ps[i] == NA_INTEGER)
				UPRET(4, "'s' slot contains NA");
			if (ps[i] < 0 || ps[i] >= n)
				UPRET(4, "'s' slot has elements not in {0,...,Dim[1]-1}");
			if (i < nc) {
				if (ps[i] != j + i)
					UPRET(4, "'s' slot is wrong within diagonal blocks (row and column indices do not coincide)");
			} else {
				if (ps[i] <= ps[i-1])
					UPRET(4, "'s' slot is not increasing within supernodes");
			}
		}
		ps += nr;
	}
	UNPROTECT(4); /* s, px, pi, super */

	return ScalarLogical(1);
}

SEXP dCHMsimpl_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	if (TYPEOF(x) != REALSXP)
		UPRET(1, "'x' slot is not of type \"double\"");

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
	if (XLENGTH(x) != pp[n])
		UPRET(2, "'x' slot does not have length p[length(p)]");

	SEXP type = PROTECT(GET_SLOT(obj, install("type")));
	if (INTEGER(type)[1]) {
		int j;
		double *px = REAL(x);

		/* Non-negative diagonal elements are necessary _and_ sufficient */
		for (j = 0; j < n; ++j)
			if (!ISNAN(px[pp[j]]) && px[pp[j]] < 0.0)
				UPRET(3, "Cholesky factor has negative diagonal elements");
	}
	UNPROTECT(3); /* type, p, x */

	return ScalarLogical(1);
}

SEXP dCHMsuper_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	if (TYPEOF(x) != REALSXP)
		UPRET(1, "'x' slot is not of type \"double\"");

	SEXP px = PROTECT(GET_SLOT(obj, install("px")));
	int *ppx = INTEGER(px), nsuper = (int) (XLENGTH(px) - 1);
	if (XLENGTH(x) != ppx[nsuper])
		UPRET(2, "'x' slot does not have length px[length(px)]");

	SEXP pi = PROTECT(GET_SLOT(obj, install("pi"))),
		super = PROTECT(GET_SLOT(obj, install("super")));
	int *ppi = INTEGER(pi), *psuper = INTEGER(super), k, j, nc;
	double *pu = REAL(x), *pv;
	R_xlen_t nr1a;

	/* Non-negative diagonal elements are necessary _and_ sufficient */
	for (k = 0; k < nsuper; ++k) {
		nc = psuper[k+1] - psuper[k];
		nr1a = (R_xlen_t) (ppi[k+1] - ppi[k]) + 1;
		pv = pu + ppx[k];
		for (j = 0; j < nc; ++j) {
			if (!ISNAN(*pv) && *pv < 0.0)
				UPRET(4, "Cholesky factor has negative diagonal elements");
			pv += nr1a;
		}
	}
	UNPROTECT(4); /* super, pi, px, x */

	return ScalarLogical(1);
}

SEXP Schur_validate(SEXP obj)
{
	/* MJ: assuming for simplicity that 'Q' and 'T' slots are formally valid */

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		UPRET(1, "Dim[1] != Dim[2] (matrix is not square)");
	UNPROTECT(1); /* dim */

	SEXP Q = PROTECT(GET_SLOT(obj, Matrix_QSym));
	PROTECT(dim = GET_SLOT(Q, Matrix_DimSym));
	pdim = INTEGER(dim);
	if (pdim[0] != n || pdim[1] != n)
		UPRET(2, "dimensions of 'Q' slot are not identical to 'Dim'");
	UNPROTECT(2); /* dim, Q */

	SEXP T = PROTECT(GET_SLOT(obj, Matrix_TSym));
	PROTECT(dim = GET_SLOT(T, Matrix_DimSym));
	pdim = INTEGER(dim);
	if (pdim[0] != n || pdim[1] != n)
		UPRET(2, "dimensions of 'T' slot are not identical to 'Dim'");
	UNPROTECT(2); /* dim, T */

	SEXP v = PROTECT(GET_SLOT(obj, install("EValues")));
	SEXPTYPE tv = TYPEOF(v);
	if (tv != REALSXP && tv != CPLXSXP)
		UPRET(1, "'EValues' slot does not have type \"double\" or type \"complex\"");
	if (XLENGTH(v) != n)
		UPRET(1, "'EValues' slot does not have length n=Dim[1]");
	UNPROTECT(1); /* v */

	return ScalarLogical(1);
}

/* where 'cl' must be an element of VALID_NONVIRTUAL_MATRIX */
void validObject(SEXP obj, const char *cl)
{
#ifndef MATRIX_DISABLE_VALIDITY

	SEXP status;

# define IS_VALID(_CLASS_) \
	do { \
		status = _CLASS_ ## _validate(obj); \
		if (TYPEOF(status) == STRSXP) \
			error(_("invalid class \"%s\" object: %s"), \
				  cl, CHAR(STRING_ELT(status, 0))); \
	} while (0)

#define IS_VALID_SPARSE(_C_) \
	do { \
		IS_VALID(_C_ ## sparseMatrix); \
		if (cl[0] == 'n') { \
			if (cl[1] == 't') \
				IS_VALID(t ## _C_ ## Matrix); \
			else if (cl[1] == 's') \
				IS_VALID(s ## _C_ ## Matrix); \
		} else { \
			if (cl[1] == 'g') \
				IS_VALID(xg ## _C_ ## Matrix); \
			else if (cl[1] == 't') \
				IS_VALID(xt ## _C_ ## Matrix); \
			else if (cl[1] == 's') \
				IS_VALID(xs ## _C_ ## Matrix); \
		} \
	} while (0)

	IS_VALID(Matrix);

	if ((cl[0] == 'i' && cl[1] == 'n' && cl[2] == 'd') ||
		(cl[0] == 'p' && cl[1] != 'C' && cl[1] != 'c')) {
		IS_VALID(indMatrix);
		if (cl[0] == 'p')
			IS_VALID(pMatrix);
		return;
	}

	const char *cl_ = cl;
	if (cl[0] == 'C')
		cl = "dtrMatrix";
	else if (cl[0] == 'p' && cl[1] == 'C')
		cl = "dtpMatrix";
	else if (cl[0] == 'c')
		cl = "dpoMatrix";
	else if (cl[0] == 'p' && cl[1] == 'c')
		cl = "dppMatrix";

	if (cl[0] == 'n' && cl[2] != 'C' && cl[2] != 'R' && cl[2] != 'T')
		IS_VALID(ndenseMatrix);
	else if (cl[0] == 'l')
		IS_VALID(lMatrix);
	else if (cl[0] == 'i')
		IS_VALID(iMatrix);
	else if (cl[0] == 'd')
		IS_VALID(dMatrix);
	else if (cl[0] == 'z')
		IS_VALID(zMatrix);

	if (cl[1] == 't')
		IS_VALID(triangularMatrix);
	else if (cl[1] == 's' || cl[1] == 'p')
		IS_VALID(symmetricMatrix);
	else if (cl[1] == 'd') {
		IS_VALID(diagonalMatrix);
		return;
	}

	if (cl[2] == 'C')
		IS_VALID_SPARSE(C);
	else if (cl[2] == 'R')
		IS_VALID_SPARSE(R);
	else if (cl[2] == 'T')
		IS_VALID_SPARSE(T);
	else if (cl[2] != 'p') {
		IS_VALID(unpackedMatrix);
		if (cl_[0] == 'C')
			IS_VALID(Cholesky);
		else if (cl[1] == 'p') {
			IS_VALID(dpoMatrix);
			if (cl_[0] == 'c')
				IS_VALID(corMatrix);
		}
	} else {
		IS_VALID(packedMatrix);
		if (cl_[0] == 'p' && cl_[1] == 'C')
			IS_VALID(pCholesky);
		else if (cl[1] == 'p') {
			IS_VALID(dppMatrix);
			if (cl_[0] == 'p' && cl_[1] == 'c')
				IS_VALID(pcorMatrix);
		}
	}

# undef IS_VALID_SPARSE
# undef IS_VALID

#endif

	return;
}
