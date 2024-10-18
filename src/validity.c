#include <math.h> /* trunc */
#include "Mdefines.h"
#include "validity.h"

#define     MK(_FORMAT_     )       mkString(_FORMAT_             )
#define     MS(_FORMAT_, ...) Matrix_sprintf(_FORMAT_, __VA_ARGS__)

#define    RMK(_FORMAT_     ) \
	return MK(   _FORMAT_              )
#define    RMS(_FORMAT_, ...) \
	return    MS(_FORMAT_, __VA_ARGS__)
#define  RMKMS(_FORMAT_, ...) \
	return MK(MS(_FORMAT_, __VA_ARGS__))

#define FRMKMS(_FORMAT_, ...) \
	do { \
		Matrix_Free(work, lwork); \
		RMKMS(_FORMAT_, __VA_ARGS__); \
	} while (0)


/* Slot validity methods ===============================================
   Called by various class validity methods (see below).
*/

char *Dim_validate(SEXP dim)
{
	if (TYPEOF(dim) != INTSXP)
		RMS(_("'%s' slot is not of type \"%s\""), "Dim", "integer");
	if (XLENGTH(dim) != 2)
		RMS(_("'%s' slot does not have length %d"), "Dim", 2);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m == NA_INTEGER || n == NA_INTEGER)
		RMS(_("'%s' slot contains NA"), "Dim");
	if (m < 0 || n < 0)
		RMS(_("'%s' slot has negative elements"), "Dim");

	return NULL;
}

SEXP R_Dim_validate(SEXP dim)
{
	char *msg = Dim_validate(dim);
	return (msg) ? mkString(msg) : ScalarLogical(1);
}

char *DimNames_validate(SEXP dimnames, int *pdim)
{
	if (TYPEOF(dimnames) != VECSXP)
		RMS(_("'%s' slot is not a list"), "Dimnames");
	if (XLENGTH(dimnames) != 2)
		RMS(_("'%s' slot does not have length %d"), "Dimnames", 2);

	/* Behave as do_matrix() from src/main/array.c:
	   Dimnames[[i]] must be NULL or _coercible to_ character
	   of length Dim[i] or 0 ... see R_Dimnames_fixup() below
	*/

	SEXP s;
	int i;
	R_xlen_t ns;

	for (i = 0; i < 2; ++i) {
		s = VECTOR_ELT(dimnames, i);
		if (s == R_NilValue)
			continue;
		if (!isVector(s))
			RMS(_("%s[[%d]] is not NULL or a vector"), "Dimnames", i + 1);
		ns = XLENGTH(s);
		if (ns != pdim[i] && ns != 0)
			RMS(_("length of %s[[%d]] (%lld) is not equal to %s[%d] (%d)"),
			    "Dimnames", i + 1, (long long) ns,
			    "Dim"     , i + 1,        pdim[i]);
	}

	return NULL;
}

SEXP R_DimNames_validate(SEXP dimnames, SEXP dim)
{
	char *msg = DimNames_validate(dimnames, INTEGER(dim));
	return (msg) ? mkString(msg) : ScalarLogical(1);
}

SEXP R_DimNames_fixup(SEXP dimnames)
{
	SEXP s;
	int i, fixup = 0;
	for (i = 0; i < 2 && !fixup; ++i)
		fixup =
			(s = VECTOR_ELT(dimnames, i)) != R_NilValue &&
			(LENGTH(s) == 0 || TYPEOF(s) != STRSXP);
	if (!fixup)
		return dimnames;
	SEXP dimnames_ = PROTECT(allocVector(VECSXP, 2));
	for (i = 0; i < 2; ++i) {
		if ((s = VECTOR_ELT(dimnames, i)) == R_NilValue || LENGTH(s) == 0)
			continue;
		if (TYPEOF(s) == STRSXP)
			PROTECT(s);
		else if (TYPEOF(s) == INTSXP && inherits(s, "factor"))
			PROTECT(s = asCharacterFactor(s));
		else {
			PROTECT(s = coerceVector(s, STRSXP));
			CLEAR_ATTRIB(s);
		}
		SET_VECTOR_ELT(dimnames_, i, s);
		UNPROTECT(1); /* s */
	}
	s = getAttrib(dimnames, R_NamesSymbol);
	if (s != R_NilValue) {
		PROTECT(s);
		setAttrib(dimnames_, R_NamesSymbol, s);
		UNPROTECT(1); /* s */
	}
	UNPROTECT(1); /* dimnames_ */
	return dimnames_;
}


/* Class validity methods ==============================================
   NB: These assume that validity methods for superclasses
   have already been called via validObject() ...
*/

SEXP Matrix_validate(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	char *msg = Dim_validate(dim);
	if (!msg) {
		SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
		msg = DimNames_validate(dimnames, INTEGER(dim));
		UNPROTECT(1); /* dimnames */
	}
	UNPROTECT(1); /* dim */
	return (msg) ? mkString(msg) : ScalarLogical(1);
}

SEXP MatrixFactorization_validate(SEXP obj)
{
	return Matrix_validate(obj);
}

#define KINDMATRIX_VALIDATE(_PREFIX_, _SEXPTYPE_) \
SEXP _PREFIX_ ## Matrix_validate(SEXP obj) \
{ \
	SEXP x = GET_SLOT(obj, Matrix_xSym); \
	if (TYPEOF(x) != _SEXPTYPE_) \
		RMKMS(_("'%s' slot is not of type \"%s\""), "x", type2char(_SEXPTYPE_)); \
	return ScalarLogical(1); \
}
KINDMATRIX_VALIDATE(n,  LGLSXP)
KINDMATRIX_VALIDATE(l,  LGLSXP)
KINDMATRIX_VALIDATE(i,  INTSXP)
KINDMATRIX_VALIDATE(d, REALSXP)
KINDMATRIX_VALIDATE(z, CPLXSXP)
#undef KINDMATRIX_VALIDATE

SEXP generalMatrix_validate(SEXP obj)
{
	SEXP factors = GET_SLOT(obj, Matrix_factorsSym);
	if (TYPEOF(factors) != VECSXP)
		RMKMS(_("'%s' slot is not a list"), "factors");
	if (XLENGTH(factors) > 0) {
		PROTECT(factors);
		SEXP nms = getAttrib(factors, R_NamesSymbol);
		UNPROTECT(1); /* factors */
		if (nms == R_NilValue)
			RMKMS(_("'%s' slot has no '%s' attribute"), "factors", "names");
	}

	return ScalarLogical(1);
}

SEXP symmetricMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

#ifdef ENFORCE_SYMMETRIC_DIMNAMES

	/* This check can be expensive when both rownames and colnames have
	   nonzero length, and even more so when coercions to character are
	   required ... Users can avoid the expense by setting at least one
	   of rownames and colnames to NULL or by ensuring that they are the
	   same object, as testing for pointer equality is fast ...
	*/

# define ANY_TO_STRING(x) \
	((TYPEOF(x) == INTSXP && inherits(x, "factor")) \
	 ? asCharacterFactor(x) \
	 : coerceVector(x, STRSXP))

	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		ndn = getAttrib(dn, R_NamesSymbol);
	UNPROTECT(1); /* dn */

	const char *ndn0, *ndn1;
	if (ndn != R_NilValue &&
	    *(ndn0 = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
	    *(ndn1 = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
	    strcmp(ndn0, ndn1) != 0)
		RMKMS(_("%s[1] differs from %s[2]"), "Dimnames", "Dimnames");
	if (n > 0) {
		/* NB: It is already known that the length of 'dn[[i]]' is 0 or 'n' */
		SEXP rn, cn;
		if ((rn = VECTOR_ELT(dn, 0)) != R_NilValue &&
		    (cn = VECTOR_ELT(dn, 1)) != R_NilValue &&
		    LENGTH(rn) == n && LENGTH(cn) == n && rn != cn) {
			PROTECT(rn);
			PROTECT(cn);
			PROTECT(rn = ANY_TO_STRING(rn));
			PROTECT(cn = ANY_TO_STRING(cn));
			UNPROTECT(4); /* cn, rn */
			if (!equal_character_vectors(rn, cn, n))
				RMKMS(_("%s[1] differs from %s[2]"), "Dimnames", "Dimnames");
		}
	}

# undef ANY_TO_STRING

#endif

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	if (TYPEOF(uplo) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "uplo", "character");
	if (XLENGTH(uplo) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "uplo", 1);
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "uplo", "U", "L");

	return generalMatrix_validate(obj);
}

SEXP triangularMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	if (TYPEOF(uplo) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "uplo", "character");
	if (XLENGTH(uplo) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "uplo", 1);
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "uplo", "U", "L");

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	if (TYPEOF(diag) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "diag", "character");
	if (XLENGTH(diag) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "diag", 1);
	const char *di = CHAR(STRING_ELT(diag, 0));
	if (di[0] == '\0' || di[1] != '\0' || (di[0] != 'N' && di[0] != 'U'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "diag", "N", "U");

	return ScalarLogical(1);
}

SEXP unpackedMatrix_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	UNPROTECT(2); /* dim, x */
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (XLENGTH(x) != (Matrix_int_fast64_t) m * n)
		RMKMS(_("'%s' slot does not have length %s"), "x", "prod(Dim)");
	return ScalarLogical(1);
}

SEXP packedMatrix_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	UNPROTECT(2); /* dim, x */
	int n = INTEGER(dim)[0];
	if (XLENGTH(x) != n + ((Matrix_int_fast64_t) n * (n - 1)) / 2)
		RMKMS(_("'%s' slot does not have length %s"), "x", "Dim[1]*(Dim[1]+1)/2");
	return ScalarLogical(1);
}

SEXP CsparseMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	UNPROTECT(2); /* i, p */

	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (XLENGTH(p) - 1 != n)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[2]+1");
	int *pp = INTEGER(p);
	if (pp[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "p");
	int j;
	for (j = 1; j <= n; ++j) {
		if (pp[j] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "p");
		if (pp[j] < pp[j - 1])
			RMKMS(_("'%s' slot is not nondecreasing"), "p");
		if (pp[j] - pp[j - 1] > m)
			RMKMS(_("first differences of '%s' slot exceed %s"), "p", "Dim[1]");
	}

	if (TYPEOF(i) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "i", "integer");
	if (XLENGTH(i) < pp[n])
		RMKMS(_("'%s' slot has length less than %s"), "i", "p[length(p)]");
	int *pi = INTEGER(i), k, kend, ik, i0;
	for (j = 1, k = 0; j <= n; ++j) {
		kend = pp[j];
		i0 = -1;
		while (k < kend) {
			ik = pi[k];
			if (ik == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (ik < 0 || ik >= m)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "i", "0,...,Dim[1]-1");
			if (ik <= i0)
				RMKMS(_("'%s' slot is not increasing within columns"), "i");
			i0 = ik;
			++k;
		}
	}

	return ScalarLogical(1);
}

SEXP RsparseMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	UNPROTECT(2); /* j, p */

	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (XLENGTH(p) - 1 != m)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[1]+1");
	int *pp = INTEGER(p);
	if (pp[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "p");
	int i;
	for (i = 1; i <= m; ++i) {
		if (pp[i] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "p");
		if (pp[i] < pp[i - 1])
			RMKMS(_("'%s' slot is not nondecreasing"), "p");
		if (pp[i] - pp[i - 1] > n)
			RMKMS(_("first differences of '%s' slot exceed %s"), "p", "Dim[2]");
	}

	if (TYPEOF(j) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "j", "integer");
	if (XLENGTH(j) < pp[m])
		RMKMS(_("'%s' slot has length less than %s"), "j", "p[length(p)]");
	int *pj = INTEGER(j), k, kend, jk, j0;
	for (i = 1, k = 0; i <= m; ++i) {
		kend = pp[i];
		j0 = -1;
		while (k < kend) {
			jk = pj[k];
			if (jk == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "j");
			if (jk < 0 || jk >= n)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "j", "0,...,Dim[2]-1");
			if (jk <= j0)
				RMKMS(_("'%s' slot is not increasing within rows"), "j");
			j0 = jk;
			++k;
		}
	}

	return ScalarLogical(1);
}

SEXP TsparseMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	UNPROTECT(2); /* j, i */

	if (TYPEOF(i) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "i", "integer");
	if (TYPEOF(j) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "j", "integer");
	R_xlen_t nnz = XLENGTH(i);
	if (XLENGTH(j) != nnz)
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "i", "j");
	if (nnz > 0) {
		if (m == 0 || n == 0)
			RMKMS(_("'%s' slot has nonzero length but %s is 0"), "i", "prod(Dim)");
		int *pi = INTEGER(i), *pj = INTEGER(j);
		while (nnz--) {
			if (*pi == NA_LOGICAL)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (*pj == NA_LOGICAL)
				RMKMS(_("'%s' slot contains NA"), "j");
			if (*pi < 0 || *pi >= m)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "i", "0,...,Dim[1]-1");
			if (*pj < 0 || *pj >= n)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "j", "0,...,Dim[2]-1");
			++pi;
			++pj;
		}
	}

	return ScalarLogical(1);
}

SEXP diagonalMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	if (TYPEOF(diag) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "diag", "character");
	if (XLENGTH(diag) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "diag", 1);
	const char *di = CHAR(STRING_ELT(diag, 0));
	if (di[0] == '\0' || di[1] != '\0' || (di[0] != 'N' && di[0] != 'U'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "diag", "N", "U");
	int nonunit = di[0] == 'N';

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (nonunit) {
		if (XLENGTH(x) != n)
			RMKMS(_("'%s' slot is \"%s\" but '%s' slot does not have length %s"),
			      "diag", "N", "x", "Dim[1]");
	} else {
		if (XLENGTH(x) != 0)
			RMKMS(_("'%s' slot is \"%s\" but '%s' slot does not have length %s"),
			      "diag", "U", "x",      "0");
	}

	return ScalarLogical(1);
}

SEXP indMatrix_validate(SEXP obj)
{
	SEXP margin = GET_SLOT(obj, Matrix_marginSym);
	if (TYPEOF(margin) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "margin", "integer");
	if (XLENGTH(margin) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "margin", 1);
	int mg = INTEGER(margin)[0] - 1;
	if (mg != 0 && mg != 1)
		RMKMS(_("'%s' slot is not %d or %d"), "margin", 1, 2);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[mg], n = pdim[!mg];
	if (m > 0 && n == 0) {
		if (mg == 0)
			RMKMS(_("%s-by-%s %s invalid for positive '%s' when %s=%d"),
			      "m", "0", "indMatrix", "m", "margin", 1);
		else
			RMKMS(_("%s-by-%s %s invalid for positive '%s' when %s=%d"),
			      "0", "n", "indMatrix", "n", "margin", 2);
	}

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != m)
		RMKMS(_("'%s' slot does not have length %s"), "perm", "Dim[margin]");
	int *pperm = INTEGER(perm);
	while (m--) {
		if (*pperm == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "perm");
		if (*pperm < 1 || *pperm > n)
			RMKMS(_("'%s' slot has elements not in {%s}"),
			      "perm", "1,...,Dim[1+margin%%2]");
		++pperm;
	}

	return ScalarLogical(1);
}

SEXP pMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	if (n > 1) {
		SEXP perm = GET_SLOT(obj, Matrix_permSym);
		char *work;
		int lwork = n;
		Matrix_Calloc(work, lwork, char);
		int j, *pperm = INTEGER(perm);
		for (j = 0; j < n; ++j) {
			if (work[*pperm - 1])
				FRMKMS(_("'%s' slot contains duplicates"), "perm");
			work[*(pperm++) - 1] = 1;
		}
		Matrix_Free(work, lwork);
	}

	return ScalarLogical(1);
}

SEXP sCMatrix_validate(SEXP obj)
{
	SEXP p = GET_SLOT(obj, Matrix_pSym);
	int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
	if (pp[n] > 0) {
		PROTECT(p);

		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));

		SEXP i = GET_SLOT(obj, Matrix_iSym);
		int *pi = INTEGER(i), j, k, kend;

		UNPROTECT(1); /* p */
		if (ul == 'U') {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j + 1];
				while (k < kend) {
					if (pi[k] > j)
						RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
						      "uplo", "U");
					++k;
				}
			}
		} else {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j + 1];
				while (k < kend) {
					if (pi[k] < j)
						RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
						      "uplo", "L");
					++k;
				}
			}
		}
	}

	return ScalarLogical(1);
}

SEXP tCMatrix_validate(SEXP obj)
{
	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char di = *CHAR(STRING_ELT(diag, 0));
	if (di == 'N')
		return sCMatrix_validate(obj);

	SEXP p = GET_SLOT(obj, Matrix_pSym);
	int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
	if (pp[n] > 0) {
		PROTECT(p);

		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));

		SEXP i = GET_SLOT(obj, Matrix_iSym);
		int *pi = INTEGER(i), j, k, kend;

		UNPROTECT(1); /* p */
		if (ul == 'U') {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j + 1];
				while (k < kend) {
					if (pi[k] >  j)
						RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
						      "uplo", "U");
					if (pi[k] == j)
						RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
						      "diag", "U");
					++k;
				}
			}
		} else {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j + 1];
				while (k < kend) {
					if (pi[k] <  j)
						RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
						      "uplo", "L");
					if (pi[k] == j)
						RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
						      "diag", "U");
					++k;
				}
			}
		}
	}

	return ScalarLogical(1);
}

SEXP sRMatrix_validate(SEXP obj)
{
	SEXP p = GET_SLOT(obj, Matrix_pSym);
	int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
	if (pp[m] > 0) {
		PROTECT(p);

		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));

		SEXP j = GET_SLOT(obj, Matrix_jSym);
		int *pj = INTEGER(j), i, k, kend;

		UNPROTECT(1); /* p */
		if (ul == 'U') {
			for (i = 0, k = 0; i < m; ++i) {
				kend = pp[i + 1];
				while (k < kend) {
					if (pj[k] < i)
						RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
						      "uplo", "U");
					++k;
				}
			}
		} else {
			for (i = 0, k = 0; i < m; ++i) {
				kend = pp[i + 1];
				while (k < kend) {
					if (pj[k] > i)
						RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
						      "uplo", "L");
					++k;
				}
			}
		}
	}

	return ScalarLogical(1);
}

SEXP tRMatrix_validate(SEXP obj)
{
	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char di = *CHAR(STRING_ELT(diag, 0));
	if (di == 'N')
		return sRMatrix_validate(obj);

	SEXP p = GET_SLOT(obj, Matrix_pSym);
	int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
	if (pp[m] > 0) {
		PROTECT(p);

		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));

		SEXP j = GET_SLOT(obj, Matrix_jSym);
		int *pj = INTEGER(j), i, k, kend;

		UNPROTECT(1); /* p */
		if (ul == 'U') {
			for (i = 0, k = 0; i < m; ++i) {
				kend = pp[i + 1];
				while (k < kend) {
					if (pj[k] <  i)
						RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
						      "uplo", "U");
					if (pj[k] == i)
						RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
						      "diag", "U");
					++k;
				}
			}
		} else {
			for (i = 0, k = 0; i < m; ++i) {
				kend = pp[i + 1];
				while (k < kend) {
					if (pj[k] >  i)
						RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
						      "uplo", "L");
					if (pj[k] == i)
						RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
						      "diag", "U");
					++k;
				}
			}
		}
	}

	return ScalarLogical(1);
}

SEXP sTMatrix_validate(SEXP obj)
{
	SEXP i = GET_SLOT(obj, Matrix_iSym);
	R_xlen_t nnz = XLENGTH(i);
	if (nnz > 0) {
		PROTECT(i);

		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));

		SEXP j = GET_SLOT(obj, Matrix_jSym);
		int *pi = INTEGER(i), *pj = INTEGER(j);

		UNPROTECT(1); /* i */
		if (ul == 'U') {
			while (nnz--)
				if (*(pi++) > *(pj++))
					RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
					      "uplo", "U");
		} else {
			while (nnz--)
				if (*(pi++) < *(pj++))
					RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
					      "uplo", "L");
		}
	}

	return ScalarLogical(1);
}

SEXP tTMatrix_validate(SEXP obj)
{
	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char di = *CHAR(STRING_ELT(diag, 0));
	if (di == 'N')
		return sTMatrix_validate(obj);

	SEXP i = GET_SLOT(obj, Matrix_iSym);
	R_xlen_t nnz = XLENGTH(i);
	if (nnz > 0) {
		PROTECT(i);

		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));

		SEXP j = GET_SLOT(obj, Matrix_jSym);
		int *pi = INTEGER(i), *pj = INTEGER(j);

		UNPROTECT(1); /* i */
		if (ul == 'U') {
			while (nnz--) {
				if (*pi >  *pj)
					RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
					      "uplo", "U");
				if (*pi == *pj)
					RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
					      "diag", "U");
				++pi;
				++pj;
			}
		} else {
			while (nnz--) {
				if (*pi <  *pj)
					RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
					      "uplo", "L");
				if (*pi == *pj)
					RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
					      "diag", "U");
				++pi;
				++pj;
			}
		}
	}

	return ScalarLogical(1);
}

SEXP xgCMatrix_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	UNPROTECT(2); /* i, x */
	if (XLENGTH(x) != XLENGTH(i))
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "i", "x");
	return ScalarLogical(1);
}

SEXP xsCMatrix_validate(SEXP obj)
{
	SEXP val = xgCMatrix_validate(obj);
	if (TYPEOF(val) != STRSXP)
		val = sCMatrix_validate(obj);
	return val;
}

SEXP xtCMatrix_validate(SEXP obj)
{
	SEXP val = xgCMatrix_validate(obj);
	if (TYPEOF(val) != STRSXP)
		val = tCMatrix_validate(obj);
	return val;
}

SEXP xgRMatrix_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	UNPROTECT(2); /* j, x */
	if (XLENGTH(x) != XLENGTH(j))
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "j", "x");
	return ScalarLogical(1);
}

SEXP xsRMatrix_validate(SEXP obj)
{
	SEXP val = xgRMatrix_validate(obj);
	if (TYPEOF(val) != STRSXP)
		val = sRMatrix_validate(obj);
	return val;
}

SEXP xtRMatrix_validate(SEXP obj)
{
	SEXP val = xgRMatrix_validate(obj);
	if (TYPEOF(val) != STRSXP)
		val = tRMatrix_validate(obj);
	return val;
}

SEXP xgTMatrix_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	UNPROTECT(2); /* i, x */
	if (XLENGTH(x) != XLENGTH(i))
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "i", "x");
	return ScalarLogical(1);
}

SEXP xsTMatrix_validate(SEXP obj)
{
	SEXP val = xgTMatrix_validate(obj);
	if (TYPEOF(val) != STRSXP)
		val = sTMatrix_validate(obj);
	return val;
}

SEXP xtTMatrix_validate(SEXP obj)
{
	SEXP val = xgTMatrix_validate(obj);
	if (TYPEOF(val) != STRSXP)
		val = tTMatrix_validate(obj);
	return val;
}

SEXP dpoMatrix_validate(SEXP obj)
{
	/* NB: Non-finite entries are "valid" because we consider
	   crossprod(x) and tcrossprod(x) to be positive semidefinite
	   even if 'x' contains non-finite entries (for speed) ...
	*/

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int j, n = INTEGER(dim)[0];
	R_xlen_t n1a = (R_xlen_t) n + 1;

	/* Non-negative diagonal elements are necessary but _not_ sufficient */
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	double *px = REAL(x);
	for (j = 0; j < n; ++j, px += n1a)
		if (!ISNAN(*px) && *px < 0.0)
			RMK(_("matrix has negative diagonal elements"));

	return ScalarLogical(1);
}

SEXP dppMatrix_validate(SEXP obj)
{
	/* NB: Non-finite entries are "valid" because we consider
	   crossprod(x) and tcrossprod(x) to be positive semidefinite
	   even if 'x' contains non-finite entries (for speed) ...
	*/

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int j, n = INTEGER(dim)[0];

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = *CHAR(STRING_ELT(uplo, 0));

	/* Non-negative diagonal elements are necessary but _not_ sufficient */
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	double *px = REAL(x);
	if (ul == 'U') {
		for (j = 0; j < n; px += (++j)+1)
			if (!ISNAN(*px) && *px < 0.0)
				RMK(_("matrix has negative diagonal elements"));
	} else {
		for (j = 0; j < n; px += n-(j++))
			if (!ISNAN(*px) && *px < 0.0)
				RMK(_("matrix has negative diagonal elements"));
	}

	return ScalarLogical(1);
}

SEXP corMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int j, n = INTEGER(dim)[0];
	R_xlen_t n1a = (R_xlen_t) n + 1;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	double *px = REAL(x);
	for (j = 0; j < n; ++j, px += n1a)
		if (ISNAN(*px) || *px != 1.0)
			RMK(_("matrix has nonunit diagonal elements"));

	SEXP sd = GET_SLOT(obj, Matrix_sdSym);
	if (TYPEOF(sd) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "sd", "double");
	if (XLENGTH(sd) != n)
		RMKMS(_("'%s' slot does not have length %s"), "sd", "Dim[1]");
	double *psd = REAL(sd);
	for (j = 0; j < n; ++j)
		if (!ISNAN(psd[j]) && psd[j] < 0.0)
			RMKMS(_("'%s' slot has negative elements"), "sd");

	return ScalarLogical(1);
}

SEXP copMatrix_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int j, n = INTEGER(dim)[0];

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = *CHAR(STRING_ELT(uplo, 0));

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	double *px = REAL(x);
	if (ul == 'U') {
		for (j = 0; j < n; px += (++j)+1)
			if (ISNAN(*px) || *px != 1.0)
				RMK(_("matrix has nonunit diagonal elements"));
	} else {
		for (j = 0; j < n; px += n-(j++))
			if (ISNAN(*px) || *px != 1.0)
				RMK(_("matrix has nonunit diagonal elements"));
	}

	SEXP sd = GET_SLOT(obj, Matrix_sdSym);
	if (TYPEOF(sd) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "sd", "double");
	if (XLENGTH(sd) != n)
		RMKMS(_("'%s' slot does not have length %s"), "sd", "Dim[1]");
	double *psd = REAL(sd);
	for (j = 0; j < n; ++j)
		if (!ISNAN(psd[j]) && psd[j] < 0.0)
			RMKMS(_("'%s' slot has negative elements"), "sd");

	return ScalarLogical(1);
}

SEXP sparseVector_validate(SEXP obj)
{
	SEXP length = GET_SLOT(obj, Matrix_lengthSym);
	if (TYPEOF(length) != INTSXP && TYPEOF(length) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "length", "integer", "double");
	if (XLENGTH(length) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "length", 1);
	Matrix_int_fast64_t n;
	if (TYPEOF(length) == INTSXP) {
		int n_ = INTEGER(length)[0];
		if (n_ == NA_INTEGER)
			RMKMS(_("'%s' slot is NA"), "length");
		if (n_ < 0)
			RMKMS(_("'%s' slot is negative"), "length");
		n = (Matrix_int_fast64_t) n_;
	} else {
		double n_ = REAL(length)[0];
		if (ISNAN(n_))
			RMKMS(_("'%s' slot is NA"), "length");
		if (n_ < 0.0)
			RMKMS(_("'%s' slot is negative"), "length");
		if (n_ > 0x1.0p+53)
			RMKMS(_("'%s' slot exceeds %s"), "length", "2^53");
		n = (Matrix_int_fast64_t) n_;
	}

	SEXP i = GET_SLOT(obj, Matrix_iSym);
	if (TYPEOF(i) != INTSXP && TYPEOF(i) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "i", "integer", "double");
	R_xlen_t nnz = XLENGTH(i);
	if (nnz > n)
		RMKMS(_("'%s' slot has length greater than '%s' slot"), "i", "length");
	if (TYPEOF(i) == INTSXP) {
		int *pi = INTEGER(i), max = (n > INT_MAX) ? INT_MAX : (int) n, last = 0;
		while (nnz--) {
			if (*pi == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (*pi < 1 || *pi > max)
				RMKMS(_("'%s' slot has elements not in {%s}"),
					  "i", "1,...,length");
			if (*pi <= last)
				RMKMS(_("'%s' slot is not increasing"), "i");
			last = *(pi++);
		}
	} else {
		double *pi = REAL(i), max = (double) n, last = 0.0, tmp;
		while (nnz--) {
			if (ISNAN(*pi))
				RMKMS(_("'%s' slot contains NA"), "i");
			tmp = trunc(*(pi++));
			if (tmp < 1.0 || tmp > max)
				RMKMS(_("'%s' slot has elements not in {%s} after truncation towards zero"),
					  "i", "1,...,length");
			if (tmp <= last)
				RMKMS(_("'%s' slot is not increasing after truncation towards zero"), "i");
			last = tmp;
		}
	}

	return ScalarLogical(1);
}

#define KINDVECTOR_VALIDATE(_PREFIX_, _SEXPTYPE_) \
SEXP _PREFIX_ ## sparseVector_validate(SEXP obj) \
{ \
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)), \
		i = PROTECT(GET_SLOT(obj, Matrix_iSym)); \
	UNPROTECT(2); /* i, x */ \
	if (TYPEOF(x) != _SEXPTYPE_) \
		RMKMS(_("'%s' slot is not of type \"%s\""), "x", type2char(_SEXPTYPE_)); \
	if (XLENGTH(x) != XLENGTH(i)) \
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "i", "x"); \
	return ScalarLogical(1); \
}
KINDVECTOR_VALIDATE(l,  LGLSXP)
KINDVECTOR_VALIDATE(i,  INTSXP)
KINDVECTOR_VALIDATE(d, REALSXP)
KINDVECTOR_VALIDATE(z, CPLXSXP)
#undef KINDVECTOR_VALIDATE

SEXP denseLU_validate(SEXP obj)
{
	/* In R, we start by checking that 'obj' would be a valid dgeMatrix */

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != r)
		RMKMS(_("'%s' slot does not have length %s"), "perm", "min(Dim)");
	int *pperm = INTEGER(perm);
	while (r--) {
		if (*pperm == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "perm");
		if (*pperm < 1 || *pperm > m)
			RMKMS(_("'%s' slot has elements not in {%s}"),
			      "perm", "1,...,Dim[1]");
		++pperm;
	}

	return ScalarLogical(1);
}

SEXP sparseLU_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym), uplo, diag;
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP L = PROTECT(GET_SLOT(obj, Matrix_LSym));
	PROTECT(dim = GET_SLOT(L, Matrix_DimSym));
	PROTECT(uplo = GET_SLOT(L, Matrix_uploSym));
	PROTECT(diag = GET_SLOT(L, Matrix_diagSym));
	UNPROTECT(4); /* diag, uplo, dim, L */

	pdim = INTEGER(dim);
	if (pdim[0] != n || pdim[1] != n)
		RMKMS(_("dimensions of '%s' slot are not identical to '%s'"), "L", "Dim");
	if (*CHAR(STRING_ELT(uplo, 0)) == 'U')
		RMKMS(_("'%s' slot is upper (not lower) triangular"), "L");
	if (*CHAR(STRING_ELT(diag, 0)) == 'N') {
		PROTECT(L);
		SEXP p = PROTECT(GET_SLOT(L, Matrix_pSym)),
			i = PROTECT(GET_SLOT(L, Matrix_iSym)),
			x = PROTECT(GET_SLOT(L, Matrix_xSym));
		UNPROTECT(4); /* x, i, p, L */

		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
		double *px = REAL(x);
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j + 1];
			if (kend == k || pi[k] != j || px[k] != 1.0)
				RMKMS(_("'%s' slot has nonunit diagonal elements"), "L");
			k = kend;
		}
	}

	SEXP U = PROTECT(GET_SLOT(obj, Matrix_USym));
	PROTECT(dim = GET_SLOT(U, Matrix_DimSym));
	PROTECT(uplo = GET_SLOT(U, Matrix_uploSym));
	UNPROTECT(3); /* uplo, dim, U */

	pdim = INTEGER(dim);
	if (pdim[0] != n || pdim[1] != n)
		RMKMS(_("dimensions of '%s' slot are not identical to '%s'"), "U", "Dim");
	if (*CHAR(STRING_ELT(uplo, 0)) != 'U')
		RMKMS(_("'%s' slot is lower (not upper) triangular"), "U");

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		q = PROTECT(GET_SLOT(obj, Matrix_qSym));
	UNPROTECT(2); /* q, p */
	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (TYPEOF(q) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "q", "integer");
	if (XLENGTH(p) != n)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[1]");
	if (XLENGTH(q) != n && XLENGTH(q) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"), "q", "Dim[2]", "0");
	char *work;
	int lwork = n;
	Matrix_Calloc(work, lwork, char);
	int j, *pp = INTEGER(p);
	for (j = 0; j < n; ++j) {
		if (*pp == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "p");
		if (*pp < 0 || *pp >= n)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			       "p", "0,...,Dim[1]-1");
		if (work[*pp])
			FRMKMS(_("'%s' slot contains duplicates"), "p");
		work[*(pp++)] = 1;
	}
	if (LENGTH(q) == n) {
	int *pq = INTEGER(q);
	for (j = 0; j < n; ++j) {
		if (*pq == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "q");
		if (*pq < 0 || *pq >= n)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			       "q", "0,...,Dim[2]-1");
		if (!work[*pq])
			FRMKMS(_("'%s' slot contains duplicates"), "q");
		work[*(pq++)] = 0;
	}
	}
	Matrix_Free(work, lwork);

	return ScalarLogical(1);
}

SEXP sparseQR_validate(SEXP obj)
{
	/* MJ: assuming for simplicity that 'V' and 'R' slots are formally valid */

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m < n)
		RMK(_("matrix has more columns than rows"));

	SEXP beta = GET_SLOT(obj, Matrix_betaSym);
	if (TYPEOF(beta) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "beta", "double");
	if (XLENGTH(beta) != n)
		RMKMS(_("'%s' slot does not have length %s"), "beta", "Dim[2]");

	SEXP p, i, q;
	int *pp, *pi, *pq, j, k, kend;

	SEXP V = PROTECT(GET_SLOT(obj, Matrix_VSym));
	PROTECT(dim = GET_SLOT(V, Matrix_DimSym));
	PROTECT(p = GET_SLOT(V, Matrix_pSym));
	PROTECT(i = GET_SLOT(V, Matrix_iSym));
	UNPROTECT(4); /* i, p, dim, V */

	pdim = INTEGER(dim);
	int m0 = pdim[0];
	if (m0 < m)
		RMKMS(_("'%s' slot has fewer than %s rows"), "V", "Dim[1]");
	if (m0 > m + n)
		RMKMS(_("'%s' slot has more than %s rows"), "V", "Dim[1]+Dim[2]");
	if (pdim[1] != n)
		RMKMS(_("'%s' slot does not have %s columns"), "V", "Dim[2]");
	pp = INTEGER(p);
	pi = INTEGER(i);
	for (j = 0, k = 0; j < n; ++j) {
		kend = pp[j + 1];
		if (k < kend) {
			if (pi[k] < j)
				RMKMS(_("'%s' slot must be lower trapezoidal but has entries above the diagonal"), "V");
		}
		k = kend;
	}

	SEXP R = PROTECT(GET_SLOT(obj, Matrix_RSym));
	PROTECT(dim = GET_SLOT(R, Matrix_DimSym));
	PROTECT(p = GET_SLOT(R, Matrix_pSym));
	PROTECT(i = GET_SLOT(R, Matrix_iSym));
	UNPROTECT(4); /* i, p, dim, R */

	pdim = INTEGER(dim);
	if (pdim[0] != m0)
		RMKMS(_("'%s' slot does not have %s row"), "R", "nrow(V)");
	if (pdim[1] != n)
		RMKMS(_("'%s' slot does not have %s columns"), "R", "Dim[2]");
	pp = INTEGER(p);
	pi = INTEGER(i);
	for (j = 0, k = 0; j < n; ++j) {
		kend = pp[j + 1];
		if (k < kend) {
			if (pi[kend - 1] > j)
				RMKMS(_("'%s' slot must be upper trapezoidal but has entries below the diagonal"), "R");
#if 0 /* cs_house imposes diag(R) >= 0 in CSparse but not in CXSparse */
			if (pi[kend - 1] == j &&
			    !ISNAN(px[kend - 1]) && px[kend - 1] < 0.0)
				RMKMS(_("'%s' slot has negative diagonal elements"), "R");
#endif
		}
		k = kend;
	}

	PROTECT(p = GET_SLOT(obj, Matrix_pSym));
	PROTECT(q = GET_SLOT(obj, Matrix_qSym));
	UNPROTECT(2); /* q, p */
	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (TYPEOF(q) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "q", "integer");
	if (XLENGTH(p) != m0)
		RMKMS(_("'%s' slot does not have length %s"), "p", "nrow(V)");
	if (XLENGTH(q) != n && XLENGTH(q) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"), "q", "Dim[2]", "0");
	char *work;
	int lwork = m0; /* n <= m <= m0 */
	Matrix_Calloc(work, lwork, char);
	pp = INTEGER(p);
	for (j = 0; j < m0; ++j) {
		if (*pp == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "p");
		if (*pp < 0 || *pp >= m0)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			       "p", "0,...,nrow(V)-1");
		if (work[*pp])
			FRMKMS(_("'%s' slot contains duplicates"), "p");
		work[*(pp++)] = 1;
	}
	if (LENGTH(q) == n) {
	pq = INTEGER(q);
	for (j = 0; j < n; ++j) {
		if (*pq == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "q");
		if (*pq < 0 || *pq >= n)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			       "q", "0,...,Dim[2]-1");
		if (!work[*pq])
			FRMKMS(_("'%s' slot contains duplicates"), "q");
		work[*(pq++)] = 0;
	}
	}
	Matrix_Free(work, lwork);

	return ScalarLogical(1);
}

SEXP BunchKaufman_validate(SEXP obj)
{
	/* In R, we start by checking that 'obj' would be a valid dtrMatrix */

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != n)
		RMKMS(_("'%s' slot does not have length %s"), "perm", "Dim[1]");
	int n_ = n, *pperm = INTEGER(perm);
	while (n_) {
		if (*pperm == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "perm");
		if (*pperm < -n || *pperm == 0 || *pperm > n)
			RMKMS(_("'%s' slot has elements not in {%s}\\{%s}"),
			      "perm", "-Dim[1],...,Dim[1]", "0");
		if (*pperm > 0) {
			pperm += 1;
			n_ -= 1;
		} else if (n_ > 1 && *(pperm + 1) == *pperm) {
			pperm += 2;
			n_ -= 2;
		} else
			RMKMS(_("'%s' slot has unpaired negative elements"), "perm");
	}

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

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int j, n = INTEGER(dim)[0];
	R_xlen_t n1a = (R_xlen_t) n + 1;

	/* Non-negative diagonal elements are necessary _and_ sufficient */
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	double *px = REAL(x);
	for (j = 0; j < n; ++j, px += n1a)
		if (!ISNAN(*px) && *px < 0.0)
			RMK(_("Cholesky factor has negative diagonal elements"));

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != n && XLENGTH(perm) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"), "perm", "Dim[1]", "0");
	if (LENGTH(perm) == n) {
		char *work;
		int lwork = n;
		Matrix_Calloc(work, lwork, char);
		int *pperm = INTEGER(perm);
		for (j = 0; j < n; ++j) {
			if (*pperm == NA_INTEGER)
				FRMKMS(_("'%s' slot contains NA"), "perm");
			if (*pperm < 0 || *pperm >= n)
				FRMKMS(_("'%s' slot has elements not in {%s}"),
				       "perm", "0,...,Dim[1]-1");
			if (work[*pperm])
				FRMKMS(_("'%s' slot contains duplicates"), "perm");
			work[*(pperm++)] = 1;
		}
		Matrix_Free(work, lwork);
	}

	return ScalarLogical(1);
}

SEXP pCholesky_validate(SEXP obj)
{
	/* In R, we start by checking that 'obj' would be a valid dtpMatrix */

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int j, n = INTEGER(dim)[0];

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = *CHAR(STRING_ELT(uplo, 0));

	/* Non-negative diagonal elements are necessary _and_ sufficient */
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	double *px = REAL(x);
	if (ul == 'U') {
		for (j = 0; j < n; px += (++j)+1)
			if (!ISNAN(*px) && *px < 0.0)
				RMK(_("Cholesky factor has negative diagonal elements"));
	} else {
		for (j = 0; j < n; px += n-(j++))
			if (!ISNAN(*px) && *px < 0.0)
				RMK(_("Cholesky factor has negative diagonal elements"));
	}

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != n && XLENGTH(perm) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"), "perm", "Dim[1]", "0");
	if (LENGTH(perm) == n) {
		char *work;
		int lwork = n;
		Matrix_Calloc(work, lwork, char);
		int *pperm = INTEGER(perm);
		for (j = 0; j < n; ++j) {
			if (*pperm == NA_INTEGER)
				FRMKMS(_("'%s' slot contains NA"), "perm");
			if (*pperm < 0 || *pperm >= n)
				FRMKMS(_("'%s' slot has elements not in {%s}"),
				       "perm", "0,...,Dim[1]-1");
			if (work[*pperm])
				FRMKMS(_("'%s' slot contains duplicates"), "perm");
			work[*(pperm++)] = 1;
		}
		Matrix_Free(work, lwork);
	}

	return ScalarLogical(1);
}

SEXP CHMfactor_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP type = GET_SLOT(obj, install("type"));
	if (TYPEOF(type) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "type", "integer");
	if (XLENGTH(type) != 6)
		RMKMS(_("'%s' slot does not have length %d"), "type", 6);
	int order = INTEGER(type)[0];
	if (order < 0 || order > 4)
		RMKMS(_("%s[%d] (%s) is not in %s"),
		      "type", 1, "cholmod_factor.ordering", "0:4");

	SEXP colcount = GET_SLOT(obj, install("colcount"));
	if (TYPEOF(colcount) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "colcount", "integer");
	if (XLENGTH(colcount) != n)
		RMKMS(_("'%s' slot does not have length %s"), "colcount", "Dim[2]");
	int j, *pcolcount = INTEGER(colcount);
	for (j = 0; j < n; ++j) {
		if (pcolcount[j] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "colcount");
		if (pcolcount[j] < 0 || pcolcount[j] > n - j)
			RMKMS(_("%s is not in {%s}"), "colcount[j]", "0,...,Dim[1]-j+1");
	}

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (order == 0) {
		if (XLENGTH(perm) != 0)
			RMKMS(_("'%s' slot does not have length %d"), "perm", 0);
	} else {
		if (XLENGTH(perm) != n)
			RMKMS(_("'%s' slot does not have length %s"), "perm", "Dim[1]");
		char *work;
		int lwork = n;
		Matrix_Calloc(work, lwork, char);
		int *pperm = INTEGER(perm);
		for (j = 0; j < n; ++j) {
			if (*pperm == NA_INTEGER)
				FRMKMS(_("'%s' slot contains NA"), "perm");
			if (*pperm < 0 || *pperm >= n)
				FRMKMS(_("'%s' slot has elements not in {%s}"),
				       "perm", "0,...,Dim[1]-1");
			if (work[*pperm])
				FRMKMS(_("'%s' slot contains duplicates"), "perm");
			work[*(pperm++)] = 1;
		}
		Matrix_Free(work, lwork);
	}

	return ScalarLogical(1);
}

SEXP CHMsimpl_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == INT_MAX)
		RMKMS(_("%s is not representable as \"%s\""), "Dim[1]+1", "integer");

	SEXP type = GET_SLOT(obj, install("type"));
	int *ptype = INTEGER(type), mono = ptype[3];
	if (ptype[1] != 0 && ptype[1] != 1)
		RMKMS(_("%s[%d] (%s) is not %d or %d"),
		      "type", 2, "cholmod_factor.is_ll", 0, 1);
	if (ptype[2] != 0)
		RMKMS(_("%s[%d] (%s) is not %d"),
		      "type", 3, "cholmod_factor.is_super", 0);
	if (ptype[3] != 0 && ptype[3] != 1)
		RMKMS(_("%s[%d] (%s) is not %d or %d"),
		      "type", 4, "cholmod_factor.is_monotonic", 0, 1);

	SEXP nxt = PROTECT(GET_SLOT(obj, install("nxt"))),
		prv = PROTECT(GET_SLOT(obj, install("prv"))),
		nz = PROTECT(GET_SLOT(obj, install("nz"))),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	UNPROTECT(5); /* i, p, nz, prv, nxt */

	if (TYPEOF(nxt) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "nxt", "integer");
	if (TYPEOF(prv) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "prv", "integer");
	if (XLENGTH(nxt) - 2 != n)
		RMKMS(_("'%s' slot does not have length %s"), "nxt", "Dim[2]+2");
	if (XLENGTH(prv) - 2 != n)
		RMKMS(_("'%s' slot does not have length %s"), "prv", "Dim[2]+2");
	int *pnxt = INTEGER(nxt), *pprv = INTEGER(prv),
		j1 = pnxt[n + 1], j2 = pprv[n], count = n + 1;
	while (count--) {
		if (j1 < 0 || j1 > n)
			RMKMS(_("%s has elements not in {%s}"),
			      "nxt[-(Dim[2]+1)]", "0,...,Dim[2]");
		if (j2 < 0 || j2 > n + 1 || j2 == n)
			RMKMS(_("%s has elements not in {%s}\\{%s}"),
			      "prv[-(Dim[2]+2)]", "0,...,Dim[2]+1", "Dim[2]");
		if ((count >  1) && mono && (pnxt[j1] != j1 + 1 || pprv[j2] != j2 - 1))
			RMKMS(_("%s is %d but columns are not stored in increasing order"),
			      "type[4]", 1);
		if ((count >= 1) ? j1 == n : j1 != n)
			RMKMS(_("traversal of '%s' slot does not complete in exactly %s steps"),
			      "nxt", "length(nxt)");
		if ((count >= 1) ? j2 == n + 1 : j2 != n + 1)
			RMKMS(_("traversal of '%s' slot does not complete in exactly %s steps"),
			      "prv", "length(prv)");
		j1 = pnxt[j1];
		j2 = pprv[j2];
	}
	if (j1 != -1)
		RMKMS(_("%s is not %d"), "nxt[Dim[2]+1]", -1);
	if (j2 != -1)
		RMKMS(_("%s is not %d"), "prv[Dim[2]+2]", -1);

	if (TYPEOF(nz) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "nz", "integer");
	if (XLENGTH(nz) != n)
		RMKMS(_("'%s' slot does not have length %s"), "nz", "Dim[2]");
	int j, *pnz = INTEGER(nz);
	for (j = 0; j < n; ++j) {
		if (pnz[j] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "nz");
		if (pnz[j] < 1 || pnz[j] > n - j)
			RMKMS(_("%s is not in {%s}"), "nz[j]", "1,...,Dim[1]-j+1");
	}

	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (XLENGTH(p) - 1 != n)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[2]+1");
	j1 = pnxt[n + 1];
	int *pp = INTEGER(p);
	if (pp[j1] != 0)
		RMKMS(_("column '%s' is stored first but %s is not 0"), "j", "p[j]");
	for (j = 0; j < n; ++j) {
		j2 = pnxt[j1];
		if (pp[j2] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "p");
		if (pp[j2] < pp[j1])
			RMKMS(_("'%s' slot is not increasing when traversed in stored column order"), "p");
		if (pp[j2] - pp[j1] < pnz[j1])
			RMKMS(_("'%s' slot allocates fewer than %s elements for column '%s'"),
			      "i", "nz[j]", "j");
		if (pp[j2] - pp[j1] > n - j1)
			RMKMS(_("'%s' slot allocates more than %s elements for column '%s'"),
			      "i", "Dim[2]-j+1", "j");
		j1 = j2;
	}

	if (TYPEOF(i) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "i", "integer");
	if (XLENGTH(i) != pp[n])
		RMKMS(_("'%s' slot does not have length %s"), "i", "p[length(p)]");
	int *pi = INTEGER(i), *pi_, k;
	j1 = pnxt[n + 1];
	for (j = 0; j < n; ++j) {
		pi_ = pi + pp[j1];
		if (pi_[0] != j1)
			RMKMS(_("first entry in column '%s' does not have row index '%s'"),
			      "j", "j");
		for (k = 1; k < pnz[j1]; ++k) {
			if (pi_[k] == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (pi_[k] < 0 || pi_[k] >= n)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "i", "0,...,Dim[1]-1");
			if (pi_[k] <= pi_[k - 1])
				RMKMS(_("'%s' slot is not increasing within columns"), "i");
		}
		j1 = pnxt[j1];
	}

	return ScalarLogical(1);
}

SEXP CHMsuper_validate(SEXP obj)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];

	SEXP type = GET_SLOT(obj, install("type"));
	int *ptype = INTEGER(type);
	if (ptype[1] != 1)
		RMKMS(_("%s[%d] (%s) is not %d"),
		      "type", 2, "cholmod_factor.is_ll", 1);
	if (ptype[2] != 1)
		RMKMS(_("%s[%d] (%s) is not %d"),
		      "type", 3, "cholmod_factor.is_super", 1);
	if (ptype[3] != 1)
		RMKMS(_("%s[%d] (%s) is not %d"),
		      "type", 4, "cholmod_factor.is_monotonic", 1);
	if (ptype[4] < 0)
		RMKMS(_("%s[%d] (%s) is negative"),
		      "type", 5, "cholmod_factor.maxcsize");
	if (ptype[5] < 0)
		RMKMS(_("%s[%d] (%s) is negative"),
		      "type", 6, "cholmod_factor.maxesize");
	if (n > 0 && ptype[5] >= n)
		RMKMS(_("%s[%d] (%s) is not less than %s"),
		      "type", 6, "cholmod_factor.maxesize", "Dim[1]");

	/* FIXME: maxcsize and maxesize are well-defined properties of the
	   factorization, so we should also test that the values are
	   _correct_ ... see CHOLMOD/Supernodal/cholmod_super_symbolic.c
	*/

	SEXP super = PROTECT(GET_SLOT(obj, install("super"))),
		pi = PROTECT(GET_SLOT(obj, install("pi"))),
		px = PROTECT(GET_SLOT(obj, install("px"))),
		s = PROTECT(GET_SLOT(obj, install("s")));
	UNPROTECT(4); /* s, px, pi, super */

	if (TYPEOF(super) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "super", "integer");
	R_xlen_t nsuper1a = XLENGTH(super);
	if (nsuper1a - 1 < ((n > 0) ? 1 : 0))
		RMKMS(_("'%s' slot has length less than %d"), "super", 2);
	if (nsuper1a - 1 > n)
		RMKMS(_("'%s' slot has length greater than %s"), "super", "Dim[2]+1");
	int k, nsuper = (int) (nsuper1a - 1), *psuper = INTEGER(super);
	if (psuper[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "super");
	if (psuper[nsuper] != n)
		RMKMS(_("last element of '%s' slot is not %s"), "super", "Dim[2]");
	for (k = 1; k <= nsuper; ++k) {
		if (psuper[k] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "super");
		if (psuper[k] <= psuper[k - 1])
			RMKMS(_("'%s' slot is not increasing"), "super");
	}

	if (TYPEOF(pi) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "pi", "integer");
	if (TYPEOF(px) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "px", "integer");
	if (XLENGTH(pi) != nsuper1a)
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "pi", "super");
	if (XLENGTH(px) != nsuper1a)
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "px", "super");
	int *ppi = INTEGER(pi), *ppx = INTEGER(px), nr, nc;
	if (ppi[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "pi");
	if (ppx[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "px");
	for (k = 1; k <= nsuper; ++k) {
		if (ppi[k] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "pi");
		if (ppx[k] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "px");
		if (ppi[k] <= ppi[k - 1])
			RMKMS(_("'%s' slot is not increasing"), "pi");
		if (ppx[k] <= ppx[k - 1])
			RMKMS(_("'%s' slot is not increasing"), "px");
		nr = ppi[k] - ppi[k - 1];
		nc = psuper[k] - psuper[k - 1];
		if (nr < nc)
			RMKMS(_("first differences of '%s' slot are less than those of '%s' slot"),
			      "pi", "super");
		if ((Matrix_int_fast64_t) nr * nc > INT_MAX)
			RMKMS(_("supernode lengths exceed %s"), "2^31-1");
		if (ppx[k] - ppx[k - 1] != nr * nc)
			RMKMS(_("first differences of '%s' slot are not equal to supernode lengths"),
			      "px");
	}

	if (TYPEOF(s) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "s", "integer");
	if (XLENGTH(s) != ppi[nsuper])
		RMKMS(_("'%s' slot does not have length %s"), "s", "pi[length(pi)]");
	int i, j, *ps = INTEGER(s);
	for (k = 1; k <= nsuper; ++k) {
		nr = ppi[k] - ppi[k-1];
		nc = psuper[k] - (j = psuper[k-1]);
		for (i = 0; i < nr; ++i) {
			if (ps[i] == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "s");
			if (ps[i] < 0 || ps[i] >= n)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "s", "0,...,Dim[1]-1");
			if (i < nc) {
				if (ps[i] != j + i)
					RMKMS(_("'%s' slot is wrong within diagonal blocks (row and column indices do not coincide)"), "s");
			} else {
				if (ps[i] <= ps[i-1])
					RMKMS(_("'%s' slot is not increasing within supernodes"), "s");
			}
		}
		ps += nr;
	}

	return ScalarLogical(1);
}

SEXP dCHMsimpl_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		type = PROTECT(GET_SLOT(obj, install("type")));
	UNPROTECT(3); /* type, p, x */

	if (TYPEOF(x) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "x", "double");

	int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
	if (XLENGTH(x) != pp[n])
		RMKMS(_("'%s' slot does not have length %s"), "x", "p[length(p)]");

	if (INTEGER(type)[1]) {
		int j;
		double *px = REAL(x);

		/* Non-negative diagonal elements are necessary _and_ sufficient */
		for (j = 0; j < n; ++j)
			if (!ISNAN(px[pp[j]]) && px[pp[j]] < 0.0)
				RMK(_("Cholesky factor has negative diagonal elements"));
	}

	return ScalarLogical(1);
}

SEXP dCHMsuper_validate(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		px = PROTECT(GET_SLOT(obj, install("px"))),
		pi = PROTECT(GET_SLOT(obj, install("pi"))),
		super = PROTECT(GET_SLOT(obj, install("super")));
	UNPROTECT(4); /* super, pi, px, x */

	if (TYPEOF(x) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "x", "double");

	int *ppx = INTEGER(px), nsuper = (int) (XLENGTH(px) - 1);
	if (XLENGTH(x) != ppx[nsuper])
		RMKMS(_("'%s' slot does not have length %s"), "x", "px[length(px)]");

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
				RMK(_("Cholesky factor has negative diagonal elements"));
			pv += nr1a;
		}
	}

	return ScalarLogical(1);
}

SEXP Schur_validate(SEXP obj)
{
	/* MJ: assuming for simplicity that 'Q' and 'T' slots are formally valid */

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP Q = PROTECT(GET_SLOT(obj, Matrix_QSym));
	dim = GET_SLOT(Q, Matrix_DimSym);
	pdim = INTEGER(dim);
	UNPROTECT(1); /* Q */
	if (pdim[0] != n || pdim[1] != n)
		RMKMS(_("dimensions of '%s' slot are not identical to '%s'"), "Q", "Dim");

	SEXP T = PROTECT(GET_SLOT(obj, Matrix_TSym));
	dim = GET_SLOT(T, Matrix_DimSym);
	pdim = INTEGER(dim);
	UNPROTECT(1); /* T */
	if (pdim[0] != n || pdim[1] != n)
		RMKMS(_("dimensions of '%s' slot are not identical to '%s'"), "T", "Dim");

	SEXP v = GET_SLOT(obj, install("EValues"));
	SEXPTYPE tv = TYPEOF(v);
	if (tv != REALSXP && tv != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "EValues", "double", "complex");
	if (XLENGTH(v) != n)
		RMKMS(_("'%s' slot does not have length %s"), "EValues", "Dim[1]");

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
			if (cl[1] == 's') \
				IS_VALID(s ## _C_ ## Matrix); \
			else if (cl[1] == 't') \
				IS_VALID(t ## _C_ ## Matrix); \
		} else { \
			if (cl[1] == 'g') \
				IS_VALID(xg ## _C_ ## Matrix); \
			else if (cl[1] == 's') \
				IS_VALID(xs ## _C_ ## Matrix); \
			else if (cl[1] == 't') \
				IS_VALID(xt ## _C_ ## Matrix); \
		} \
	} while (0)

	IS_VALID(Matrix);

	if ((cl[0] == 'i' && cl[1] == 'n' && cl[2] == 'd') ||
		(cl[0] == 'p' && cl[1] != 'c')) {
		IS_VALID(indMatrix);
		if (cl[0] == 'p')
			IS_VALID(pMatrix);
		return;
	}

	const char *cl_ = cl;
	if (cl[0] == 'c')
		cl = "dpoMatrix";
	else if (cl[0] == 'p' && cl[1] == 'c')
		cl = "dppMatrix";

	if (cl[0] == 'n' && cl[2] != 'C' && cl[2] != 'R' && cl[2] != 'T')
		IS_VALID(nMatrix);
	else if (cl[0] == 'l')
		IS_VALID(lMatrix);
	else if (cl[0] == 'i')
		IS_VALID(iMatrix);
	else if (cl[0] == 'd')
		IS_VALID(dMatrix);
	else if (cl[0] == 'z')
		IS_VALID(zMatrix);

	if (cl[1] == 's' || cl[1] == 'p')
		IS_VALID(symmetricMatrix);
	else if (cl[1] == 't')
		IS_VALID(triangularMatrix);
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
		if (cl[1] == 'p') {
			IS_VALID(dpoMatrix);
			if (cl_[0] == 'c')
				IS_VALID(corMatrix);
		}
	} else {
		IS_VALID(packedMatrix);
		if (cl[1] == 'p') {
			IS_VALID(dppMatrix);
			if (cl_[0] == 'c')
				IS_VALID(copMatrix);
		}
	}

# undef IS_VALID_SPARSE
# undef IS_VALID

#endif

	return;
}
