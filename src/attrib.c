#include "Mdefines.h"
#include "attrib.h"

/* .... Dimnames .................................................... */

int DimNames_is_trivial(SEXP dn)
{
	return
		isNull(VECTOR_ELT(dn, 0)) &&
		isNull(VECTOR_ELT(dn, 1)) &&
		isNull(getAttrib(dn, R_NamesSymbol));
}

int DimNames_is_symmetric(SEXP dn)
{
	SEXP rn, cn, ndn;
	const char *nrn, *ncn;
	int n;

	return
		!((!isNull(rn = VECTOR_ELT(dn, 0)) &&
		   !isNull(cn = VECTOR_ELT(dn, 1)) &&
		   rn != cn &&
		   ((n = LENGTH(rn)) != LENGTH(cn) ||
		    !equal_character_vectors(rn, cn, n))) ||
		  ((!isNull(ndn = getAttrib(dn, R_NamesSymbol)) &&
		    *(nrn = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
		    *(ncn = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
		    strcmp(nrn, ncn) != 0)));
}

SEXP R_DimNames_is_symmetric(SEXP dn)
{
	return ScalarLogical(DimNames_is_symmetric(dn));
}

void symDN(SEXP dest, SEXP src, int J /* -1|0|1 */)
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

SEXP R_symDN(SEXP dn)
{
	if (DimNames_is_trivial(dn))
		return dn;
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	symDN(newdn, dn, -1);
	UNPROTECT(1);
	return newdn;
}

SEXP R_revDN(SEXP dn)
{
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
	symDN(newdn, dn, J);
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
		symDN(newdn, dn, J);
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


/* .... factors ..................................................... */

static
int strmatch(const char *x, SEXP valid)
{
	int i, n = LENGTH(valid);
	for (i = 0; i < n; ++i)
		if (strcmp(x, CHAR(STRING_ELT(valid, i))) == 0)
			return i;
	return -1;
}

static
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

SEXP get_factor(SEXP obj, const char *nm)
{
	SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorsSym)), val = R_NilValue;
	if (LENGTH(factors) > 0) {
		SEXP valid = PROTECT(getAttrib(factors, R_NamesSymbol));
		int i = strmatch(nm, valid);
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
	PROTECT_WITH_INDEX(factors = GET_SLOT(obj, Matrix_factorsSym), &pid);
	if (LENGTH(factors) > 0) {
		SEXP valid = PROTECT(getAttrib(factors, R_NamesSymbol));
		int i = strmatch(nm, valid);
		UNPROTECT(1);
		if (i >= 0) {
			SET_VECTOR_ELT(factors, i, val);
			UNPROTECT(2);
			return;
		}
	}
	REPROTECT(factors = append_to_named_list(factors, nm, val), pid);
	SET_SLOT(obj, Matrix_factorsSym, factors);
	UNPROTECT(2);
	return;
}

SEXP R_set_factor(SEXP obj, SEXP nm, SEXP val, SEXP warn)
{
	if (TYPEOF(nm) != STRSXP || LENGTH(nm) < 1 ||
	    (nm = STRING_ELT(nm, 0)) == NA_STRING)
		error(_("invalid factor name"));
	else if (TYPEOF(getAttrib(obj, Matrix_factorsSym)) == VECSXP)
		set_factor(obj, CHAR(nm), val);
	else if (asLogical(warn) != 0)
		warning(_("attempt to set factor on %s without '%s' slot"),
		        "Matrix", "factors");
	return val;
}
