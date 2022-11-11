#include <ctype.h> /* toupper */
#include "Mutils.h"

/**
 * A safe version of `NEW_OBJECT(MAKE_CLASS(what))`, protecting the
 * intermediate R object.  The caller must protect the return value
 * of this function.
 *
 * @param A string specifying the name of a defined S4 class.
 */
SEXP NEW_OBJECT_OF_CLASS(const char* what)
{
    SEXP class = PROTECT(MAKE_CLASS(what)), obj = NEW_OBJECT(class);
    UNPROTECT(1);
    return obj;
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
    SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorSym));
    if (LENGTH(factors) > 0) {
	SEXP valid = PROTECT(getAttrib(factors, R_NamesSymbol));
	int i = strmatch2(nm, valid);
	UNPROTECT(1);
	if (i >= 0) {
	    factors = VECTOR_ELT(factors, i);
	    UNPROTECT(1);
	    return factors;
	}
    }
    UNPROTECT(1);
    return R_NilValue;
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
	warning(_("attempt to set factor on Matrix "
		  "without 'factors' slot"));
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
	warning(_("attempt to discard factors from Matrix "
		  "without 'factors' slot"));
    return ScalarLogical(0); /* no-op */
}


/* For inheritance ================================================== */

char type2kind(SEXPTYPE type)
{
    switch (type) {
    case LGLSXP:
	return 'l';
    case INTSXP:
#ifdef HAVE_PROPER_IMATRIX
	return 'i';
#endif
    case REALSXP:
	return 'd';
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
	return 'z';
#endif
    default:
	error(_("unexpected type \"%s\" in 'type2kind()'"), type2char(type)); 
	return '\0';
    }
}

SEXPTYPE kind2type(char kind)
{
    switch (kind) {
    case 'n':
    case 'l':
	return LGLSXP;
#ifdef HAVE_PROPER_IMATRIX
    case 'i':
	return INTSXP;
#endif
    case 'd':
	return REALSXP;
#ifdef HAVE_PROPER_ZMATRIX
    case 'z':
	return CPLXSXP;
#endif
    default:
	error(_("unexpected kind \"%c\" in 'kind2type()'"), kind);
	return NILSXP;
    }
}

size_t kind2size(char kind)
{
    switch (kind) {
    case 'n':
    case 'l':
#ifdef HAVE_PROPER_IMATRIX
    case 'i':
#endif
	return sizeof(int);
    case 'd':
	return sizeof(double);
#ifdef HAVE_PROPER_ZMATRIX
    case 'z':
	return sizeof(Rcomplex);
#endif
    default:
	error(_("unexpected kind \"%c\" in 'kind2size()'"), kind);
	return 0;
    }
}

/* TODO: compare with macros in ./Mdefines.h */

#define VALID_NONVIRTUAL						\
/*  0 */ "dgCMatrix", "dgRMatrix", "dgTMatrix", "dgeMatrix",		\
/*  4 */ "dsCMatrix", "dsRMatrix", "dsTMatrix", "dspMatrix", "dsyMatrix", \
/*  9 */ "dtCMatrix", "dtRMatrix", "dtTMatrix", "dtpMatrix", "dtrMatrix", \
/* 14 */ "ddiMatrix", "dsparseVector",					\
/* 16 */ "lgCMatrix", "lgRMatrix", "lgTMatrix", "lgeMatrix",		\
/* 20 */ "lsCMatrix", "lsRMatrix", "lsTMatrix", "lspMatrix", "lsyMatrix", \
/* 25 */ "ltCMatrix", "ltRMatrix", "ltTMatrix", "ltpMatrix", "ltrMatrix", \
/* 30 */ "ldiMatrix", "lsparseVector",					\
/* 32 */ "ngCMatrix", "ngRMatrix", "ngTMatrix", "ngeMatrix",		\
/* 36 */ "nsCMatrix", "nsRMatrix", "nsTMatrix", "nspMatrix", "nsyMatrix", \
/* 41 */ "ntCMatrix", "ntRMatrix", "ntTMatrix", "ntpMatrix", "ntrMatrix", \
/* 46 */ /* "ndiMatrix", */ "nsparseVector",				\
/* 47 */ "igCMatrix", "igRMatrix", "igTMatrix", "igeMatrix",		\
/* 51 */ "isCMatrix", "isRMatrix", "isTMatrix", "ispMatrix", "isyMatrix", \
/* 56 */ "itCMatrix", "itRMatrix", "itTMatrix", "itpMatrix", "itrMatrix", \
/* 61 */ "idiMatrix", "isparseVector",					\
/* 63 */ "zgCMatrix", "zgRMatrix", "zgTMatrix", "zgeMatrix",		\
/* 67 */ "zsCMatrix", "zsRMatrix", "zsTMatrix", "zspMatrix", "zsyMatrix", \
/* 72 */ "ztCMatrix", "ztRMatrix", "ztTMatrix", "ztpMatrix", "ztrMatrix", \
/* 77 */ "zdiMatrix", "zsparseVector",					\
/* 79 */ "indMatrix", "pMatrix"

char Matrix_kind(SEXP obj, int i2d)
{
    if (IS_S4_OBJECT(obj)) {
	static const char *valid[] = { VALID_NONVIRTUAL, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
	    error(_("\"kind\" not yet defined for objects of class \"%s\""),
		  class_P(obj));
	if (ivalid >= 79)
	    return 'n'; /* indMatrix, pMatrix */
	return valid[ivalid][0];
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
	    error(_("\"kind\" not yet defined for objects of type \"%s\""),
		  type2char(TYPEOF(obj)));
	    return '\0';
	}
    }
}

SEXP R_Matrix_kind(SEXP obj, SEXP i2d)
{
    char k = Matrix_kind(obj, asLogical(i2d));
    char s[] = { k, '\0' };
    return mkString(s);
}

char Matrix_shape(SEXP obj)
{
    if (!IS_S4_OBJECT(obj))
	error(_("\"shape\" not yet defined for objects of type \"%s\""),
		  type2char(TYPEOF(obj)));
    static const char *valid[] = { VALID_NONVIRTUAL, "" };
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	error(_("\"shape\" not yet defined for objects of class \"%s\""),
	      class_P(obj));
    if (ivalid >= 79 || valid[ivalid][3] != 'M')
	return 'g'; /* indMatrix, pMatrix, .sparseVector */
    return valid[ivalid][1];
}

SEXP R_Matrix_shape(SEXP obj)
{
    char k = Matrix_shape(obj);
    char s[] = { k, '\0' };
    return mkString(s);
}

/* TODO: compare with macros in ./Mdefines.h */

#define VALID_CRTSPARSE				\
"dgCMatrix", "dsCMatrix", "dtCMatrix",		\
"lgCMatrix", "lsCMatrix", "ltCMatrix",		\
"ngCMatrix", "nsCMatrix", "ntCMatrix",		\
"dgRMatrix", "dsRMatrix", "dtRMatrix",		\
"lgRMatrix", "lsRMatrix", "ltRMatrix",		\
"ngRMatrix", "nsRMatrix", "ntRMatrix",		\
"dgTMatrix", "dsTMatrix", "dtTMatrix",		\
"lgTMatrix", "lsTMatrix", "ltTMatrix",		\
"ngTMatrix", "nsTMatrix", "ntTMatrix"

char Matrix_repr(SEXP obj)
{
    if (!IS_S4_OBJECT(obj))
	error(_("\"repr\" not yet defined for objects of type \"%s\""),
	      type2char(TYPEOF(obj)));
    static const char *valid[] = { VALID_CRTSPARSE, "" };
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	return '\0'; /* useful to _not_ throw an error, for inheritance tests */
    return valid[ivalid][2];
}

SEXP R_Matrix_repr(SEXP obj)
{
    char k = Matrix_repr(obj);
    if (k == '\0')
	return mkString("");
    char s[] = { k, '\0' };
    return mkString(s);
}


/* For indexing ===================================================== */

SEXP R_index_triangle(SEXP n_, SEXP upper_, SEXP diag_, SEXP packed_)
{
    int n = asInteger(n_), packed = asLogical(packed_);
    double nn = (double) n * n, nx = (packed) ? nn : 0.5 * (nn + n);
    if (nx > R_XLEN_T_MAX)
	error(_("cannot index a vector of length exceeding R_XLEN_T_MAX"));
    SEXP r;
    int i, j, upper = asLogical(upper_), diag = asLogical(diag_);
    double nr = (diag) ? 0.5 * (nn + n) : 0.5 * (nn - n);
    if (nx > INT_MAX) {
	
	PROTECT(r = allocVector(REALSXP, (R_xlen_t) nr));
	double k = 1.0, *pr = REAL(r);

#define DO_INDEX(_ONE_, _NR_)				\
	do {						\
	    if (packed) {				\
		if (diag) {				\
		    while (k <= _NR_) {			\
			*(pr++) = k;			\
			k += _ONE_;			\
		    }					\
		} else if (upper) {			\
		    for (j = 0; j < n; ++j) {		\
			for (i = 0; i < j; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
			k += _ONE_;			\
		    }					\
		} else {				\
		    for (j = 0; j < n; ++j) {		\
			k += 1.0;			\
			for (i = j+1; i < n; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
		    }					\
		}					\
	    } else if (diag) {				\
		if (upper) {				\
		    for (j = 0; j < n; ++j) {		\
			for (i = 0; i <= j; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
			k += n-j-1;			\
		    }					\
		} else {				\
		    for (j = 0; j < n; ++j) {		\
			k += j;				\
			for (i = j; i < n; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
		    }					\
		}					\
	    } else {					\
		if (upper) {				\
		    for (j = 0; j < n; ++j) {		\
			for (i = 0; i < j; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
			k += n-j;			\
		    }					\
		} else {				\
		    for (j = 0; j < n; ++j) {		\
			k += j+1;			\
			for (i = j+1; i < n; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
		    }					\
		}					\
	    }						\
	} while (0)

	DO_INDEX(1.0, nr);
	
    } else {
	
	PROTECT(r = allocVector(INTSXP, (R_xlen_t) nr));
	int k = 1, nr_ = (int) nr, *pr = INTEGER(r);

	DO_INDEX(1, nr_);

#undef DO_INDEX
	
    }
    
    UNPROTECT(1);
    return r;
}

SEXP R_index_diagonal(SEXP n_, SEXP upper_, SEXP packed_)
{
    int n = asInteger(n_), packed = asLogical(packed_);
    double nn = (double) n * n, nx = (packed) ? nn : 0.5 * (nn + n);
    if (nx > R_XLEN_T_MAX)
	error(_("cannot index a vector of length exceeding R_XLEN_T_MAX"));
    SEXP r;
    int j, upper = (packed) ? asLogical(upper_) : NA_LOGICAL;
    if (nx > INT_MAX) {
	
	PROTECT(r = allocVector(REALSXP, n));
	double k = 1.0, *pr = REAL(r);

#define DO_INDEX				\
	do {					\
	    if (!packed) {			\
		for (j = 0; j < n; ++j) {	\
		    *(pr++) = k;		\
		    k += n+1;			\
		}				\
	    } else if (upper) {			\
		for (j = 0; j < n; ++j) {	\
		    *(pr++) = k;		\
		    k += j+2;			\
		}				\
	    } else {				\
		for (j = 0; j < n; ++j) {	\
		    *(pr++) = k;		\
		    k += n-j;			\
		}				\
	    }					\
	} while (0)

	DO_INDEX;
	
    } else {

	PROTECT(r = allocVector(INTSXP, n));
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

#define DO_NNZ(_CTYPE_, _PTR_, _NA_, _NZ_, _STRICTLY_NZ_)	\
    do {							\
	_CTYPE_ *px = _PTR_(x);					\
	if (do_countNA == NA_LOGICAL) {				\
	    while (n-- > 0) {					\
		if (_NA_(*px))					\
		    return ScalarInteger(NA_INTEGER);		\
		if (_NZ_(*px))					\
		    ++nnz;					\
		++px;						\
	    }							\
	} else if (do_countNA != 0) {				\
	    while (n-- > 0) {					\
		if (_NZ_(*px))					\
		    ++nnz;					\
		++px;						\
	    }							\
	} else {						\
	    while (n-- > 0) {					\
		if (_STRICTLY_NZ_(*px))				\
		    ++nnz;					\
		++px;						\
	    }							\
	}							\
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
	ERROR_INVALID_TYPE("'x'", TYPEOF(x), "R_nnz");
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
	ERROR_INVALID_TYPE("'x'", TYPEOF(x), "na2one");
	break;
    }
    return;
}

SEXP v2spV(SEXP from)
{
    SEXPTYPE tx = TYPEOF(from);
    SEXP to = NULL, length = NULL, i = NULL, x = NULL;
    R_xlen_t n_ = XLENGTH(from);

#define V2SPV(_KIND_, _NZ_,						\
	      _CTYPE1_, _SEXPTYPE1_, _PTR1_,				\
	      _CTYPE2_, _SEXPTYPE2_, _PTR2_)				\
    do {								\
	PROTECT(to = NEW_OBJECT_OF_CLASS(#_KIND_ "sparseVector"));	\
	_CTYPE1_ *py = _PTR1_(from);					\
	for (k = 0; k < n; ++k)						\
	    if (_NZ_(py[k]))						\
		++nnz;							\
	PROTECT(i = allocVector(_SEXPTYPE2_, nnz));			\
	PROTECT(x = allocVector(_SEXPTYPE1_, nnz));			\
	_CTYPE2_ *pi = _PTR2_(i);					\
	_CTYPE1_ *px = _PTR1_(x);					\
	for (k = 0; k < n; ++k) {					\
	    if (_NZ_(py[k])) {						\
		*(pi++) = (_CTYPE2_) (k + 1);				\
		*(px++) = py[k];					\
	    }								\
	}								\
    } while (0)

#define V2SPV_CASES(_CTYPE2_, _SEXPTYPE2_, _PTR2_)			\
    do {								\
	switch (tx) {							\
	case LGLSXP:							\
	    V2SPV(l, ISNZ_LOGICAL, int, LGLSXP, LOGICAL,		\
		  _CTYPE2_, _SEXPTYPE2_, _PTR2_);			\
	    break;							\
	case INTSXP:							\
	    V2SPV(i, ISNZ_INTEGER, int, INTSXP, INTEGER,		\
		  _CTYPE2_, _SEXPTYPE2_, _PTR2_);			\
	    break;							\
	case REALSXP:							\
	    V2SPV(d, ISNZ_REAL, double, REALSXP, REAL,			\
		  _CTYPE2_, _SEXPTYPE2_, _PTR2_);			\
	    break;							\
	case CPLXSXP:							\
	    V2SPV(z, ISNZ_COMPLEX, Rcomplex, CPLXSXP, COMPLEX,		\
		  _CTYPE2_, _SEXPTYPE2_, _PTR2_);			\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE("object", tx, "v2spV");			\
	    break;							\
	}								\
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
    
/* La_norm_type() and La_rcond_type() have been in src/include/R_ext/Lapack.h
   and later in src/modules/lapack/Lapack.c but have still not been available 
   to package writers ...
*/
char La_norm_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = (char) toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* aliases */
    else if (typup == 'E')
	typup = 'F';
    else if (typup != 'M' && typup != 'O' && typup != 'I' && typup != 'F')
	error(_("argument type[1]='%s' must be one of 'M','1','O','I','F', or 'E'"),
	      typstr);
    return typup;
}

char La_rcond_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = (char) toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* alias */
    else if (typup != 'O' && typup != 'I')
	error(_("argument type[1]='%s' must be one of '1','O', or 'I'"),
	      typstr);
    return typup; /* 'O' or 'I' */
}

SEXP as_det_obj(double mod, int log, int sign)
{
    SEXP det = PROTECT(allocVector(VECSXP, 2)),
	nms = PROTECT(allocVector(STRSXP, 2)),
	val = PROTECT(ScalarReal(mod));

    setAttrib(det, R_NamesSymbol, nms);
    SET_STRING_ELT(nms, 0, mkChar("modulus"));
    SET_STRING_ELT(nms, 1, mkChar("sign"));
    setAttrib(val, install("logarithm"), ScalarLogical(log));
    SET_VECTOR_ELT(det, 0, val);
    SET_VECTOR_ELT(det, 1, ScalarInteger(sign));
    setAttrib(det, R_ClassSymbol, mkString("det"));
    UNPROTECT(3);
    return det;
}

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

    if(TYPEOF(di) != INTSXP) {di = PROTECT(coerceVector(di, INTSXP)); nprot++; }
    if(TYPEOF(ij) != INTSXP) {ij = PROTECT(coerceVector(ij, INTSXP)); nprot++; }
    if(!isMatrix(ij) ||
       (ij_di = INTEGER(getAttrib(ij, R_DimSymbol)))[1] != 2)
	error(_("Argument ij must be 2-column integer matrix"));
    n = ij_di[0];
    int *Di = INTEGER(di), *IJ = INTEGER(ij),
	*j_ = IJ+n;/* pointer offset! */

    if((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
	ans = PROTECT(allocVector(REALSXP, n));
	double *ii = REAL(ans), nr = (double) Di[0];
#define do_ii_FILL(_i_, _j_)						\
	int i;								\
	if(check_bounds) {						\
	    for(i=0; i < n; i++) {					\
		if(_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER)	\
		    ii[i] = NA_INTEGER;					\
		else {							\
		    register int i_i, j_i;				\
	            if(one_ind) { i_i = _i_[i]-1; j_i = _j_[i]-1; }	\
	            else        { i_i = _i_[i]  ; j_i = _j_[i]  ; }	\
		    if(i_i < 0 || i_i >= Di[0])				\
			error(_("subscript 'i' out of bounds in M[ij]")); \
		    if(j_i < 0 || j_i >= Di[1])				\
			error(_("subscript 'j' out of bounds in M[ij]")); \
		    ii[i] = i_i + j_i * nr;				\
		}							\
	    }								\
	} else {							\
	    for(i=0; i < n; i++)					\
		ii[i] = (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER)	\
		    ? NA_INTEGER					\
 	            : (one_ind ? ((_i_[i]-1) + (_j_[i]-1)*nr)		\
	                       :   _i_[i]    +  _j_[i]   *nr);		\
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

    if(TYPEOF(di)!= INTSXP) {di = PROTECT(coerceVector(di,INTSXP)); nprot++; }
    if(TYPEOF(i) != INTSXP) { i = PROTECT(coerceVector(i, INTSXP)); nprot++; }
    if(TYPEOF(j) != INTSXP) { j = PROTECT(coerceVector(j, INTSXP)); nprot++; }
    if(LENGTH(j) != n)
	error(_("i and j must be integer vectors of the same length"));
    int *Di = INTEGER(di), *i_ = INTEGER(i), *j_ = INTEGER(j);

    if((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
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
//	     data, nrow, ncol, byrow, dimnames,
//	     missing(nrow), missing(ncol))
SEXP Mmatrix(SEXP args)
{
    SEXP vals, ans, snr, snc, dimnames;
    int nr = 1, nc = 1, byrow, miss_nr, miss_nc;
    R_xlen_t lendat;

    args = CDR(args); /* skip 'name' */
    vals = CAR(args); args = CDR(args);
    /* Supposedly as.vector() gave a vector type, but we check */
    switch(TYPEOF(vals)) {
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

    if(lendat > 0) {
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
	}
	else if ((lendat > 1) && (nrc == 0)){
	    warning(_("data length exceeds size of matrix"));
	}
    }

#ifndef LONG_VECTOR_SUPPORT
   if ((double)nr * (double)nc > INT_MAX)
	error(_("too many elements specified"));
#endif

    PROTECT(ans = allocMatrix(TYPEOF(vals), nr, nc));
    if(lendat) {
	if (isVector(vals))
	    copyMatrix(ans, vals, byrow);
	else
	    copyListMatrix(ans, vals, byrow);
    } else if (isVector(vals)) { /* fill with NAs */
	R_xlen_t N = (R_xlen_t) nr * nc, i;
	switch(TYPEOF(vals)) {
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
	    Rcomplex zna = { NA_REAL, 0.0 };
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
    if(!isNull(dimnames)&& length(dimnames) > 0)
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
 */SEXP R_rbind2_vector(SEXP a, SEXP b) {
    int *d_a = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*d_b = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	n1 = d_a[0], m = d_a[1],
	n2 = d_b[0];
    if(d_b[1] != m)
	error(_("the number of columns differ in R_rbind2_vector: %d != %d"),
	      m, d_b[1]);
    SEXP
	a_x = GET_SLOT(a, Matrix_xSym),
	b_x = GET_SLOT(b, Matrix_xSym);
    int nprot = 1;
    // Care: can have "ddenseMatrix" "ldenseMatrix" or "ndenseMatrix"
    if(TYPEOF(a_x) != TYPEOF(b_x)) { // choose the "common type"
	// Now know: either LGLSXP or REALSXP. FIXME for iMatrix, zMatrix,..
	if(TYPEOF(a_x) != REALSXP) {
	    a_x = PROTECT(duplicate(coerceVector(a_x, REALSXP))); nprot++;
	} else if(TYPEOF(b_x) != REALSXP) {
	    b_x = PROTECT(duplicate(coerceVector(b_x, REALSXP))); nprot++;
	}
    }

    SEXP ans = PROTECT(allocVector(TYPEOF(a_x), m * (n1 + n2)));
    int ii = 0;
    switch(TYPEOF(a_x)) {
    case LGLSXP: {
	int
	    *r = LOGICAL(ans),
	    *ax= LOGICAL(a_x),
	    *bx= LOGICAL(b_x);

#define COPY_a_AND_b_j					\
	for(int j=0; j < m; j++) {			\
	    Memcpy(r+ii, ax+ j*n1, n1); ii += n1;	\
	    Memcpy(r+ii, bx+ j*n2, n2); ii += n2;	\
	} ; break

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
// all0     <- function(x) !any(is.na(x)) && all(!x) ## ~= allFalse
// allFalse <- function(x) !any(x) && !any(is.na(x)) ## ~= all0
SEXP R_all0(SEXP x) {
    if (!isVectorAtomic(x)) {
	if(length(x) == 0) return TRUE_;
	// Typically S4.  TODO: Call the R code above, instead!
	error(_("Argument must be numeric-like atomic vector"));
    }
    R_xlen_t i, n = XLENGTH(x);
    if(n == 0) return TRUE_;

    switch(TYPEOF(x)) {
    case LGLSXP: {
	int *xx = LOGICAL(x);
	for(i=0; i < n; i++)
	    if(xx[i] == NA_LOGICAL || xx[i] != 0) return FALSE_;
	return TRUE_;
    }
    case INTSXP: {
	int *xx = INTEGER(x);
	for(i=0; i < n; i++)
	    if(xx[i] == NA_INTEGER || xx[i] != 0) return FALSE_;
	return TRUE_;
    }
    case REALSXP: {
	double *xx = REAL(x);
	for(i=0; i < n; i++)
	    if(ISNAN(xx[i]) || xx[i] != 0.) return FALSE_;
	return TRUE_;
    }
    case RAWSXP: {
	unsigned char *xx = RAW(x);
	for(i=0; i < n; i++)
	    if(xx[i] != 0) return FALSE_;
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
	if(length(x) == 0) return FALSE_;
	// Typically S4.  TODO: Call the R code above, instead!
	error(_("Argument must be numeric-like atomic vector"));
    }
    R_xlen_t i, n = XLENGTH(x);
    if(n == 0) return FALSE_;

    switch(TYPEOF(x)) {
    case LGLSXP: {
	int *xx = LOGICAL(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    case INTSXP: {
	int *xx = INTEGER(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    case REALSXP: {
	double *xx = REAL(x);
	for(i=0; i < n; i++) if(xx[i] == 0.) return TRUE_;
	return FALSE_;
    }
    case RAWSXP: {
	unsigned char *xx = RAW(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    }
    error(_("Argument must be numeric-like atomic vector"));
    return R_NilValue; // -Wall
}

#undef TRUE_
#undef FALSE_
