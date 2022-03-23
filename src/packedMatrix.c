#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "packedMatrix.h"

/* int i, j, n; R_xlen_t n2; */
#define PM_AR21_UP(i, j) (i) + ((j) * (((R_xlen_t) (j)) + 1)) / 2
#define PM_AR21_LO(i, j, n2) (i) + ((j) * ((n2) - (j) - 1)) / 2
#define PM_LENGTH(n) ((n) * (((R_xlen_t) (n)) + 1)) / 2

/* An alternative to the existing utility 'symmetric_DimNames' 
   enabling users to avoid a copy in cases where it is avoidable ...
   perhaps the existing utility can be rebuilt around this one 
   so that they remain in sync ...
*/
void fast_symmetric_DimNames(SEXP dn, SEXP *vec, SEXP *nm)
{
    *vec = VECTOR_ELT(dn, 0);
    if (isNull(*vec)) {
	*vec = VECTOR_ELT(dn, 1);
    }
    *nm = getAttrib(dn, R_NamesSymbol);
    if (!isNull(*nm)) {
	*nm = STRING_ELT(*nm, 0);
	if (*nm == NA_STRING) {
	    *nm = STRING_ELT(*nm, 1);
	}
    }
}

/* FIXME: could avoid some duplication by calling 'symmetricMatrix_validate'
   or 'triangularMatrix_validate' conditional on existence of 'diag' slot ...
   would still need to check length(.@x) though
*/
SEXP packedMatrix_validate(SEXP obj)
{
    SEXP val = GET_SLOT(obj, Matrix_DimSym);
    if (LENGTH(val) != 2) {
        return mkString(_("'Dim' slot does not have length 2"));
    }
    int n = INTEGER(val)[0];
    if (INTEGER(val)[1] != n) {
        return mkString(_("matrix is not square"));
    }
    val = check_scalar_string(GET_SLOT(obj, Matrix_uploSym), "LU", "uplo");
    if (isString(val)) {
        return val;
    }
    if (R_has_slot(obj, Matrix_diagSym)) {
	val = check_scalar_string(GET_SLOT(obj, Matrix_diagSym), "NU", "diag");
	if (isString(val)) {
	    return val;
	}
    }
    val = GET_SLOT(obj, Matrix_xSym);
    if (XLENGTH(val) != PM_LENGTH(n)) {
        return mkString(_("'x' slot does not have length 'n*(n+1)/2', n=Dim[1]"));
    }
    return ScalarLogical(1);
}

#define PM_T_LOOP							\
    do {								\
	if (up) {							\
	    for (int j = 0; j < n; ++j) {				\
		for (int i = j; i < n; ++i) {				\
		    *(px1++) = px0[PM_AR21_UP(j, i)];			\
		}							\
	    }								\
	} else {							\
	    R_xlen_t n2 = ((R_xlen_t) n) * 2;				\
	    for (int j = 0; j < n; ++j) {				\
		for (int i = 0; i <= j; ++i) {				\
		    *(px1++) = px0[PM_AR21_LO(j, i, n2)];		\
		}							\
	    }								\
	}								\
    } while (0)

#define PM_T(_datatype_, _sexptype_, _accessor_)			\
    do {								\
	SEXP x1 = PROTECT(allocVector(_sexptype_, XLENGTH(x0)));	\
	_datatype_ *px0 = _accessor_(x0);				\
	_datatype_ *px1 = _accessor_(x1);				\
	PM_T_LOOP;							\
	SET_SLOT(res, Matrix_xSym, x1);					\
	UNPROTECT(1);							\
    } while (0)

/* t(x) */
SEXP packedMatrix_t(SEXP obj)
{
    /* Initialize result of same class */
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(class_P(obj)));
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    /* ? upper -> lower : lower -> upper */
    Rboolean up = uplo_P(obj)[0] == 'U';
    
    SEXP x0 = GET_SLOT(obj, Matrix_xSym);
    if (n > 1) {
	/* Permute 'x' slot */
	if (isReal(x0)) { /* "d" */
	    PM_T(double, REALSXP, REAL);
	} else { /* "[ln]" */
	    PM_T(int, LGLSXP, LOGICAL);
	}
    } else {
	/* Preserve 'x' slot */
	SET_SLOT(res, Matrix_xSym, x0);
    }

    /* Toggle 'uplo' slot */
    SET_SLOT(res, Matrix_uploSym, mkString(up ? "L" : "U"));
    /* Preserve 'Dim' slot */
    SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
    /* Reverse 'Dimnames' slot and 'names(Dimnames)' (if not absent) */
    SEXP dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    SEXP dn1 = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn1, 0, VECTOR_ELT(dn0, 1));
    SET_VECTOR_ELT(dn1, 1, VECTOR_ELT(dn0, 0));
    SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
    if (!isNull(ndn0)) {
	SEXP ndn1 = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(ndn1, 0, STRING_ELT(ndn0, 1));
	SET_STRING_ELT(ndn1, 1, STRING_ELT(ndn0, 0));
	setAttrib(dn1, R_NamesSymbol, ndn1);
	UNPROTECT(1);
    }
    SET_SLOT(res, Matrix_DimNamesSym, dn1);
    UNPROTECT(2);
    return res;
}

#undef PM_T
#undef PM_T_LOOP

#define PM_D_G_UDIAG(_one_)			\
    do {					\
	for (int j = 0; j < n; ++j) {		\
	    *(pres++) = _one_;			\
	}					\
    } while(0)

#define PM_D_G_NDIAG							\
    do {								\
	for (int j = 0; j < n; px += (up ? 1 + (++j) : n - (j++))) {	\
	    *(pres++) = *px;						\
	}								\
    } while (0)

#define PM_D_G(_datatype_, _sexptype_, _accessor_, _one_)	\
    do {							\
	res = PROTECT(allocVector(_sexptype_, n));		\
	_datatype_ *pres = _accessor_(res);			\
	if (utr) {						\
	    PM_D_G_UDIAG(_one_);				\
	} else {						\
	    _datatype_ *px = _accessor_(x);			\
	    PM_D_G_NDIAG;					\
	}							\
    } while (0)

/* diag(x, names) */
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL) {
	error(_("'names' must be TRUE or FALSE"));
    }

    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    /* ? triangular : symmetric */
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);
    /* ? unit triangular : other */
    Rboolean utr = tr && diag_P(obj)[0] == 'U';
    /* ? upper : lower */
    Rboolean up = uplo_P(obj)[0] == 'U';

    SEXP res;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (isReal(x)) { /* "d" */
	PM_D_G(double, REALSXP, REAL, 1.0);
    } else { /* "[ln]" */
	PM_D_G(int, LGLSXP, LOGICAL, 1);
    }

    if (do_nms) {
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym);
	SEXP rn = VECTOR_ELT(dn, 0);
	SEXP cn = VECTOR_ELT(dn, 1);
	if (isNull(rn)) {
	    if (!tr && !isNull(cn)) {
		setAttrib(res, R_NamesSymbol, cn);
	    }
	} else {
	    if (!tr || R_compute_identical(rn, cn, 16)) {
		setAttrib(res, R_NamesSymbol, rn);
	    }
	}
    }
    UNPROTECT(1);
    return res;
}

#undef PM_D_G
#undef PM_D_G_NDIAG
#undef PM_D_G_UDIAG

#define PM_D_S_ONE							\
    do {								\
	for (int j = 0; j < n; px += (up ? 1 + (++j) : n - (j++))) {	\
	    *px = d;							\
	}								\
    } while (0)

#define PM_D_S_FULL							\
    do {								\
	for (int j = 0; j < n; px += (up ? 1 + (++j) : n - (j++))) {	\
	    *px = *(pval++);						\
	}								\
    } while (0)

#define PM_D_S(_datatype_, _accessor_)		\
    do {					\
	_datatype_ *px = _accessor_(x);		\
	_datatype_ *pval = _accessor_(val);	\
	if (nv1) {				\
	    _datatype_ d = pval[0];		\
	    PM_D_S_ONE;				\
	} else {				\
	    PM_D_S_FULL;			\
	}					\
    } while (0)

/* diag(x) <- value */
SEXP packedMatrix_diag_set(SEXP obj, SEXP val)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    int nv = LENGTH(val);
    Rboolean nv1 = (nv == 1);
    if (!(nv1 || nv == n)) {
	error(_("replacement diagonal has wrong length"));
    }

    /* Initialize result object of same class */
    SEXP res;
    int nprotect = 0;
    if (MAYBE_REFERENCED(obj)) {
	res = PROTECT(NEW_OBJECT_OF_CLASS(class_P(obj))); ++nprotect;
	SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
	SET_SLOT(res, Matrix_DimNamesSym, GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));
	slot_dup(res, obj, Matrix_xSym);
    } else {
	res = obj;
    }

    /* Toggle 'diag' slot when assigning to unit diagonal "[dln]tpMatrix",
       _even if_ RHS of assignment is unity
     */
    if (Diag_P(res)[0] == 'U') {
	SET_SLOT(res, Matrix_diagSym, mkString("N"));
    }
    /* Reset 'factors' slot when assigning to "[dln]spMatrix" */
    if (R_has_slot(res, Matrix_factorSym) &&
	LENGTH(GET_SLOT(res, Matrix_factorSym)) > 0) {
	SET_SLOT(res, Matrix_factorSym, allocVector(VECSXP, 0));
    }

    /* ? double result : logical result */
    Rboolean dbl = TRUE;
    /* ? upper : lower */
    Rboolean up = uplo_P(res)[0] == 'U';

    /* Test that LHS and RHS of assignment have compatible types,
       coercing one or the other from logical to double if necessary
     */
    SEXP x = GET_SLOT(res, Matrix_xSym);
    switch (TYPEOF(x)) {
    case LGLSXP:
	switch (TYPEOF(val)) {
	case LGLSXP:
	    dbl = FALSE;
	    break;
	case INTSXP:
	    val = PROTECT(coerceVector(val, REALSXP)); ++nprotect;
	case REALSXP:
	{
	    /* [ln][st]pMatrix -> d[st]pMatrix */
	    SEXP strcl = getAttrib(res, R_ClassSymbol);
	    char *cl = strdup(CHAR(STRING_ELT(strcl, 0)));
	    cl[0] = 'd';
	    SET_STRING_ELT(strcl, 0, mkChar(cl));
	    free(cl);
	    x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
	    SET_SLOT(res, Matrix_xSym, x);
	    break;
	}
	default:
	    error(_("replacement diagonal has incompatible type '%s'"),
		  type2char(TYPEOF(val)));
	}
	break;
    case REALSXP:
	switch (TYPEOF(val)) {
	case LGLSXP:
	case INTSXP:
	    /* logical, integer -> double */
	    val = PROTECT(coerceVector(val, REALSXP)); ++nprotect;
	case REALSXP:
	    break;
	default:
	    error(_("replacement diagonal has incompatible type '%s'"),
		  type2char(TYPEOF(val)));
	}
	break;
    default:
	error(_("'x' slot is not of type '%s' or '%s', which should never happen; please report"),
	      type2char(LGLSXP), type2char(REALSXP));
    }

    if (dbl) {
	PM_D_S(double, REAL);
    } else {
	PM_D_S(int, LOGICAL);
    }
    UNPROTECT(nprotect);
    return res;
}

#undef PM_D_S
#undef PM_D_S_FULL
#undef PM_D_S_ONE

#define PM_IJ2K(_px_, _zero_, _one_)				\
    (utr && i == j						\
     ? _one_							\
     : (up							\
	? (i <= j						\
	   ? _px_[PM_AR21_UP(i, j)]				\
	   : (tr						\
	      ? _zero_						\
	      : _px_[PM_AR21_UP(j, i)]))			\
	: (i >= j						\
	   ? _px_[PM_AR21_LO(i, j, n2)]				\
	   : (tr						\
	      ? _zero_						\
	      : _px_[PM_AR21_LO(j, i, n2)]))))

#define PM_SUB1_LOOP(_na_, _zero_, _one_)				\
    do {								\
	R_xlen_t n2 = ((R_xlen_t) n) * 2;				\
	R_xlen_t nn = ((R_xlen_t) n) * n;				\
	if (TYPEOF(index) == INTSXP) {					\
	    int pos, *pindex = INTEGER(index);				\
	    int i, j;							\
	    for (R_xlen_t k = 0; k < nindex; ++k) {			\
		pos = *(pindex++);					\
		if (pos == NA_INTEGER || pos > nn) {			\
		    *(pres++) = _na_;					\
		} else {						\
		    pos -= 1; /* 1-index -> 0-index */			\
		    i = pos % n;					\
		    j = pos / n;					\
		    *(pres++) = PM_IJ2K(px, _zero_, _one_);		\
		}							\
	    }								\
	} else {							\
	    double pos, *pindex = REAL(index);				\
	    int i, j;							\
	    for (R_xlen_t k = 0, truncpos; k < nindex; ++k) {		\
		pos = *(pindex++);					\
		if (!R_FINITE(pos) || (truncpos = (R_xlen_t) pos) > nn) { \
		    *(pres++) = _na_;					\
		} else {						\
		    truncpos -= 1; /* 1-index -> 0-index */		\
		    i = truncpos % n;					\
		    j = truncpos / n;					\
		    *(pres++) = PM_IJ2K(px, _zero_, _one_);		\
		}							\
	    }								\
	}								\
    } while (0)

#define PM_SUB1_END(_datatype_, _sexptype_, _accessor_, _na_, _zero_, _one_) \
    do {								\
	SEXP res = PROTECT(allocVector(_sexptype_, nindex));		\
	_datatype_ *pres = _accessor_(res);				\
	_datatype_ *px = _accessor_(x);					\
	PM_SUB1_LOOP(_na_, _zero_, _one_);				\
	UNPROTECT(1);							\
	return res;							\
    } while (0)

#define PM_SUB1								\
    do {								\
	SEXP x = GET_SLOT(obj, Matrix_xSym);				\
	int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];		\
	/* ? triangular : symmetric */					\
	Rboolean tr = R_has_slot(obj, Matrix_diagSym);			\
	/* ? unit triangular : other */					\
	Rboolean utr = tr && diag_P(obj)[0] == 'U';			\
	/* ? upper : lower */						\
	Rboolean up = uplo_P(obj)[0] == 'U';				\
									\
	if (isReal(x)) { /* "d" */					\
	    PM_SUB1_END(double, REALSXP, REAL, NA_REAL, 0.0, 1.0);	\
	} else { /* "[ln]" */						\
	    PM_SUB1_END(int, LGLSXP, LOGICAL, NA_LOGICAL, 0, 1);	\
	}								\
    } while (0)

/* 'x[index]' where 'index' is an integer _or_ double vector
   with elements greater than or equal to 1
*/
SEXP packedMatrix_sub1(SEXP obj, SEXP index)
{
    R_xlen_t nindex = XLENGTH(index);
    PM_SUB1;
}

#undef PM_SUB1_LOOP

#define PM_SUB1_LOOP(_na_, _zero_, _one_)				\
    do {								\
	R_xlen_t n2 = ((R_xlen_t) n) * 2;				\
	int *pi = INTEGER(index);					\
	int *pj = pi + nindex;						\
	for (int k = 0, i, j; k < nindex; ++k) {			\
	    i = *(pi++);						\
	    j = *(pj++);						\
	    if (i == NA_INTEGER || j == NA_INTEGER) {			\
		*(pres++) = _na_;					\
	    } else {							\
		i -= 1; /* 1-index -> 0-index */			\
		j -= 1;							\
		*(pres++) = PM_IJ2K(px, _zero_, _one_);			\
	    }								\
	}								\
    } while (0)

/* 'x[index]' where 'index' is a 2-column integer matrix supplying 
   integers in 'c(1:n, NA)' _only_
*/
SEXP packedMatrix_sub1_mat(SEXP obj, SEXP index)
{
    int nindex = INTEGER(getAttrib(index, R_DimSymbol))[0];
    PM_SUB1;
}

#undef PM_SUB1
#undef PM_SUB1_LOOP

#define PM_SUB2_LOOP(_na_, _zero_, _one_, _FOR_)			\
    do {								\
        R_xlen_t n2 = ((R_xlen_t) n) * 2;				\
	int i, j;							\
	for (int kj = 0; kj < nj; ++kj) {				\
	    if (mj) {							\
		j = kj;							\
	    } else {							\
		j = pj[kj];						\
		if (j == NA_INTEGER) {					\
		    _FOR_ {						\
			*(px1++) = _na_;				\
		    }							\
		    if (do_cn) {					\
			SET_STRING_ELT(cn1, kj, NA_STRING);		\
		    }							\
		    continue;						\
		} else {						\
		    j -= 1;						\
		}							\
	    }								\
	    _FOR_ {							\
		if (mi) {						\
		    i = ki;						\
		} else {						\
		    i = pi[ki];						\
		    if (i == NA_INTEGER) {				\
			*(px1++) = _na_;				\
			continue;					\
		    } else {						\
			i -= 1;						\
		    }							\
		}							\
		*(px1++) = PM_IJ2K(px0, _zero_, _one_);			\
	    }								\
	    if (do_cn) {						\
		SET_STRING_ELT(cn1, kj, STRING_ELT(cn0, j));		\
	    }								\
	}								\
	if (do_rn) {							\
	    for (int ki = 0; ki < ni; ++ki) {				\
		if (mi) {						\
		    i = ki;						\
		} else {						\
		    i = pi[ki];						\
		    if (i == NA_INTEGER) {				\
			SET_STRING_ELT(rn1, ki, NA_STRING);		\
			continue;					\
		    } else {						\
			i -= 1;						\
		    }							\
		}							\
		SET_STRING_ELT(rn1, ki, STRING_ELT(rn0, i));		\
	    }								\
	}								\
    } while (0)

#define PM_SUB2(_datatype_, _sexptype_, _accessor_, _na_, _zero_, _one_) \
    do {								\
	R_xlen_t len = (do_sp ? PM_LENGTH(ni) : ni * nj);		\
	SEXP x1 = PROTECT(allocVector(_sexptype_, len));		\
	_datatype_ *px0 = _accessor_(x0);				\
	_datatype_ *px1 = _accessor_(x1);				\
	if (do_sp) {							\
	    if (up) {							\
		PM_SUB2_LOOP(_na_, _zero_, _one_,			\
			     for (int ki = 0; ki <= kj; ++ki));		\
	    } else {							\
	    	PM_SUB2_LOOP(_na_, _zero_, _one_,			\
			     for (int ki = kj; ki < ni; ++ki));		\
	    }								\
	} else {							\
	    PM_SUB2_LOOP(_na_, _zero_, _one_,				\
			 for (int ki = 0; ki < ni; ++ki));		\
	}								\
	SET_SLOT(res, Matrix_xSym, x1);					\
	UNPROTECT(1);							\
    } while (0)

/* 'x[index1, ]', 'x[, index2]', and 'x[index1, index2]'
   where 'index1' and 'index2' are integer vectors supplying
   integers in 'c(1:n, NA)' _only_ ... NULL indicates missingness
*/
SEXP packedMatrix_sub2(SEXP obj, SEXP index1, SEXP index2, SEXP drop)
{
    Rboolean do_drop = asLogical(drop) != 0;
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    int nprotect = 0;
    /* ? triangular : symmetric */
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);
    /* ? unit triangular : other */
    Rboolean utr = tr && diag_P(obj)[0] == 'U';
    /* ? upper : lower */
    Rboolean up = uplo_P(obj)[0] == 'U';

    Rboolean mi, mj;
    R_xlen_t ni, nj;
    int *pi, *pj;
    
    mi = isNull(index1);
    mj = isNull(index2);
    if (mi) {
	ni = n;
    } else {
	ni = XLENGTH(index1);
	if (ni > INT_MAX) {
	    error(_("dimensions of result cannot exceed 2^31-1"));
	}
	pi = INTEGER(index1);
    }
    if (mj) {
	nj = n;
    } else {
	nj = XLENGTH(index2);
	if (nj > INT_MAX) {
	    error(_("dimensions of result cannot exceed 2^31-1"));
	}
	pj = INTEGER(index2);
    }
    
    /* Initialize result of same type but "general" class, except
       for symmetric indexing of symmetric matrix, when class is
       retained also
    */
    static const char *valid[] = {
	"dspMatrix", "lspMatrix", "nspMatrix",
	"dtpMatrix", "ltpMatrix", "ntpMatrix", ""};
    char *cl = strdup(valid[R_check_class_etc(obj, valid)]);
    Rboolean do_sp = (cl[1] == 's' && !mi && !mj && ni == nj &&
		      !memcmp(pi, pj, ni * sizeof(int))); 
    if (!do_sp) {
	cl[1] = 'g';
	cl[2] = 'e';
    }
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl)); ++nprotect;
    free(cl);

    /* Set 'uplo' slot (if retained) */
    if (do_sp) {
	SET_SLOT(res, Matrix_uploSym, mkString(up ? "U" : "L"));
    }
    
    /* Set 'Dim' slot */
    SEXP d1 = PROTECT(GET_SLOT(res, Matrix_DimSym));
    INTEGER(d1)[0] = ni;
    INTEGER(d1)[1] = nj;
    UNPROTECT(1);

    /* Set 'Dimnames' slot and 'names(Dimnames)' (if not absent) */
    SEXP dn0, dn1, rn0, rn1, cn0, cn1;
    dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    dn1 = PROTECT(GET_SLOT(res, Matrix_DimNamesSym)); ++nprotect;
    if (tr) {
	rn0 = VECTOR_ELT(dn0, 0);
	cn0 = VECTOR_ELT(dn0, 1);
	SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
	if (!isNull(ndn0)) {
	    setAttrib(dn1, R_NamesSymbol, ndn0);
	}
    } else {
	SEXP s;
	fast_symmetric_DimNames(dn0, &rn0, &s);
	cn0 = rn0;
	if (!isNull(s)) {
	    SEXP ndn1 = PROTECT(allocVector(STRSXP, 2));
	    SET_STRING_ELT(ndn1, 0, s);
	    SET_STRING_ELT(ndn1, 1, s);
	    setAttrib(dn1, R_NamesSymbol, ndn1);
	    UNPROTECT(1);
	}
    }
    Rboolean has_rn, has_cn, do_rn, do_cn;
    has_rn = !isNull(rn0) && ni;
    has_cn = !isNull(cn0) && nj;
    do_rn = do_cn = FALSE;
    if (has_rn) {
	if (mi) {
	    SET_VECTOR_ELT(dn1, 0, rn0);
	} else {
	    rn1 = PROTECT(allocVector(STRSXP, ni)); ++nprotect;
	    SET_VECTOR_ELT(dn1, 0, rn1);
	    do_rn = TRUE;
	}
    }
    if (has_cn && !do_sp) {
	if (mj) {
	    SET_VECTOR_ELT(dn1, 1, cn0);
	} else {
	    cn1 = PROTECT(allocVector(STRSXP, nj)); ++nprotect;
	    SET_VECTOR_ELT(dn1, 1, cn1);
	    do_cn = TRUE;
	}
    }

    /* Set 'x' slot */
    SEXP x0 = GET_SLOT(obj, Matrix_xSym);
    if (isReal(x0)) { /* "d" */
	PM_SUB2(double, REALSXP, REAL, NA_REAL, 0.0, 1.0);
    } else { /* "[ln]" */
	PM_SUB2(int, LGLSXP, LOGICAL, NA_LOGICAL, 0, 1);
    }

    /* Drop dimensions in this special case */
    if (do_drop && (ni == 1 || nj == 1)) {
	res = PROTECT(GET_SLOT(res, Matrix_xSym)); ++nprotect;
	if (has_rn && nj == 1 && ni != 1) {
	    setAttrib(res, R_NamesSymbol, VECTOR_ELT(dn1, 0));
	} else if (has_cn && ni == 1 && nj != 1) {
	    setAttrib(res, R_NamesSymbol, VECTOR_ELT(dn1, 1));
	}
    }
    UNPROTECT(nprotect);
    return res;
}

#undef PM_SUB2
#undef PM_SUB2_LOOP
#undef PM_IJ2K
