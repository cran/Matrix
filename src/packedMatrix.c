#include "packedMatrix.h"

/* unpack(x), returning unpackedMatrix */
SEXP packedMatrix_unpack(SEXP from, SEXP strict)
{
    static const char *valid_from[] = {
	/* 0 */ "pCholesky", "pBunchKaufman", /* must match before dtpMatrix */
	/* 2 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 5 */ "dppMatrix", /* must match before dspMatrix */
	/* 6 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    static const char *valid_to[] = {
	/* 0 */ "Cholesky", "BunchKaufman",
	/* 2 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 5 */ "dpoMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
    int ivalid = R_check_class_etc(from, valid_from);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "packedMatrix_unpack");
    if (asLogical(strict) == 0) {
	if (ivalid < 2)
	    ivalid = 2; /* pCholesky,pBunchKaufman->dtrMatrix */
	else if (ivalid == 5)
	    ivalid = 6; /* dppMatrix->dsyMatrix */
    }

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(valid_to[ivalid])),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = GET_SLOT(from, Matrix_uploSym),
	x_from = GET_SLOT(from, Matrix_xSym),
	x_to;

    int n = INTEGER(dim)[0];
    if ((double) n * n > R_XLEN_T_MAX)
	error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));

    char ul = *CHAR(STRING_ELT(uplo, 0));
    SEXPTYPE tx = TYPEOF(x_from);
    R_xlen_t nx = (R_xlen_t) n * n;
    PROTECT(x_to = allocVector(tx, nx));

#define UNPACK(_PREFIX_, _PTR_)						\
    do {								\
	Memzero(_PTR_(x_to), nx);					\
	_PREFIX_ ## dense_unpack(_PTR_(x_to), _PTR_(x_from), n, ul, 'N'); \
    } while (0)
    
    switch (tx) {
    case REALSXP: /* d..Matrix */
	UNPACK(d, REAL);
	break; 
    case LGLSXP: /* [ln]..Matrix */
	UNPACK(i, LOGICAL);
	break;
    case INTSXP: /* i..Matrix */
	UNPACK(i, INTEGER);
	break;
    case CPLXSXP: /* z..Matrix */
	UNPACK(z, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_unpack");
	break;
    }

#undef UNPACK
    
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    SET_SLOT(to, Matrix_uploSym, uplo);
    SET_SLOT(to, Matrix_xSym, x_to);
    if (ivalid < 5) {
	/* .tpMatrix */
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
	if (ivalid == 1)
	    /* pBunchKaufman */
	    SET_SLOT(to, Matrix_permSym, GET_SLOT(from, Matrix_permSym));
    } else {
	/* .spMatrix */
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    }
    UNPROTECT(2);
    return to;
}

/* forceSymmetric(x, uplo), returning .spMatrix */
SEXP packedMatrix_force_symmetric(SEXP from, SEXP uplo_to)
{
    static const char *valid[] = {
	"dspMatrix", "lspMatrix", "nspMatrix", /* be fast */
	"dtpMatrix", "ltpMatrix", "ntpMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "packedMatrix_force_symmetric");
    const char *clf = valid[ivalid];
    
    char ulf = *uplo_P(from),
	ult = (isNull(uplo_to)) ? ulf : *CHAR(asChar(uplo_to));
    if (clf[1] == 's') {
	/* .spMatrix */
	if (ulf == ult)
	    return from;
	SEXP to = PROTECT(packedMatrix_transpose(from));
	if (clf[0] == 'z') {
	    /* Need _conjugate_ transpose */
	    SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	    conjugate(x);
	    UNPROTECT(1);
	}
	UNPROTECT(1);
	return to;
    }
    
    /* Now handling just .tpMatrix ... */
    
    char clt[] = ".spMatrix";
    clt[0] = clf[0];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x_from = GET_SLOT(from, Matrix_xSym);
    
    if (ulf == ult) {
	/* .tpMatrix with correct uplo */
	SET_SLOT(to, Matrix_xSym, x_from);
    } else {
	/* .tpMatrix with incorrect uplo */
	int n = INTEGER(dim)[0];
	char di = *diag_P(from);
	SEXPTYPE tx = TYPEOF(x_from);
	R_xlen_t nx = XLENGTH(x_from);
	SEXP x_to = PROTECT(allocVector(tx, nx));

#define COPY_DIAGONAL(_PREFIX_, _PTR_)					\
	do {								\
	    Memzero(_PTR_(x_to), nx);					\
	    _PREFIX_ ## dense_packed_copy_diagonal(			\
		_PTR_(x_to), _PTR_(x_from), n, nx, ult, ulf, di);	\
	} while (0)
	
	switch (tx) {
	case REALSXP: /* d..Matrix */
	    COPY_DIAGONAL(d, REAL);
	    break;
	case LGLSXP: /* [ln]..Matrix */
	    COPY_DIAGONAL(i, LOGICAL);
	    break;
	case INTSXP: /* i..Matrix */
	    COPY_DIAGONAL(i, INTEGER);
	    break;
	case CPLXSXP: /* z..Matrix */
	    COPY_DIAGONAL(z, COMPLEX);
	    break;
	default:
	    ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_force_symmetric");
	    break;
	}

#undef COPY_DIAGONAL
	
	SET_SLOT(to, Matrix_xSym, x_to);
	UNPROTECT(1);
    }
    
    SET_SLOT(to, Matrix_DimSym, dim);
    set_symmetrized_DimNames(to, dimnames, -1);
    SET_SLOT(to, Matrix_uploSym, mkString((ult == 'U') ? "U" : "L"));
    UNPROTECT(1);
    return to;
}

#define PM_IS_DI(_RES_, _X_, _N_, _UPLO_, _METHOD_)			\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_packed_is_diagonal(REAL(_X_), _N_, _UPLO_);	\
	    break;							\
	case LGLSXP:							\
	    _RES_ = idense_packed_is_diagonal(LOGICAL(_X_), _N_, _UPLO_); \
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_packed_is_diagonal(INTEGER(_X_), _N_, _UPLO_); \
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_packed_is_diagonal(COMPLEX(_X_), _N_, _UPLO_); \
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", TYPEOF(_X_), _METHOD_);	\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

/* isSymmetric(x, tol = 0, checkDN) */
/* FIXME: not checking for real diagonal in complex case */
SEXP packedMatrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    static const char *valid[] = {
	"dspMatrix", "lspMatrix", "nspMatrix", /* be fast */
	"dtpMatrix", "ltpMatrix", "ntpMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0) {
	ERROR_INVALID_CLASS(class_P(obj), "packedMatrix_is_symmetric");
	return R_NilValue;
    } else if (ivalid < 3) {
	/* .spMatrix: symmetric by definition */
	return ScalarLogical(1);
    } else {
	/* .tpMatrix: symmetric iff diagonal */
	if (asLogical(checkDN) != 0 &&
	    !DimNames_is_symmetric(GET_SLOT(obj, Matrix_DimNamesSym)))
	    return ScalarLogical(0);
	Rboolean res = FALSE;
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
	char ul = *uplo_P(obj);
	PM_IS_DI(res, x, n, ul, "packedMatrix_is_symmetric");
	return ScalarLogical(res);
    }
}

#define RETURN_TRUE_OF_KIND(_KIND_)					\
    do {								\
	SEXP ans = PROTECT(allocVector(LGLSXP, 1));			\
	LOGICAL(ans)[0] = 1;						\
	setAttrib(ans, install("kind"), _KIND_);			\
	UNPROTECT(1);							\
	return ans;							\
    } while (0)

/* isTriangular(x, upper) */
SEXP packedMatrix_is_triangular(SEXP obj, SEXP upper)
{
    static const char *valid[] = {
	"dtpMatrix", "ltpMatrix", "ntpMatrix", /* be fast */
	"dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "packedMatrix_is_triangular");

    SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
    char ul = *CHAR(STRING_ELT(uplo, 0));
    int need_upper = asLogical(upper);

#define IF_DIAGONAL							\
    Rboolean res = FALSE;						\
    SEXP x = GET_SLOT(obj, Matrix_xSym);				\
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];			\
    PM_IS_DI(res, x, n, ul, "packedMatrix_is_triangular");		\
    if (res)
    
    if (ivalid < 3) {
	/* .tpMatrix: be fast if 'upper', 'uplo' agree; else need diagonal */
	if (need_upper == NA_LOGICAL) {
	    RETURN_TRUE_OF_KIND(uplo);
	} else if ((need_upper != 0 && ul == 'U') ||
		   (need_upper == 0 && ul != 'U')) {
	    return ScalarLogical(1);
	} else {
	    IF_DIAGONAL {
		return ScalarLogical(1);
	    }
	}
    } else {
	/* .spMatrix: triangular iff diagonal */
	IF_DIAGONAL {
	    if (need_upper == NA_LOGICAL) {
		RETURN_TRUE_OF_KIND(mkString("U"));
	    } else {
		return ScalarLogical(1);
	    }
	}
    }

#undef IF_DIAGONAL

    return ScalarLogical(0);
}

#undef RETURN_TRUE_OF_KIND

/* isDiagonal(x) */
SEXP packedMatrix_is_diagonal(SEXP obj)
{
    /* _Not_ checking class of 'obj' */
    Rboolean res = FALSE;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    char ul = *uplo_P(obj);
    PM_IS_DI(res, x, n, ul, "packedMatrix_is_diagonal");
    return ScalarLogical(res);
}

#undef PM_IS_DI

/* t(x), typically preserving class */
SEXP packedMatrix_transpose(SEXP from)
{
    static const char *valid[] = {
	/* 0 */ "pCholesky", "pBunchKaufman",
	/* 2 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 5 */ "dppMatrix",
	/* 6 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "packedMatrix_transpose");
    if (ivalid == 1)
	ivalid = 2; /* pBunchKaufman->dtpMatrix */

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(valid[ivalid])),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x = GET_SLOT(from, Matrix_xSym);
    int n = INTEGER(dim)[0];
    char ul = *uplo_P(from);
    
    /* Preserve 'Dim' slot */
    SET_SLOT(to, Matrix_DimSym, dim);
    /* Reverse or preserve 'Dimnames' slot (preserving if symmetric) */
    if (ivalid < 5)
	set_reversed_DimNames(to, dimnames);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    /* Toggle 'uplo' slot */
    SET_SLOT(to, Matrix_uploSym, mkString((ul == 'U') ? "L" : "U"));
    if (ivalid < 5)
	/* Preserve 'diag' slot */
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    else
	/* Preserve 'factors' slot */
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    
    /* Permute 'x' slot */
    SET_SLOT(to, Matrix_xSym, packed_transpose(x, n, ul));
    
    UNPROTECT(1);
    return to;
}

/* diag(x, names) */
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL)
	error(_("'names' must be TRUE or FALSE"));

    SEXPTYPE tx;
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    SEXP x = GET_SLOT(obj, Matrix_xSym),
	res = PROTECT(allocVector(tx = TYPEOF(x), n));
    char ul = *uplo_P(obj), di = *Diag_P(obj);
    
#define PM_D_G(_CTYPE_, _PTR_, _ONE_)					\
    do {								\
	_CTYPE_ *pres = _PTR_(res);					\
	int j;								\
	if (di == 'U') {						\
	    for (j = 0; j < n; ++j)					\
		*(pres++) = _ONE_;					\
	} else {							\
	    _CTYPE_ *px = _PTR_(x);					\
	    if (ul == 'U')						\
		for (j = 0; j < n; px += (++j)+1)			\
		    *(pres++) = *px;					\
	    else							\
		for (j = 0; j < n; px += n-(j++))			\
		    *(pres++) = *px;					\
	}								\
    } while (0)

    switch (tx) {
    case REALSXP: /* d..Matrix */
	PM_D_G(double, REAL, 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PM_D_G(int, LOGICAL, 1);
	break;
    case INTSXP: /* i..Matrix */
	PM_D_G(int, INTEGER, 1);
	break;
    case CPLXSXP: /* z..Matrix */
	PM_D_G(Rcomplex, COMPLEX, Matrix_zone);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_diag_get");
	break;
    }

#undef PM_D_G
    
    if (do_nms) {
	/* NB: The logic here must be adjusted once the validity method 
	       for 'symmetricMatrix' enforces symmetric 'Dimnames' */
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym),
	    rn = VECTOR_ELT(dn, 0),
	    cn = VECTOR_ELT(dn, 1);
	if (isNull(cn)) {
	    if (di == ' ' && !isNull(rn))
		setAttrib(res, R_NamesSymbol, rn);
	} else {
	    if (di == ' ' || (!isNull(rn) &&
			      (rn == cn || equal_string_vectors(rn, cn, n))))
		setAttrib(res, R_NamesSymbol, cn);
	}
    }
    UNPROTECT(1);
    return res;
}

/* diag(x) <- value */
SEXP packedMatrix_diag_set(SEXP obj, SEXP val)
{
    SEXP dim = GET_SLOT(obj, Matrix_DimSym);
    int n = INTEGER(dim)[0];
    R_xlen_t nv = XLENGTH(val);
    if (nv != 1 && nv != n)
	error(_("replacement diagonal has wrong length"));
    
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    SEXPTYPE tx = TYPEOF(x), tv = TYPEOF(val);
    if (tx < LGLSXP || tx > CPLXSXP)
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_diag_set");
    if (tv < LGLSXP || tv > REALSXP)
	/* Upper bound can become CPLXSXP once we have proper zMatrix */
	error(_("replacement diagonal has incompatible type \"%s\""),
	      type2char(tv));
    
    static const char *valid[] = {
	"dtpMatrix", "dspMatrix",
	"ltpMatrix", "lspMatrix",
	"ntpMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "packedMatrix_diag_set");

    SEXP res,
	dimnames = GET_SLOT(obj, Matrix_DimNamesSym),
	uplo = GET_SLOT(obj, Matrix_uploSym);
    char ul = *CHAR(STRING_ELT(uplo, 0));
    int nprotect = 2;

    /* Allocate and coerce as necessary */
    if (tv <= tx) {
	if (tv < tx) {
	    PROTECT(val = coerceVector(val, tv = tx));
	    ++nprotect;
	}
	PROTECT(res = NEW_OBJECT_OF_CLASS(valid[ivalid]));
	PROTECT(x = duplicate(x));
    } else { /* tv > tx */
	/* dMatrix result is only possibility until we have proper [iz]Matrix */
	if (tv < REALSXP) {
	    PROTECT(val = coerceVector(val, tv = REALSXP));
	    ++nprotect;
	}
	char cl[] = "d.pMatrix";
	cl[1] = valid[ivalid][1];
	PROTECT(res = NEW_OBJECT_OF_CLASS(cl));
	PROTECT(x = coerceVector(x, tx = tv));
    }
    SET_SLOT(res, Matrix_xSym, x);
    
    /* Transfer slots other than 'x', 'diag', and 'factors';
       latter two should keep their prototypes "N" and list() */
    SET_SLOT(res, Matrix_DimSym, dim);
    SET_SLOT(res, Matrix_DimNamesSym, dimnames);
    SET_SLOT(res, Matrix_uploSym, uplo);
    
#define PM_D_S(_CTYPE_, _PTR_)						\
    do {								\
	_CTYPE_ *px = _PTR_(x), *pval = _PTR_(val);			\
	int j;								\
	if (nv == 1) {							\
	    if (ul == 'U')						\
		for (j = 0; j < n; px += (++j)+1)			\
		    *px = *pval;					\
	    else							\
		for (j = 0; j < n; px += n-(j++))			\
		    *px = *pval;					\
	} else {							\
	    if (ul == 'U')						\
		for (j = 0; j < n; px += (++j)+1)			\
		    *px = *(pval++);					\
	    else							\
		for (j = 0; j < n; px += n-(j++))			\
		    *px = *(pval++);					\
	}								\
    } while (0)

    switch (tx) {
    case REALSXP:
	PM_D_S(double, REAL);
	break;
    case LGLSXP:
	PM_D_S(int, LOGICAL);
	break;
    case INTSXP:
	PM_D_S(int, INTEGER);
	break;
    case CPLXSXP:
	PM_D_S(Rcomplex, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_diag_set");
	break;
    }
    
#undef PM_D_S
    
    UNPROTECT(nprotect);
    return res;
}

/* symmpart(x) */
SEXP packedMatrix_symmpart(SEXP from)
{
    static const char *valid[] = {
	"dtpMatrix", "dspMatrix",
	"ltpMatrix", "lspMatrix",
	"ntpMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "packedMatrix_symmpart");
    const char *clf = valid[ivalid];
    if (clf[1] == 's' && clf[0] == 'd')
	return from;

    char clt[] = ".spMatrix";
    clt[0] = (clf[0] != 'z') ? 'd' : 'z';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = GET_SLOT(from, Matrix_uploSym),
	x = GET_SLOT(from, Matrix_xSym);
    
    PROTECT(x = (clf[0] == clt[0]) ? duplicate(x) : coerceVector(x, REALSXP));
    if (clf[0] == 'n')
	na2one(x);

    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_uploSym, uplo);
    SET_SLOT(to, Matrix_xSym, x);
    
    if (clf[1] == 't') {

	int i, j, n = INTEGER(dim)[0];
	char ul = *CHAR(STRING_ELT(uplo, 0)), di = *diag_P(from);

#define PM_SYMMPART_TP(_CTYPE_, _PTR_, _ASSIGN_OFFDIAG_, _ASSIGN_ONDIAG_) \
	do {								\
	    _CTYPE_ *px = _PTR_(x);					\
	    if (ul == 'U') {						\
		for (j = 0; j < n; ++j) {				\
		    for (i = 0; i < j; ++i, ++px)			\
			_ASSIGN_OFFDIAG_;				\
		    ++px;						\
		}							\
		if (di != 'N') {					\
		    px = _PTR_(x);					\
		    for (j = 0; j < n; px += (++j)+1)			\
			_ASSIGN_ONDIAG_;				\
		}							\
	    } else {							\
		for (j = 0; j < n; ++j) {				\
		    ++px;						\
		    for (i = j+1; i < n; ++i, ++px)			\
			_ASSIGN_OFFDIAG_;				\
		}							\
		if (di != 'N') {					\
		    px = _PTR_(x);					\
		    for (j = 0; j < n; px += n-(j++))			\
			_ASSIGN_ONDIAG_;				\
		}							\
	    }								\
	} while (0)

	if (clt[0] != 'z')
	    PM_SYMMPART_TP(double, REAL,
			   *px *= 0.5,
			   *px  = 1.0);
	else
	    PM_SYMMPART_TP(Rcomplex, COMPLEX,
			   do { (*px).r *= 0.5; (*px).i *= 0.5; } while (0),
			   do { (*px).r  = 1.0; (*px).i  = 0.0; } while (0));

#undef PM_SYMMPART_TP
	
	set_symmetrized_DimNames(to, dimnames, -1);
	
    } else { /* clf[1] == 's' */
	
	if (clt[0] == 'z')
	    /* Symmetric part of Hermitian matrix is real part */
	    zeroIm(x);

	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	
    }

    UNPROTECT(2);
    return to;
}

/* skewpart(x) */
SEXP packedMatrix_skewpart(SEXP from)
{
    static const char *valid[] = {
	"dtpMatrix", "dspMatrix",
	"ltpMatrix", "lspMatrix",
	"ntpMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "packedMatrix_skewpart");
    const char *clf = valid[ivalid];

    SEXP to,
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = GET_SLOT(from, Matrix_uploSym),
	x = GET_SLOT(from, Matrix_xSym);
    int n = INTEGER(dim)[0];
    
    if (clf[1] == 't') {
	
	if ((double) n * n > R_XLEN_T_MAX)
	    error(_("attempt to allocate vector of length exceeding "
		    "R_XLEN_T_MAX"));
	SEXP y;
	char ul = *CHAR(STRING_ELT(uplo, 0));
	R_xlen_t upos = 0, lpos = 0;
	int i, j;

#define PM_SKEWPART(_CTYPE_, _PTR_, _ASSIGN_OFFDIAG_, _ASSIGN_ONDIAG_)	\
	do {								\
	    _CTYPE_ *px = _PTR_(x), *py = _PTR_(y);			\
	    if (ul == 'U') {						\
		for (j = 0; j < n; ++j) {				\
		    lpos = j;						\
		    for (i = 0; i < j; ++i) {				\
			_ASSIGN_OFFDIAG_(upos, lpos);			\
			++px; ++upos; lpos += n;			\
		    }							\
		    _ASSIGN_ONDIAG_(upos);				\
		    ++px; upos += n-j;					\
		}							\
	    } else {							\
		for (j = 0; j < n; ++j) {				\
		    upos = lpos;					\
		    _ASSIGN_ONDIAG_(lpos);				\
		    for (i = j+1; i < n; ++i) {				\
			++px; upos += n; ++lpos;			\
			_ASSIGN_OFFDIAG_(lpos, upos);			\
		    }							\
		    ++px; lpos += j+2;					\
		}							\
	    }								\
	} while (0)	    
	
	if (clf[0] != 'z') {
	    
	    PROTECT(to = NEW_OBJECT_OF_CLASS("dgeMatrix"));
	    PROTECT(y = allocVector(REALSXP, (R_xlen_t) n * n));
	    PROTECT(x = coerceVector(x, REALSXP));
	    if (clf[0] == 'n')
		na2one(x);

#define ASSIGN_OFFDIAG_DTP(_UPOS_, _LPOS_)			\
	    do {						\
		py[_UPOS_] = 0.5 * *px;				\
		py[_LPOS_] = -py[_UPOS_];			\
	    } while (0)

#define ASSIGN_ONDIAG_DTP(_DPOS_) py[_DPOS_] = 0.0

	    PM_SKEWPART(double, REAL,
			ASSIGN_OFFDIAG_DTP, ASSIGN_ONDIAG_DTP);
	    UNPROTECT(1);
	    
#undef ASSIGN_OFFDIAG_DTP
#undef ASSIGN_ONDIAG_DTP
	    
	} else { /* clf[0] == 'z' */

	    PROTECT(to = NEW_OBJECT_OF_CLASS("zgeMatrix"));
	    PROTECT(y = allocVector(CPLXSXP, (R_xlen_t) n * n));
	    
#define ASSIGN_OFFDIAG_ZTP(_UPOS_, _LPOS_)			\
	    do {						\
		py[_UPOS_].r = 0.5 * (*px).r;			\
		py[_UPOS_].i = 0.5 * (*px).i;			\
		py[_LPOS_].r = -py[upos].r;			\
		py[_LPOS_].i = -py[upos].i;			\
	    } while (0)

#define ASSIGN_ONDIAG_ZTP(_DPOS_) py[_DPOS_].r = py[_DPOS_].i = 0.0
		
	    PM_SKEWPART(Rcomplex, COMPLEX,
			ASSIGN_OFFDIAG_ZTP, ASSIGN_ONDIAG_ZTP);

#undef ASSIGN_OFFDIAG_ZTP
#undef ASSIGN_ONDIAG_ZTP

	}

#undef PM_SKEWPART
	
	SET_SLOT(to, Matrix_DimSym, dim);
	set_symmetrized_DimNames(to, dimnames, -1);
	SET_SLOT(to, Matrix_xSym, y);
	
    } else { /* clf[1] == 's' */
	
	if (clf[0] != 'z') {
	    /* Skew-symmetric part of symmetric matrix is zero matrix */
	    PROTECT(to = NEW_OBJECT_OF_CLASS("dsCMatrix"));
	    R_xlen_t n1a = (R_xlen_t) n + 1;
	    SEXP p = PROTECT(allocVector(INTSXP, n1a));
	    int *pp = INTEGER(p);
	    Memzero(pp, n1a);
	    SET_SLOT(to, Matrix_pSym, p);
	} else {
	    /* Skew-symmetric part of Hermitian matrix is imaginary part */
	    PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
	    PROTECT(x = duplicate(x));
	    zeroRe(x);
	    SET_SLOT(to, Matrix_xSym, x);
	}

	SET_SLOT(to, Matrix_DimSym, dim);
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	SET_SLOT(to, Matrix_uploSym, uplo);

    }
    
    UNPROTECT(2);
    return to;
}

#define PM_XIJ_TR_UP_NUN(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ <= _J_							\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : _ZERO_)

#define PM_XIJ_TR_UP_UNT(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ < _J_							\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : (_I_ == _J_ ? _ONE_ : _ZERO_))

#define PM_XIJ_TR_LO_NUN(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ >= _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : _ZERO_)

#define PM_XIJ_TR_LO_UNT(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ > _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : (_I_ == _J_ ? _ONE_ : _ZERO_))

#define PM_XIJ_SY_UP(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ <= _J_				       			\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))	       			\
     : *(_X_ + PM_AR21_UP(_J_, _I_)))

#define PM_XIJ_SY_LO(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ >= _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : *(_X_ + PM_AR21_LO(_J_, _I_, _N2_)))

#define PM_SUB1_LOOP(_XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	int i, j;							\
	R_xlen_t k;							\
	if (TYPEOF(index) == INTSXP) {					\
	    int pos, *pindex = INTEGER(index);				\
	    for (k = 0; k < nindex; ++k) {				\
		pos = *(pindex++);					\
		if (pos == NA_INTEGER || pos > nn) {			\
		    *(pres++) = _NA_;					\
		} else {						\
		    pos -= 1; /* 1-index -> 0-index */			\
		    i = pos % n;					\
		    j = pos / n;					\
		    *(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);	\
		}							\
	    }								\
	} else {							\
	    double pos, *pindex = REAL(index);				\
	    R_xlen_t truncpos;						\
	    for (k = 0; k < nindex; ++k) {				\
		pos = *(pindex++);					\
		if (!R_FINITE(pos) || (truncpos = (R_xlen_t) pos) > nn) { \
		    *(pres++) = _NA_;					\
		} else {						\
		    truncpos -= 1; /* 1-index -> 0-index */		\
		    i = truncpos % n;					\
		    j = truncpos / n;					\
		    *(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);	\
		}							\
	    }								\
	}								\
    } while (0)

#define PM_SUB1_END(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_)		\
    do {								\
	_CTYPE_ *pres = _PTR_(res), *px = _PTR_(x);			\
	if (diag == ' ') { /* symmetric */				\
	    if (uplo == 'U') {						\
		PM_SUB1_LOOP(PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB1_LOOP(PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else if (diag == 'N') { /* non-unit triangular */		\
	    if (uplo == 'U') {						\
		PM_SUB1_LOOP(PM_XIJ_TR_UP_NUN, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB1_LOOP(PM_XIJ_TR_LO_NUN, _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else { /* unit triangular */					\
	    if (uplo == 'U') {						\
		PM_SUB1_LOOP(PM_XIJ_TR_UP_UNT, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB1_LOOP(PM_XIJ_TR_LO_UNT, _NA_, _ZERO_, _ONE_);	\
	    }								\
	}    								\
    } while (0)

#define PM_SUB1								\
    do {								\
	SEXPTYPE tx;							\
	SEXP x = GET_SLOT(obj, Matrix_xSym),				\
	    res = PROTECT(allocVector(tx = TYPEOF(x), nindex));		\
	char uplo = *uplo_P(obj), diag = *Diag_P(obj);			\
									\
	switch (tx) {							\
	case REALSXP: /* d..Matrix */					\
	    PM_SUB1_END(double, REAL, NA_REAL, 0.0, 1.0);		\
	    break;							\
	case LGLSXP: /* [ln]..Matrix */					\
	    PM_SUB1_END(int, LOGICAL, NA_LOGICAL, 0, 1);		\
	    break;							\
	case INTSXP: /* i..Matrix */					\
	    PM_SUB1_END(int, INTEGER, NA_INTEGER, 0, 1);		\
	    break;							\
	case CPLXSXP: /* z..Matrix */					\
	{								\
	    Rcomplex na, zero, one;					\
	    na.r = NA_REAL; zero.r = 0.0; one.r = 1.0;			\
	    na.i = NA_REAL; zero.i = 0.0; one.i = 0.0;			\
	    PM_SUB1_END(Rcomplex, COMPLEX, na, zero, one);		\
	    break;							\
	}								\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_sub1");	\
	    break;							\
	}								\
	UNPROTECT(1);							\
	return res;							\
    } while (0)
	  
/* 'x[i]' where 'i' is an integer or double vector with elements 
   greater than or equal to 1
*/
SEXP packedMatrix_sub1(SEXP obj, SEXP index)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    if ((double) n * n > R_XLEN_T_MAX)
	error(_("indexing n-by-n packedMatrix is not supported "
		"for n*n exceeding R_XLEN_T_MAX"));
    R_xlen_t n2 = (R_xlen_t) n * 2, nn = (R_xlen_t) n * n,
	nindex = XLENGTH(index);
    PM_SUB1;
}

#undef PM_SUB1_LOOP

#define PM_SUB1_LOOP(_XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	int i, j, k, *pi = INTEGER(index), *pj = pi + nindex;		\
	for (k = 0; k < nindex; ++k) {					\
	    i = *(pi++);						\
	    j = *(pj++);						\
	    if (i == NA_INTEGER || j == NA_INTEGER) {			\
		*(pres++) = _NA_;					\
	    } else {							\
		i -= 1; /* 1-index -> 0-index */			\
		j -= 1;							\
		*(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);		\
	    }								\
	}								\
    } while (0)

/* 'x[i]' where 'i' is a 2-column integer matrix supplying 
   integers in 'c(1:n, NA)' _only_
*/
SEXP packedMatrix_sub1_mat(SEXP obj, SEXP index)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0],
	nindex = INTEGER(getAttrib(index, R_DimSymbol))[0];
    if ((double) n * n > R_XLEN_T_MAX)
	error(_("indexing n-by-n packedMatrix is not supported "
		"for n*n exceeding R_XLEN_T_MAX"));
    R_xlen_t n2 = (R_xlen_t) n * 2;
    PM_SUB1;
}

#undef PM_SUB1
#undef PM_SUB1_LOOP

#define PM_SUB2_LOOP(_FOR_, _XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
        int i, j, ki, kj;						\
	for (kj = 0; kj < nj; ++kj) {					\
	    if (mj) {							\
		j = kj;							\
	    } else {							\
		j = pj[kj];						\
		if (j == NA_INTEGER) {					\
		    _FOR_ {						\
			*(px1++) = _NA_;				\
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
			*(px1++) = _NA_;				\
			continue;					\
		    } else {						\
			i -= 1;						\
		    }							\
		}							\
		*(px1++) = _XIJ_(px0, i, j, n2, _ZERO_, _ONE_);		\
	    }								\
	    if (do_cn) {						\
		SET_STRING_ELT(cn1, kj, STRING_ELT(cn0, j));		\
	    }								\
	}								\
	if (do_rn) {							\
	    for (ki = 0; ki < ni; ++ki) {				\
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

#define PM_SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);			\
	if (diag == ' ') { /* symmetric */				\
	    if (uplo == 'U') {						\
		if (do_sp) {						\
		    PM_SUB2_LOOP(for (ki = 0; ki <= kj; ++ki),		\
				 PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
		} else {						\
		    PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				 PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
		}							\
	    } else {							\
		if (do_sp) {						\
		    PM_SUB2_LOOP(for (ki = kj; ki < ni; ++ki),		\
				 PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
		} else { 						\
		    PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				 PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
		}							\
	    }								\
	} else if (diag == 'N') { /* non-unit triangular */		\
	    if (uplo == 'U') {						\
		PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_UP_NUN, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_LO_NUN, _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else { /* unit triangular */					\
	    if (uplo == 'U') {						\
		PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_UP_UNT, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_LO_UNT, _NA_, _ZERO_, _ONE_);	\
	    }								\
	}								\
	SET_SLOT(res, Matrix_xSym, x1);					\
    } while (0)

/* 'x[i, ]', 'x[, j]', and 'x[i, j]' where 'i' and 'j' are integer vectors 
   supplying integers in 'c(1:n, NA)' _only_ ... NULL indicates missingness
*/
SEXP packedMatrix_sub2(SEXP obj, SEXP index1, SEXP index2, SEXP drop)
{
    Rboolean mi = isNull(index1), mj = isNull(index2);
    int *pi = NULL, *pj = NULL, n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    R_xlen_t ni, nj, n2 = (R_xlen_t) n * 2;
    
    if ((double) n * n > R_XLEN_T_MAX)
	error(_("indexing of n-by-n packedMatrix with n*n "
		"exceeding R_XLEN_T_MAX is not supported"));
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

    static const char *valid[] = {
	"dspMatrix", "lspMatrix", "nspMatrix",
	"dtpMatrix", "ltpMatrix", "ntpMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid), nprotect = 0;
    if (ivalid < 0) {
	ERROR_INVALID_CLASS(class_P(obj), "packedMatrix_sub2");
    }
    char uplo = *uplo_P(obj), diag = *Diag_P(obj);

    /* Initialize result of same type but "general" class,
       except for symmetric indexing of symmetric matrix, 
       when class is retained also
    */
    SEXP res;
    Rboolean do_sp = (ivalid < 3 && !mi && !mj && ni == nj &&
		      memcmp(pi, pj, ni * sizeof(int)) == 0);
    if (do_sp) {
	PROTECT(res = NEW_OBJECT_OF_CLASS(valid[ivalid]));
	SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));
	++nprotect;
    } else {
	char cl[] = ".geMatrix";
	cl[0] = valid[ivalid][0];
	PROTECT(res = NEW_OBJECT_OF_CLASS(cl));
	++nprotect;
    }

    /* Set 'Dim' slot */
    int *pdim = INTEGER(GET_SLOT(res, Matrix_DimSym));
    pdim[0] = (int) ni;
    pdim[1] = (int) nj;

    /* Set 'Dimnames' slot and 'names(Dimnames)' */
    SEXP rn0, rn1 = R_NilValue, cn0, cn1 = R_NilValue,
	dn0 = GET_SLOT(obj, Matrix_DimNamesSym),
	dn1 = PROTECT(GET_SLOT(res, Matrix_DimNamesSym));
    ++nprotect;
    if (diag == ' ') { /* symmetric */
    	int J;
	if (isNull(rn0 = cn0 = VECTOR_ELT(dn0, J = 1))) {
	    rn0 = cn0 = VECTOR_ELT(dn0, J = 0);
	}
	SEXP s;
	if (!isNull(s = getAttrib(dn0, R_NamesSymbol))) {
	    SEXP ndn1 = PROTECT(allocVector(STRSXP, 2));
	    SET_STRING_ELT(ndn1, 0, s = STRING_ELT(s, J));
	    SET_STRING_ELT(ndn1, 1, s);
	    setAttrib(dn1, R_NamesSymbol, ndn1);
	    UNPROTECT(1);
	}
    } else { /* triangular */
	rn0 = VECTOR_ELT(dn0, 0);
	cn0 = VECTOR_ELT(dn0, 1);
	SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
	if (!isNull(ndn0)) {
	    setAttrib(dn1, R_NamesSymbol, ndn0);
	}
    }
    Rboolean has_rn, has_cn, do_rn, do_cn;
    has_rn = !isNull(rn0) && ni > 0;
    has_cn = !isNull(cn0) && nj > 0;
    do_rn = do_cn = FALSE;
    if (has_rn) {
	if (mi) {
	    SET_VECTOR_ELT(dn1, 0, rn0);
	} else {
	    PROTECT(rn1 = allocVector(STRSXP, ni));
	    SET_VECTOR_ELT(dn1, 0, rn1);
	    ++nprotect;
	    do_rn = TRUE;
	}
    }
    if (has_cn && !do_sp) {
	if (mj) {
	    SET_VECTOR_ELT(dn1, 1, cn0);
	} else {
	    PROTECT(cn1 = allocVector(STRSXP, nj));
	    SET_VECTOR_ELT(dn1, 1, cn1);
	    ++nprotect;
	    do_cn = TRUE;
	}
    }

    /* Set 'x' slot */
    SEXPTYPE tx;
    R_xlen_t nx = (do_sp ? ni + (ni * (ni - 1)) / 2 : ni * nj);
    SEXP x0 = GET_SLOT(obj, Matrix_xSym),
	x1 = PROTECT(allocVector(tx = TYPEOF(x0), nx));
    ++nprotect;
    switch (tx) {
    case REALSXP: /* d..Matrix */
	PM_SUB2(double, REAL, NA_REAL, 0.0, 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PM_SUB2(int, LOGICAL, NA_LOGICAL, 0, 1);
	break;
    case INTSXP: /* i..Matrix */
	PM_SUB2(int, INTEGER, NA_INTEGER, 0, 1);
	break;
    case CPLXSXP: /* z..Matrix */
    {
	Rcomplex na, zero, one;
	na.r = NA_REAL; zero.r = 0.0; one.r = 1.0;
	na.i = NA_REAL; zero.i = 0.0; one.i = 0.0;
	PM_SUB2(Rcomplex, COMPLEX, na, zero, one);
	break;
    }
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_sub2");
	break;
    }

    /* Drop dimensions in this special case */
    if (asLogical(drop) != 0 && (ni == 1 || nj == 1)) {
	PROTECT(res = GET_SLOT(res, Matrix_xSym));
	++nprotect;
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
