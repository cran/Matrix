#include "unpackedMatrix.h"

/* pack(x), returning packedMatrix */
SEXP unpackedMatrix_pack(SEXP from, SEXP strict, SEXP tr_if_ge, SEXP up_if_ge)
{
    static const char *valid_from[] = {
	/*  0 */ "Cholesky", "BunchKaufman", /* must match before dtrMatrix */
	/*  2 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/*  5 */ "corMatrix", "dpoMatrix", /* must match before dsyMatrix */
	/*  7 */ "dsyMatrix", "lsyMatrix", "nsyMatrix",
	/* 10 */ "dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    static const char *valid_to[] = {
	/*  0 */ "pCholesky", "pBunchKaufman",
	/*  2 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/*  5 */ "dppMatrix", "dppMatrix", /* no pcorMatrix _yet_ */
	/*  7 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid_from);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "unpackedMatrix_pack");
    if (asLogical(strict) == 0) {
	if (ivalid < 2)
	    ivalid = 2; /* Cholesky,BunchKaufman->dtpMatrix */
	else if (ivalid == 5 || ivalid == 6)
	    ivalid = 7; /* corMatrix,dpoMatrix->dspMatrix */
    }

    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int *pdim = INTEGER(dim), n = pdim[0], shift = 0;
    if (ivalid >= 10) {
	/* .geMatrix */
	if (pdim[1] != n)
	    error(_("attempt to pack non-square matrix"));
	shift = (asLogical(tr_if_ge) != 0) ?  3+2+3 :      3;
	/*                                 ? ge->tp : ge->sp */
    }

    SEXPTYPE tx;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(valid_to[ivalid - shift])),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x_from = GET_SLOT(from, Matrix_xSym),
	x_to = PROTECT(allocVector(tx = TYPEOF(x_from), PM_LENGTH(n))),
	uplo;

    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    SET_SLOT(to, Matrix_xSym, x_to);
    if (ivalid >= 10) {
	/* .geMatrix */
	uplo = mkString((asLogical(up_if_ge) != 0) ? "U" : "L");
    } else {
	/* .(sy|tr)Matrix */
	uplo = GET_SLOT(from, Matrix_uploSym);
	if (ivalid < 5) {
	    /* .trMatrix */
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
	    if (ivalid == 1)
		/* BunchKaufman */
		SET_SLOT(to, Matrix_permSym, GET_SLOT(from, Matrix_permSym));
	} else {
	    /* .syMatrix */
	    SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
	}
    }
    SET_SLOT(to, Matrix_uploSym, uplo);
    char ul = *CHAR(STRING_ELT(uplo, 0));

#define PACK(_PREFIX_, _PTR_)						\
    _PREFIX_ ## dense_pack(_PTR_(x_to), _PTR_(x_from), n, ul, 'N')
    
    switch (tx) {
    case REALSXP: /* d..Matrix */
	PACK(d, REAL);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PACK(i, LOGICAL);
	break;
    case INTSXP: /* i..Matrix */
	PACK(i, INTEGER);
	break;
    case CPLXSXP: /* z..Matrix */
    	PACK(z, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_pack");
	break;
    }

#undef PACK

    UNPROTECT(2);
    return to;
}

/* forceSymmetric(x, uplo), returning .syMatrix */
SEXP unpackedMatrix_force_symmetric(SEXP from, SEXP uplo_to)
{
    static const char *valid[] = {
	"dsyMatrix", "lsyMatrix", "nsyMatrix", /* be fast */
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "unpackedMatrix_force_symmetric");
    const char *clf = valid[ivalid];
    
    SEXP uplo_from;
    char ulf = 'U', ult;
    if (clf[1] != 'g') {
	/* .(sy|tr)Matrix */
	uplo_from = GET_SLOT(from, Matrix_uploSym);
	ulf = *CHAR(STRING_ELT(uplo_from, 0));
    }
    ult = (isNull(uplo_to)) ? ulf : *CHAR(asChar(uplo_to));
    if (clf[1] == 's') {
	/* .syMatrix */
	if (ulf == ult)
	    return from;
	SEXP to = PROTECT(unpackedMatrix_transpose(from));
	if (clf[0] == 'z') {
	    /* Need _conjugate_ transpose */
	    SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	    conjugate(x);
	    UNPROTECT(1);
	}
	UNPROTECT(1);
	return to;
    }
    
    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to symmetrize a non-square matrix"));

    /* Now handling just square .(tr|ge)Matrix ... */
    
    char clt[] = ".syMatrix";
    clt[0] = clf[0];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x_from = GET_SLOT(from, Matrix_xSym);
    
    if (clf[1] == 'g' || ulf == ult) {
	/* .geMatrix or .trMatrix with correct uplo */
	SET_SLOT(to, Matrix_xSym, x_from);
    } else {
	/* .trMatrix with incorrect uplo */
	char di = *diag_P(from);
	SEXPTYPE tx = TYPEOF(x_from);
	R_xlen_t nx = XLENGTH(x_from);
	SEXP x_to = PROTECT(allocVector(tx, nx));

#define COPY_DIAGONAL(_PREFIX_, _PTR_)					\
	do {								\
	    Memzero(_PTR_(x_to), nx);					\
	    _PREFIX_ ## dense_unpacked_copy_diagonal(			\
		_PTR_(x_to), _PTR_(x_from), n, nx, 'U' /* unused */, di); \
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
	   ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_force_symmetric");
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

#define UPM_IS_SY(_RES_, _X_, _N_, _LDENSE_, _WHAT_, _METHOD_)		\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_unpacked_is_symmetric(REAL(_X_), _N_);	\
	    break;							\
	case LGLSXP:							\
	    _RES_ = (_LDENSE_						\
		     ? ldense_unpacked_is_symmetric(LOGICAL(_X_), _N_)	\
		     : ndense_unpacked_is_symmetric(LOGICAL(_X_), _N_)); \
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_unpacked_is_symmetric(INTEGER(_X_), _N_);	\
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_unpacked_is_symmetric(COMPLEX(_X_), _N_);	\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE(_WHAT_, TYPEOF(_X_), _METHOD_);		\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

#define UPM_IS_TR(_RES_, _X_, _N_, _UPLO_, _WHAT_, _METHOD_)		\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_unpacked_is_triangular(REAL(_X_), _N_, _UPLO_); \
	    break;							\
	case LGLSXP:							\
	    _RES_ = idense_unpacked_is_triangular(LOGICAL(_X_), _N_, _UPLO_); \
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_unpacked_is_triangular(INTEGER(_X_), _N_, _UPLO_); \
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_unpacked_is_triangular(COMPLEX(_X_), _N_, _UPLO_); \
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE(_WHAT_, TYPEOF(_X_), _METHOD_);		\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

#define UPM_IS_DI(_RES_, _X_, _N_, _WHAT_, _METHOD_)			\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_unpacked_is_diagonal(REAL(_X_), _N_);	\
	    break;							\
	case LGLSXP:							\
	    _RES_ = idense_unpacked_is_diagonal(LOGICAL(_X_), _N_);	\
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_unpacked_is_diagonal(INTEGER(_X_), _N_);	\
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_unpacked_is_diagonal(COMPLEX(_X_), _N_);	\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE(_WHAT_, TYPEOF(_X_), _METHOD_);		\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

/* isSymmetric(x, tol = 0, checkDN) */
/* FIXME: not checking for real diagonal in complex case */
SEXP unpackedMatrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    static const char *valid[] = {
	"dsyMatrix", "lsyMatrix", "nsyMatrix", /* be fast */
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0) {
	ERROR_INVALID_CLASS(class_P(obj), "unpackedMatrix_is_symmetric");
	return R_NilValue;
    } else if (ivalid < 3) {
	/* .syMatrix: symmetric by definition */
	return ScalarLogical(1);
    } else {
	/* .(ge|tr)Matrix */
	int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
	if (ivalid >= 6 && pdim[1] != n)
	    return ScalarLogical(0);
	if (asLogical(checkDN) != 0 &&
	    !DimNames_is_symmetric(GET_SLOT(obj, Matrix_DimNamesSym)))
	    return ScalarLogical(0);
	Rboolean res = FALSE;
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (ivalid < 6) {
	    /* .trMatrix: symmetric iff diagonal (upper _and_ lower tri.) */
	    char ul = (*uplo_P(obj) == 'U') ? 'L' : 'U';
	    UPM_IS_TR(res, x, n, ul,
		      "'x' slot", "unpackedMatrix_is_symmetric");
	} else {
	    /* .geMatrix: need to do a complete symmetry check */
	    UPM_IS_SY(res, x, n, ivalid == 7,
		      "'x' slot", "unpackedMatrix_is_symmetric");
	}
	return ScalarLogical(res);
    }
}

/* isSymmetric(x, tol = 0) */
SEXP matrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    if (asLogical(checkDN) != 0) {
	SEXP dn = getAttrib(obj, R_DimNamesSymbol);
	if (!isNull(dn) && !DimNames_is_symmetric(dn))
	    return ScalarLogical(0);
    }
    Rboolean res = FALSE;
    UPM_IS_SY(res, obj, n, 1, "matrix", "matrix_is_symmetric");
    return ScalarLogical(res);
}

#define RETURN_TRUE_OF_KIND(_KIND_)					\
    do {								\
	SEXP ans = PROTECT(allocVector(LGLSXP, 1));			\
	LOGICAL(ans)[0] = 1;						\
	setAttrib(ans, install("kind"), _KIND_);			\
	UNPROTECT(1);							\
	return ans;							\
    } while (0)

#define RETURN_GE_IS_TR(_X_, _N_, _UPPER_, _WHAT_, _METHOD_)		\
    do {								\
	Rboolean res = FALSE;						\
	if (_UPPER_ == NA_LOGICAL) {					\
	    UPM_IS_TR(res, _X_, _N_, 'U', _WHAT_, _METHOD_);		\
	    if (res)							\
	        RETURN_TRUE_OF_KIND(mkString("U"));			\
	    UPM_IS_TR(res, _X_, _N_, 'L', _WHAT_, _METHOD_);		\
	    if (res)							\
		RETURN_TRUE_OF_KIND(mkString("L"));			\
	} else {							\
	    UPM_IS_TR(res, _X_, _N_, (_UPPER_ != 0) ? 'U' : 'L',	\
		      _WHAT_, _METHOD_);				\
	    if (res)							\
		return ScalarLogical(1);				\
	}								\
	return ScalarLogical(0);					\
    } while (0)

/* isTriangular(x, upper) */
SEXP unpackedMatrix_is_triangular(SEXP obj, SEXP upper)
{
    static const char *valid[] = {
	"dtrMatrix", "ltrMatrix", "ntrMatrix", /* be fast */
	"dsyMatrix", "lsyMatrix", "nsyMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "unpackedMatrix_is_triangular");
    
    SEXP uplo;
    char ul;
    int need_upper = asLogical(upper);
    if (ivalid < 6) {
	uplo = GET_SLOT(obj, Matrix_uploSym);
	ul = *CHAR(STRING_ELT(uplo, 0));
    }

#define IF_DIAGONAL							\
    Rboolean res = FALSE;						\
    SEXP x = GET_SLOT(obj, Matrix_xSym);				\
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];			\
    UPM_IS_TR(res, x, n, (ul == 'U') ? 'L' : 'U',			\
	      "'x' slot", "unpackedMatrix_is_triangular");		\
    if (res)
    
    if (ivalid < 3) {
	/* .trMatrix: be fast if 'upper', 'uplo' agree; else need diagonal */
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
	return ScalarLogical(0);
    } else if (ivalid < 6) {
	/* .syMatrix: triangular iff diagonal (upper _and_ lower tri.) */
	IF_DIAGONAL {
	    if (need_upper == NA_LOGICAL) {
		RETURN_TRUE_OF_KIND(mkString("U"));
	    } else {
		return ScalarLogical(1);
	    }
	}
	return ScalarLogical(0);

#undef IF_DIAGONAL
	
    } else {
	/* .geMatrix: need to do a complete triangularity check */
	int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
	if (pdim[1] != n)
	    return ScalarLogical(0);
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	RETURN_GE_IS_TR(x, n, need_upper,
			"'x' slot", "unpackedMatrix_is_triangular");
    }
}

/* isTriangular(x, upper) */
SEXP matrix_is_triangular(SEXP obj, SEXP upper)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    int need_upper = asLogical(upper);
    RETURN_GE_IS_TR(obj, n, need_upper,
		    "matrix", "matrix_is_triangular");
}

#undef RETURN_GE_IS_TR
#undef RETURN_TRUE

/* isDiagonal(x) */
SEXP unpackedMatrix_is_diagonal(SEXP obj)
{
    static const char *valid[] = {
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dsyMatrix", "lsyMatrix", "nsyMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "unpackedMatrix_is_diagonal");
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    Rboolean res = FALSE;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (ivalid < 6) {
	/* .(sy|tr)Matrix: diagonal iff stored triangle is zero off diagonal */
	char ul = (*uplo_P(obj) == 'U') ? 'L' : 'U';
	UPM_IS_TR(res, x, n, ul, "'x' slot", "unpackedMatrix_is_diagonal");
    } else {
	/* .geMatrix: need to do a complete diagonality check */
	UPM_IS_DI(res, x, n, "'x' slot", "unpackedMatrix_is_diagonal");
    }
    return ScalarLogical(res);
}

/* isDiagonal(x) */
SEXP matrix_is_diagonal(SEXP obj)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    Rboolean res = FALSE;
    UPM_IS_DI(res, obj, n, "matrix", "matrix_is_diagonal");
    return ScalarLogical(res);
}

#undef UPM_IS_DI
#undef UPM_IS_TR
#undef UPM_IS_SY

/* t(x), typically preserving class */
/* MJ: Technically no need to do full transpose of .(sy|tr)Matrix ...  */
/*     but then identical(.@x, t(t(.))@x) can be FALSE ...             */
SEXP unpackedMatrix_transpose(SEXP from)
{
    static const char *valid[] = {
	/*  0 */ "Cholesky", "BunchKaufman",
	/*  2 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/*  5 */ "corMatrix", "dpoMatrix",
	/*  7 */ "dsyMatrix", "lsyMatrix", "nsyMatrix",
	/* 10 */ "dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "unpackedMatrix_transpose");
    if (ivalid == 1)
	ivalid = 2; /* BunchKaufman->dtrMatrix */
    const char *cl = valid[ivalid];
    
    SEXPTYPE tx;
    R_xlen_t nx;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x0 = GET_SLOT(from, Matrix_xSym),
	x1 = PROTECT(allocVector(tx = TYPEOF(x0), nx = XLENGTH(x0)));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    char ul = (ivalid < 10) ? *uplo_P(from) : '\0';
    
    /* Preserve or reverse 'Dim' slot (preserving if square) */
    if (m == n) {
	SET_SLOT(to, Matrix_DimSym, dim);
    } else {
	pdim = INTEGER(GET_SLOT(to, Matrix_DimSym));
	pdim[0] = n;
	pdim[1] = m;
    }
    /* Reverse or preserve 'Dimnames' slot (preserving if symmetric) */
    if (ivalid < 5 || ivalid >= 10)
	set_reversed_DimNames(to, dimnames);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (ivalid < 10) {
	/* Toggle 'uplo' slot */
	SET_SLOT(to, Matrix_uploSym, mkString((ul == 'U') ? "L" : "U"));
	if (ivalid < 5)
	    /* Preserve 'diag' slot */
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
	else
	    /* Preserve 'factors' slot */
	    SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
	if (ivalid == 5)
	    /* Preserve 'sd' slot */
	    SET_SLOT(to, install("sd"), GET_SLOT(from, install("sd")));
    }
    /* NB: Nothing to do for 'factors' slot: prototype is already list() ...
       FIXME: However, it would be much better to also "transpose" each 
       factorization ... */

#define UPM_T(_CTYPE_, _PTR_)						\
    do {								\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);			\
	int i, j;							\
	/* if (ivalid < 10) { */					\
	/*     Memzero(px1, (size_t) n * n); */				\
	/*     R_xlen_t upos = 0, lpos = 0; */				\
	/*     if (ul == 'U') */					\
	/* 	for (j = 0; j < n; upos = (lpos += (++j))) */		\
	/* 	    for (i = j; i < n; ++i, ++lpos, upos += n) */	\
	/* 		px1[lpos] = px0[upos]; */			\
	/*     else */							\
	/* 	for (j = 0; j < n; upos = (lpos += (++j))) */		\
	/* 	    for (i = j; i < n; ++i, ++lpos, upos += n) */	\
	/* 		px1[upos] = px0[lpos]; */			\
	/* } else { */							\
	R_xlen_t nx1s = nx - 1;						\
	for (j = 0; j < m; ++j, px0 -= nx1s)				\
	    for (i = 0; i < n; ++i, px0 += m)				\
		*(px1++) = *px0;					\
	/* } */								\
    } while (0)

    /* Permute 'x' slot */
    switch (tx) {
    case REALSXP: /* d..Matrix */
	UPM_T(double, REAL);
	break;
    case LGLSXP: /* [ln]..Matrix */
	UPM_T(int, LOGICAL);
	break;
    case INTSXP: /* i..Matrix */
	UPM_T(int, INTEGER);
	break;
    case CPLXSXP: /* z..Matrix */
	UPM_T(Rcomplex, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_transpose");
	break;
    }
    SET_SLOT(to, Matrix_xSym, x1);

#undef UPM_T

    UNPROTECT(2);
    return to;
}

/* diag(x, names) */
SEXP unpackedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL)
	error(_("'names' must be TRUE or FALSE"));

    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
    SEXPTYPE tx;
    SEXP x = GET_SLOT(obj, Matrix_xSym),
	res = PROTECT(allocVector(tx = TYPEOF(x), r));
    char ul = *Uplo_P(obj), di = *Diag_P(obj);
    
#define UPM_D_G(_CTYPE_, _PTR_, _ONE_)					\
    do {								\
	_CTYPE_ *pres = _PTR_(res);					\
	int j;								\
	if (di == 'U') {						\
	    for (j = 0; j < r; ++j)					\
		*(pres++) = _ONE_;					\
	} else {							\
	    _CTYPE_ *px = _PTR_(x);					\
	    R_xlen_t m1a = ((R_xlen_t) m) + 1;				\
	    for (j = 0; j < r; ++j, px += m1a)				\
		*(pres++) = *px;					\
	}								\
    } while (0)
    
    switch (tx) {
    case REALSXP: /* d..Matrix */
	UPM_D_G(double, REAL, 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	UPM_D_G(int, LOGICAL, 1);
	break;
    case INTSXP: /* i..Matrix */
	UPM_D_G(int, INTEGER, 1);
	break;
    case CPLXSXP: /* z..Matrix */
	UPM_D_G(Rcomplex, COMPLEX, Matrix_zone);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_diag_get");
	break;
    }

#undef UPM_D_G

    if (do_nms) {
	/* NB: The logic here must be adjusted once the validity method 
	       for 'symmetricMatrix' enforces symmetric 'Dimnames' */
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym),
	    rn = VECTOR_ELT(dn, 0),
	    cn = VECTOR_ELT(dn, 1);
	if (isNull(cn)) {
	    if (ul != ' ' && di == ' ' && !isNull(rn))
		setAttrib(res, R_NamesSymbol, rn);
	} else {
	    if (ul != ' ' && di == ' ')
		setAttrib(res, R_NamesSymbol, cn);
	    else if (!isNull(rn) &&
		     (rn == cn || equal_string_vectors(rn, cn, r)))
		setAttrib(res, R_NamesSymbol, (r == m) ? rn : cn);
	}
    }
    UNPROTECT(1);
    return res;
}

/* diag(x) <- value */
SEXP unpackedMatrix_diag_set(SEXP obj, SEXP val)
{
    SEXP dim = GET_SLOT(obj, Matrix_DimSym);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
    R_xlen_t nv = XLENGTH(val);
    if (nv != 1 && nv != r)
	error(_("replacement diagonal has wrong length"));

    SEXP x = GET_SLOT(obj, Matrix_xSym);
    SEXPTYPE tx = TYPEOF(x), tv = TYPEOF(val);
    if (tx < LGLSXP || tx > CPLXSXP)
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_diag_set");
    if (tv < LGLSXP || tv > REALSXP)
	/* Upper bound can become CPLXSXP once we have proper zMatrix */
	error(_("replacement diagonal has incompatible type \"%s\""),
	      type2char(tv));
    
    static const char *valid[] = {
	"dgeMatrix", "dsyMatrix", "dtrMatrix",
	"lgeMatrix", "lsyMatrix", "ltrMatrix",
	"ngeMatrix", "nsyMatrix", "ntrMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "unpackedMatrix_diag_set");

    SEXP res;
    int nprotect = 0;
    
    /* Allocate and coerce as necessary */
    if (tv <= tx) {
	if (tv < tx) {
	    PROTECT(val = coerceVector(val, tv = tx));
	    ++nprotect;
	}
	PROTECT(res = NEW_OBJECT_OF_CLASS(valid[ivalid]));
	PROTECT(x = duplicate(x));
	nprotect += 2;
    } else { /* tv > tx */
	/* dMatrix result is only possibility until we have proper [iz]Matrix */
	if (tv < REALSXP) {
	    PROTECT(val = coerceVector(val, tv = REALSXP));
	    ++nprotect;
	}
	char cl[] = "d..Matrix";
	cl[1] = valid[ivalid][1];
	cl[2] = valid[ivalid][2];
	PROTECT(res = NEW_OBJECT_OF_CLASS(cl));
	PROTECT(x = coerceVector(x, tx = tv));
	++nprotect;
    }
    SET_SLOT(res, Matrix_xSym, x);
    
    /* Transfer slots other than 'x', 'diag', and 'factors';
       latter two should keep their prototypes "N" and list() */
    SET_SLOT(res, Matrix_DimSym, dim);
    SET_SLOT(res, Matrix_DimNamesSym, GET_SLOT(obj, Matrix_DimNamesSym));
    if (R_has_slot(res, Matrix_uploSym))
	SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));

#define UPM_D_S(_CTYPE_, _PTR_)						\
    do {								\
	_CTYPE_ *px = _PTR_(x), *pval = _PTR_(val);			\
	R_xlen_t m1a = (R_xlen_t) m + 1;				\
	int j;								\
	if (nv == 1)							\
	    for (j = 0; j < r; ++j, px += m1a)				\
		*px = *pval;						\
	else								\
	    for (j = 0; j < r; ++j, px += m1a)				\
		*px = *(pval++);					\
    } while (0)
    
    switch (tx) {
    case REALSXP:
	UPM_D_S(double, REAL);
	break;
    case LGLSXP:
	UPM_D_S(int, LOGICAL);
	break;
    case INTSXP:
	UPM_D_S(int, INTEGER);
	break;
    case CPLXSXP:
	UPM_D_S(Rcomplex, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_diag_set");
	break;
    }

#undef UPM_D_S
    
    UNPROTECT(nprotect);
    return res;
}

/* symmpart(x) */
SEXP unpackedMatrix_symmpart(SEXP from)
{
    static const char *valid[] = {
	"dgeMatrix", "dtrMatrix", "dsyMatrix",
	"lgeMatrix", "ltrMatrix", "lsyMatrix",
	"ngeMatrix", "ntrMatrix", "nsyMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "unpackedMatrix_symmpart");
    const char *clf = valid[ivalid];
    if (clf[1] == 's' && clf[0] == 'd')
	return from;

    char clt[] = ".syMatrix";
    clt[0] = (clf[0] != 'z') ? 'd' : 'z';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = (clf[1] != 'g') ? GET_SLOT(from, Matrix_uploSym) : R_NilValue,
	x = GET_SLOT(from, Matrix_xSym);
    
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)	
	error(_("attempt to get symmetric part of non-square matrix"));
    
    PROTECT(x = (clf[0] == clt[0]) ? duplicate(x) : coerceVector(x, REALSXP));
    if (clf[0] == 'n')
	na2one(x);

    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_xSym, x);

    if (clf[1] == 'g') {

	int i, j;
	R_xlen_t upos = 0, lpos = 0;
	
#define UPM_SYMMPART_GE(_CTYPE_, _PTR_, _ASSIGN_OFFDIAG_)	\
	do {							\
	    _CTYPE_ *px = _PTR_(x);				\
	    for (j = 0; j < n; ++j) {				\
		for (i = j+1; i < n; ++i) {			\
		    upos += n; ++lpos;				\
		    _ASSIGN_OFFDIAG_(upos, lpos);		\
		}						\
		upos = (lpos += j+2);				\
	    }							\
	} while (0)
	
#define ASSIGN_OFFDIAG_DGE(_UPOS_, _LPOS_)	\
	do {					\
	    px[_UPOS_] += px[_LPOS_];		\
	    px[_UPOS_] *= 0.5;			\
	} while (0)

#define ASSIGN_OFFDIAG_ZGE(_UPOS_, _LPOS_)	\
	do {					\
	    px[_UPOS_].r += px[_LPOS_].r;	\
	    px[_UPOS_].i += px[_LPOS_].i;	\
	    px[_UPOS_].r *= 0.5;		\
	    px[_UPOS_].i *= 0.5;		\
	} while (0)
	
	if (clf[0] != 'z')
	    UPM_SYMMPART_GE(double, REAL, ASSIGN_OFFDIAG_DGE);
	else
	    UPM_SYMMPART_GE(Rcomplex, COMPLEX, ASSIGN_OFFDIAG_ZGE);
	
	set_symmetrized_DimNames(to, dimnames, -1);
	
    } else if (clf[1] == 't') {

	int i, j;
	char ul = *CHAR(STRING_ELT(uplo, 0)), di = *diag_P(from);

#define UPM_SYMMPART_TR(_CTYPE_, _PTR_, _ASSIGN_OFFDIAG_, _ASSIGN_ONDIAG_) \
	do {								\
	    _CTYPE_ *px = _PTR_(x);					\
	    if (ul == 'U') {						\
		for (j = 0; j < n; ++j) {				\
		    for (i = 0; i < j; ++i, ++px)			\
			_ASSIGN_OFFDIAG_;				\
		    px += n-j;						\
		}							\
	    } else {							\
		for (j = 0; j < n; ++j) {				\
		    px += j+1;						\
		    for (i = j+1; i < n; ++i, ++px)			\
			_ASSIGN_OFFDIAG_;				\
		}							\
	    }								\
	    if (di != 'N') {						\
		R_xlen_t n1a = (R_xlen_t) n + 1;			\
		px = _PTR_(x);						\
		for (j = 0; j < n; ++j, px += n1a)			\
		    _ASSIGN_ONDIAG_;					\
	    }								\
	} while (0)
	
	if (clt[0] != 'z')
	    UPM_SYMMPART_TR(double, REAL,
			    *px *= 0.5,
			    *px  = 1.0);
	else
	    UPM_SYMMPART_TR(Rcomplex, COMPLEX,
			    do { (*px).r *= 0.5; (*px).i *= 0.5; } while (0),
			    do { (*px).r  = 1.0; (*px).i  = 0.0; } while (0));
	
	set_symmetrized_DimNames(to, dimnames, -1);
	SET_SLOT(to, Matrix_uploSym, uplo);
	
    } else { /* clf[1] == 's' */

	if (clt[0] == 'z')
	    /* Symmetric part of Hermitian matrix is real part */
	    zeroIm(x);

	SET_SLOT(to, Matrix_DimNamesSym, dimnames);	
	SET_SLOT(to, Matrix_uploSym, uplo);
	
    }

    UNPROTECT(2);
    return to;
}

/* symmpart(x) */
SEXP matrix_symmpart(SEXP from)
{
    SEXP to,
	dim = getAttrib(from, R_DimSymbol),
	dimnames = getAttrib(from, R_DimNamesSymbol),
	x = from;
    int *pdim = INTEGER(dim), n = pdim[0], i, j, nprotect = 1;
    R_xlen_t upos = 0, lpos = 0, nn = (R_xlen_t) n * n;

    if (pdim[1] != n)	
	error(_("attempt to get symmetric part of non-square matrix"));
    
    switch (TYPEOF(x)) {
    case LGLSXP:
    case INTSXP:
	PROTECT(x = coerceVector(x, REALSXP));
	++nprotect;
    case REALSXP:
	PROTECT(to = NEW_OBJECT_OF_CLASS("dsyMatrix"));
	if (MAYBE_REFERENCED(x)) {
	    PROTECT(x = allocVector(REALSXP, nn));
	    Memcpy(REAL(x), REAL(from), nn);
	    ++nprotect;
	} else {
	    SET_ATTRIB(x, R_NilValue);
	}
	UPM_SYMMPART_GE(double, REAL, ASSIGN_OFFDIAG_DGE);
	break;
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
	PROTECT(to = NEW_OBJECT_OF_CLASS("zsyMatrix"));
	if (MAYBE_REFERENCED(x)) {
	    PROTECT(x = allocVector(CPLXSXP, nn));
	    Memcpy(COMPLEX(x), COMPLEX(from), nn);
	    ++nprotect;
	} else {
	    SET_ATTRIB(x, R_NilValue);
	}
	UPM_SYMMPART_GE(Rcomplex, COMPLEX, ASSIGN_OFFDIAG_ZGE);
	break;
#endif
    default:
	ERROR_INVALID_TYPE("matrix", TYPEOF(x), "matrix_symmpart");
	break;
    }
    
    SET_SLOT(to, Matrix_DimSym, dim);
    if (!isNull(dimnames))
	set_symmetrized_DimNames(to, dimnames, -1);
    SET_SLOT(to, Matrix_xSym, x);

    UNPROTECT(nprotect);
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
	"dgeMatrix", "dtrMatrix", "dsyMatrix",
	"lgeMatrix", "ltrMatrix", "lsyMatrix",
	"ngeMatrix", "ntrMatrix", "nsyMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "unpackedMatrix_skewpart");
    const char *clf = valid[ivalid];

    SEXP to,
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = (clf[1] != 'g') ? GET_SLOT(from, Matrix_uploSym) : R_NilValue,
	x = GET_SLOT(from, Matrix_xSym);
    
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)	
	error(_("attempt to get skew-symmetric part of non-square matrix"));

    if (clf[1] != 's') {

	SEXP y;
	char ul = (clf[1] != 'g') ? *CHAR(STRING_ELT(uplo, 0)) : 'U';
	R_xlen_t upos = 0, lpos = 0;
	int i, j;

#define UPM_SKEWPART(_CTYPE_, _PTR_, _X_, _Y_,				\
		     _ASSIGN_OFFDIAG_, _ASSIGN_ONDIAG_)			\
	do {								\
	    _CTYPE_ *px = _PTR_(_X_), *py = _PTR_(_Y_);			\
	    if (ul == 'U') {						\
		for (j = 0; j < n; ++j) {				\
		    lpos = j;						\
		    for (i = 0; i < j; ++i) {				\
			_ASSIGN_OFFDIAG_(upos, lpos);			\
			++upos; lpos += n;				\
		    }							\
		    _ASSIGN_ONDIAG_(upos);				\
		    upos += n-j;					\
		}							\
	    } else {							\
		for (j = 0; j < n; ++j) {				\
		    upos = lpos;					\
		    _ASSIGN_ONDIAG_(lpos);				\
		    for (i = j+1; i < n; ++i) {				\
			upos += n; ++lpos;				\
			_ASSIGN_OFFDIAG_(lpos, upos);			\
		    }							\
		    lpos += j+2;					\
		}							\
	    }								\
	} while (0)

#define ASSIGN_ONDIAG_DGE(_DPOS_) py[_DPOS_] = 0.0
	
#define ASSIGN_OFFDIAG_DGE(_UPOS_, _LPOS_)				\
	    do {							\
		py[_UPOS_] = 0.5 * (px[_UPOS_] - px[_LPOS_]);		\
		py[_LPOS_] = -py[_UPOS_];				\
	    } while (0)

#define ASSIGN_OFFDIAG_DTR(_UPOS_, _LPOS_)			\
	    do {						\
		py[_UPOS_] = 0.5 * px[_UPOS_];			\
		py[_LPOS_] = -py[_UPOS_];			\
	    } while (0)

#define ASSIGN_ONDIAG_ZGE(_DPOS_) py[_DPOS_].r = py[_DPOS_].i = 0.0
	
#define ASSIGN_OFFDIAG_ZGE(_UPOS_, _LPOS_)				\
	    do {							\
		py[_UPOS_].r = 0.5 * (px[_UPOS_].r - px[_LPOS_].r);	\
		py[_UPOS_].i = 0.5 * (px[_UPOS_].i - px[_LPOS_].i);	\
		py[_LPOS_].r = -py[upos].r;				\
		py[_LPOS_].i = -py[upos].i;				\
	    } while (0)

#define ASSIGN_OFFDIAG_ZTR(_UPOS_, _LPOS_)				\
	    do {							\
		py[_UPOS_].r = 0.5 * px[_UPOS_].r;			\
		py[_UPOS_].i = 0.5 * px[_UPOS_].i;			\
		py[_LPOS_].r = -py[upos].r;				\
		py[_LPOS_].i = -py[upos].i;				\
	    } while (0)
	
	if (clf[0] != 'z') {
	    PROTECT(to = NEW_OBJECT_OF_CLASS("dgeMatrix"));
	    if (clf[0] == 'd') {
		PROTECT(y = allocVector(REALSXP, (R_xlen_t) n * n));
		if (clf[1] == 'g')
		    UPM_SKEWPART(double, REAL, x, y,
				 ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
		else
		    UPM_SKEWPART(double, REAL, x, y,
				 ASSIGN_OFFDIAG_DTR, ASSIGN_ONDIAG_DGE);
	    } else {
		PROTECT(y = coerceVector(x, REALSXP));
		if (clf[1] == 'g')
		    UPM_SKEWPART(double, REAL, y, y,
				 ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
		else
		    UPM_SKEWPART(double, REAL, y, y,
				 ASSIGN_OFFDIAG_DTR, ASSIGN_ONDIAG_DGE);
	    }
	} else { /* clf[0] == 'z' */
	    PROTECT(to = NEW_OBJECT_OF_CLASS("zgeMatrix"));
	    PROTECT(y = allocVector(CPLXSXP, (R_xlen_t) n * n));
	    if (clf[1] == 'g')
		UPM_SKEWPART(Rcomplex, COMPLEX, x, y,
			     ASSIGN_OFFDIAG_ZGE, ASSIGN_ONDIAG_ZGE);
	    else
		UPM_SKEWPART(Rcomplex, COMPLEX, x, y,
			     ASSIGN_OFFDIAG_ZTR, ASSIGN_ONDIAG_ZGE);
	}

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
	    PROTECT(x = duplicate(GET_SLOT(from, Matrix_xSym)));
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

/* skewpart(x) */
SEXP matrix_skewpart(SEXP from)
{
    SEXP to,
	dim = getAttrib(from, R_DimSymbol),
	dimnames = getAttrib(from, R_DimNamesSymbol),
	x = from;
    int *pdim = INTEGER(dim), n = pdim[0], i, j, nprotect = 1;
    R_xlen_t upos = 0, lpos = 0;
    char ul = 'U';
    
    if (pdim[1] != n)	
	error(_("attempt to get skew-symmetric part of non-square matrix"));
    
    switch (TYPEOF(x)) {
    case LGLSXP:
    case INTSXP:
	PROTECT(x = coerceVector(x, REALSXP));
	++nprotect;
    case REALSXP:
	PROTECT(to = NEW_OBJECT_OF_CLASS("dgeMatrix"));
	if (MAYBE_REFERENCED(x)) {
	    PROTECT(x = allocVector(REALSXP, (R_xlen_t) n * n));
	    UPM_SKEWPART(double, REAL, from, x,
			 ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
	    ++nprotect;
	} else {
	    SET_ATTRIB(x, R_NilValue);
	    UPM_SKEWPART(double, REAL, x, x,
			 ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
	}
	break;
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
	PROTECT(to = NEW_OBJECT_OF_CLASS("zgeMatrix"));
	if (MAYBE_REFERENCED(from)) {
	    PROTECT(x = allocVector(CPLXSXP, (R_xlen_t) n * n));
	    UPM_SKEWPART(Rcomplex, COMPLEX, from, x,
			 ASSIGN_OFFDIAG_ZGE, ASSIGN_ONDIAG_ZGE);
	    ++nprotect;
	} else {
	    SET_ATTRIB(x, R_NilValue);
	    UPM_SKEWPART(Rcomplex, COMPLEX, x, x,
			 ASSIGN_OFFDIAG_DGE, ASSIGN_ONDIAG_DGE);
	}
	break;
#endif
    default:
	ERROR_INVALID_TYPE("matrix", TYPEOF(x), "matrix_skewpart");
	break;
    }
    
    SET_SLOT(to, Matrix_DimSym, dim);
    if (!isNull(dimnames))
	set_symmetrized_DimNames(to, dimnames, -1);
    SET_SLOT(to, Matrix_xSym, x);
    
    UNPROTECT(nprotect);
    return to;
}

#undef ASSIGN_ONDIAG_DGE
#undef ASSIGN_ONDIAG_ZGE
#undef ASSIGN_OFFDIAG_DGE
#undef ASSIGN_OFFDIAG_DTR
#undef ASSIGN_OFFDIAG_ZGE
#undef ASSIGN_OFFDIAG_ZTR
#undef UPM_SKEWPART
