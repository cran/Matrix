#include "dense.h"

SEXP matrix_as_dense(SEXP from, const char *code, char uplo, char diag,
		     int new, int transpose_if_vector)
{
    /* NB: also handing vectors 'from' _without_ a 'dim' attribute */
    SEXPTYPE tf = TYPEOF(from);
    switch (tf) {
    case LGLSXP:
#ifdef HAVE_PROPER_IMATRIX
    case INTSXP:
#endif
    case REALSXP:
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
#endif
	break;
#ifndef HAVE_PROPER_IMATRIX
    case INTSXP:
	if (!inherits(from, "factor"))
	    break;
#endif
    default:
	if (OBJECT(from))
	    ERROR_INVALID_CLASS(from, "matrix_as_dense");
	else
	    ERROR_INVALID_TYPE("object", tf, "matrix_as_dense");
	break;
    }

    char cl[] = "...Matrix";
    cl[0] = (code[0] == '.' || code[0] == ',') ? type2kind(tf) : code[0];
    cl[1] = code[1];
    cl[2] = code[2];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));
    
    SEXP dim, dimnames;
    int *pdim, m, n, doDN, isM = (int) isMatrix(from), nprotect = 1;
    R_xlen_t len = XLENGTH(from);
    
    if (isM) {
	
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(dim = getAttrib(from, R_DimSymbol), &pid);
	pdim = INTEGER(dim);
	if ((m = pdim[0]) != (n = pdim[1]) || n > 0) {
	    if (new > 1)
		REPROTECT(dim = duplicate(dim), pid);
	    SET_SLOT(to, Matrix_DimSym, dim);
	}
	UNPROTECT(1); /* dim */
	
	PROTECT(dimnames = getAttrib(from, R_DimNamesSymbol));
	++nprotect;
	doDN = !isNull(dimnames);
	
    } else {
	
	if (len > INT_MAX)
	    error(_("dimensions cannot exceed 2^31-1"));
	PROTECT(dim = allocVector(INTSXP, 2));
	pdim = INTEGER(dim);
	if (transpose_if_vector) {
	    pdim[0] = m = 1;
	    pdim[1] = n = (int) len;
	} else {
	    pdim[0] = m = (int) len;
	    pdim[1] = n = 1;
	}
	SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP nms = PROTECT(getAttrib(from, R_NamesSymbol));
	++nprotect;
	doDN = !isNull(nms);
	if (doDN) {
	    PROTECT(dimnames = allocVector(VECSXP, 2));
	    ++nprotect;
	    SET_VECTOR_ELT(dimnames, transpose_if_vector ? 1 : 0, nms);
	}
	
    }
    
    if (cl[1] != 'g' && m != n)
	error(_("attempt to construct triangular or symmetric "
		"denseMatrix from non-square matrix"));

    if (doDN) {
	if (cl[1] == 's')
	    set_symmetrized_DimNames(to, dimnames, -1);
	else if (isM && new > 1)
	    set_DimNames(to, dimnames);
	else
	    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    }
    
    if (cl[1] != 'g' && uplo != 'U') {
	SEXP val = PROTECT(mkString("L"));
	SET_SLOT(to, Matrix_uploSym, val);
	UNPROTECT(1); /* val */
    }

    if (cl[1] == 't' && diag != 'N') {
	SEXP val = PROTECT(mkString("U"));
	SET_SLOT(to, Matrix_diagSym, val);
	UNPROTECT(1); /* val */
    }
    
    SEXPTYPE tt = kind2type(cl[0]);
    if (tf != tt) {
	PROTECT(from = coerceVector(from, tt));
	++nprotect;
    }

    SEXP x;
    
    if (cl[2] != 'p') {
	
	if (tf != tt || new < 0 || (new == 0 && !MAYBE_REFERENCED(from)))
	    x = from;
	else {
	    PROTECT(x = allocVector(tt, len));
	    ++nprotect;
	    switch (tt) {
	    case LGLSXP:
		Matrix_memcpy(LOGICAL(x), LOGICAL(from), len, sizeof(int));
		break;
	    case INTSXP:
		Matrix_memcpy(INTEGER(x), INTEGER(from), len, sizeof(int));
		break;
	    case REALSXP:
		Matrix_memcpy(REAL(x), REAL(from), len, sizeof(double));
		break;
	    case CPLXSXP:
		Matrix_memcpy(COMPLEX(x), COMPLEX(from), len, sizeof(Rcomplex));
		break;
	    default:
		break;
	    }
	}
	if (!isNull(ATTRIB(x))) {
	    SET_ATTRIB(x, R_NilValue);
	    if (OBJECT(x))
		SET_OBJECT(x, 0);
	}
	
    } else {

	PROTECT(x = allocVector(tt, PM_LENGTH(n)));
	++nprotect;
	
#define PACK(_PREFIX_, _PTR_)						\
	_PREFIX_ ## dense_pack(_PTR_(x), _PTR_(from), n, uplo, diag)
	
	switch (tt) {
	case LGLSXP:
	    PACK(i, LOGICAL);
	    break;
#ifdef HAVE_PROPER_IMATRIX
	case INTSXP:
	    PACK(i, INTEGER);
	    break;
#endif
	case REALSXP:
	    PACK(d, REAL);
	    break;
#ifdef HAVE_PROPER_ZMATRIX
	case CPLXSXP:
	    PACK(z, COMPLEX);
	    break;
#endif
	default:
	    break;
	}

#undef PACK
	
    }
    
    SET_SLOT(to, Matrix_xSym, x);

    UNPROTECT(nprotect);
    return to;
}

/* as(<matrix>, ".(ge|tr|sy|tp|sp)Matrix") */
SEXP R_matrix_as_dense(SEXP from, SEXP code, SEXP uplo, SEXP diag)
{
    const char *zzz;
    if (TYPEOF(code) != STRSXP || LENGTH(code) < 1 ||
	(code = STRING_ELT(code, 0)) == NA_STRING ||
	(zzz = CHAR(code))[0] == '\0' ||
	(zzz             )[1] == '\0' ||
	!((zzz[1] == 'g' && (zzz[2] == 'e'                 )) ||
	  (zzz[1] == 't' && (zzz[2] == 'r' || zzz[2] == 'p')) ||
	  (zzz[1] == 's' && (zzz[2] == 'y' || zzz[2] == 'p'))))
	error(_("invalid 'code' to 'R_matrix_as_dense()'"));
    char ul = 'U', di = 'N';
    if (zzz[1] != 'g') {
	if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
	    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
	    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
	    error(_("invalid 'uplo' to 'R_matrix_as_dense()'"));
	if (zzz[1] == 't') {
	    if (TYPEOF(diag) != STRSXP || LENGTH(diag) < 1 ||
		(diag = STRING_ELT(diag, 0)) == NA_STRING ||
		((di = *CHAR(diag)) != 'N' && di != 'U'))
		error(_("invalid 'diag' to 'R_matrix_as_dense()'"));
	}
    }
    return matrix_as_dense(from, zzz, ul, di, 0, 0);
}

/* as(<denseMatrix>, "[CRT]sparseMatrix") */
SEXP R_dense_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP diag)
{
    static const char *valid[] = {
	VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
    int ivalid = R_check_class_etc(from, valid), nprotect = 0;
    
    const char *zzz;
    char z0, z1, z2;
    if (TYPEOF(code) != STRSXP || LENGTH(code) < 1 ||
	(code = STRING_ELT(code, 0)) == NA_STRING ||
	(z0 = (zzz = CHAR(code))[0]) == '\0' ||
	((z1 = zzz[1]) != '.' && z1 != 'g' && z1 != 't' && z1 != 's') ||
	((z2 = zzz[2]) != 'C' && z2 != 'R' && z2 != 'T'))
	error(_("invalid 'code' to 'R_dense_as_sparse()'"));
    
    SEXP dim, dimnames, x_from;
    SEXPTYPE txf, txt = (z0 == '.') ? NILSXP : kind2type(z0);
    PROTECT_INDEX pidA, pidB;
    int *pdim = NULL, doDN = 1, packed = 0;
    char clt[] = "...Matrix", ul = 'U', di = 'N';
    clt[2] = z2;
    
    if (ivalid >= 0) {
	
	const char *clf = valid[ivalid];
	packed = (clf[2] == 'p');
	
	PROTECT(dim = GET_SLOT(from, Matrix_DimSym));
	++nprotect;
	pdim = INTEGER(dim);
	
	PROTECT(dimnames = GET_SLOT(from, Matrix_DimNamesSym));
	++nprotect;

	PROTECT(x_from = GET_SLOT(from, Matrix_xSym));
	++nprotect;
	txf = TYPEOF(x_from);
	
	if (clf[1] != 'g') {
	    PROTECT(uplo = GET_SLOT(from, Matrix_uploSym));
	    ++nprotect;
	    ul = *CHAR(STRING_ELT(uplo, 0));

	    if (clf[1] == 't') {
		PROTECT(diag = GET_SLOT(from, Matrix_diagSym));
		++nprotect;
	    	di = *CHAR(STRING_ELT(diag, 0));
	    }
	}
	
	clt[0] = (z0 != '.') ? z0 : clf[0]; 
	clt[1] = clf[1];

    } else {

	/* 'from' is a base matrix or base vector, but we behave as though
	   it is the 'x' slot of a .(ge|tr|sy)Matrix (depending on 'code')
	   for efficiency, relying on the user to specify 'uplo' and 'diag'
	   as necessary ...
	*/
	
	if (isMatrix(from)) {

	    PROTECT(dim = getAttrib(from, R_DimSymbol));
	    ++nprotect;
	    pdim = INTEGER(dim);
	    
	    PROTECT(dimnames = getAttrib(from, R_DimNamesSymbol));
	    ++nprotect;
	    doDN = !isNull(dimnames);
	    
	} else {

	    R_xlen_t len = XLENGTH(from);
	    if (len > INT_MAX)
		error(_("vector of length exceeding 2^31-1 "
			"to 'R_dense_as_sparse()'"));

	    PROTECT(dim = allocVector(INTSXP, 2));
	    ++nprotect;
	    pdim = INTEGER(dim);
	    pdim[0] = (int) len;
	    pdim[1] = 1;

	    SEXP nms = PROTECT(getAttrib(from, R_NamesSymbol));
	    ++nprotect;
	    doDN = !isNull(nms);
	    if (doDN) {
		PROTECT(dimnames = allocVector(VECSXP, 2));
		++nprotect;
		SET_VECTOR_ELT(dimnames, 0, nms);
	    }

	}
	
	if (z1 == 't' || z1 == 's') {
	    
	    if (pdim[0] != pdim[1])
		error(_("attempt to construct triangular or symmetric "
			"%csparseMatrix from non-square matrix"), z2);

	    if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		(uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		(ul = *CHAR(uplo)) == '\0')
		error(_("invalid 'uplo' to 'R_dense_as_sparse()'"));
	    PROTECT(uplo = mkString((ul == 'U') ? "U" : "L"));
	    ++nprotect;

	    if (z1 == 't') {
		if (TYPEOF(diag) != STRSXP || LENGTH(diag) < 1 ||
		    (diag = STRING_ELT(diag, 0)) == NA_STRING ||
		    (di = *CHAR(diag)) == '\0')
		    error(_("invalid 'diag' to 'R_dense_as_sparse()'"));
		PROTECT(diag = mkString((di == 'N') ? "N" : "U"));
		++nprotect;
	    }

	}
	
	PROTECT_WITH_INDEX(x_from = from, &pidA);
	++nprotect;
	txf = TYPEOF(x_from);
	
#ifndef HAVE_PROPER_IMATRIX
	if (z0 == '.' && txf == INTSXP)
	    REPROTECT(x_from = coerceVector(x_from, txf = REALSXP), pidA);
#endif
	
	clt[0] = (z0 != '.') ? z0 : type2kind(txf);
	clt[1] = (z1 != '.') ? z1 : 'g';

    }
    
    if (z0 == '.')
	txt = txf;

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)), p_to, i_to, j_to, x_to;
    ++nprotect;
    p_to = i_to = j_to = x_to = NULL;
    
    R_xlen_t nnz = 0;
    int m = pdim[0], n = pdim[1], i, j, *pp, *pi, *pj;
    pp = pi = pj = NULL;
    
    if (m != n || n > 0)
	SET_SLOT(to, Matrix_DimSym, dim);
    if (doDN)
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (clt[1] != 'g' && ul != 'U')
	SET_SLOT(to, Matrix_uploSym, uplo);
    if (clt[1] == 't' && di != 'N')
	SET_SLOT(to, Matrix_diagSym, diag);
    if (clt[2] != 'T') {
	PROTECT(p_to = allocVector(
		    INTSXP, (R_xlen_t) ((clt[2] == 'C') ? n : m) + 1));
	++nprotect;
	SET_SLOT(to, Matrix_pSym, p_to);
	pp = INTEGER(p_to);
	*(pp++) = 0;
	if (n > 0 && di != 'N' && ul == ((clt[2] == 'C') ? 'U' : 'L'))
	    *(pp++) = 0; /* first row or column skipped in these loops */
    }
    
#define DAS_LOOP_GE2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	for (j = 0; j < n; ++j) {				\
	    for (i = 0; i < m; ++i, ++_X_)			\
		if (_NZ_(*_X_)) _DO_INNER_;			\
	    _DO_OUTER_;						\
	}							\
    } while (0)
    
#define DAS_LOOP_GE2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t mn1s = (R_xlen_t) m * n - 1;			\
	for (i = 0; i < m; ++i, _X_ -= mn1s) {			\
	    for (j = 0; j < n; ++j, _X_ += m)			\
		if (_NZ_(*_X_)) _DO_INNER_;			\
	    _DO_OUTER_;						\
	}							\
    } while (0)
    
#define DAS_LOOP_TRN2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	if (ul == 'U') {					\
	    for (j = 0; j < n; _X_ += n-(++j)) {		\
		for (i = 0; i <= j; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    for (j = 0; j < n; _X_ += (++j)) {			\
		for (i = j; i < n; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)
    
#define DAS_LOOP_TRN2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t d;						\
	if (ul == 'U') {					\
	    d = (R_xlen_t) n * n - 1;				\
	    for (i = 0; i < n; ++i, _X_ -= (d -= n)) {		\
		for (j = i; j < n; ++j, _X_ += n)		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    d = -1;						\
	    for (i = 0; i < n; ++i, _X_ -= (d += n)) {		\
		for (j = 0; j <= i; ++j, _X_ += n)		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)

#define DAS_LOOP_TRU2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	if (ul == 'U') {					\
	    _X_ += n;						\
	    for (j = 1; j < n; ++j) {				\
		for (i = 0; i < j; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
		_X_ += n-j;					\
	    }							\
	} else {						\
	    for (j = 0; j < n; ++j) {				\
		_X_ += j+1;					\
		for (i = j+1; i < n; ++i, ++_X_)		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)
    
#define DAS_LOOP_TRU2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t d;						\
	if (ul == 'U') {					\
	    d = (R_xlen_t) n * (n - 1) - 1;			\
	    for (i = 0; i < n; ++i) {				\
		for (j = i+1; j < n; ++j) {			\
		    _X_ += n;					\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		}						\
		_DO_OUTER_;					\
		_X_ -= (d -= n);				\
	    }							\
	} else {						\
	    ++_X_;						\
	    d = -1;						\
	    for (i = 1; i < n; ++i) {				\
		for (j = 0; j < i; ++j) {			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		    _X_ += n;					\
		}						\
		_DO_OUTER_;					\
		_X_ -= (d += n);				\
	    }							\
	}							\
    } while (0)
    
#define DAS_LOOP_TPN2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	if (ul == 'U') {					\
	    for (j = 0; j < n; ++j) {				\
		for (i = 0; i <= j; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    for (j = 0; j < n; ++j) {				\
		for (i = j; i < n; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)
    
#define DAS_LOOP_TPN2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t d;						\
	if (ul == 'U') {					\
	    d = PM_LENGTH(n) - 1;				\
	    for (i = 0; i < n; _X_ -= (d -= (++i))) {		\
		for (j = i; j < n; _X_ += (++j))		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    d = -1;						\
	    for (i = 0; i < n; _X_ -= (d += n-(++i))) {		\
		for (j = 0; j <= i; _X_ += n-(++j))		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)

#define DAS_LOOP_TPU2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	if (ul == 'U') {					\
	    for (j = 1; j < n; ++j) {				\
		++_X_;						\
		for (i = 0; i < j; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    for (j = 0; j < n; ++j) {				\
		++_X_;						\
		for (i = j+1; i < n; ++i, ++_X_)		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)

#define DAS_LOOP_TPU2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t d;						\
	if (ul == 'U') {					\
	    d = PM_LENGTH(n-1) - 1;				\
	    for (i = 0; i < n; ++i) {				\
		for (j = i+1; j < n; ++j) {			\
		    _X_ += j;					\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		}						\
		_DO_OUTER_;					\
		_X_ -= (d -= i+1);				\
	    }							\
	} else {						\
	    ++_X_;						\
	    d = -1;						\
	    for (i = 1; i < n; ++i) {				\
		for (j = 0; j < i; ++j) {			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		    _X_ += n-j-1;				\
		}						\
		_DO_OUTER_;					\
		_X_ -= (d += n-i);				\
	    }							\
	}							\
    } while (0)
    
#define DAS_VALID2T						\
    if (nnz > INT_MAX)						\
	error(_("attempt to construct sparse matrix with "	\
		"more than 2^31-1 nonzero elements"))
    
#define DAS_VALID2CR						\
    do { DAS_VALID2T; else *(pp++) = (int) nnz; } while (0)
    
#define DAS_SUBSUBCASES(_X_, _NZ_, _LOOP_2CT_, _LOOP_2R_)	\
    do {							\
	switch (clt[2]) {					\
	case 'C':						\
	    _LOOP_2CT_(_X_, _NZ_, ++nnz, DAS_VALID2CR);		\
	    break;						\
	case 'R':						\
	    _LOOP_2R_ (_X_, _NZ_, ++nnz, DAS_VALID2CR);		\
	    break;						\
	case 'T':						\
	    _LOOP_2CT_(_X_, _NZ_, ++nnz, DAS_VALID2T);		\
	    break;						\
	default:						\
	    break;						\
	}							\
    } while (0)
    
#define DAS_SUBCASES(_CTYPE_, _PTR_, _NZ_)				\
    do {								\
	_CTYPE_ *px = _PTR_(x_from);					\
	if (clt[1] == 'g')						\
	    /* .geMatrix */						\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_GE2C,  DAS_LOOP_GE2R);	\
	else if (!packed && di == 'N')					\
	    /* .syMatrix, non-unit diagonal .trMatrix */		\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_TRN2C, DAS_LOOP_TRN2R);	\
	else if (!packed)						\
	    /* unit diagonal .trMatrix */				\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_TRU2C, DAS_LOOP_TRU2R);	\
	else if (di == 'N')						\
	    /* .spMatrix, non-unit diagonal .tpMatrix */		\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_TPN2C, DAS_LOOP_TPN2R);	\
	else								\
	    /* unit diagonal .tpMatrix */				\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_TPU2C, DAS_LOOP_TPU2R);	\
    } while (0)
    
#define DAS_CASES(_SEXPTYPE_)						\
    do {								\
	switch (_SEXPTYPE_) {						\
	case LGLSXP:							\
	    DAS_SUBCASES(int, LOGICAL, ISNZ_LOGICAL);			\
	    break;							\
	case INTSXP:							\
	    DAS_SUBCASES(int, INTEGER, ISNZ_INTEGER);			\
	    break;							\
	case REALSXP:							\
	    DAS_SUBCASES(double, REAL, ISNZ_REAL);			\
	    break;							\
	case CPLXSXP:							\
	    DAS_SUBCASES(Rcomplex, COMPLEX, ISNZ_COMPLEX);		\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", _SEXPTYPE_, "R_dense_as_sparse"); \
	    break;							\
	}								\
    } while (0)

    /* First we loop over the _nontrivial part_ of the denseMatrix 'from',
       by row ('R' case) or by column ('C' and 'T' cases), counting the
       nonzero elements and filling the 'p' slot of the result accordingly 
       ('C' and 'R' cases) ... */
    DAS_CASES(txf);

#undef DAS_SUBCASES
#undef DAS_SUBSUBCASES
#undef DAS_VALID2CR
#undef DAS_VALID2T

    /* Then we allocate ... */
    if (clt[2] != 'R') {
	PROTECT(i_to = allocVector(INTSXP, nnz));
	++nprotect;
	SET_SLOT(to, Matrix_iSym, i_to);
	pi = INTEGER(i_to);
    }
    if (clt[2] != 'C') {
	PROTECT(j_to = allocVector(INTSXP, nnz));
	++nprotect;
	SET_SLOT(to, Matrix_jSym, j_to);
	pj = INTEGER(j_to);
    }
    if (clt[0] != 'n') {
	PROTECT_WITH_INDEX(x_to = allocVector(txf, nnz), &pidB);
	++nprotect;
    }

#define DAS_SUBSUBCASES(_X_, _Y_, _NZ_, _LOOP_2CT_, _LOOP_2R_)		\
    do {								\
	switch (clt[2]) {						\
	case 'C':							\
	    if (clt[0] == 'n')						\
		_LOOP_2CT_(_X_, _NZ_,					\
			   *(pi++) = i, );				\
	    else							\
		_LOOP_2CT_(_X_, _NZ_,					\
			   do {						\
			       *(pi++) = i;				\
			       *(_Y_++) = *_X_;				\
			   } while (0), );				\
	    break;							\
	case 'R':							\
	    if (clt[0] == 'n')						\
		_LOOP_2R_ (_X_, _NZ_,					\
			  *(pj++) = j, );				\
	    else							\
		_LOOP_2R_ (_X_, _NZ_,					\
			   do {						\
			       *(pj++) = j;				\
			       *(_Y_++) = *_X_;				\
			   } while (0), );				\
	    break;							\
	case 'T':							\
	    if (clt[0] == 'n')						\
		_LOOP_2CT_(_X_, _NZ_,					\
			   do {						\
			       *(pi++) = i;				\
			       *(pj++) = j;				\
			   } while (0), );				\
	    else							\
		_LOOP_2CT_(_X_, _NZ_,					\
			   do {						\
			       *(pi++) = i;				\
			       *(pj++) = j;				\
			       *(_Y_++) = *_X_;				\
			   } while (0), );				\
	    break;							\
	default:							\
	    break;							\
	}								\
    } while (0)
    
#define DAS_SUBCASES(_CTYPE_, _PTR_, _NZ_)				\
    do {								\
	_CTYPE_ *px = _PTR_(x_from), *py = NULL;			\
	if (clt[0] != 'n')						\
	    py = _PTR_(x_to);						\
	if (clt[1] == 'g')						\
	    /* .geMatrix */						\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_GE2C,  DAS_LOOP_GE2R); \
	else if (!packed && di == 'N')					\
	    /* .syMatrix, non-unit diagonal .trMatrix */		\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_TRN2C, DAS_LOOP_TRN2R); \
	else if (!packed)						\
	    /* unit diagonal .trMatrix */				\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_TRU2C, DAS_LOOP_TRU2R); \
	else if (di == 'N')						\
	    /* .spMatrix, non-unit diagonal .tpMatrix */		\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_TPN2C, DAS_LOOP_TPN2R); \
	else								\
	    /* unit diagonal .tpMatrix */				\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_TPU2C, DAS_LOOP_TPU2R); \
    } while (0)
    
    /* Then we loop back over the same elements in order to fill
       the 'i', 'j', and 'x' slots of the result (whichever exist) ... */
    DAS_CASES(txf);

#undef DAS_CASES
#undef DAS_SUBCASES
#undef DAS_SUBSUBCASES
#undef DAS_LOOP_GE2C
#undef DAS_LOOP_GE2R
#undef DAS_LOOP_TRN2C
#undef DAS_LOOP_TRN2R
#undef DAS_LOOP_TRU2C
#undef DAS_LOOP_TRU2R
#undef DAS_LOOP_TPN2C
#undef DAS_LOOP_TPN2R
#undef DAS_LOOP_TPU2C
#undef DAS_LOOP_TPU2R

    if (clt[0] != 'n') {
	REPROTECT(x_to = coerceVector(x_to, txt), pidB);
	SET_SLOT(to, Matrix_xSym, x_to);
    }
    
    UNPROTECT(nprotect);
    return to;
}

/* as(<denseMatrix>, "matrix") */
SEXP R_dense_as_matrix(SEXP from)
{
    /* Result must be newly allocated because we add attributes */
    PROTECT(from = dense_as_general(from, ',', 1, 0));
    SEXP to = PROTECT(GET_SLOT(from, Matrix_xSym)),
	dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    setAttrib(to, R_DimSymbol, dim);
    if (!DimNames_is_trivial(dimnames))
	setAttrib(to, R_DimNamesSymbol, dimnames);
    UNPROTECT(4); /* dimnames, dim, to, from */
    return to;
}

/* as(<.geMatrix>, "matrix") */
SEXP R_geMatrix_as_matrix(SEXP from, SEXP pattern)
{
    /* Result must be newly allocated because we add attributes */
    SEXP to,
	dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(to = GET_SLOT(from, Matrix_xSym), &pid);
    REPROTECT(to = duplicate(to), pid);
    if (asLogical(pattern) != 0)
	na2one(to);
    setAttrib(to, R_DimSymbol, dim);
    if (!DimNames_is_trivial(dimnames))
	setAttrib(to, R_DimNamesSymbol, dimnames);
    UNPROTECT(3); /* dimnames, dim, to */
    return to;
}

/* as(<denseMatrix>, "vector") */
SEXP R_dense_as_vector(SEXP from)
{
    /* Result must be newly allocated if and only if different from 'x' slot */
    PROTECT(from = dense_as_general(from, ',', 0, 0));
    from = GET_SLOT(from, Matrix_xSym);
    UNPROTECT(1); /* from */
    return from;
}

/* as(<.geMatrix>, "vector") */
SEXP R_geMatrix_as_vector(SEXP from, SEXP pattern)
{
    /* Result must be newly allocated if and only if different from 'x' slot */
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(from = GET_SLOT(from, Matrix_xSym), &pid);
    if (asLogical(pattern) != 0) {
	int *px = LOGICAL(from);
	R_xlen_t nx = XLENGTH(from);
	while (nx--) {
	    if (*(px++) == NA_LOGICAL) {
		REPROTECT(from = duplicate(from), pid);
		na2one(from);
		break;
	    }
	}
    }
    UNPROTECT(1); /* from */
    return from;
}

/* as(<denseMatrix>, "[nlidz](dense)?Matrix") */
SEXP R_dense_as_kind(SEXP from, SEXP kind)
{
    static const char *valid[] = {
	VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(from, "R_dense_as_kind");
    const char *clf = valid[ivalid];
    
    char k;
    if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	(kind = STRING_ELT(kind, 0)) == NA_STRING ||
	(k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_dense_as_kind()'"));
    if (k == '.' || k == clf[0])
	return from;
    SEXPTYPE tt = kind2type(k); /* validating before doing more */
    
    char clt[] = "...Matrix";
    clt[0] = k;
    clt[1] = clf[1];
    clt[2] = clf[2];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

    SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    if (m != n || n > 0)
	SET_SLOT(to, Matrix_DimSym, dim);
    UNPROTECT(1); /* dim */

    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */
    
    if (clf[1] != 'g') {
	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	if (ul != 'U')
	    SET_SLOT(to, Matrix_uploSym, uplo);
	UNPROTECT(1); /* uplo */
	if (clf[1] == 't') {
	    SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	    char di = *CHAR(STRING_ELT(diag, 0));
	    if (di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	    UNPROTECT(1); /* diag */
	}
    }

    SEXP x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);
    SEXPTYPE tf = TYPEOF(x);
    
    if (clf[0] != 'n')
	REPROTECT(x = coerceVector(x, tt), pid);
    else {
	R_xlen_t ix, nx = XLENGTH(x);
	if (tf == tt) {
	    /* n->l ... allocate iff 'x' contains NA */
	    int *px = LOGICAL(x);
	    for (ix = 0; ix < nx; ++ix, ++px)
		if (*px == NA_LOGICAL)
		    break;
	    if (ix != nx) {
		REPROTECT(x = duplicate(x), pid);
		px = LOGICAL(x);
		for (ix = 0; ix < nx; ++ix, ++px)
		    if (*px == NA_LOGICAL)
			*px = 1;
	    }
	} else {
	    /* n->[idz] */
	    REPROTECT(x = coerceVector(x, tt), pid);
	    switch (tt) {
	    case INTSXP:
	    {
		int *px = INTEGER(x);
		for (ix = 0; ix < nx; ++ix, ++px)
		    if (*px == NA_INTEGER)
			*px = 1;
		break;
	    }
	    case REALSXP:
	    {
		double *px = REAL(x);
		for (ix = 0; ix < nx; ++ix, ++px)
		    if (ISNAN(*px))
			*px = 1.0;
		break;
	    }
	    case CPLXSXP:
	    {
		Rcomplex *px = COMPLEX(x);
		for (ix = 0; ix < nx; ++ix, ++px) {
		    if (ISNAN((*px).r) || ISNAN((*px).i)) {
			(*px).r = 1.0;
			(*px).i = 0.0;
		    }
		}
		break;
	    }
	    default:
		break;
	    }
	}
    }

    SET_SLOT(to, Matrix_xSym, x);
    
    UNPROTECT(2); /* x, to */
    return to;
}

/** @brief Coerce `denseMatrix` (and others) to `.geMatrix`.
 *
 *  This utility supports the many `*_{prod,crossprod,tcrossprod,...}`
 *  functions that should work with both classed and unclassed matrices.
 *  It is used in many places for `.geMatrix` ("generalized") dispatch.
 *
 * @param from A `denseMatrix`, a `diagonalMatrix`, a numeric or logical 
 *     `matrix`, or a numeric or logical vector.
 * @param kind A `char` flag, one of `'.'`, `','`, `'d'`, `'l'`, and `'n'`, 
 *     indicating the "kind" of `.geMatrix` desired.  A dot `'.'` means 
 *     to preserve the kind of `from`.  A comma `','` is equivalent to
 *     `'.'` with the additional step of replacing `NA` with 1 for `from`
 *     inheriting from `nMatrix`.
 * @param new An `int` flag allowing the user to control allocation.
 *     If 0, then usual copy-on-modify rules are employed.  
 *     If less than 0, then the `x` slot of the result is the result 
 *     of modifying in place the `x` slot of `from` (or `from` itself 
 *     if it is a `matrix` or vector), unless the two have different 
 *     lengths or types.
 *     If greater than 0, then the `x` slot of the result is always 
 *     newly allocated.
 *     If greater than 1, then the `Dim` and `Dimnames` slots are also 
 *     always newly allocated (but never the _elements_ of `Dimnames`).
 * @param transpose_if_vector Should length-`n` vectors without a `dim` 
 *     attribute become 1-by-`n` (rather than `n`-by-1) matrices?
 *
 * @return A `.geMatrix`.
 */
SEXP dense_as_general(SEXP from, char kind, int new, int transpose_if_vector)
{
    /* NB: diagonalMatrix is no longer a subclass of denseMatrix 
           but .diMatrix are nonetheless dealt with here ... 
    */
    
    static const char *valid[] = {
	"dgeMatrix", "lgeMatrix", "ngeMatrix",
	"dtrMatrix", "dsyMatrix", "dtpMatrix", "dspMatrix", "ddiMatrix",
	"ltrMatrix", "lsyMatrix", "ltpMatrix", "lspMatrix", "ldiMatrix",
	"ntrMatrix", "nsyMatrix", "ntpMatrix", "nspMatrix", /* "ndiMatrix", */
	""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0) {
	/* dispatch to method for base matrix and base vector */
	char zzz[] = ".ge";
	zzz[0] = kind;
	return matrix_as_dense(from, zzz, '\0', '\0', new, transpose_if_vector);
    }
    
    /* Now handling just denseMatrix and diagonalMatrix ...
       We want to be fast if 'from' is already a .geMatrix,
       and _especially_ fast if it is already a .geMatrix
       of the right "kind" ... 
    */
    
    const char *clf = valid[ivalid];
    int do_na2one = clf[0] == 'n' && kind != 'n' && kind != '.';
    if (kind == '.' || kind == ',')
	kind = clf[0];
    int ge = clf[1] == 'g', ge0 = ge && kind == clf[0];
    if (ge0 && new <= 0 && !do_na2one)
	return from;
    
    SEXP to;
    if (ge0)
	PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
    else {
	char clt[] = ".geMatrix";
	clt[0] = kind;
	PROTECT(to = NEW_OBJECT_OF_CLASS(clt));
    }

    SEXP dim;
    PROTECT_INDEX pidA;
    PROTECT_WITH_INDEX(dim = GET_SLOT(from, Matrix_DimSym), &pidA);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    if (m != n || n > 0) {
	if (new > 1)
	    REPROTECT(dim = duplicate(dim), pidA);
	SET_SLOT(to, Matrix_DimSym, dim);
    }
    UNPROTECT(1); /* dim */
    
    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    if (clf[1] == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else if (new > 1)
	set_DimNames(to, dimnames);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */
    
    SEXP x0;
    PROTECT_INDEX pidB;
    PROTECT_WITH_INDEX(x0 = GET_SLOT(from, Matrix_xSym), &pidB);
    
    if (ge0) {
	REPROTECT(x0 = duplicate(x0), pidB);
	if (do_na2one)
	    na2one(x0);
	SET_SLOT(to, Matrix_xSym, x0);
	UNPROTECT(2); /* x0, to */
	return to;
    }
    
    SEXPTYPE tf = TYPEOF(x0), tt = kind2type(kind);
    if (ge) {
	if (new == 0 && do_na2one) {
	    /* Try to avoid an allocation ... */
	    R_xlen_t ix, nx = XLENGTH(x0);
	    int *px0 = LOGICAL(x0);
	    for (ix = 0; ix < nx; ++ix)
		if (*(px0++) == NA_LOGICAL)
		    break;
	    do_na2one = (ix < nx);
	}
	if (tf != tt)
	    REPROTECT(x0 = coerceVector(x0, tt), pidB);
	else if (new > 0 || (new == 0 && do_na2one))
	    REPROTECT(x0 = duplicate(x0), pidB);
	if (do_na2one)
	    na2one(x0);
	SET_SLOT(to, Matrix_xSym, x0);
	UNPROTECT(2); /* x0, to */
	return to;
    }
    
    /* Now handling 'from' inheriting from .(tr|sy|tp|sp|di)Matrix ... */

    char ul = 'U', di = 'N';
    if (clf[1] != 'd') {
	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */
    }
    if (clf[1] != 's') {
	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
    } else if (kind == clf[0]) {
	SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
	if (LENGTH(factors) > 0)
	    SET_SLOT(to, Matrix_factorSym, factors);
	UNPROTECT(1); /* factors */
    }

    SEXP x1;
    int packed = clf[1] == 'd' || clf[2] == 'p';
    if (!packed) {
	/* (tr|sy)->ge */
	if (tf != tt)
	    REPROTECT(x1 = coerceVector(x0, tt), pidB);
	else if (new >= 0)
	    REPROTECT(x1 = duplicate(x0), pidB);
	else
	    x1 = x0;
    } else {
	/* (tp|sp|di)->ge */
	if ((double) n * n > R_XLEN_T_MAX)
	    error(_("attempt to allocate vector of length exceeding "
		    "R_XLEN_T_MAX"));
	if (tf != tt)
	    REPROTECT(x0 = coerceVector(x0, tt), pidB);
	PROTECT(x1 = allocVector(tt, (R_xlen_t) n * n));
    }

#define AS_GE(_PREFIX_, _CTYPE_, _PTR_)					\
    do {								\
	_CTYPE_ *px1 = _PTR_(x1);					\
	if (clf[1] == 'd') {						\
	    /* di->ge */						\
	    Matrix_memset(px1, 0, (R_xlen_t) n * n, sizeof(_CTYPE_));	\
	    _PREFIX_ ## dense_unpacked_copy_diagonal(			\
		px1, _PTR_(x0), n, n, ul /* unused */, di);		\
	} else {							\
	    /* (tr|sy|tp|sp)->ge */					\
	    if (clf[2] == 'p')						\
		_PREFIX_ ## dense_unpack(px1, _PTR_(x0), n, ul, di);	\
	    if (clf[1] == 't')						\
		_PREFIX_ ## dense_unpacked_make_triangular(px1, n, n, ul, di); \
	    else							\
		_PREFIX_ ## dense_unpacked_make_symmetric(px1, n, ul);	\
	}								\
    } while (0)

    switch (kind) {
    case 'n':
    case 'l':
	AS_GE(i, int, LOGICAL);
	break;
    case 'i':
	AS_GE(i, int, INTEGER);
	break;
    case 'd': 
	AS_GE(d, double, REAL);
	break;
    case 'z':
	AS_GE(z, Rcomplex, COMPLEX);
	break;
    default:
	break;
    }

#undef AS_GE

    if (do_na2one)
	na2one(x1);
    SET_SLOT(to, Matrix_xSym, x1);

    if (packed)
	UNPROTECT(3); /* x1, x0, to */
    else
	UNPROTECT(2); /* x1,     to */
    return to;
}

/* as(<denseMatrix>, "generalMatrix") */
SEXP R_dense_as_general(SEXP from, SEXP kind)
{
    char k;
    if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	(kind = STRING_ELT(kind, 0)) == NA_STRING ||
	(k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_dense_as_general()'"));
    return dense_as_general(from, k, 0, 0);
}

/* band(<denseMatrix>, k1, k2), tri[ul](<denseMatrix>, k)
   band(     <matrix>, k1, k2), tri[ul](     <matrix>, k) */
SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2)
{
    static const char *valid[] = {
	VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
    int ivalid = R_check_class_etc(from, valid), nprotect = 0;
    const char *clf;
    if (ivalid >= 0)
	clf = valid[ivalid];
    else {
	/* matrix->.geMatrix with unreferenced 'x' slot ... modify directly! */
	PROTECT(from = matrix_as_dense(from, ".ge", '\0', '\0', 0, 0));
	++nprotect;
	clf = valid[R_check_class_etc(from, valid)];
    }
    
    SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
    ++nprotect;
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], a, b;
    if (isNull(k1))
	a = (m > 0) ? 1-m : 0;
    else if ((a = asInteger(k1)) == NA_INTEGER || a < -m || a > n)
	error(_("'k1' must be an integer from -Dim[1] to Dim[2]"));
    if (isNull(k2))
	b = (n > 0) ? n-1 : 0;
    else if ((b = asInteger(k2)) == NA_INTEGER || b < -m || b > n)
	error(_("'k2' must be an integer from -Dim[1] to Dim[2]"));
    else if (b < a)
	error(_("'k1' must be less than or equal to 'k2'"));
    /* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
    if (a <= 1-m && b >= n-1 && (clf[1] == 't' || m != n || m > 1 || n > 1)) {
	UNPROTECT(nprotect);
	return from;
    }
    
    char ulf = 'U', ult = 'U', di = 'N';
    if (clf[1] != 'g') {
	SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
	ulf = *CHAR(STRING_ELT(uplo_from, 0));
	UNPROTECT(1); /* uplo_from */
	if (clf[1] == 't') {
	    /* Be fast if band contains entire triangle */
	    if ((ulf == 'U') ? (a <= 0 && b >= n-1) : (b >= 0 && a <= 1-m)) {
		UNPROTECT(nprotect);
		return from;
	    } else if (a <= 0 && b >= 0) {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */
	    }
	}
    }
    
    SEXP x_from = NULL, x_to = NULL;

#define UNPACKED_MAKE_BANDED(_PREFIX_, _CTYPE_, _PTR_)			\
    _PREFIX_ ## dense_unpacked_make_banded(_PTR_(x_to), m, n, a, b, di)
    
#define PACKED_MAKE_BANDED(_PREFIX_, _CTYPE_, _PTR_)			\
    _PREFIX_ ## dense_packed_make_banded(_PTR_(x_to), n, a, b, ult, di)
    
#define DENSE_BAND(_MAKE_BANDED_)					\
    do {								\
	switch (TYPEOF(x_to)) {						\
	case LGLSXP: /* [nl]..Matrix */					\
	    _MAKE_BANDED_(i, int, LOGICAL);				\
	    break;							\
	case INTSXP: /* i..Matrix */					\
	    _MAKE_BANDED_(i, int, INTEGER);				\
	    break;							\
	case REALSXP: /* d..Matrix */					\
	    _MAKE_BANDED_(d, double, REAL);				\
	    break;							\
	case CPLXSXP: /* z..Matrix */					\
	    _MAKE_BANDED_(z, Rcomplex, COMPLEX);			\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", TYPEOF(x_to), "R_dense_band"); \
	    break;							\
	}								\
    } while (0)

#define IF_UNPACKED_ELSE(_IF_, _ELSE_)		\
    do {					\
	if (clf[2] != 'p')			\
	    DENSE_BAND(_IF_);			\
	else					\
	    DENSE_BAND(_ELSE_);			\
    } while (0)

    int ge = 0, tr = 0, sy = 0;
    ge = m != n || !((tr = a >= 0 || b <= 0 || clf[1] == 't') ||
		     (sy = a == -b && clf[1] == 's'));
    
    if (ge) {
	/* Returning .geMatrix */
	if (ivalid >= 0) {
	    PROTECT(from = dense_as_general(from, '.', 1, 0));
	    ++nprotect;
	}
	PROTECT(x_to = GET_SLOT(from, Matrix_xSym));
	++nprotect;
	DENSE_BAND(UNPACKED_MAKE_BANDED);
	UNPROTECT(nprotect);
	return from;
    }
    
    /* Returning .(tr|sy|tp|sp)Matrix ... */

    SEXP to;
    if (tr) {
	char clt[] = ".t.Matrix";
	clt[0] = clf[0];
	clt[2] = (clf[2] != 'p') ? 'r' : 'p';
	PROTECT(to = NEW_OBJECT_OF_CLASS(clt));
	++nprotect;
    } else {
	PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
	++nprotect;
    }

    if (m != n || n > 0)
	SET_SLOT(to, Matrix_DimSym, dim);

    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    if (tr && clf[1] == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */

    ult = (tr && clf[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ulf;
    if (ult != 'U') {
	SEXP uplo_to = PROTECT(mkString("L"));
	SET_SLOT(to, Matrix_uploSym, uplo_to);
	UNPROTECT(1); /* uplo_to */
    }
    
    if (di != 'N') {
	SEXP diag = PROTECT(mkString("U"));
	SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */
    }
    
    /* It remains to set 'x' ... */

    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x_from = GET_SLOT(from, Matrix_xSym), &pid);
    ++nprotect;
    
    if (tr) {
	/* Returning .t[rp]Matrix */
	if (clf[1] == 't') {
	    if (ulf == ult) {
		REPROTECT(x_to = duplicate(x_from), pid);
		IF_UNPACKED_ELSE(UNPACKED_MAKE_BANDED, PACKED_MAKE_BANDED);
	    } else {
		/* Result is either a diagonal matrix or a zero matrix */
		R_xlen_t nx = XLENGTH(x_from);
		PROTECT(x_to = allocVector(TYPEOF(x_from), nx));
		++nprotect;
		
#define UNPACKED_COPY_DIAGONAL(_PREFIX_, _CTYPE_, _PTR_)		\
		do {							\
		    Matrix_memset(_PTR_(x_to), 0, nx, sizeof(_CTYPE_));	\
		    if (a <= 0 && b >= 0)				\
			_PREFIX_ ## dense_unpacked_copy_diagonal(	\
			    _PTR_(x_to), _PTR_(x_from),			\
			    n, nx, 'U' /* unused */, di);		\
		} while (0)

#define PACKED_COPY_DIAGONAL(_PREFIX_, _CTYPE_, _PTR_)			\
		do {							\
		    Matrix_memset(_PTR_(x_to), 0, nx, sizeof(_CTYPE_));	\
		    if (a <= 0 && b >= 0)				\
			_PREFIX_ ## dense_packed_copy_diagonal(		\
			    _PTR_(x_to), _PTR_(x_from),			\
			    n, nx, ult, ulf, di);			\
		} while (0)

		IF_UNPACKED_ELSE(UNPACKED_COPY_DIAGONAL, PACKED_COPY_DIAGONAL);
	    }
	} else { /* clf[1] == 'g' || clf[1] == 's' */
	    if (ivalid >= 0) {
		x_to = (clf[1] == 'g' || ulf == ult || n <= 1
			? duplicate(x_from)
			/* band is "opposite" the stored triangle: */
			: (clf[2] == 'p'
			   ? packed_transpose(x_from, n, ulf)
			   : unpacked_force(x_from, n, ulf, '\0')));
		REPROTECT(x_to, pid);
	    } else
		x_to = x_from;
	    IF_UNPACKED_ELSE(UNPACKED_MAKE_BANDED, PACKED_MAKE_BANDED);
	}
    } else {
	/* Returning .s[yp]Matrix */
	REPROTECT(x_to = duplicate(x_from), pid);
	IF_UNPACKED_ELSE(UNPACKED_MAKE_BANDED, PACKED_MAKE_BANDED);
    }

#undef IF_UNPACKED_ELSE
#undef UNPACKED_MAKE_BANDED
#undef UNPACKED_COPY_DIAGONAL
#undef PACKED_MAKE_BANDED
#undef PACKED_COPY_DIAGONAL
#undef DENSE_BAND

    SET_SLOT(to, Matrix_xSym, x_to);
    
    UNPROTECT(nprotect);
    return to;
}

/* colSums(<denseMatrix>) */
SEXP R_dense_colSums(SEXP obj, SEXP narm, SEXP mean)
{
    static const char *valid[] = {
	VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(obj, "R_dense_colSums");
    const char *cl = valid[ivalid];
    if (cl[1] == 's')
	return R_dense_rowSums(obj, narm, mean);

    int doNaRm = asLogical(narm) != 0,
	doMean = asLogical(mean) != 0,
	doCount = doNaRm && doMean;
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1); /* dim */
    
    char ul = 'U', di = 'N';
    if (cl[1] == 't') {
	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */
	
	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
    }

    SEXP res = PROTECT(allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, n)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int i, j, count = m;
    
#define DENSE_COLSUMS_LOOP						\
    do {								\
	if (cl[1] == 'g') { /* general */				\
	    for (j = 0; j < n; ++j, ++pres) {				\
		DO_INIT(ZERO);						\
		for (i = 0; i < m; ++i, ++px)				\
		    DO_INCR;						\
		DO_SCALE;						\
	    }								\
	} else if (cl[2] != 'p') {					\
	    if (ul == 'U') { /* unpacked upper triangular */		\
		if (di == 'N') {					\
		    for (j = 0; j < n; ++j, ++pres) {			\
			DO_INIT(ZERO);					\
			for (i = 0; i <= j; ++i, ++px)			\
			    DO_INCR;					\
			DO_SCALE;					\
			px += n-j-1;					\
		    }							\
		} else {						\
		    for (j = 0; j < n; ++j, ++pres) {			\
			DO_INIT(ONE);					\
			for (i = 0; i < j; ++i, ++px)			\
			    DO_INCR;					\
			DO_SCALE;					\
			px += n-j;					\
		    }							\
		}							\
	    } else { /* unpacked lower triangular */			\
		if (di == 'N') {					\
		    for (j = 0; j < n; ++j, ++pres) {			\
			px += j;					\
			DO_INIT(ZERO);					\
			for (i = j; i < n; ++i, ++px)			\
			    DO_INCR;					\
			DO_SCALE;					\
		    }							\
		} else {						\
		    for (i = j = 0; j < n; ++j, i = j, ++pres) {	\
			px += j+1;					\
			DO_INIT(ONE);					\
			for (i = j+1; i < n; ++i, ++px)			\
			    DO_INCR;					\
			DO_SCALE;					\
		    }							\
		}							\
	    }								\
	} else {							\
	    if (ul == 'U') { /* packed upper triangular */		\
		if (di == 'N') {					\
		    for (j = 0; j < n; ++j, ++pres) {			\
			DO_INIT(ZERO);					\
			for (i = 0; i <= j; ++i, ++px)			\
			    DO_INCR;					\
			DO_SCALE;					\
		    }							\
		} else {						\
		    for (j = 0; j < n; ++j, ++pres) {			\
			DO_INIT(ONE);					\
			for (i = 0; i < j; ++i, ++px)			\
			    DO_INCR;					\
			DO_SCALE;					\
			++px;						\
		    }							\
		}							\
	    } else { /* packed lower triangular */			\
		if (di == 'N') {					\
		    for (j = 0; j < n; ++j, ++pres) {			\
			DO_INIT(ZERO);					\
			for (i = j; i < n; ++i, ++px)			\
			    DO_INCR;					\
			DO_SCALE;					\
		    }							\
		} else {						\
		    for (i = j = 0; j < n; ++j, i = j, ++pres) {	\
			++px;						\
			DO_INIT(ONE);					\
			for (i = j+1; i < n; ++i, ++px)			\
			    DO_INCR;					\
			DO_SCALE;					\
		    }							\
		}							\
	    }								\
	}								\
    } while (0)

#define DENSE_COLSUMS(_CTYPE1_, _PTR1_, _CTYPE2_, _PTR2_)		\
    do {								\
	_CTYPE1_ *pres = _PTR1_(res);					\
	_CTYPE2_ *px   = _PTR2_(x);					\
	DENSE_COLSUMS_LOOP;						\
    } while (0)
    
    switch (cl[0]) {
    case 'n':

#define ZERO         0.0
#define ONE          1.0
#define DO_INIT(_U_) *pres = _U_
#define DO_INCR	     if (*px) *pres += 1.0
#define DO_SCALE     if (doMean) *pres /= count

	DENSE_COLSUMS(double, REAL, int, LOGICAL);
	break;

#undef DO_INIT
#undef DO_INCR
	
    case 'l':
	
#define DO_INIT(_U_)				\
	do {					\
	    *pres = _U_;			\
	    if (doCount)			\
		count = m;			\
	} while (0)
#define DO_INCR					\
	do {					\
	    if (*px != NA_LOGICAL) {		\
		if (*px) *pres += 1.0;		\
	    } else if (!doNaRm)			\
		*pres = NA_REAL;		\
	    else if (doMean)			\
		--count;			\
	} while (0)
	
	DENSE_COLSUMS(double, REAL, int, LOGICAL);
	break;
	
#undef DO_INCR
	
    case 'i':
	    
#define DO_INCR					\
	do {					\
	    if (*px != NA_INTEGER)		\
		*pres += *px;			\
	    else if (!doNaRm)			\
		*pres = NA_REAL;		\
	    else if (doMean)			\
		--count;			\
	} while (0)
	
	DENSE_COLSUMS(double, REAL, int, INTEGER);
	break;

#undef DO_INCR

    case 'd':

#define DO_INCR					\
	do {					\
	    if (!(doNaRm && ISNAN(*px)))	\
		*pres += *px;			\
	    else if (doMean)			\
		--count;			\
	} while (0)

	DENSE_COLSUMS(double, REAL, double, REAL);
	break;

#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_SCALE
	
    case 'z':

#define ZERO         Matrix_zzero
#define ONE          Matrix_zone
#define DO_INCR								\
	do {								\
	    if (!(doNaRm && (ISNAN((*px).r) || ISNAN((*px).i)))) {	\
		(*pres).r += (*px).r;					\
		(*pres).i += (*px).i;					\
	    } else if (doMean)						\
		--count;						\
	} while (0)
#define DO_SCALE				\
	do {					\
	    if (doMean) {			\
		(*pres).r /= count;		\
		(*pres).i /= count;		\
	    }					\
	} while (0)

	DENSE_COLSUMS(Rcomplex, COMPLEX, Rcomplex, COMPLEX);
	break;

#undef ZERO
#undef ONE
#undef DO_INIT
#undef DO_INCR
#undef DO_SCALE

    default:
	break;
    }

#undef DENSE_COLSUMS
#undef DENSE_COLSUMS_LOOP
    
    SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	nms = VECTOR_ELT(dimnames, 1);
    if (!isNull(nms))
	setAttrib(res, R_NamesSymbol, nms);
    
    UNPROTECT(3); /* dimnames, x, res */
    return res;
}

/* rowSums(<denseMatrix>) */
SEXP R_dense_rowSums(SEXP obj, SEXP narm, SEXP mean)
{
    static const char *valid[] = {
	VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(obj, "R_dense_rowSums");
    const char *cl = valid[ivalid];

    int doNaRm = asLogical(narm) != 0,
	doMean = asLogical(mean) != 0;
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1); /* dim */
    
    char ul = 'U', di = 'N';
    if (cl[1] != 'g') {
	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */

	if (cl[1] == 't') {
	    SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	    di = *CHAR(STRING_ELT(diag, 0));
	    UNPROTECT(1); /* diag */
	}
    }

    SEXP res = PROTECT(allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, m)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int i, j, *pcount = NULL;
    
#define DENSE_ROWSUMS_LOOP						\
    do {								\
	if (cl[1] == 'g') { /* general */				\
	    for (j = 0; j < n; ++j)					\
		for (i = 0; i < m; ++i, ++px)				\
		    DO_INCR;						\
	} else if (cl[1] == 't') {					\
	    if (cl[2] != 'p') {						\
		if (ul == 'U') { /* unpacked upper triangular */	\
		    if (di == 'N') {					\
			for (j = 0; j < n; ++j) {			\
			    for (i = 0; i <= j; ++i, ++px)		\
				DO_INCR;				\
			    px += n-j-1;				\
			}						\
		    } else {						\
			for (j = 0; j < n; ++j) {			\
			    for (i = 0; i < j; ++i, ++px)		\
				DO_INCR;				\
			    px += n-j;					\
			}						\
		    }							\
		} else { /* unpacked lower triangular */		\
		    if (di == 'N') {					\
			for (j = 0; j < n; ++j) {			\
			    px += j;					\
			    for (i = j; i < n; ++i, ++px)		\
				DO_INCR;				\
			}						\
		    } else {						\
			for (i = j = 0; j < n; ++j, i = j) {		\
			    px += j+1;					\
			    for (i = j+1; i < n; ++i, ++px)		\
				DO_INCR;				\
			}						\
		    }							\
		}							\
	    } else {							\
		if (ul == 'U') { /* packed upper triangular */		\
		    if (di == 'N') {					\
			for (j = 0; j < n; ++j)				\
			    for (i = 0; i <= j; ++i, ++px)		\
				DO_INCR;				\
		    } else {						\
			for (j = 0; j < n; ++j) {			\
			    for (i = 0; i < j; ++i, ++px)		\
				DO_INCR;				\
			    ++px;					\
			}						\
		    }							\
		} else { /* packed lower triangular */			\
		    if (di == 'N') {					\
			for (j = 0; j < n; ++j)				\
			    for (i = j; i < n; ++i, ++px)		\
				DO_INCR;				\
		    } else {						\
			for (i = j = 0; j < n; ++j, i = j) {		\
			    ++px;					\
			    for (i = j+1; i < n; ++i, ++px)		\
				DO_INCR;				\
			}						\
		    }							\
		}							\
	    }								\
	} else {							\
	    if (cl[2] != 'p') {						\
		if (ul == 'U') { /* unpacked upper symmetric */		\
		    for (j = 0; j < n; ++j) {				\
			for (i = 0; i < j; ++i, ++px)			\
			    DO_INCR_SYMM;				\
			DO_INCR;					\
			px += n-j;					\
		    }							\
		} else { /* unpacked lower symmetric */			\
		    for (i = j = 0; j < n; ++j, i = j) {		\
			px += j;					\
			DO_INCR;					\
			++px;						\
			for (i = j+1; i < n; ++i, ++px)			\
			    DO_INCR_SYMM;				\
		    }							\
		}							\
	    } else {							\
		if (ul == 'U') { /* packed upper symmetric */		\
		    for (j = 0; j < n; ++j) {				\
			for (i = 0; i < j; ++i, ++px)			\
			    DO_INCR_SYMM;				\
			DO_INCR;					\
			++px;						\
		    }							\
		} else { /* packed lower symmetric */			\
		    for (i = j = 0; j < n; ++j, i = j) {		\
			DO_INCR;					\
			++px;						\
			for (i = j+1; i < n; ++i, ++px)			\
			    DO_INCR_SYMM;				\
		    }							\
		}							\
	    }								\
	}								\
    } while (0)
    
#define DENSE_ROWSUMS(_CTYPE1_, _PTR1_, _CTYPE2_, _PTR2_)		\
    do {								\
	_CTYPE1_ *pres = _PTR1_(res), u = (di == 'N') ? ZERO : ONE;	\
	_CTYPE2_ *px   = _PTR2_(x);					\
	if (doNaRm && doMean && cl[0] != 'n') {				\
	    Calloc_or_Alloca_TO(pcount, m, int);			\
	    for (i = 0; i < m; ++i) {					\
		pres[i] = u;						\
		pcount[i] = n;						\
	    }								\
	} else {							\
	    for (i = 0; i < m; ++i)					\
		pres[i] = u;						\
	}								\
	DENSE_ROWSUMS_LOOP;						\
    } while (0)
    
    switch (cl[0]) {
    case 'n':
	
#define ZERO         0.0
#define ONE          1.0
#define DO_INCR      if (*px) pres[i] += 1.0
#define DO_INCR_SYMM				\
	do {					\
	    if (*px) {				\
		pres[i] += 1.0;			\
		pres[j] += 1.0;			\
	    }					\
	} while (0)
	
	DENSE_ROWSUMS(double, REAL, int, LOGICAL);
	break;

#undef DO_INCR
#undef DO_INCR_SYMM
	
    case 'l':

#define DO_INCR					\
	do {					\
	    if (*px != NA_LOGICAL) {		\
		if (*px)			\
		    pres[i] += 1.0;		\
	    } else if (!doNaRm)			\
		pres[i] = NA_REAL;		\
	    else if (doMean)			\
		--pcount[i];			\
	} while (0)
#define DO_INCR_SYMM				\
	do {					\
	    if (*px != NA_LOGICAL) {		\
		if (*px) {			\
		    pres[i] += 1.0;		\
		    pres[j] += 1.0;		\
		}				\
	    } else if (!doNaRm) {		\
		pres[i] = NA_REAL;		\
		pres[j] = NA_REAL;		\
	    } else if (doMean) {		\
		--pcount[i];			\
		--pcount[j];			\
	    }					\
	} while (0)

	DENSE_ROWSUMS(double, REAL, int, LOGICAL);
	break;

#undef DO_INCR
#undef DO_INCR_SYMM
	
    case 'i':

#define DO_INCR					\
	do {					\
	    if (*px != NA_INTEGER)		\
		pres[i] += *px;			\
	    else if (!doNaRm)			\
		pres[i] = NA_REAL;		\
	    else if (doMean)			\
		--pcount[i];			\
	} while (0)
#define DO_INCR_SYMM				\
	do {					\
	    if (*px != NA_INTEGER) {		\
		pres[i] += *px;			\
		pres[j] += *px;			\
	    } else if (!doNaRm) {		\
		pres[i] = NA_REAL;		\
		pres[j] = NA_REAL;		\
	    } else if (doMean) {		\
		--pcount[i];			\
		--pcount[j];			\
	    }					\
	} while (0)
	
	DENSE_ROWSUMS(double, REAL, int, INTEGER);
	break;
	
#undef DO_INCR
#undef DO_INCR_SYMM
	
    case 'd':

#define DO_INCR					\
	do {					\
	    if (!(doNaRm && ISNAN(*px)))	\
		pres[i] += *px;			\
	    else if (doMean)			\
		--pcount[i];			\
	} while (0)
#define DO_INCR_SYMM				\
	do {					\
	    if (!(doNaRm && ISNAN(*px))) {	\
		pres[i] += *px;			\
		pres[j] += *px;			\
	    } else if (doMean) {		\
		--pcount[i];			\
		--pcount[j];			\
	    }					\
	} while (0)
	
	DENSE_ROWSUMS(double, REAL, double, REAL);
	break;
	
#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_INCR_SYMM

    case 'z':

#define ZERO         Matrix_zzero
#define ONE          Matrix_zone
#define DO_INCR								\
	do {								\
	    if (!(doNaRm && (ISNAN((*px).r) || ISNAN((*px).i)))) {	\
		pres[i].r += (*px).r;					\
		pres[i].i += (*px).i;					\
	    } else if (doMean)						\
		--pcount[i];						\
	} while (0)	
#define DO_INCR_SYMM							\
	do {								\
	    if (!(doNaRm && (ISNAN((*px).r) || ISNAN((*px).i)))) {	\
		pres[i].r += (*px).r;					\
		pres[i].i += (*px).i;					\
		pres[j].r += (*px).r;					\
		pres[j].i += (*px).i;					\
	    } else if (doMean) {					\
		--pcount[i];						\
		--pcount[j];						\
	    }								\
	} while (0)
	
	DENSE_ROWSUMS(Rcomplex, COMPLEX, Rcomplex, COMPLEX);
	break;
	
#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_INCR_SYMM
	
    default:
	break;
    }

#undef DENSE_ROWSUMS
#undef DENSE_ROWSUMS_LOOP

    if (doMean) {
	if (cl[0] != 'z') {
	    double *pres = REAL(res);
	    if (doNaRm && cl[0] != 'n') {
		for (i = 0; i < m; ++i)
		    pres[i] /= pcount[i];
		Free_FROM(pcount, m);
	    } else {
		for (i = 0; i < m; ++i)
		    pres[i] /= n;
	    }
	} else {
	    Rcomplex *pres = COMPLEX(res);
	    if (doNaRm) {
		for (i = 0; i < m; ++i) {
		    pres[i].r /= pcount[i];
		    pres[i].i /= pcount[i];
		}
		Free_FROM(pcount, m);
	    } else {
		for (i = 0; i < m; ++i) {
		    pres[i].r /= n;
		    pres[i].i /= n;
		}
	    }
	}
    }

    SEXP dimnames;
    if (cl[1] != 's')
	PROTECT(dimnames = GET_SLOT(obj, Matrix_DimNamesSym));
    else
	PROTECT(dimnames = get_symmetrized_DimNames(obj, -1));
    SEXP nms = VECTOR_ELT(dimnames, 0);
    if (!isNull(nms))
	setAttrib(res, R_NamesSymbol, nms);
    
    UNPROTECT(3); /* dimnames, x, res */
    return res;
}

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
int left_cyclic(double x[], int ldx, int j, int k,
		double cosines[], double sines[])
{
    if (j >= k)
	error(_("incorrect left cyclic shift, j (%d) >= k (%d)"), j, k);
    if (j < 0)
	error(_("incorrect left cyclic shift, j (%d) < 0"), j, k);
    if (ldx < k)
	error(_("incorrect left cyclic shift, k (%d) > ldx (%d)"), k, ldx);

    double *lastcol = (double *) R_alloc((size_t) k + 1, sizeof(double));
    int i;
				/* keep a copy of column j */
    for(i = 0; i <= j; i++) lastcol[i] = x[i + j*ldx];
				/* For safety, zero the rest */
    for(i = j+1; i <= k; i++) lastcol[i] = 0.;
    for(int jj = j+1, ind = 0; jj <= k; jj++, ind++) { /* columns to be shifted */
	int diagind = jj*(ldx+1); //  ind == (jj-j) - 1
	double tmp = x[diagind], cc, ss;
	/* Calculate the Givens rotation. */
				/* This modified the super-diagonal element */
	F77_CALL(drotg)(x + diagind-1, &tmp, cosines + ind, sines + ind);
	cc = cosines[ind]; ss = sines[ind];
				/* Copy column jj+1 to column jj. */
	for(i = 0; i < jj; i++) x[i + (jj-1)*ldx] = x[i+jj*ldx];
				/* Apply rotation to columns up to k */
	for(i = jj; i < k; i++) {
	    tmp = cc*x[(jj-1)+i*ldx] + ss*x[jj+i*ldx];
	    x[jj+i*ldx] = cc*x[jj+i*ldx] - ss*x[(jj-1)+i*ldx];
	    x[(jj-1)+i*ldx] = tmp;
	}
				/* Apply rotation to lastcol */
	lastcol[jj] = -ss*lastcol[jj-1]; lastcol[jj-1] *= cc;
    }
				/* Copy lastcol to column k */
    for(i = 0; i <= k; i++) x[i+k*ldx] = lastcol[i];
    return 0;
}

static
SEXP getGivens(double x[], int ldx, int jmin, int rank)
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
	error(_("Unknown error in getGivens"));
    UNPROTECT(1);
    return ans;
}

SEXP checkGivens(SEXP X, SEXP jmin, SEXP rank)
{
    SEXP ans = PROTECT(allocVector(VECSXP, 2)),
	Xcp = PROTECT(duplicate(X));
    int  *Xdims;

    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP));
    SET_VECTOR_ELT(ans, 1, getGivens(REAL(Xcp), Xdims[0],
				     asInteger(jmin), asInteger(rank)));
    SET_VECTOR_ELT(ans, 0, Xcp);
    UNPROTECT(2);
    return ans;
}

SEXP lsq_dense_Chol(SEXP X, SEXP y)
{
    SEXP ans;
    double d_one = 1., d_zero = 0.;

    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    int *Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP)),
	n = Xdims[0],
	p = Xdims[1];
    if (!(isReal(y) & isMatrix(y)))
	error(_("y must be a numeric (double precision) matrix"));
    int *ydims = INTEGER(coerceVector(getAttrib(y, R_DimSymbol), INTSXP));
    if (ydims[0] != n)
	error(_(
	    "number of rows in y (%d) does not match number of rows in X (%d)"),
	    ydims[0], n);
    int k = ydims[1];
    if (k < 1 || p < 1) return allocMatrix(REALSXP, p, k);
    ans = PROTECT(allocMatrix(REALSXP, p, k));
    F77_CALL(dgemm)("T", "N", &p, &k, &n, &d_one, REAL(X), &n, REAL(y), &n,
		    &d_zero, REAL(ans), &p FCONE FCONE);
    double *xpx = (double *) R_alloc((size_t) p * p, sizeof(double));
    F77_CALL(dsyrk)("U", "T", &p, &n, &d_one, REAL(X), &n, &d_zero,
		    xpx, &p FCONE FCONE);
    int info;
    F77_CALL(dposv)("U", &p, &k, xpx, &p, REAL(ans), &p, &info FCONE);
    if (info) error(_("Lapack routine dposv returned error code %d"), info);
    UNPROTECT(1);
    return ans;
}

SEXP lsq_dense_QR(SEXP X, SEXP y)
{
    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    int *Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP)),
	n = Xdims[0],
	p = Xdims[1];
    if (!(isReal(y) & isMatrix(y)))
	error(_("y must be a numeric (double precision) matrix"));
    int *ydims = INTEGER(coerceVector(getAttrib(y, R_DimSymbol), INTSXP));
    if (ydims[0] != n)
	error(_(
	    "number of rows in y (%d) does not match number of rows in X (%d)"),
	    ydims[0], n);
    int k = ydims[1];
    if (k < 1 || p < 1) return allocMatrix(REALSXP, p, k);
    double tmp,
	*xvals = (double *) Memcpy(R_alloc((size_t) n * p, sizeof(double)),
				   REAL(X),
				   (size_t) n * p);
    SEXP ans = PROTECT(duplicate(y));
    int lwork = -1, info;
    F77_CALL(dgels)("N", &n, &p, &k, xvals, &n, REAL(ans), &n,
		    &tmp, &lwork, &info FCONE);
    if (info)
	error(_("First call to Lapack routine dgels returned error code %d"),
	      info);
    lwork = (int) tmp;
    double *work = (double *) R_alloc((size_t) lwork, sizeof(double));
    F77_CALL(dgels)("N", &n, &p, &k, xvals, &n, REAL(ans), &n,
		    work, &lwork, &info FCONE);
    if (info)
	error(_("Second call to Lapack routine dgels returned error code %d"),
	      info);
    UNPROTECT(1);
    return ans;
}

/* Rank-Correcting/Adapting LAPACK  QR Decomposition
 * From Doug Bates' initial import; __unused__

 * Provides a qr() with 'rcond' and rank reduction while(rcond < tol),
 * possibly via Givens rotations but WITHOUT PIVOTING

 * .Call(Matrix:::lapack_qr, A, 1e-17) --> ~/R/MM/Pkg-ex/Matrix/qr-rank-deficient.R

 * TODO: export as Matrix::qrNoPiv() or qr1()  or similar
 */
SEXP lapack_qr(SEXP Xin, SEXP tl)
{
    if (!(isReal(Xin) & isMatrix(Xin)))
	error(_("X must be a real (numeric) matrix"));
    double tol = asReal(tl);
    if (tol < 0.) error(_("tol, given as %g, must be non-negative"), tol);
    if (tol > 1.) error(_("tol, given as %g, must be <= 1"), tol);
    SEXP ans = PROTECT(allocVector(VECSXP,5)), X, qraux, pivot;
    SET_VECTOR_ELT(ans, 0, X = duplicate(Xin));
    int *Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP)),
	n = Xdims[0], i,
	p = Xdims[1],
	trsz = (n < p) ? n : p ; /* size of triangular part of decomposition */

    SET_VECTOR_ELT(ans, 2, qraux = allocVector(REALSXP, trsz));
    SET_VECTOR_ELT(ans, 3, pivot = allocVector(INTSXP, p));
    for (i = 0; i < p; i++) INTEGER(pivot)[i] = i + 1;
    SEXP nms,
	Givens = PROTECT(allocVector(VECSXP, trsz - 1));
    setAttrib(ans, R_NamesSymbol, nms = allocVector(STRSXP, 5));
    SET_STRING_ELT(nms, 0, mkChar("qr"));
    SET_STRING_ELT(nms, 1, mkChar("rank"));
    SET_STRING_ELT(nms, 2, mkChar("qraux"));
    SET_STRING_ELT(nms, 3, mkChar("pivot"));
    SET_STRING_ELT(nms, 4, mkChar("Givens"));
    int rank = trsz,
	nGivens = 0;
    double rcond = 0.;
    if (n > 0 && p > 0) {
	int  info, *iwork, lwork;
	double *xpt = REAL(X), *work, tmp;

	lwork = -1;
	F77_CALL(dgeqrf)(&n, &p, xpt, &n, REAL(qraux), &tmp, &lwork, &info);
	if (info)
	    error(_("First call to dgeqrf returned error code %d"), info);
	lwork = (int) tmp;
	work = (double *) R_alloc(((size_t) lwork < (size_t) 3 * trsz)
				  ? (size_t) 3 * trsz : (size_t) lwork,
				  sizeof(double));
	F77_CALL(dgeqrf)(&n, &p, xpt, &n, REAL(qraux), work, &lwork, &info);
	if (info)
	    error(_("Second call to dgeqrf returned error code %d"), info);
	iwork = (int *) R_alloc((size_t) trsz, sizeof(int));
	F77_CALL(dtrcon)("1", "U", "N", &rank, xpt, &n, &rcond,
			 work, iwork, &info FCONE FCONE FCONE);
	if (info)
	    error(_("Lapack routine dtrcon returned error code %d"), info);
	while (rcond < tol) {	/* check diagonal elements */
	    double minabs = (xpt[0] < 0.) ? -xpt[0]: xpt[0];
	    int jmin = 0;
	    for (i = 1; i < rank; i++) {
		double el = xpt[i*n]; // had  i*(n+1)  which looks wrong to MM
		if(el < 0.) el = -el;
		if (el < minabs) {
		    jmin = i;
		    minabs = el;
		}
	    }
	    if (jmin < (rank - 1)) {
		SET_VECTOR_ELT(Givens, nGivens, getGivens(xpt, n, jmin, rank));
		nGivens++;
	    } // otherwise jmin == (rank - 1) , so just "drop that column"
	    rank--;
	    // new  rcond := ... for reduced rank
	    F77_CALL(dtrcon)("1", "U", "N", &rank, xpt, &n, &rcond,
			     work, iwork, &info FCONE FCONE FCONE);
	    if (info)
		error(_("Lapack routine dtrcon returned error code %d"), info);
	}
    }
    SEXP Gcpy, sym;
    SET_VECTOR_ELT(ans, 4, Gcpy = allocVector(VECSXP, nGivens));
    for (i = 0; i < nGivens; i++)
	SET_VECTOR_ELT(Gcpy, i, VECTOR_ELT(Givens, i));
    SET_VECTOR_ELT(ans, 1, ScalarInteger(rank));
    sym = PROTECT(install("useLAPACK")); setAttrib(ans, sym, ScalarLogical(1)); UNPROTECT(1);
    sym = PROTECT(install("rcond"));     setAttrib(ans, sym, ScalarReal(rcond));UNPROTECT(1);
    UNPROTECT(2);
    return ans;
}

/* MJ: no longer needed ... prefer R_dense_as_sparse() above */
#if 0

SEXP dense_to_Csparse(SEXP x)
{
    SEXP ge_x = PROTECT(dense_as_general(x, '.', 2, 0)),
	Dim = GET_SLOT(ge_x, Matrix_DimSym);
    int *dims = INTEGER(Dim);
    Rboolean longi = (dims[0] * (double)dims[1] > INT_MAX);
    // int itype = longi ? CHOLMOD_LONG : CHOLMOD_INT;
    CHM_DN chxd = AS_CHM_xDN(ge_x); // cholmod_dense (has no itype)
    CHM_SP chxs;
    /* cholmod_dense_to_sparse() in CHOLMOD/Core/ below does only work for
       "REAL" 'xtypes', i.e. *not* for "nMatrix".
       ===> need "_x" in above AS_CHM_xDN() call.

       Also it cannot keep symmetric / triangular, hence the
       dense_as_general() above.  Note that this is already a *waste*
       for symmetric matrices; However, we could conceivably use an
       enhanced cholmod_dense_to_sparse(), with an extra boolean
       argument for symmetry.
    */
#define DLONG
/* You can try defining DLONG -- then just get a seg.fault :
 *  I think it is because of this in  ./CHOLMOD/Include/cholmod_core.h :
 *
 The itype of all parameters for all CHOLMOD routines must match.
 --- ^^^^^ ------------------------------------------------------
 but then as_cholmod_dense should *not* make a difference: cholmod_dense has *no* itype   (????)
*/
    if(longi) { // calling cholmod_dense_to_sparse() gives wrong matrix
#ifdef DLONG
	chxs = cholmod_l_dense_to_sparse(chxd, 1, &cl);
	// in gdb, I found that 'chxs' seems "basically empty": all
	// p chxs->foo   give ''Cannot access memory at address 0x....''
	// for now rather give error:
	if(cl.status)
	    error(_("dense_to_Csparse(<LARGE>): cholmod_l_dense_to_sparse failure status=%d"),
		  cl.status);
#else
        error(_("Matrix dimension %d x %d (= %g) too large [FIXME calling cholmod_l_dense_to_sparse]"),
	      m,n, m * (double)n);
#endif
    } else { // fits, using integer (instead of long int) 'itype'
	chxs = cholmod_dense_to_sparse(chxd, 1, &c);
    }

    int Rkind = (chxd->xtype == CHOLMOD_REAL) ? Real_KIND2(x) : 0;
    /* Note: when 'x' was integer Matrix, Real_KIND(x) = -1, but *_KIND2(.) = 0 */
    R_CheckStack();

    UNPROTECT(1);
    /* chm_sparse_to_SEXP() *could* deal with symmetric
     * if chxs had such an stype; and we should be able to use uplo below */
    return chm_sparse_to_SEXP(chxs, 1, 0/*TODO: uplo_P(x) if x has an uplo slot*/,
			      Rkind, "",
			      isMatrix(x) ? getAttrib(x, R_DimNamesSymbol)
			      : GET_SLOT(x, Matrix_DimNamesSym));
}

#endif /* MJ */

/* MJ: no longer needed ... prefer (un)?packedMatrix_(symm|skew)part() */
#if 0

SEXP ddense_symmpart(SEXP x)
/* Class of the value will be dsyMatrix */
{
    SEXP dx = PROTECT(dense_as_general(x, 'd', 2, 0));
    int *adims = INTEGER(GET_SLOT(dx, Matrix_DimSym)), n = adims[0];

    if(n != adims[1]) {
	error(_("matrix is not square! (symmetric part)"));
	return R_NilValue; /* -Wall */
    } else {
	SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS("dsyMatrix"));
	double *xx = REAL(GET_SLOT(dx, Matrix_xSym));

	/* only need to assign the *upper* triangle (uplo = "U");
	 * noting that diagonal remains unchanged */
	R_xlen_t n_ = n;
	for (int j = 0; j < n; j++) {
	    for (int i = 0; i < j; i++) {
		xx[j * n_ + i] = (xx[j * n_ + i] + xx[i * n_ + j]) / 2.;
	    }
	}

/* Copy dx to ans:
   Because slots of dx are freshly allocated and dx will not
   be used, we use the slots themselves and don't duplicate */	
#define MK_SYMMETRIC_DIMNAMES_AND_RETURN(_J_ /* -1|0|1 */)		\
	SET_SLOT(ans, Matrix_DimSym,  GET_SLOT(dx, Matrix_DimSym));	\
	set_symmetrized_DimNames(ans, GET_SLOT(dx, Matrix_DimNamesSym), _J_); \
	SET_SLOT(ans, Matrix_uploSym, mkString((_J_) ? "U" : "L"));	\
	SET_SLOT(ans, Matrix_xSym,    GET_SLOT(dx, Matrix_xSym));	\
	UNPROTECT(2);							\
	return ans

        MK_SYMMETRIC_DIMNAMES_AND_RETURN(-1);
    }
}

SEXP ddense_skewpart(SEXP x)
/* Class of the value will be dgeMatrix */
{
    SEXP dx = PROTECT(dense_as_general(x, 'd', 2, 0));
    int *adims = INTEGER(GET_SLOT(dx, Matrix_DimSym)), n = adims[0];

    if(n != adims[1]) {
	error(_("matrix is not square! (skew-symmetric part)"));
	return R_NilValue; /* -Wall */
    } else {
	SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix"));
	double *xx = REAL(GET_SLOT(dx, Matrix_xSym));
	R_xlen_t n_ = n;

	for (int j = 0; j < n_; j++) {
	    xx[j * n_ + j] = 0.;
	    for (int i = 0; i < j; i++) {
		double s = (xx[j * n_ + i] - xx[i * n_ + j]) / 2.;
		xx[j * n_ + i] =  s;
		xx[i * n_ + j] = -s;
	    }
	}

        MK_SYMMETRIC_DIMNAMES_AND_RETURN(-1);
    }
}

#endif /* MJ */

/* MJ: no longer needed ... prefer (un)?packedMatrix_force_symmetric() */
#if 0

SEXP dense_to_symmetric(SEXP x, SEXP uplo, SEXP symm_test)
/* Class of result will be [dln]syMatrix */
{
/*== FIXME: allow  uplo = NA   and then behave a bit like symmpart():
 *== -----  would use the *dimnames* to determine U or L   (??)
 */
    
    int symm_tst = asLogical(symm_test);
    SEXP dx = PROTECT(dense_as_general(x, '.', 2, 0));
    SEXP ans;
    const char *cl = class_P(dx);
    /* same as in ..._geMatrix() above:*/
    enum dense_enum M_type = ( (cl[0] == 'd') ? ddense :
			       ((cl[0] == 'l') ? ldense : ndense) );
    int *adims = INTEGER(GET_SLOT(dx, Matrix_DimSym)), n = adims[0];
    if(n != adims[1]) {
	UNPROTECT(1);
	error(_("dense_to_symmetric(): matrix is not square!"));
	return R_NilValue; /* -Wall */
    }

    if(symm_tst) {
	int i, j;
	R_xlen_t n_ = n;
	
#define CHECK_SYMMETRIC							\
	do {								\
	    for (j = 0; j < n; j++) {					\
		for (i = 0; i < j; i++)	{				\
		    if(xx[j * n_ + i] != xx[i * n_ + j]) {		\
			UNPROTECT(1);					\
			error(_("matrix is not symmetric [%d,%d]"), i+1, j+1); \
			return R_NilValue; /* -Wall */			\
		    }							\
		}							\
	    }								\
	} while (0)
	
	if(M_type == ddense) {
	    double *xx = REAL(GET_SLOT(dx, Matrix_xSym));
	    CHECK_SYMMETRIC;
	} else { /* (M_type == ldense || M_type == ndense) */
	    int *xx = LOGICAL(GET_SLOT(dx, Matrix_xSym));
	    CHECK_SYMMETRIC;
	}
    }
#undef CHECK_SYMMETRIC
    
    ans = PROTECT(NEW_OBJECT_OF_CLASS(M_type == ddense
				      ? "dsyMatrix"
				      : (M_type == ldense
					 ? "lsyMatrix"
					 : "nsyMatrix")));
    int uploT = (*CHAR(asChar(uplo)) == 'U');
    MK_SYMMETRIC_DIMNAMES_AND_RETURN(uploT);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_dense_band() above */
#if 0

SEXP dense_band(SEXP x, SEXP k1P, SEXP k2P)
/* Always returns a full matrix with entries outside the band zeroed
 * Class of the value can be [dln]trMatrix or [dln]geMatrix
 */
{
    int k1 = asInteger(k1P), k2 = asInteger(k2P);

    if (k1 > k2) {
	error(_("Lower band %d > upper band %d"), k1, k2);
	return R_NilValue; /* -Wall */
    }
    else {
	SEXP ans = PROTECT(dense_as_general(x, '.', 2, 0));
	int *adims = INTEGER(GET_SLOT(ans, Matrix_DimSym)),
	    j, m = adims[0], n = adims[1],
	    sqr = (adims[0] == adims[1]),
	    tru = (k1 >= 0), trl = (k2 <= 0);
	const char *cl = class_P(ans);
	enum dense_enum M_type = ( (cl[0] == 'd') ? ddense :
			      ((cl[0] == 'l') ? ldense : ndense));


#define SET_ZERO_OUTSIDE				\
	for (j = 0; j < n; j++) {			\
	    int i, i1 = j - k2, i2 = j + 1 - k1;	\
	    R_xlen_t jm = j * (R_xlen_t) m;		\
	    if(i1 > m) i1 = m;				\
	    if(i2 < 0) i2 = 0;				\
	    for (i = 0; i < i1; i++) xx[i + jm] = 0;	\
	    for (i = i2; i < m; i++) xx[i + jm] = 0;	\
	}

	if(M_type == ddense) {
	    double *xx = REAL(GET_SLOT(ans, Matrix_xSym));
	    SET_ZERO_OUTSIDE
	}
	else { /* (M_type == ldense || M_type == ndense) */
	    int *xx = LOGICAL(GET_SLOT(ans, Matrix_xSym));
	    SET_ZERO_OUTSIDE
	}

	if (!sqr || (!tru && !trl)) { /* return the *geMatrix */
	    UNPROTECT(1);
	    return ans;
	}
	else {
	    /* Copy ans to a *trMatrix object (must be square) */
	    SEXP aa = PROTECT(NEW_OBJECT_OF_CLASS(M_type == ddense
						  ? "dtrMatrix"
						  : (M_type == ldense
						     ? "ltrMatrix"
						     : "ntrMatrix")));
	    /* Because slots of ans are freshly allocated and ans will not be
	     * used, we use the slots themselves and don't duplicate */
	    SET_SLOT(aa, Matrix_xSym,        GET_SLOT(ans, Matrix_xSym));
	    SET_SLOT(aa, Matrix_DimSym,      GET_SLOT(ans, Matrix_DimSym));
	    SET_SLOT(aa, Matrix_DimNamesSym, GET_SLOT(ans, Matrix_DimNamesSym));
	    SET_SLOT(aa, Matrix_diagSym,     mkString("N"));
	    SET_SLOT(aa, Matrix_uploSym,     mkString(tru ? "U" : "L"));
	    UNPROTECT(2);
	    return aa;
	}
    }
}

#endif /* MJ */
