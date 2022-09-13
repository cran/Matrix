#include "sparse.h"

SEXP sparse_as_dense(SEXP from, int packed)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_as_dense");
    const char *clf = valid[ivalid];
    packed = packed && clf[1] != 'g';

    char clt[] = "...Matrix";
    clt[0] = clf[0];
    clt[1] = clf[1];
    clt[2] = (clf[1] == 'g')
	? 'e' : ((packed) ? 'p' : ((clf[1] == 't') ? 'r' : 'y'));
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x0 = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue;
    
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    double mn = (double) m * n, dnx = (packed) ? 0.5 * (mn + n) : mn,
	dsize = dnx * kind2size(clt[0]);
    if (dnx > R_XLEN_T_MAX)
	error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
    if (packed && clf[2] != 'C' && mn > R_XLEN_T_MAX)
	error(_("coercing n-by-n [RT]sparseMatrix to packedMatrix "
		"is not supported for n*n exceeding R_XLEN_T_MAX"));
    if (dsize > 0x1p+30 /* 1 GiB */)
	warning(_("sparse->dense coercion: "
		  "allocating vector of size %0.1f GiB"),
		0x1p-30 * dsize);
    
    char ul = 'U', di = 'N';
    SEXPTYPE tx = kind2type(clf[0]);
    R_xlen_t nx = (R_xlen_t) dnx;
    SEXP x1 = PROTECT(allocVector(tx, nx));
    
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    SET_SLOT(to, Matrix_xSym, x1);
    if (clf[1] != 'g') {
	SEXP uplo = GET_SLOT(from, Matrix_uploSym);
	SET_SLOT(to, Matrix_uploSym, uplo);
	ul = *CHAR(STRING_ELT(uplo, 0));
	if (clf[1] == 't') {
	    SEXP diag = GET_SLOT(from, Matrix_diagSym);
	    SET_SLOT(to, Matrix_diagSym, diag);
	    di = *CHAR(STRING_ELT(diag, 0));
	}
    }

    /* It remains to fill 'x' ... */
    
    if (clf[2] == 'C') {
	
	SEXP p = GET_SLOT(from, Matrix_pSym),
	    i = GET_SLOT(from, Matrix_iSym);
	int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
	++pp;
	
#define SAD_C(_VAL_)						\
	do {							\
	    if (!packed) {					\
		for (j = 0, k = 0; j < n; ++j, px1 += m) {	\
		    kend = pp[j];				\
		    while (k < kend) {				\
			px1[*pi] = _VAL_;			\
			++pi; ++k;				\
		    }						\
		}						\
	    } else if (ul == 'U') {				\
		for (j = 0, k = 0; j < n; px1 += (++j)) {	\
		    kend = pp[j];				\
		    while (k < kend) {				\
			px1[*pi] = _VAL_;			\
			++pi; ++k;				\
		    }						\
		}						\
	    } else {						\
		for (j = 0, k = 0; j < n; px1 += n-(j++)) {	\
		    kend = pp[j];				\
		    while (k < kend) {				\
			px1[*pi - j] = _VAL_;			\
			++pi; ++k;				\
		    }						\
		}						\
	    }							\
	} while (0)
	
#define SAD_C_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)		\
	do {						\
	    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	    Memzero(px1, nx);				\
	    SAD_C(*(px0++));				\
	} while (0)
	    
	if (clf[0] == 'n') {
	    int *px1 = LOGICAL(x1);
	    Memzero(px1, nx);
	    SAD_C(1);
	} else {
	    SPARSE_CASES(tx, SAD_C_X);
	}
	
#undef SAD_C_X
#undef SAD_C
	
    } else if (clf[2] == 'R') {
	
	SEXP p = GET_SLOT(from, Matrix_pSym),
	    j = GET_SLOT(from, Matrix_jSym);
	int *pp = INTEGER(p), *pj = INTEGER(j), i, k, kend;
	++pp;
	
#define SAD_R(_VAL_)						\
	do {							\
	    if (!packed) {					\
		R_xlen_t m1 = (R_xlen_t) m;			\
		for (i = 0, k = 0; i < m; ++i, ++px1) {		\
		    kend = pp[i];				\
		    while (k < kend) {				\
			px1[m1 * *pj] = _VAL_;			\
			++pj; ++k;				\
		    }						\
		}						\
	    } else if (ul == 'U') {				\
		for (i = 0, k = 0; i < n; ++i) {		\
		    kend = pp[i];				\
		    while (k < kend) {				\
			px1[PM_AR21_UP(i, *pj)] = _VAL_;	\
			++pj; ++k;				\
		    }						\
		}						\
	    } else {						\
		R_xlen_t n2 = (R_xlen_t) n * 2;			\
		for (i = 0, k = 0; i < n; ++i) {		\
		    kend = pp[i];				\
		    while (k < kend) {				\
			px1[PM_AR21_LO(i, *pj, n2)] = _VAL_;	\
			++pj; ++k;				\
		    }						\
		}						\
	    }							\
	} while (0)

#define SAD_R_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)		\
	do {						\
	    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	    Memzero(px1, nx);				\
	    SAD_R(*(px0++));				\
	} while (0)

	if (clf[0] == 'n') {
	    int *px1 = LOGICAL(x1);
	    Memzero(px1, nx);
	    SAD_R(1);
	} else {
	    SPARSE_CASES(tx, SAD_R_X);
	}
	
#undef SAD_R_X
#undef SAD_R
	
    } else { /* clf[2] == 'T' */
	
	SEXP i = GET_SLOT(from, Matrix_iSym),
	    j = GET_SLOT(from, Matrix_jSym);
	int *pi = INTEGER(i), *pj = INTEGER(j);
	R_xlen_t index, k, nnz = XLENGTH(i);

#define SAD_T(_DO_FOR_)					\
	do {						\
	    if (!packed) {				\
		R_xlen_t m1 = (R_xlen_t) m;		\
		_DO_FOR_(m1 * *pj + *pi);		\
	    } else if (ul == 'U') {			\
		_DO_FOR_(PM_AR21_UP(*pi, *pj));		\
	    } else {					\
		R_xlen_t n2 = (R_xlen_t) n * 2;		\
		_DO_FOR_(PM_AR21_LO(*pi, *pj, n2));	\
	    }						\
	} while (0)
	
	switch (clf[0]) {
	case 'd':
	{
	    double *px0 = REAL(x0), *px1 = REAL(x1);
	    Memzero(px1, nx);

#define SAD_T_D(_INDEX_)					\
	    do {						\
		for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0)	\
		    px1[_INDEX_] += *px0;			\
	    } while (0)
	    
	    SAD_T(SAD_T_D);
	    break;
	}
	case 'l':
	{
	    int *px0 = LOGICAL(x0), *px1 = LOGICAL(x1);
	    Memzero(px1, nx);

#define SAD_T_L(_INDEX_)					\
	    do {						\
		for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0) {	\
		    if (*px0 != 0) {				\
			index = _INDEX_;			\
			if (*px0 != NA_LOGICAL)			\
			    px1[index] = 1;			\
			else if (px1[index] == 0)		\
			    px1[index] = NA_LOGICAL;		\
		    }						\
		}						\
	    } while (0)

	    SAD_T(SAD_T_L);
	    break;
	}
	case 'n':
	{
	    int *px1 = LOGICAL(x1);
	    Memzero(px1, nx);
	    
#define SAD_T_N(_INDEX_)				\
	    do {					\
		for (k = 0; k < nnz; ++k, ++pi, ++pj)	\
		    px1[_INDEX_] = 1;			\
	    } while (0)

	    SAD_T(SAD_T_N);
	    break;
	}
	case 'i':
	{
	    int *px0 = INTEGER(x0), *px1 = INTEGER(x1);
	    Memzero(px1, nx);

/* FIXME: not detecting integer overflow here */
#define SAD_T_I(_INDEX_) SAD_T_D(_INDEX_)
	    
	    SAD_T(SAD_T_I);
	    break;
	}
	case 'z':
	{
	    Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
	    Memzero(px1, nx);
	    R_xlen_t index;

#define SAD_T_Z(_INDEX_)					\
	    do {						\
		for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0) {	\
		    index = _INDEX_;				\
		    px1[index].r += (*px0).r;			\
		    px1[index].i += (*px0).i;			\
		}						\
	    } while (0)
	    
	    SAD_T(SAD_T_Z);
	    break;
	}
	default:
	    ERROR_INVALID_TYPE("'x' slot", tx, "sparse_as_dense");
	    break;
	}

#undef SAD_T_D
#undef SAD_T_L
#undef SAD_T_N
#undef SAD_T_I
#undef SAD_T_Z
#undef SAD_T
    
    }
    
    if (di != 'N') {

	int j;
	
#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_)		\
	do {						\
	    _CTYPE_ *px1 = _PTR_(x1);			\
	    if (!packed) {				\
		R_xlen_t n1a = (R_xlen_t) n + 1;	\
		for (j = 0; j < n; ++j, px1 += n1a)	\
		    *px1 = _ONE_;			\
	    } else if (ul == 'U') {			\
		for (j = 0; j < n; px1 += (++j)+1)	\
		    *px1 = _ONE_;			\
	    } else {					\
		for (j = 0; j < n; px1 += n-(j++))	\
		    *px1 = _ONE_;			\
	    }						\
	} while (0)
	    
	SPARSE_CASES(tx, SET1);

#undef SET1
	
    }

    UNPROTECT(2);
    return to;
}

/* as(<[CRT]sparseMatrix>,       "denseMatrix") */
/* as(<[CRT]sparseMatrix>, "(un)?packedMatrix") */
SEXP R_sparse_as_dense(SEXP from, SEXP packed)
{
    return sparse_as_dense(from, asLogical(packed));
}

/* as(<[CRT]sparseMatrix>, "matrix") */
SEXP R_sparse_as_matrix(SEXP from)
{
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(from = sparse_as_dense(from, 0), &pid);
    REPROTECT(from = dense_as_general(from, '.', -1, 0), pid); /* in-place */
    SEXP to = PROTECT(GET_SLOT(from, Matrix_xSym)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    setAttrib(to, R_DimSymbol, dim);
    if (!TRIVIAL_DIMNAMES(dimnames))
	setAttrib(to, R_DimNamesSymbol, dimnames);
    UNPROTECT(2);
    return to;
}

/* as(<[CRT]sparseMatrix>, "vector") */
SEXP R_sparse_as_vector(SEXP from)
{
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(from = sparse_as_dense(from, 0), &pid);
    REPROTECT(from = dense_as_general(from, '.', -1, 0), pid); /* in-place */
    from = GET_SLOT(from, Matrix_xSym);
    UNPROTECT(1);
    return from;
}

/* as(<[CRT]sparseMatrix>, "[nlidz]Matrix") */
SEXP R_sparse_as_kind(SEXP from, SEXP kind, SEXP drop0)
{
    char k;
    if ((kind = asChar(kind)) == NA_STRING || (k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_sparse_as_kind()'"));
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_as_kind");
    const char *clf = valid[ivalid]; // clf := class~(from)
    if (k == '.' || k == clf[0])
	return from;
    SEXPTYPE tx = kind2type(k); /* validating 'k' before doing more */

    int do_drop0 =
#if 0 /* MJ: Matrix <= 1.4-1 forces for (only??) d.[CR]Matrix->lMatrix */
	(clf[0] == 'd' && k != 'n' && clf[2] != 'T') ||
#endif
	(clf[0] != 'n' && asLogical(drop0) != 0);
    if (do_drop0)
	PROTECT(from = R_sparse_drop0(from));

    int do_aggr = (clf[2] == 'T' &&
		   (clf[0] == 'n' || clf[0] == 'l') && k != 'n' && k != 'l');
    if (do_aggr)
	PROTECT(from = Tsparse_aggregate(from));
    
    char clt[] = "...Matrix";
    clt[0] = k;
    clt[1] = clf[1];
    clt[2] = clf[2];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	p = R_NilValue, i = R_NilValue;
    
    SET_SLOT(to, Matrix_DimSym, GET_SLOT(from, Matrix_DimSym));
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    if (clf[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (clf[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    if (clf[2] != 'T')
	SET_SLOT(to, Matrix_pSym, p = GET_SLOT(from, Matrix_pSym));
    if (clf[2] != 'R')
	SET_SLOT(to, Matrix_iSym, i = GET_SLOT(from, Matrix_iSym));
    if (clf[2] != 'C')
	SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_jSym));
    if (clf[0] != 'n' && k != 'n') {
	SET_SLOT(to, Matrix_xSym,
		 coerceVector(GET_SLOT(from, Matrix_xSym), tx));
    } else if (clf[0] == 'n') {
	R_xlen_t nx = (clf[2] == 'T')
	    ? XLENGTH(i) : INTEGER(p)[XLENGTH(p) - 1];
	SEXP x = PROTECT(allocVector(tx, nx));

#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_)		\
	do {						\
	    _CTYPE_ *px = _PTR_(x);			\
	    while (nx--)				\
		*(px++) = _ONE_;			\
	} while (0)
	
	SPARSE_CASES(tx, SET1);

#undef SET1
	
	SET_SLOT(to, Matrix_xSym, x);
	UNPROTECT(1);
    }

    UNPROTECT((do_drop0 || do_aggr) ? 2 : 1);
    return to;
}

/* as(<[CRT]sparseMatrix>, "generalMatrix") */
SEXP R_sparse_as_general(SEXP from)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_as_general");
    const char *clf = valid[ivalid];
    if (clf[1] == 'g')
	return from;
    
    char clt[] = ".g.Matrix";
    clt[0] = clf[0];
    clt[2] = clf[2];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    
    SET_SLOT(to, Matrix_DimSym, dim);
    if (clf[1] == 's') {
	set_symmetrized_DimNames(to, dimnames, -1);
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    } else {
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	if (*diag_P(from) == 'N') {
	    if (clf[2] != 'T')
		SET_SLOT(to, Matrix_pSym, GET_SLOT(from, Matrix_pSym));
	    if (clf[2] != 'R')
		SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_iSym));
	    if (clf[2] != 'C')
		SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_jSym));
	    if (clf[0] != 'n')
		SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
	    UNPROTECT(1);
	    return to;
	}
    }

    int n = INTEGER(dim)[0];
    char ul = *uplo_P(from);

    if (clf[2] == 'T') {
	
	SEXP i0 = GET_SLOT(from, Matrix_iSym),
	    j0 = GET_SLOT(from, Matrix_jSym);
	int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
	R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = nnz0;
	
	if (clf[1] == 's') {
	    for (k = 0; k < nnz0; ++k)
		if (pi0[k] == pj0[k])
		    --nnz1;
	    nnz1 += nnz0;
	} else {
	    nnz1 += n;
	}

	SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
	    j1 = PROTECT(allocVector(INTSXP, nnz1));
	int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
	SET_SLOT(to, Matrix_iSym, i1);
	SET_SLOT(to, Matrix_jSym, j1);
	Memcpy(pi1, pi0, nnz0);
	Memcpy(pj1, pj0, nnz0);
	pi1 += nnz0;
	pj1 += nnz0;

	SEXP x0, x1;
	SEXPTYPE tx = NILSXP;
	if (clf[0] != 'n') {
	    x0 = GET_SLOT(from, Matrix_xSym);
	    PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
	    SET_SLOT(to, Matrix_xSym, x1);
	}
	
	if (clf[1] == 's') {

#define SAG_T(_XASSIGN_, _XINCR_)				\
	    do {						\
		for (k = 0; k < nnz0; ++k) {			\
		    if (*pi0 != *pj0) {				\
			*(pi1++) = *pj0;			\
			*(pj1++) = *pi0;			\
			_XASSIGN_; /* *(px1++) = *px0; */	\
		    }						\
		    ++pi0; ++pj0; _XINCR_; /* ++px0; */		\
		}						\
	    } while (0)
	    
#define SAG_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		Memcpy(px1, px0, nnz0);				\
		px1 += nnz0;					\
		SAG_T(*(px1++) = *px0, ++px0);			\
	    } while (0)

	    if (clf[0] == 'n')
		SAG_T(, );
	    else
		SPARSE_CASES(tx, SAG_T_X);

#undef SAG_T_X
#undef SAG_T

	} else {

#define SAG_T(_XASSIGN_ONE_)					\
	    do {						\
		int j;						\
		for (j = 0; j < n; ++j) {			\
		    *(pi1++) = *(pj1++) = j;			\
		    _XASSIGN_ONE_; /* *(px1++) = _ONE_; */	\
		}						\
	    } while (0)

#define SAG_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		Memcpy(px1, px0, nnz0);				\
		px1 += nnz0;					\
		SAG_T(*(px1++) = _ONE_);			\
	    } while (0)
	    
	    if (clf[0] == 'n')
		SAG_T();
	    else
		SPARSE_CASES(tx, SAG_T_X);

#undef SAG_T_X
#undef SAG_T
	    
	}
	
    } else { /* clf[2] != 'T' */

	SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
	    p0 = GET_SLOT(from, Matrix_pSym),
	    p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
	    i0 = GET_SLOT(from, iSym);
	int j, k, kend, nnz1,
	    *pp0 = INTEGER(p0), *pp1 = INTEGER(p1), *pi0 = INTEGER(i0);
	pp0++;
	*(pp1++) = 0;
	
	if (clf[1] == 's') {
	    Memzero(pp1, n);
	    for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if (pi0[k] != j)
			++pp1[pi0[k]];
		    ++k;
		}
	    }
	    for (j = 1; j < n; ++j)
		pp1[j] += pp1[j-1];
	    for (j = 0; j < n; ++j)
		pp1[j] += pp0[j];
	} else {
	    for (j = 0; j < n; ++j)
		pp1[j] = pp0[j] + j + 1;
	}

	SEXP i1 = PROTECT(allocVector(INTSXP, nnz1 = pp1[n-1]));
	int *pi1 = INTEGER(i1);
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to, iSym, i1);

	SEXP x0 = R_NilValue, x1 = R_NilValue;
	SEXPTYPE tx = NILSXP;
	if (clf[0] != 'n') {
	    x0 = GET_SLOT(from, Matrix_xSym);
	    PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
	    SET_SLOT(to, Matrix_xSym, x1);
	}

	if (clf[1] == 's') {

	    int *pp1_;
	    Calloc_or_Alloca_TO(pp1_, n, int);
	    Memcpy(pp1_, pp1 - 1, n);
	    
#define SAG_CR(_XASSIGN_IJ_, _XASSIGN_JI_)				\
	    do {							\
		for (j = 0, k = 0; j < n; ++j) {			\
		    kend = pp0[j];					\
		    while (k < kend) {					\
			pi1[pp1_[j]] = pi0[k];				\
			_XASSIGN_IJ_; /* px1[pp1_[j]] = px0[k]; */	\
			++pp1_[j];					\
			if (pi0[k] != j) {				\
			    pi1[pp1_[pi0[k]]] = j;			\
			    _XASSIGN_JI_; /* px1[pp1_[pi0[k]]] = px0[k]; */ \
			    ++pp1_[pi0[k]];				\
			}						\
			++k;						\
		    }							\
		}							\
	    } while (0)
	    
#define SAG_CR_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		SAG_CR(px1[pp1_[j]] = px0[k],			\
		       px1[pp1_[pi0[k]]] = px0[k]);		\
	    } while (0)

	    if (clf[0] == 'n')
		SAG_CR(, );
	    else
		SPARSE_CASES(tx, SAG_CR_X);

#undef SAG_CR_X
#undef SAG_CR

	    Free_FROM(pp1_, n);
	    
	} else {

#define SAG_CR(_XASSIGN_, _XASSIGN_ONE_)				\
	    do {							\
		if (ul == ((clf[2] == 'C') ? 'U' : 'L')) {		\
		    for (j = 0, k = 0; j < n; ++j) {			\
			kend = pp0[j];					\
			while (k < kend) {				\
			    *(pi1++) = *(pi0++);			\
			    _XASSIGN_; /* *(px1++) = *(px0++); */	\
			    ++k;					\
			}						\
			*(pi1++) = j;					\
			_XASSIGN_ONE_; /* *(px1++) = _ONE_; */		\
		    }							\
		} else {						\
		    for (j = 0, k = 0; j < n; ++j) {			\
			kend = pp0[j];					\
			*(pi1++) = j;					\
			_XASSIGN_ONE_; /* *(px1++) = _ONE_; */		\
			while (k < kend) {				\
			    *(pi1++) = *(pi0++);			\
			    _XASSIGN_; /* *(px1++) = *(px0++); */	\
			    ++k;					\
			}						\
		    }							\
		}							\
	    } while (0)

#define SAG_CR_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		SAG_CR(*(px1++) = *(px0++), *(px1++) = _ONE_);	\
	    } while (0)

	    if (clf[0] == 'n')
		SAG_CR(, );
	    else
		SPARSE_CASES(tx, SAG_CR_X);

#undef SAG_CR_X
#undef SAG_CR

	}
	
    }

    UNPROTECT((clf[0] == 'n') ? 3 : 4);
    return to;
}

/* as(<diagonalMatrix>, "[nlidz][gts][CRT]Matrix") */
SEXP R_diagonal_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP drop0)
{
    const char *zzz;
    char z0, z1, z2, ul = '\0';
    if ((code = asChar(code)) == NA_STRING ||
	(z0 = (zzz = CHAR(code))[0]) == '\0' ||
	(z1 = zzz[1]) == '\0' ||
	(z2 = zzz[2]) == '\0')
	error(_("invalid 'code' to 'R_diagonal_as_sparse()'"));
    if ((z1 == 't' || z1 == 's') &&
	(((uplo = asChar(uplo)) == NA_STRING || (ul = *CHAR(uplo)) == '\0')))
	error(_("invalid 'uplo' to 'R_diagonal_as_sparse()'"));
    
    static const char *valid[] = { VALID_DIAGONAL, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_diagonal_as_sparse");
    const char *clf = valid[ivalid];
    if (z0 == '.')
	z0 = clf[0];
    SEXPTYPE tx = kind2type(z0); /* validating 'z0' before doing more */

    char clt[] = "...Matrix";
    clt[0] = z0;
    clt[1] = z1;
    clt[2] = z2;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	diag = GET_SLOT(from, Matrix_diagSym);
    char di = *CHAR(STRING_ELT(diag, 0));

    SET_SLOT(to, Matrix_DimSym, dim);
    if (z1 == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (z1 != 'g')
	SET_SLOT(to, Matrix_uploSym, mkString((ul == 'U') ? "U" : "L"));
    if (z1 == 't')
	SET_SLOT(to, Matrix_diagSym, diag);
    
    SEXP p = R_NilValue, i = R_NilValue, x = R_NilValue;
    int k, n = INTEGER(dim)[0], nprotect = 1;
    R_xlen_t n1a = (R_xlen_t) n + 1;
    
    if (di != 'N') { /* unit diagonal */
	
	if (z2 != 'T') {
	    PROTECT(p = allocVector(INTSXP, n1a));
	    ++nprotect;
	    int *pp = INTEGER(p);
	    if (z1 == 't')
		Memzero(pp, n1a);
	    else
		for (k = 0; k <= n; ++k)
		    *(pp++) = k;
	}
	if (z1 == 't') {
	    PROTECT(i = allocVector(INTSXP, 0));
	    ++nprotect;
	    if (z0 != 'n') {
		PROTECT(x = allocVector(tx, 0));
		++nprotect;
	    }
	} else {
	    PROTECT(i = allocVector(INTSXP, n));
	    ++nprotect;
	    int *pi = INTEGER(i);
	    if (z0 == 'n') {
		for (k = 0; k < n; ++k)
		    *(pi++) = k;
	    } else {
		PROTECT(x = allocVector(tx, n));
		++nprotect;

#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
		do {				\
		    _CTYPE_ *px = _PTR_(x);	\
		    for (k = 0; k < n; ++k) {	\
			*(pi++) = k;		\
			*(px++) = _ONE_;	\
		    }				\
		} while (0)
		
		SPARSE_CASES(tx, SET1);
		
#undef SET1

	    }
	}
	
    } else { /* non-unit diagonal */
	
	if (TYPEOF(x = GET_SLOT(from, Matrix_xSym)) != tx) {
	    PROTECT(x = coerceVector(x, tx));
	    ++nprotect;
	}
	    
	if (z0 == 'n' || asLogical(drop0) != 0) { /* _do_ drop zeros (if any) */
	    
	    int nnz = 0;
	    
#define DROP0_DIAGONAL(_CTYPE_, _PTR_, _NZ_)				\
	    do {							\
		_CTYPE_ *px = _PTR_(x);					\
		if (z2 == 'T') {					\
		    for (k = 0; k < n; ++k)				\
			if (_NZ_(px[k]))				\
			    ++nnz;					\
		} else {						\
		    PROTECT(p = allocVector(INTSXP, n1a));		\
		    ++nprotect;						\
		    int *pp = INTEGER(p);				\
		    *(pp++) = 0;					\
		    for (k = 0; k < n; ++k)				\
			*(pp++) = (_NZ_(px[k])) ? ++nnz : nnz;		\
		}							\
		PROTECT(i = allocVector(INTSXP, nnz));			\
		++nprotect;						\
		if (nnz != n && z0 != 'n') {				\
		    PROTECT(x = allocVector(tx, nnz));			\
		    ++nprotect;						\
		}							\
		if (nnz == 0)						\
		    continue;						\
		int *pi = INTEGER(i);					\
		if (nnz == n) {						\
		    for (k = 0; k < n; ++k)				\
			*(pi++) = k;					\
		} else if (z0 == 'n') {					\
		    for (k = 0; k < n; ++k)				\
			if (_NZ_(px[k]))				\
			    *(pi++) = k;				\
		} else {						\
		    _CTYPE_ *newpx = _PTR_(x);				\
		    for (k = 0; k < n; ++k) {				\
			if (_NZ_(px[k])) {				\
			    *(pi++) = k;				\
			    *(newpx++) = px[k];				\
			}						\
		    }							\
		}							\
	    } while (0)

	    DROP0_CASES(tx, DROP0_DIAGONAL);

#undef DROP0_DIAGONAL
	    
	} else { /* _do not_ drop zeros */
	    
	    PROTECT(i = allocVector(INTSXP, n));
	    ++nprotect;
	    int *pi = INTEGER(i);
	    if (z2 == 'T') {
		PROTECT(p = allocVector(INTSXP, n1a));
		++nprotect;
		int *pp = INTEGER(p);
		*(pp++) = 0;
		for (k = 0; k < n; ++k)
		    *(pi++) = *(pp++) = k;
	    } else {
		for (k = 0; k < n; ++k)
		    *(pi++) = k;
	    }
	    
	}
    }

    if (z2 != 'T')
	SET_SLOT(to, Matrix_pSym, p);
    if (z2 != 'R')
	SET_SLOT(to, Matrix_iSym, i);
    if (z2 != 'C')
	SET_SLOT(to, Matrix_jSym, i);
    if (z0 != 'n')
	SET_SLOT(to, Matrix_xSym, x);
    UNPROTECT(nprotect);
    return to;
}

/* as(<diagonalMatrix>, ".(ge|tr|sy|tp|sp)Matrix") */
SEXP R_diagonal_as_dense(SEXP from, SEXP code, SEXP uplo)
{
    const char *zzz;
    char z0, z1, z2, ul = '\0';
    if ((code = asChar(code)) == NA_STRING ||
	(z0 = (zzz = CHAR(code))[0]) == '\0' ||
	(z1 = zzz[1]) == '\0' ||
	(z2 = zzz[2]) == '\0')
	error(_("invalid 'code' to 'R_diagonal_as_dense()'"));
    if ((z1 == 't' || z1 == 's') &&
	(((uplo = asChar(uplo)) == NA_STRING || (ul = *CHAR(uplo)) == '\0')))
	error(_("invalid 'uplo' to 'R_diagonal_as_dense()'"));
	
    static const char *valid[] = { VALID_DIAGONAL, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_diagonal_as_dense");
    const char *clf = valid[ivalid];
    if (z0 == '.')
	z0 = clf[0];
    
    SEXPTYPE tx = kind2type(z0); /* validating 'z0' before doing more */
    char clt[] = "...Matrix";
    clt[0] = z0;
    clt[1] = z1;
    clt[2] = z2;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	diag = GET_SLOT(from, Matrix_diagSym);

    char di = *CHAR(STRING_ELT(diag, 0));
    int n = INTEGER(dim)[0];
    double nn = (double) n * n, dnx = (z2 == 'p') ? 0.5 * (nn + n) : nn,
	dsize = dnx * kind2size(z0);
    if (dnx > R_XLEN_T_MAX)
	error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
    if (dsize > 0x1p+30 /* 1 GiB */)
	warning(_("sparse->dense coercion: "
		  "allocating vector of size %0.1f GiB"),
		0x1p-30 * dsize);

    R_xlen_t nx = (R_xlen_t) dnx;
    SEXP x = PROTECT(allocVector(tx, nx));
    
    SET_SLOT(to, Matrix_DimSym, dim);
    if (z1 == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (z1 == 't' || z1 == 's')
	SET_SLOT(to, Matrix_uploSym, mkString((ul == 'U') ? "U" : "L"));
    if (z1 == 't')
	SET_SLOT(to, Matrix_diagSym, diag);
    SET_SLOT(to, Matrix_xSym, x);
    
#define DAD_COPY_DIAGONAL(_PREFIX_, _PTR_)				\
    do {								\
	Memzero(_PTR_(x), nx);						\
	if (z1 != 't' || di == 'N') {					\
	    SEXP y = PROTECT(coerceVector(GET_SLOT(from, Matrix_xSym), tx)); \
	    if (z2 == 'p')						\
		_PREFIX_ ## dense_packed_copy_diagonal(			\
		    _PTR_(x), _PTR_(y), n, n, ul, ul /* unused */, di);	\
	    else							\
		_PREFIX_ ## dense_unpacked_copy_diagonal(		\
		    _PTR_(x), _PTR_(y), n, n,     ul /* unused */, di);	\
	    UNPROTECT(1);						\
	}								\
    } while (0)
	
    switch (tx) {
    case REALSXP:
	DAD_COPY_DIAGONAL(d, REAL);
	break;
    case LGLSXP:
	DAD_COPY_DIAGONAL(i, LOGICAL);
	break;
    case INTSXP:
	DAD_COPY_DIAGONAL(i, INTEGER);
	break;
    case CPLXSXP:
	DAD_COPY_DIAGONAL(z, COMPLEX);
	break;
    default:
	break;
    }

    UNPROTECT(2);
    return to;
}

/* as(<diagonalMatrix>, "[lidz]Matrix") */
SEXP R_diagonal_as_kind(SEXP from, SEXP kind)
{
    char k;
    if ((kind = asChar(kind)) == NA_STRING || (k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_sparse_as_kind()'"));
    if (k == 'n')
	error(_("class ndiMatrix is unimplemented"));
    
    static const char *valid[] = { VALID_DIAGONAL, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_diagonal_as_kind");
    const char *clf = valid[ivalid];
    if (k == '.' || k == clf[0])
	return from;
    
    char clt[] = ".diMatrix";
    clt[0] = k;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
    SET_SLOT(to, Matrix_DimSym, GET_SLOT(from, Matrix_DimSym));
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    SET_SLOT(to, Matrix_xSym,
	     coerceVector(GET_SLOT(from, Matrix_xSym), kind2type(k)));
    UNPROTECT(1);
    return to;
}

/* drop0(<[CRT]sparseMatrix>) 
   TODO: support 'tol' argument, to be interpreted as modulus for zMatrix
*/
SEXP R_sparse_drop0(SEXP from)
{
    static const char *valid[] = { VALID_DSPARSE, VALID_LSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_drop0");
    const char *cl = valid[ivalid];
    
    SEXP p0, x0 = GET_SLOT(from, Matrix_xSym);
    SEXPTYPE tx = TYPEOF(x0);
    int *pp0 = NULL, n = 0;
    R_xlen_t k, kend, nnz_, nnz0, nnz1 = 0, n1a = 0;

    if (cl[2] == 'T') {
	nnz0 = XLENGTH(x0);
    } else {
	p0 = GET_SLOT(from, Matrix_pSym);
	pp0 = INTEGER(p0);
	n1a = XLENGTH(p0);
	n = (int) (n1a - 1);
	nnz0 = pp0[n];
    }
    
#define DROP0_START(_CTYPE_, _PTR_, _NZ_)			\
    do {							\
	_CTYPE_ *px0 = _PTR_(x0);				\
	while (nnz1 < nnz0 && _NZ_(*px0)) {			\
	    ++nnz1;						\
	    ++px0;						\
	}							\
	if (nnz1 == nnz0)					\
	    return from;					\
	nnz_ = nnz1;						\
	for (k = nnz_; k < nnz0; ++k, ++px0)			\
	    if (_NZ_(*px0))					\
		++nnz1;						\
    } while (0)
    
    DROP0_CASES(tx, DROP0_START);

#undef DROP0_START

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));
    SET_SLOT(to, Matrix_DimSym, GET_SLOT(from, Matrix_DimSym));
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    if (cl[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (cl[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    else
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));

    /* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

    SEXP iSym = (cl[2] == 'R') ? Matrix_jSym : Matrix_iSym,
	i0 = GET_SLOT(from, iSym),
	i1 = PROTECT(allocVector(INTSXP, nnz1)),
	x1 = PROTECT(allocVector(tx, nnz1));
    int *pi0 = INTEGER(i0),
	*pi1 = INTEGER(i1);
    SET_SLOT(to, iSym, i1);
    SET_SLOT(to, Matrix_xSym, x1);
    
    if (cl[2] == 'T') {

	SEXP j0 = GET_SLOT(from, Matrix_jSym),
	    j1 = PROTECT(allocVector(INTSXP, nnz1));
	int *pj0 = INTEGER(j0),
	    *pj1 = INTEGER(j1);
	SET_SLOT(to, Matrix_jSym, j1);

#define DROP0_END(_CTYPE_, _PTR_, _NZ_)				\
	do {							\
	    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
	    Memcpy(pi1, pi0, nnz_);				\
	    Memcpy(pj1, pj0, nnz_);				\
	    Memcpy(px1, px0, nnz_);				\
	    for (k = nnz_; k < nnz0; ++k) {			\
		if (_NZ_(px0[k])) {				\
		    pi1[nnz_] = pi0[k];				\
		    pj1[nnz_] = pj0[k];				\
		    px1[nnz_] = px0[k];				\
		    ++nnz_;					\
		}						\
	    }							\
	} while (0)

	DROP0_CASES(tx, DROP0_END);

#undef DROP0_END
	
    } else { /* cl[2] == 'C' || cl[2] == 'R' */

	SEXP p1 = PROTECT(allocVector(INTSXP, n1a));
	int *pp1 = INTEGER(p1), j;
	SET_SLOT(to, Matrix_pSym, p1);

#define DROP0_END(_CTYPE_, _PTR_, _NZ_)			\
	do {						\
	    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	    Memcpy(pi1, pi0, nnz_);			\
	    Memcpy(px1, px0, nnz_);			\
	    j = 0;					\
	    while ((kend = pp0[j]) <= nnz_)		\
		pp1[j++] = kend;			\
	    for (k = nnz_; k < kend; ++k) {		\
		if (_NZ_(px0[k])) {			\
		    pi1[nnz_] = pi0[k];			\
		    px1[nnz_] = px0[k];			\
		    ++nnz_;				\
		}					\
	    }						\
	    pp1[j] = nnz_;				\
	    while (++j <= n) {				\
		kend = pp0[j];				\
		while (k < kend) {			\
		    if (_NZ_(px0[k])) {			\
			pi1[nnz_] = pi0[k];		\
			px1[nnz_] = px0[k];		\
			++nnz_;				\
		    }					\
		    ++k;				\
		}					\
		pp1[j] = nnz_;				\
	    }						\
	} while (0)
    
	DROP0_CASES(tx, DROP0_END);

#undef DROP0_END
	
    }

    UNPROTECT(4);
    return to;
}

/* band(<[CRT]sparseMatrix>, k1, k2), tri[ul](<[CRT]sparseMatrix>, k) */
/* NB: argument validation more or less copied from R_dense_band() */
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_band");
    const char *clf = valid[ivalid];
    
    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], a, b;
    if (isNull(k1))
	a = (m > 0) ? 1 - m : 0;
    else if ((a = asInteger(k1)) == NA_INTEGER || a < -m || a > n)
	error(_("'k1' must be an integer from -Dim[1] to Dim[2]"));
    if (isNull(k2))
	b = (n > 0) ? n - 1 : 0;
    else if ((b = asInteger(k2)) == NA_INTEGER || b < -m || b > n)
	error(_("'k2' must be an integer from -Dim[1] to Dim[2]"));
    else if (b < a)
	error(_("'k1' must be less than or equal to 'k2'"));
    /* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
    if (a <= 1 - m && b >= n - 1 && (clf[1] == 't' || m != n || m > 1 || n > 1))
	return from;

    char ulf = 'U', ult = 'U', di = 'N';
    if (clf[1] != 'g') {
	ulf = *uplo_P(from);
	if (clf[1] == 't') {
	    /* Be fast if band contains entire triangle */
	    if ((ulf == 'U') ? (a <= 0 && b >= n - 1) : (b >= 0 && a <= 1 - m))
		return from;
	    if (a <= 0 && b >= 0)
		di = *diag_P(from);
	}
    }

    int ge = 0, tr = 0, sy = 0, nprotect = 0;
    ge = m != n || (!(tr = a >= 0 || b <= 0 || clf[1] == 't') &&
		    !(sy = a == -b && clf[1] == 's'));   

    char clt[] = "...Matrix";
    clt[0] = clf[0];
    clt[1] = (ge) ? 'g' : ((tr) ? 't' : 's');
    clt[2] = clf[2];

    /* band(<R>, a, b) is equivalent to t(band(t(<R>), -b, -a)) ! */
    
    if (clf[2] == 'R') {
	int tmp;
	tmp = m; m =  n; n =  tmp;
	tmp = a; a = -b; b = -tmp;
	ulf = (ulf == 'U') ? 'L' : 'U';
	clt[2] = 'C';
	PROTECT(from = tCRsparse_as_RCsparse(from));
	++nprotect;
	if (m != n)
	    dim = GET_SLOT(from, Matrix_DimSym);
    }

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    ++nprotect;
    
    SET_SLOT(to, Matrix_DimSym, dim);
    if (!sy && clf[1] == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);

    if (!ge) {
	if (sy || clf[1] == 't') {
	    /* We can tighten the band ... but does it help?? FIXME */
	    if (ulf == 'U') {
		if (b >= 0 && a < 0)
		    a = 0;
	    } else {
		if (a <= 0 && b > 0)
		    b = 0;
	    }
	}
	
	ult = (tr && clf[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ulf;
	if (ult != 'U')
	    SET_SLOT(to, Matrix_uploSym, mkString("L"));
	if (di != 'N')
	    SET_SLOT(to, Matrix_diagSym, mkString("U"));
    }

    /* It remains to set some subset of 'p', 'i', 'j', 'x' ... */
    
    SEXP p0, p1, i0, j0 = R_NilValue;
    int *pp0 = NULL, *pp1 = NULL, *pi0 = NULL, *pj0 = NULL, d, j;
    R_xlen_t k, kend, nnz0, nnz1;
    i0 = GET_SLOT(from, Matrix_iSym);
    pi0 = INTEGER(i0);
    if (clf[2] == 'T') {
	j0 = GET_SLOT(from, Matrix_jSym);
	pj0 = INTEGER(j0);
	nnz0 = XLENGTH(j0);
    } else {
	p0 = GET_SLOT(from, Matrix_pSym);
	pp0 = INTEGER(p0);
	nnz0 = pp0[n];
	pp0++;
    }
    
    /* Counting number of nonzero elements in band ... */
    
    nnz1 = 0;
    if (clf[2] == 'T') {

	if (!sy && clf[1] == 's') {
	    for (k = 0; k < nnz0; ++k) {
		if ((d = pj0[k] - pi0[k]) >= a && d <= b)
		    ++nnz1;
		if (d != 0 && -d >= a && -d <= b)
		    ++nnz1;
	    }
	} else {
	    for (k = 0; k < nnz0; ++k) {
		if ((d = pj0[k] - pi0[k]) >= a && d <= b)
		    ++nnz1;
	    }
	}

    } else {

	PROTECT(p1 = allocVector(INTSXP, (R_xlen_t) n + 1));
	SET_SLOT(to, Matrix_pSym, p1);
	pp1 = INTEGER(p1);
	++nprotect;
	*(pp1++) = 0;
	
	if (!sy && clf[1] == 's') {
	    Memzero(pp1, n);
	    for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if ((d = j - pi0[k]) >= a && d <= b)
			++pp1[j];
		    if (d != 0 && -d >= a && -d <= b)
			++pp1[pi0[k]];
		    ++k;
		}
	    }
	    for (j = 0; j < n; ++j) {
		nnz1 += pp1[j];
		pp1[j] = nnz1;
	    }
	} else {
	    for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if ((d = j - pi0[k]) >= a && d <= b)
			++nnz1;
		    ++k;
		}
		pp1[j] = nnz1;
	    }
	}
	
    }

    if (nnz1 == nnz0 && (sy || clf[1] != 's')) {
	/* No need to allocate in this case: band has all nonzero elements */
	SET_SLOT(to, Matrix_iSym, i0);
	if (clf[2] == 'T')
	    SET_SLOT(to, Matrix_jSym, j0);
	if (clf[0] != 'n')
	    SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));	
	if (clf[2] == 'R')
	    to = tCRsparse_as_RCsparse(to);
	UNPROTECT(nprotect);
	return to;
    }
    
    /* Now allocating and filling out slots ... */

    SEXP i1, j1;
    int *pi1 = NULL, *pj1 = NULL;
    
    PROTECT(i1 = allocVector(INTSXP, nnz1));
    SET_SLOT(to, Matrix_iSym, i1);
    pi1 = INTEGER(i1);
    ++nprotect;
    if (clf[2] == 'T') {
	PROTECT(j1 = allocVector(INTSXP, nnz1));
	SET_SLOT(to, Matrix_jSym, j1);
	pj1 = INTEGER(j1);
	++nprotect;
    }

#define SPARSE_BAND(_XASSIGN_, _XASSIGN_IJ_, _XASSIGN_JI_)		\
    do {								\
	if (clf[2] == 'T') {						\
	    if (!sy && clf[1] == 's') {					\
		for (k = 0; k < nnz0; ++k) {				\
		    if ((d = pj0[k] - pi0[k]) >= a && d <= b) {		\
			*(pi1++) = pi0[k];				\
			*(pj1++) = pj0[k];				\
			_XASSIGN_; /* *(px1++) = px0[k]; */		\
		    }							\
		    if (d != 0 && -d >= a && -d <= b) {			\
			*(pi1++) = pj0[k];				\
			*(pj1++) = pi0[k];				\
			_XASSIGN_; /* *(px1++) = px0[k]; */		\
		    }							\
		}							\
	    } else {							\
		for (k = 0; k < nnz0; ++k) {				\
		    if ((d = pj0[k] - pi0[k]) >= a && d <= b) {		\
			*(pi1++) = pi0[k];				\
			*(pj1++) = pj0[k];				\
			_XASSIGN_; /* *(px1++) = px0[k]; */		\
		    }							\
		}							\
	    }								\
	} else {							\
	    if (!sy && clf[1] == 's') {					\
		int *pp1_;						\
		Calloc_or_Alloca_TO(pp1_, n, int);			\
		Memcpy(pp1_, pp1 - 1, n);				\
		for (j = 0, k = 0; j < n; ++j) {			\
		    kend = pp0[j];					\
		    while (k < kend) {					\
			if ((d = j - pi0[k]) >= a && d <= b) {		\
			    pi1[pp1_[j]] = pi0[k];			\
			    _XASSIGN_IJ_; /* px1[pp1_[j]] = px0[k]; */	\
			    ++pp1_[j];					\
			}						\
			if (d != 0 && -d >= a && -d <= b) {		\
			    pi1[pp1_[pi0[k]]] = j;			\
			    _XASSIGN_JI_; /* px1[pp1_[pi0[k]]] = px0[k]; */ \
			    ++pp1_[pi0[k]];				\
			}						\
			++k;						\
		    }							\
		}							\
		Free_FROM(pp1_, n);					\
	    } else {							\
		for (j = 0, k = 0; j < n; ++j) {			\
		    kend = pp0[j];					\
		    while (k < kend) {					\
			if ((d = j - pi0[k]) >= a && d <= b) {		\
			    *(pi1++) = pi0[k];				\
			    _XASSIGN_; /* *(px1++) = px0[k]; */		\
			}						\
			++k;						\
		    }							\
		}							\
	    }								\
	}								\
    } while (0)

#define SPARSE_BAND_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	SPARSE_BAND(*(px1++) = px0[k],			\
		    px1[pp1_[j]] = px0[k],		\
		    px1[pp1_[pi0[k]]] = px0[k]);	\
    } while (0)

    if (clf[0] == 'n') {
	SPARSE_BAND(, , );
    } else {
	SEXPTYPE tx;
	SEXP x0 = GET_SLOT(from, Matrix_xSym),
	    x1 = PROTECT(allocVector(tx = TYPEOF(x0), nnz1));
	SPARSE_CASES(tx, SPARSE_BAND_X);
	SET_SLOT(to, Matrix_xSym, x1);
	UNPROTECT(1);
    }

#undef SPARSE_BAND_X
#undef SPARSE_BAND
    
    if (clf[2] == 'R')
	to = tCRsparse_as_RCsparse(to);
    UNPROTECT(nprotect);
    return to;
}

/* diag(<[CRT]sparseMatrix>, names) */
SEXP R_sparse_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL)
	error(_("'names' must be TRUE or FALSE"));

    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "R_sparse_diag_get");
    const char *cl = valid[ivalid];

    char kind = cl[0];
    SEXPTYPE type = kind2type(kind);
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
    SEXP res = PROTECT(allocVector(type, r));

    if (cl[1] == 't' && *diag_P(obj) != 'N') {

	/* .t[CRT]Matrix with unit diagonal */

	int i;
	
#define DO_ONES(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
	do {					\
	    _CTYPE_ *pres = _PTR_(res);		\
	    for (i = 0; i < r; ++i)		\
		*(pres++) = _ONE_;		\
	} while (0)

	SPARSE_CASES(type, DO_ONES);

#undef DO_ONES
	
    } else if (cl[2] == 'T') {

	/* ..TMatrix with non-unit diagonal */
	
	SEXP x = (kind != 'n') ? GET_SLOT(obj, Matrix_xSym) : R_NilValue,
	    i = GET_SLOT(obj, Matrix_iSym),
	    j = GET_SLOT(obj, Matrix_jSym);
	int *pi = INTEGER(i), *pj = INTEGER(j);
	R_xlen_t k, nnz = XLENGTH(i);

	switch (kind) {
	case 'd':
	{
	    double *px = REAL(x), *pres = REAL(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px)
		if (*pi == *pj)
		    pres[*pi] += *px;
	    break;
	}
	case 'l':
	{
	    int *px = LOGICAL(x), *pres = LOGICAL(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px) {
		if (*pi == *pj && *px != 0) {
		    if (*px != NA_LOGICAL)
			pres[*pi] = 1;
		    else if (pres[*pi] == 0)
			pres[*pi] = NA_LOGICAL;
		}
	    }
	    break;
	}
	case 'n':
	{
	    int *pres = LOGICAL(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj)
		if (*pi == *pj)
		    pres[*pi] = 1;
	    break;
	}
	case 'i':
	{
	    /* FIXME: not detecting integer overflow here */
	    int *px = INTEGER(x), *pres = INTEGER(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px)
		if (*pi == *pj)
		    pres[*pi] += *px;
	    break;
	}
	case 'z':
	{
	    Rcomplex *px = COMPLEX(x), *pres = COMPLEX(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px) {
		if (*pi == *pj) {
		    pres[*pi].r += (*px).r;
		    pres[*pi].i += (*px).i;
		}
	    }
	    break;
	}
	default:
	    break;
	}
	
    } else {

	/* ..[CR]Matrix with non-unit diagonal */
	
	SEXP x = (kind != 'n') ? GET_SLOT(obj, Matrix_xSym) : R_NilValue,
	    iSym = (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym;
	int j, k, kend,
	    *pp = INTEGER(GET_SLOT(obj, Matrix_pSym)) + 1,
	    *pi = INTEGER(GET_SLOT(obj, iSym));
	char ul = (cl[1] != 'g') ? *uplo_P(obj) : '\0';
	
#define DO_DIAG(_RES_, _VAL_GENERAL_, _VAL_TRAILING_, _VAL_LEADING_, _ZERO_) \
	do {								\
	    if (ul == '\0') {						\
		/* .g[CR]Matrix */					\
	    	for (j = 0, k = 0; j < r; ++j) {			\
		    _RES_[j] = _ZERO_;					\
		    kend = pp[j];					\
		    while (k < kend) {					\
			if (pi[k] == j) {				\
			    _RES_[j] = _VAL_GENERAL_ /* px[k] */;	\
			    k = kend;					\
			    break;					\
			}						\
			++k;						\
		    }							\
		}							\
	    } else if (ul == ((cl[2] == 'C') ? 'U' : 'L')) {		\
		/* .[ts][CR]Matrix with "trailing" diagonal */		\
	    	for (j = 0, k = 0; j < r; ++j) {			\
		    kend = pp[j];					\
		    _RES_[j] = (kend - k > 0 && pi[kend-1] == j		\
				? _VAL_TRAILING_ /* px[kend-1] */	\
				: _ZERO_);				\
		    k = kend;						\
		}							\
	    } else {							\
		/* .[ts][CR]Matrix with "leading" diagonal */		\
		for (j = 0, k = 0; j < r; ++j) {			\
		    kend = pp[j];					\
		    _RES_[j] = (kend - k > 0 && pi[k] == j		\
				? _VAL_LEADING_	/* px[k] */		\
				: _ZERO_);				\
		    k = kend;						\
		}							\
	    }								\
	} while (0)
	
#define DO_DIAG_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)		\
	do {							\
	    _CTYPE_ *px = _PTR_(x), *pres = _PTR_(res);		\
	    DO_DIAG(pres, px[k], px[kend-1], px[k], _ZERO_);	\
	} while (0)
	    
	if (kind == 'n') {
	    int *pres = LOGICAL(res);
	    DO_DIAG(pres, 1, 1, 1, 0);
	} else {
	    SPARSE_CASES(type, DO_DIAG_X);
	}

#undef DO_DIAG_X
#undef DO_DIAG
	
    }

    if (do_nms) {
	/* NB: The logic here must be adjusted once the validity method 
	       for 'symmetricMatrix' enforces symmetric 'Dimnames' */
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym),
	    rn = VECTOR_ELT(dn, 0),
	    cn = VECTOR_ELT(dn, 1);
	if (isNull(cn)) {
	    if (cl[1] == 's' && !isNull(rn))
		setAttrib(res, R_NamesSymbol, rn);
	} else {
	    if (cl[1] == 's')
		setAttrib(res, R_NamesSymbol, cn);
	    else if (!isNull(rn) &&
		     (rn == cn || equal_string_vectors(rn, cn, r)))
		setAttrib(res, R_NamesSymbol, (r == m) ? rn : cn);
	}
    }

    UNPROTECT(1);
    return res;
}

/* t(<[CRT]sparseMatrix>) */
SEXP R_sparse_transpose(SEXP from)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_transpose");
    const char *cl = valid[ivalid];

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

    if (m == n) {
	SET_SLOT(to, Matrix_DimSym, dim);
    } else {
	pdim = INTEGER(GET_SLOT(to, Matrix_DimSym));
	pdim[0] = n;
	pdim[1] = m;
    }
    if (cl[1] == 's')
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    else
	set_reversed_DimNames(to, dimnames);
    if (cl[1] != 'g') {
	SET_SLOT(to, Matrix_uploSym,
		 mkString((*uplo_P(from) == 'U') ? "L" : "U"));
	if (cl[1] == 't')
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
	else
	    SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    }

    /* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

    if (cl[2] == 'T') {
	/* No need to allocate in this case: need only reverse 'i' and 'j' */
	SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_jSym));
	SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_iSym));
	if (cl[0] != 'n')
	    SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
	UNPROTECT(1);
	return to;
    }

    /* Now dealing only with [CR]sparseMatrix ... */

    int m_, n_;
    SEXP iSym;
    if (cl[2] == 'C') {
	m_ = m;
	n_ = n;
	iSym = Matrix_iSym;
    } else {
	m_ = n;
	n_ = m;
	iSym = Matrix_jSym;
    }
    
    R_xlen_t m1a = (R_xlen_t) m_ + 1;
    SEXP p0 = GET_SLOT(from, Matrix_pSym),
	p1 = PROTECT(allocVector(INTSXP, m1a));
    int i, j, k, kend,
	*pp0 = INTEGER(p0),
	*pp1 = INTEGER(p1),
	nnz = pp0[n_];
    SEXP i0 = GET_SLOT(from, iSym),
	i1 = PROTECT(allocVector(INTSXP, nnz));
    int *pi0 = INTEGER(i0),
	*pi1 = INTEGER(i1);
    pp0++;
    
    SET_SLOT(to, Matrix_pSym, p1);
    SET_SLOT(to, iSym, i1);
    
    /* Counting number of nonzero elements, by "row" */
    Memzero(pp1, m1a);
    ++pp1;
    for (k = 0; k < nnz; ++k)
	++pp1[pi0[k]];

    /* Computing cumulative sum, in place */
    for (i = 1; i < m_; ++i)
	pp1[i] += pp1[i-1];

    /* Allocating work space */
    int *pp1_;
    Calloc_or_Alloca_TO(pp1_, m_, int);
    Memcpy(pp1_, pp1 - 1, m_);

#define SPARSE_T(_XASSIGN_)				\
    do {						\
	for (j = 0, k = 0; j < n_; ++j) {		\
	    kend = pp0[j];				\
	    while (k < kend) {				\
		i = pi0[k];				\
		pi1[pp1_[i]] = j;			\
		_XASSIGN_; /* px1[pp1_[i]] = px0[k] */	\
		++pp1_[i];				\
		++k;					\
	    }						\
	}						\
    } while (0)

#define SPARSE_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	SPARSE_T(px1[pp1_[i]] = px0[k]);		\
    } while (0);

    if (cl[0] == 'n') {
	SPARSE_T();
    } else {
	SEXPTYPE tx;
	SEXP x0 = GET_SLOT(from, Matrix_xSym),
	    x1 = PROTECT(allocVector(tx = TYPEOF(x0), nnz));
	SET_SLOT(to, Matrix_xSym, x1);
	SPARSE_CASES(tx, SPARSE_T_X);
	UNPROTECT(1);
    }
			 
#undef SPARSE_T_X
#undef SPARSE_T

    Free_FROM(pp1_, m_);
    UNPROTECT(3);
    return to;
}

/* forceSymmetric(<[CRT]sparseMatrix>, uplo) */
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo_to)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_force_symmetric");
    const char *clf = valid[ivalid];
    
    SEXP uplo_from;
    char ulf = 'U', ult = *CHAR(asChar(uplo_to)), di = 'N';
    if (clf[1] != 'g') {
	uplo_from = GET_SLOT(from, Matrix_uploSym);
	ulf = *CHAR(STRING_ELT(uplo_from, 0));
    }
    if (ult == '\0') /* to handle missing(uplo) */
	ult = ulf;
    if (clf[1] == 's') {
	/* .s[CRT]Matrix */
	if (ulf == ult)
	    return from;
	SEXP to = PROTECT(R_sparse_transpose(from));
	if (clf[0] == 'z') {
	    /* Need _conjugate_ transpose */
	    SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	    conjugate(x);
	    UNPROTECT(1);
	}
	UNPROTECT(1);
	return to;
    }
    if (clf[1] == 't')
	di = *diag_P(from);
    
    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to symmetrize a non-square matrix"));

    /* Now handling just square .[gt][CRT]Matrix ... */

    char clt[] = ".s.Matrix";
    clt[0] = clf[0];
    clt[2] = clf[2];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);

    SET_SLOT(to, Matrix_DimSym, dim);
    set_symmetrized_DimNames(to, dimnames, -1);
    SET_SLOT(to, Matrix_uploSym, mkString((ult == 'U') ? "U" : "L"));

    /* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

    if (clf[1] == 't' && di == 'N' && ulf == ult) {
	
	/* No need to allocate in this case: we have the triangle we want */
	if (clf[2] != 'T')
	    SET_SLOT(to, Matrix_pSym, GET_SLOT(from, Matrix_pSym));
	if (clf[2] != 'R')
	    SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_iSym));
	if (clf[2] != 'C')
	    SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_jSym));
	if (clf[0] != 'n')
	    SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
	UNPROTECT(1);
	return to;
	
    } else if (clf[2] == 'T') {

	/* Symmetrizing square .[gt]TMatrix ... */

	SEXP i0 = GET_SLOT(from, Matrix_iSym),
	    j0 = GET_SLOT(from, Matrix_jSym);
	int *pi0 = INTEGER(i0),
	    *pj0 = INTEGER(j0);
	R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = 0;
	
	/* Counting number of nonzero elements in triangle ... */

	if (clf[1] == 't' && di != 'N') {
	    nnz1 = (ulf == ult) ? n + nnz0 : n;
	} else {
	    if (ult == 'U') {
		for (k = 0; k < nnz0; ++k)
		    if (pi0[k] <= pj0[k])
			++nnz1;
	    } else {
		for (k = 0; k < nnz0; ++k)
		    if (pi0[k] >= pj0[k])
			++nnz1;
	    }
	}

	/* Now allocating and filling out slots ... */

	SEXP x0 = R_NilValue, x1 = R_NilValue,
	    i1 = PROTECT(allocVector(INTSXP, nnz1)),
	    j1 = PROTECT(allocVector(INTSXP, nnz1));
	SEXPTYPE tx = NILSXP;
	int *pi1 = INTEGER(i1),
	    *pj1 = INTEGER(j1);
	SET_SLOT(to, Matrix_iSym, i1);
	SET_SLOT(to, Matrix_jSym, j1);
	if (clf[0] != 'n') {
	    x0 = GET_SLOT(from, Matrix_xSym);
	    PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
	    SET_SLOT(to, Matrix_xSym, x1);
	}

	if (clf[1] == 't' && di != 'N') {
	    if (ulf == ult) {
		Memcpy(pi1, pi0, nnz0);
		Memcpy(pj1, pj0, nnz0);
		pi1 += nnz0;
		pj1 += nnz0;
	    }

#define SPARSE_FS(_XASSIGN_)				\
	    do {					\
		int j;					\
		for (j = 0; j < n; ++j) {		\
		    *(pi1++) = *(pj1++) = j;		\
		    _XASSIGN_; /* *(px1++) = _ONE_; */	\
		}					\
	    } while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {							\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
		if (ulf == ult) {					\
		    Memcpy(px1, px0, nnz0);				\
		    px1 += nnz0;					\
		}							\
		SPARSE_FS(*(px1++) = _ONE_);				\
	    } while (0)

	    if (clf[0] == 'n')
		SPARSE_FS();
	    else
		SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS
		
	} else {

#define SPARSE_FS(_XASSIGN_)					\
	    do {						\
		if (ult == 'U') {				\
		    for (k = 0; k < nnz0; ++k) {		\
			if (pi0[k] <= pj0[k]) {			\
			    *(pi1++) = pi0[k];			\
			    *(pj1++) = pj0[k];			\
			    _XASSIGN_; /* *(px1++) = px0[k]; */	\
			}					\
		    }						\
		} else {					\
		    for (k = 0; k < nnz0; ++k) {		\
			if (pi0[k] <= pj0[k]) {			\
			    *(pi1++) = pi0[k];			\
			    *(pj1++) = pj0[k];			\
			    _XASSIGN_; /* *(px1++) = px0[k]; */	\
			}					\
		    }						\
		}						\
	    } while (0)

#define SPARSE_FS_X_BASIC(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		SPARSE_FS(*(px1++) = px0[k]);			\
	    } while (0)

	    if (clf[0] == 'n')
		SPARSE_FS();
	    else
		SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS
	    
	}
	    
    } else {

	/* Symmetrizing square .[gt][CR]Matrix ... */

	SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
	    p0 = GET_SLOT(from, Matrix_pSym),
	    p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
	    i0 = GET_SLOT(from, iSym);
	int j, k, kend,
	    *pp0 = INTEGER(p0),
	    *pp1 = INTEGER(p1),
	    *pi0 = INTEGER(i0),
	    nnz0 = pp0[n],
	    nnz1 = 0;
	pp0++;
	*(pp1++) = 0;
	
	/* Counting number of nonzero elements in triangle, by "column" ... */

	if (clf[1] == 't') {
	    if (di != 'N') {
		/* Have triangular matrix with unit diagonal */
		if (ulf != ult) {
		    /* Returning identity matrix */
		    for (j = 0; j < n; ++j)
			pp1[j] = ++nnz1;
		} else {
		    /* Returning symmetric matrix with unit diagonal */
		    for (j = 0; j < n; ++j)
			pp1[j] = ++nnz1 + pp0[j];
		    nnz1 += nnz0;
		}
	    } else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
		/* Have triangular matrix with non-unit "trailing" diagonal
		   and returning diagonal part */
		for (j = 0; j < n; ++j) {
		    if (pp0[j-1] < pp0[j] && pi0[pp0[j]-1] == j)
			++nnz1;
		    pp1[j] = nnz1;
		}
	    } else {
		/* Have triangular matrix with non-unit "leading" diagonal
		   and returning diagonal part */
		for (j = 0; j < n; ++j) {
		    if (pp0[j-1] < pp0[j] && pi0[pp0[j-1]] == j)
			++nnz1;
		    pp1[j] = nnz1;
		}
	    }
	} else if (ult == ((clf[2] == 'C') ? 'U' : 'L')) {
	    /* Have general matrix and returning upper triangle */
	    for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if (pi0[k] <= j)
			++nnz1;
		    ++k;
		}
		pp1[j] = nnz1;
	    }
	} else {
	    /* Have general matrix and returning lower triangle */
	    for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if (pi0[k] >= j)
			++nnz1;
		    ++k;
		}
		pp1[j] = nnz1;
	    }
	}

	/* Now allocating and filling out slots ... */

	SEXP x0 = R_NilValue, x1 = R_NilValue,
	    i1 = PROTECT(allocVector(INTSXP, nnz1));
	SEXPTYPE tx = NILSXP;
	int *pi1 = INTEGER(i1);
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to, iSym, i1);
	if (clf[0] != 'n') {
	    x0 = GET_SLOT(from, Matrix_xSym);
	    PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
	    SET_SLOT(to, Matrix_xSym, x1);
	}
	
	if (clf[1] == 't') {
	    if (di != 'N') {
		/* Have triangular matrix with unit diagonal */
		if (ulf != ult) {
		    /* Returning identity matrix */

#define SPARSE_FS(_XASSIGN_)					\
		    do {					\
			for (j = 0; j < n; ++j) {		\
			    *(pi1++) = j;			\
			    _XASSIGN_; /* *(px1++) = _ONE_; */	\
			}					\
		    } while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
		    do {				\
			_CTYPE_ *px1 = _PTR_(x1);	\
			SPARSE_FS(*(px1++) = _ONE_);	\
		    } while (0)

		    if (clf[0] == 'n')
			SPARSE_FS();
		    else
			SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS
		    
		} else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
		    /* Returning symmetric matrix
		       with unit "trailing" diagonal */

#define SPARSE_FS(_XASSIGN_, _XASSIGN_ONE_)				\
		    do {						\
			for (j = 0, k = 0; j < n; ++j) {		\
			    kend = pp0[j];				\
			    while (k < kend) {				\
				*(pi1++) = pi0[k];			\
				_XASSIGN_; /* *(px1++) = px0[k]; */	\
				++k;					\
			    }						\
			    *(pi1++) = j;				\
			    _XASSIGN_ONE_; /* *(px1++) = _ONE_; */	\
			}						\
		    } while (0)
		    
#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
		    do {						\
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
			SPARSE_FS(*(px1++) = px0[k], *(px1++) = _ONE_);	\
		    } while (0)

		    if (clf[0] == 'n')
			SPARSE_FS(, );
		    else
			SPARSE_CASES(tx, SPARSE_FS_X);
		    
#undef SPARSE_FS
		    
		} else {
		    /* Returning symmetric matrix
		       with unit "leading" diagonal */

#define SPARSE_FS(_XASSIGN_, _XASSIGN_ONE_)				\
		    do {						\
			for (j = 0, k = 0; j < n; ++j) {		\
			    *(pi1++) = j;				\
			    _XASSIGN_ONE_; /* *(px1++) = _ONE_; */	\
			    kend = pp0[j];				\
			    while (k < kend) {				\
				*(pi1++) = pi0[k];			\
				_XASSIGN_; /* *(px1++) = px0[k]; */	\
				++k;					\
			    }						\
			}						\
		    } while (0)
		    		    
		    if (clf[0] == 'n')
			SPARSE_FS(, );
		    else
			SPARSE_CASES(tx, SPARSE_FS_X);
		    
#undef SPARSE_FS_X
#undef SPARSE_FS
		    
		}
	    } else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
		/* Have triangular matrix with non-unit "trailing" diagonal
		   and returning diagonal part */

#define SPARSE_FS(_XASSIGN_)						\
		do {							\
		    for (j = 0; j < n; ++j) {				\
			if (pp0[j-1] < pp0[j] && pi0[pp0[j]-1] == j) {	\
			    *(pi1++) = j;				\
			    _XASSIGN_; /* *(px1++) = px0[pp0[j]-1]; */	\
			}						\
		    }							\
		} while (0)
		
#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
		do {							\
		    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
		    SPARSE_FS(*(px1++) = px0[pp0[j]-1]);		\
		} while (0)

		if (clf[0] == 'n')
		    SPARSE_FS();
		else
		    SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS
		
	    } else {
		/* Have triangular matrix with non-unit "leading" diagonal 
		   and returning diagonal part */

#define SPARSE_FS(_XASSIGN_)						\
		do {							\
		    for (j = 0; j < n; ++j) {				\
			if (pp0[j-1] < pp0[j] && pi0[pp0[j-1]] == j) {	\
			    *(pi1++) = j;				\
			    _XASSIGN_; /* *(px1++) = px0[pp0[j-1]]; */	\
			}						\
		    }							\
		} while (0)
		
#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
		do {							\
		    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
		    SPARSE_FS(*(px1++) = px0[pp0[j-1]]);		\
		} while (0)

		if (clf[0] == 'n')
		    SPARSE_FS();
		else
		    SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS
		
	    }
	} else if (ult == ((clf[2] == 'C') ? 'U' : 'L')) {
	    /* Have general matrix and returning upper triangle */

#define SPARSE_FS(_XASSIGN_)					\
	    do {						\
		for (j = 0, k = 0; j < n; ++j) {		\
		    kend = pp0[j];				\
		    while (k < kend) {				\
			if (pi0[k] <= j) {			\
			    *(pi1++) = pi0[k];			\
			    _XASSIGN_; /* *(px1++) = px0[k]; */	\
			}					\
			++k;					\
		    }						\
		}						\
	    } while (0)
	    
	    if (clf[0] == 'n')
		SPARSE_FS();
	    else
		SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS
	    
	} else {
	    /* Have general matrix and returning lower triangle */

#define SPARSE_FS(_XASSIGN_)					\
	    do {						\
		for (j = 0, k = 0; j < n; ++j) {		\
		    kend = pp0[j];				\
		    while (k < kend) {				\
			if (pi0[k] >= j) {			\
			    *(pi1++) = pi0[k];			\
			    _XASSIGN_; /* *(px1++) = px0[k]; */	\
			}					\
			++k;					\
		    }						\
		}						\
	    } while (0)
	    
	    if (clf[0] == 'n')
		SPARSE_FS();
	    else
		SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS_X_BASIC
#undef SPARSE_FS

	}
    }
    
    UNPROTECT((clf[0] == 'n') ? 3 : 4);
    return to;
}

/* symmpart(<TsparseMatrix>) */
SEXP Tsparse_symmpart(SEXP from)
{
    static const char *valid[] = { VALID_TSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "Tsparse_symmpart");
    const char *clf = valid[ivalid];
    if (clf[1] == 's' && clf[0] == 'd')
	return from;

    char clt[] = ".sTMatrix";
    clt[0] = (clf[0] != 'z') ? 'd' : 'z';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = (clf[1] != 'g') ? GET_SLOT(from, Matrix_uploSym) : R_NilValue;
    
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to get symmetric part of non-square matrix"));

    int diagU = clf[1] == 't' && *diag_P(from) != 'N';
    if (diagU)
	PROTECT(from = R_sparse_as_general(from)); /* U->N */
    
    SEXP i0 = GET_SLOT(from, Matrix_iSym),
	j0 = GET_SLOT(from, Matrix_jSym),
	x0 = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue,
	x1 = R_NilValue;
    R_xlen_t k, nnz = XLENGTH(i0);
    int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);

#define SPARSE_SYMMPART_CASES(_SYMMPART_)				\
    do {								\
	switch (clf[0]) {						\
	case 'n':							\
	{								\
	    double *px1 = REAL(x1);					\
	    _SYMMPART_(px1[k] = 1.0, px1[k] *= 0.5);			\
	    break;							\
	}								\
	case 'l':							\
	{								\
	    int *px0 = LOGICAL(x0);					\
	    double *px1 = REAL(x1);					\
	    _SYMMPART_(							\
		px1[k]  = (px0[k] == NA_LOGICAL				\
			   ? NA_REAL : ((px0[k] != 0) ? 1.0 : 0.0)),	\
		px1[k] *= 0.5);						\
	    break;							\
	}								\
	case 'i':							\
	{								\
	    int *px0 = INTEGER(x0);					\
	    double *px1 = REAL(x1);					\
	    _SYMMPART_(							\
		px1[k]  = (px0[k] == NA_INTEGER				\
			   ? NA_REAL : (double) px0[k]),		\
		px1[k] *= 0.5);						\
	    break;							\
	}								\
	case 'd':							\
	{								\
	    double *px0 = REAL(x0), *px1 = REAL(x1);			\
	    _SYMMPART_(px1[k] = px0[k], px1[k] *= 0.5);			\
	    break;							\
	}								\
	case 'z':							\
	{								\
	    Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);		\
	    _SYMMPART_(							\
		px1[k] = px0[k],					\
		do { px1[k].r *= 0.5; px1[k].i *= 0.5; } while (0));	\
	    break;							\
	}								\
	default:							\
	    break;							\
	}								\
    } while (0)

    if (clf[1] == 'g') {

	SEXP i1 = PROTECT(allocVector(INTSXP, nnz)),
	    j1 = PROTECT(allocVector(INTSXP, nnz));
	int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
	PROTECT(x1 = allocVector(kind2type(clt[0]), nnz));
	
#define TSPARSE_SYMMPART_GENERAL(_DO_ASSIGN_, _DO_HALF_)		\
	do {								\
	    for (k = 0; k < nnz; ++k) {					\
		if (pi0[k] <= pj0[k]) {					\
		    pi1[k] = pi0[k];					\
		    pj1[k] = pj0[k];					\
		} else {						\
		    pi1[k] = pj0[k];					\
		    pj1[k] = pi0[k];					\
		}							\
		_DO_ASSIGN_;						\
		if (pi0[k] != pj0[k])					\
		    _DO_HALF_;						\
	    }								\
	} while (0)
	
	SPARSE_SYMMPART_CASES(TSPARSE_SYMMPART_GENERAL);

#undef TSPARSE_SYMMPART_GENERAL
	
	SET_SLOT(to, Matrix_DimSym, dim);
	set_symmetrized_DimNames(to, dimnames, -1);
	SET_SLOT(to, Matrix_iSym, i1);
	SET_SLOT(to, Matrix_jSym, j1);
	SET_SLOT(to, Matrix_xSym, x1);
	UNPROTECT(4);
	
    } else if (clf[1] == 't') {

	PROTECT(x1 = (diagU && clf[0] == clt[0]
		      ? x0 : allocVector(kind2type(clt[0]), nnz)));
	
#define TSPARSE_SYMMPART_TRIANGULAR(_DO_ASSIGN_, _DO_HALF_)		\
	do {								\
	    for (k = 0; k < nnz; ++k) {					\
		_DO_ASSIGN_;						\
		if (pi0[k] != pj0[k])					\
		    _DO_HALF_;						\
	    }								\
	} while (0)

	SPARSE_SYMMPART_CASES(TSPARSE_SYMMPART_TRIANGULAR);

#undef TSPARSE_SYMMPART_TRIANGULAR
	
	SET_SLOT(to, Matrix_DimSym, dim);
	set_symmetrized_DimNames(to, dimnames, -1);
	SET_SLOT(to, Matrix_uploSym, uplo);
	SET_SLOT(to, Matrix_iSym, i0);
	SET_SLOT(to, Matrix_jSym, j0);
	SET_SLOT(to, Matrix_xSym, x1);
	UNPROTECT((diagU) ? 3 : 2);
	
    } else { /* clf[1] == 's' */
	
#define SPARSE_SYMMPART_CASES_TRIVIAL				\
	do {							\
	    switch (clf[0]) {					\
	    case 'n':						\
	    {							\
		PROTECT(x1 = allocVector(REALSXP, nnz));	\
		double *px1 = REAL(x1);				\
		while (nnz--)					\
		    *(px1++) = 1.0;				\
		break;						\
	    }							\
	    case 'l':						\
	    case 'i':						\
	    case 'd':						\
		PROTECT(x1 = coerceVector(x0, REALSXP));	\
		break;						\
	    case 'z':						\
		PROTECT(x1 = duplicate(x0));			\
		zeroIm(x1);					\
		break;						\
	    default:						\
		break;						\
	    }							\
	} while (0)
	
	SPARSE_SYMMPART_CASES_TRIVIAL;
	SET_SLOT(to, Matrix_DimSym, dim);
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	SET_SLOT(to, Matrix_uploSym, uplo);
	SET_SLOT(to, Matrix_iSym, i0);
	SET_SLOT(to, Matrix_jSym, j0);
	SET_SLOT(to, Matrix_xSym, x1);
	UNPROTECT(2);
	
    }

    return to;
}

/* symmpart(<[CR]sparseMatrix>) */
SEXP CRsparse_symmpart(SEXP from)
{
    static const char *valid[] = { VALID_CRSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "CRsparse_symmpart");
    const char *clf = valid[ivalid];
    if (clf[1] == 's' && clf[0] == 'd')
	return from;

    char clt[] = ".s.Matrix";
    clt[0] = (clf[0] != 'z') ? 'd' : 'z';
    clt[2] = clf[2];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = (clf[1] != 'g') ? GET_SLOT(from, Matrix_uploSym) : R_NilValue;

    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to get symmetric part of non-square matrix"));

    int diagU = clf[1] == 't' && *diag_P(from) != 'N';
    if (diagU)
	PROTECT(from = R_sparse_as_general(from)); /* U->N */
    
    SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
	p0 = GET_SLOT(from, Matrix_pSym),
	i0 = GET_SLOT(from, iSym),
	x0 = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue,
	x1 = R_NilValue;
    int j, k, kend, *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), nnz = pp0[n];
    pp0++;
    
    if (clf[1] == 'g') {

	PROTECT(from = R_sparse_transpose(from));
	SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
	    i1,
	    p0_ = GET_SLOT(from, Matrix_pSym),
	    i0_ = GET_SLOT(from, iSym),
	    x0_ = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue;
	int k_, kend_,
	    *pp1 = INTEGER(p1),
	    *pi1,
	    *pp0_ = INTEGER(p0_) + 1,
	    *pi0_ = INTEGER(i0_);
	*(pp1++) = 0;

	/* Counting number of nonzero elements in each "column" of result ... */
	
	for (j = 0, k = 0, k_ = 0; j < n; ++j) {
	    pp1[j] = pp1[j-1];
	    kend = pp0[j];
	    kend_ = pp0_[j];
	    while (k < kend) {
		if (pi0[k] > j) {
		    k = kend;
		    break;
		}
		while (k_ < kend_ && pi0_[k_] < pi0[k]) {
		    ++pp1[j];
		    ++k_;
		}
		++pp1[j];
		if (k_ < kend_ && pi0_[k_] == pi0[k])
		    ++k_;
		++k;
	    }
	    while (k_ < kend_) {
		if (pi0_[k_] > j) {
		    k_ = kend_;
		    break;
		}
		++pp1[j];
		++k_;
	    }
	}
	
	PROTECT(x1 = allocVector(kind2type(clt[0]), pp1[n-1]));
	PROTECT(i1 = allocVector(INTSXP, pp1[n-1]));
	pi1 = INTEGER(i1);

#define CRSPARSE_SYMMPART_GENERAL(_DO_ASSIGN_,				\
				  _DO_ASSIGN_FROM_TRANSPOSE_,		\
				  _DO_INCR_FROM_TRANSPOSE_)		\
	do {								\
	    for (j = 0, k = 0, k_ = 0; j < n; ++j) {			\
		kend = pp0[j];						\
		kend_ = pp0_[j];					\
		while (k < kend) {					\
		    if (pi0[k] > j) {					\
			k = kend;					\
			break;						\
		    }							\
		    while (k_ < kend_ && pi0_[k_] < pi0[k]) {		\
			*pi1 = pi0_[k_];				\
			_DO_ASSIGN_FROM_TRANSPOSE_;			\
			++pi1; ++px1; ++k_;				\
		    }							\
		    *pi1 = pi0[k];					\
		    _DO_ASSIGN_;					\
		    if (k_ < kend_ && pi0_[k_] == pi0[k]) {		\
			_DO_INCR_FROM_TRANSPOSE_;			\
			++k_;						\
		    }							\
		    ++pi1; ++px1; ++k;					\
		}							\
		while (k_ < kend_) {					\
		    if (pi0_[k_] > j) {					\
			k_ = kend_;					\
			break;						\
		    }							\
		    *pi1 = pi0_[k_];					\
		    _DO_ASSIGN_FROM_TRANSPOSE_;				\
		    ++pi1; ++px1; ++k_;					\
		}							\
	    }								\
	} while (0)
	
	switch (clf[0]) {
	case 'n':
	{
	    double *px1 = REAL(x1);
	    CRSPARSE_SYMMPART_GENERAL(
		*px1  = 0.5,
		*px1  = 0.5,
		*px1 += 0.5);
	    break;
	}
	case 'l':
	{
	    int *px0 = LOGICAL(x0), *px0_ = LOGICAL(x0_);
	    double *px1 = REAL(x1);
	    CRSPARSE_SYMMPART_GENERAL(
		*px1  = (px0[k] == NA_LOGICAL
			 ? NA_REAL : ((px0[k] != 0) ? 0.5 : 0.0)),
		*px1  = (px0_[k_] == NA_LOGICAL
			 ? NA_REAL : ((px0_[k_] != 0) ? 0.5 : 0.0)),
		*px1 += (px0_[k_] == NA_LOGICAL
			 ? NA_REAL : ((px0_[k_] != 0) ? 0.5 : 0.0)));
	    break;
	}
	case 'i':
	{
	    int *px0 = INTEGER(x0), *px0_ = INTEGER(x0_);
	    double *px1 = REAL(x1);
	    CRSPARSE_SYMMPART_GENERAL(
		*px1  = (px0[k] == NA_INTEGER
			 ? NA_REAL : 0.5 * (double) px0[k]),
		*px1  = (px0_[k_] == NA_INTEGER
			 ? NA_REAL : 0.5 * (double) px0_[k_]),
		*px1 += (px0_[k_] == NA_INTEGER
			 ? NA_REAL : 0.5 * (double) px0_[k_]));
	    break;
	}
	case 'd':
	{
	    double *px0 = REAL(x0), *px0_ = REAL(x0_),
		*px1 = REAL(x1);
	    CRSPARSE_SYMMPART_GENERAL(
		*px1  = 0.5 * px0[k],
		*px1  = 0.5 * px0_[k_],
		*px1 += 0.5 * px0_[k_]);
	    break;
	}
	case 'z':
	{
	    Rcomplex *px0 = COMPLEX(x0), *px0_ = COMPLEX(x0_),
		*px1 = COMPLEX(x1);
	    CRSPARSE_SYMMPART_GENERAL(
		do {
		    (*px1).r  = 0.5 * px0[k].r;
		    (*px1).i  = 0.5 * px0[k].i;
		} while (0),
		do {
		    (*px1).r  = 0.5 * px0_[k_].r;
		    (*px1).i  = 0.5 * px0_[k_].i;
		} while (0),
		do {
		    (*px1).r += 0.5 * px0_[k_].r;
		    (*px1).i += 0.5 * px0_[k_].i;
		} while (0));
	    break;
	}
	default:
	    break;
	}
	
#undef CRSPARSE_SYMMPART_GENERAL
	
	SET_SLOT(to, Matrix_DimSym, dim);
	set_symmetrized_DimNames(to, dimnames, -1);
	if (clf[2] == 'R')
	    SET_SLOT(to, Matrix_uploSym, mkString("L"));
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to, iSym, i1);
	SET_SLOT(to, Matrix_xSym, x1);
	UNPROTECT(5);
	
    } else if (clf[1] == 't') {

	PROTECT(x1 = (diagU && clf[0] == clt[0]
		      ? x0 : allocVector(kind2type(clt[0]), nnz)));
	
#define CRSPARSE_SYMMPART_TRIANGULAR(_DO_ASSIGN_, _DO_HALF_)	\
	do {							\
	    for (j = 0, k = 0; j < n; ++j) {			\
		kend = pp0[j];					\
		while (k < kend) {				\
		    _DO_ASSIGN_;				\
		    if (pi0[k] != j)				\
			_DO_HALF_;				\
		    ++k;					\
		}						\
	    }							\
	} while (0)

	SPARSE_SYMMPART_CASES(CRSPARSE_SYMMPART_TRIANGULAR);
	
#undef CRSPARSE_SYMMPART_TRIANGULAR
#undef SPARSE_SYMMPART_CASES
	
	SET_SLOT(to, Matrix_DimSym, dim);
	set_symmetrized_DimNames(to, dimnames, -1);
	SET_SLOT(to, Matrix_uploSym, uplo);
	SET_SLOT(to, Matrix_pSym, p0);
	SET_SLOT(to, iSym, i0);
	SET_SLOT(to, Matrix_xSym, x1);
	UNPROTECT((diagU) ? 3 : 2);
	
    } else {

	SPARSE_SYMMPART_CASES_TRIVIAL;

#undef SPARSE_SYMMPART_CASES_TRIVIAL
	
	SET_SLOT(to, Matrix_DimSym, dim);
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	SET_SLOT(to, Matrix_uploSym, uplo);
	SET_SLOT(to, Matrix_pSym, p0);
	SET_SLOT(to, iSym, i0);
	SET_SLOT(to, Matrix_xSym, x1);
	UNPROTECT(2);
	
    }

    return to;
}

/* skewpart(<TsparseMatrix>) */
SEXP Tsparse_skewpart(SEXP from)
{
    static const char *valid[] = { VALID_TSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "Tsparse_skewpart");
    const char *clf = valid[ivalid];

    SEXP to,
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = (clf[1] != 'g') ? GET_SLOT(from, Matrix_uploSym) : R_NilValue,
	i0 = GET_SLOT(from, Matrix_iSym),
	j0 = GET_SLOT(from, Matrix_jSym),
	x0 = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue,
	x1;

    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to get skew-symmetric part of non-square matrix"));

    if (clf[1] == 's') {

	if (clf[0] != 'z') {
	    /* Skew-symmetric part of symmetric matrix is zero matrix */
	    PROTECT(to = NEW_OBJECT_OF_CLASS("dsTMatrix"));
	} else {
	    /* Skew-symmetric part of Hermitian matrix is imaginary part */
	    PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
	    PROTECT(x1 = duplicate(x0));
	    zeroRe(x1);
	    SET_SLOT(to, Matrix_iSym, i0);
	    SET_SLOT(to, Matrix_jSym, j0);
	    SET_SLOT(to, Matrix_xSym, x1);
	    UNPROTECT(1);
	}
	
	SET_SLOT(to, Matrix_DimSym, dim);
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	SET_SLOT(to, Matrix_uploSym, uplo);
	UNPROTECT(1);
	return to;
	
    }

    char clt[] = ".gTMatrix";
    clt[0] = (clf[0] != 'z') ? 'd' : 'z';
    PROTECT(to = NEW_OBJECT_OF_CLASS(clt));
    SET_SLOT(to, Matrix_DimSym, dim);
    set_symmetrized_DimNames(to, dimnames, -1);

    R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = nnz0;
    int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
    
    if (clf[1] != 't' || *diag_P(from) == 'N')
	for (k = 0; k < nnz0; ++k)
	    if (pi0[k] == pj0[k])
		--nnz1;
    nnz1 *= 2;

    SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
	j1 = PROTECT(allocVector(INTSXP, nnz1));
    int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
    PROTECT(x1 = allocVector(kind2type(clt[0]), nnz1));

#define TSPARSE_SKEWPART_GENERAL(_DO_ASSIGN_)				\
    do {								\
	for (k = 0; k < nnz0; ++k) {					\
	    if (pi0[k] != pj0[k]) {					\
		*(pi1++) = pi0[k];					\
		*(pj1++) = pj0[k];					\
		*(pi1++) = pj0[k];					\
		*(pj1++) = pi0[k];					\
		_DO_ASSIGN_;						\
	    }								\
	}								\
    } while (0)

    switch (clf[0]) {
    case 'n':
    {
	double *px1 = REAL(x1);
	TSPARSE_SKEWPART_GENERAL(
	    do {
		*(px1++) =  0.5;
		*(px1++) = -0.5;
	    } while (0));
	break;
    }
    case 'l':
    {
	int *px0 = LOGICAL(x0);
	double *px1 = REAL(x1);
	TSPARSE_SKEWPART_GENERAL(
	    do {
		if (px0[k] == NA_LOGICAL) {
		    *(px1++) = NA_REAL;
		    *(px1++) = NA_REAL;
		} else if (px0[k] != 0) {
		    *(px1++) =  0.5;
		    *(px1++) = -0.5;
		} else {
		    *(px1++) = 0.0;
		    *(px1++) = 0.0;
		}
	    } while (0));
	break;
    }
    case 'i':
    {
	int *px0 = INTEGER(x0);
	double *px1 = REAL(x1);
	TSPARSE_SKEWPART_GENERAL(
	    do {
		if (px0[k] == NA_INTEGER) {
		    *(px1++) = NA_REAL;
		    *(px1++) = NA_REAL;
		} else {
		    *(px1++) =  0.5 * (double) px0[k];
		    *(px1++) = -0.5 * (double) px0[k];
		}
	    } while (0));
	break;
    }
    case 'd':
    {
	double *px0 = REAL(x0), *px1 = REAL(x1);
	TSPARSE_SKEWPART_GENERAL(
	    do {
		*(px1++) =  0.5 * px0[k];
		*(px1++) = -0.5 * px0[k];
	    } while (0));
	break;
    }
    case 'z':
    {
	Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
	TSPARSE_SKEWPART_GENERAL(
	    do {
		(*(px1  )).r =  0.5 * px0[k].r;
		(*(px1++)).i =  0.5 * px0[k].i;
		(*(px1  )).r = -0.5 * px0[k].r;
		(*(px1++)).i = -0.5 * px0[k].i;
	    } while (0));
	break;
    }
    default:
	break;
    }
    
#undef TSPARSE_SKEWPART_GENERAL
    
    SET_SLOT(to, Matrix_iSym, i1);
    SET_SLOT(to, Matrix_jSym, j1);
    SET_SLOT(to, Matrix_xSym, x1);
    UNPROTECT(4);
    return to;
}

/* skewpart(<[CR]sparseMatrix>) */
SEXP CRsparse_skewpart(SEXP from)
{
    static const char *valid[] = { VALID_CRSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "CRsparse_skewpart");
    const char *clf = valid[ivalid];

    SEXP to,
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	uplo = (clf[1] != 'g') ? GET_SLOT(from, Matrix_uploSym) : R_NilValue,
	iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
	p0 = GET_SLOT(from, Matrix_pSym),
	i0 = GET_SLOT(from, iSym),
	x0 = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue,
	x1 = R_NilValue;

    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to get skew-symmetric part of non-square matrix"));
    
    if (clf[1] == 's') {

	if (clf[0] != 'z') {
	    /* Skew-symmetric part of symmetric matrix is zero matrix */
	    char clt[] = "ds.Matrix";
	    clt[2] = clf[2];
	    PROTECT(to = NEW_OBJECT_OF_CLASS(clt));
	    R_xlen_t n1a = (R_xlen_t) n + 1;
	    SEXP p1 = PROTECT(allocVector(INTSXP, n1a));
	    int *pp1 = INTEGER(p1);
	    Memzero(pp1, n1a);
	    SET_SLOT(to, Matrix_pSym, p1);
	} else {
	    /* Skew-symmetric part of Hermitian matrix is imaginary part */
	    PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
	    PROTECT(x1 = duplicate(x0));
	    zeroRe(x1);
	    SET_SLOT(to, Matrix_iSym, p0);
	    SET_SLOT(to, Matrix_iSym, i0);
	    SET_SLOT(to, Matrix_xSym, x1);
	}
	
	SET_SLOT(to, Matrix_DimSym, dim);
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	SET_SLOT(to, Matrix_uploSym, uplo);
	UNPROTECT(2);
	return to;
	
    }
    
    char clt[] = ".g.Matrix";
    clt[0] = (clf[0] != 'z') ? 'd' : 'z';
    clt[2] = clf[2];
    PROTECT(to = NEW_OBJECT_OF_CLASS(clt));
    SET_SLOT(to, Matrix_DimSym, dim);
    set_symmetrized_DimNames(to, dimnames, -1);

    PROTECT(from = R_sparse_transpose(from));
    SEXP i1,
	p1  = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
	p0_ = GET_SLOT(from, Matrix_pSym),
	i0_ = GET_SLOT(from, iSym),
	x0_ = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue;
    int j, k, kend, k_, kend_,
	*pp0  = INTEGER(p0 ) + 1, *pi0  = INTEGER(i0 ),
	*pp0_ = INTEGER(p0_) + 1, *pi0_ = INTEGER(i0_),
	*pp1  = INTEGER(p1 )    , *pi1                ,
	*pp1_;
    Calloc_or_Alloca_TO(pp1_, n, int);
    Memzero(pp1_, n);

    /* Counting number of nonzero elements in each "column" of result ... */

    for (j = 0, k = 0, k_ = 0; j < n; ++j) {
	kend = pp0[j];
	kend_ = pp0_[j];
	while (k < kend) {
	    if (pi0[k] >= j) {
		k = kend;
		break;
	    }
	    while (k_ < kend_ && pi0_[k_] < pi0[k]) {
		++pp1_[j];
		++pp1_[pi0_[k_]];
		++k_;
	    }
	    ++pp1_[j];
	    ++pp1_[pi0[k]];
	    if (k_ < kend_ && pi0_[k_] == pi0[k])
		++k_;
	    ++k;
	}
	while (k_ < kend_) {
	    if (pi0_[k_] >= j) {
		k_ = kend_;
		break;
	    }
	    ++pp1_[j];
	    ++pp1_[pi0_[k_]];
	    ++k_;
	}
    }

    *(pp1++) = 0;
    for (j = 0; j < n; ++j) {
	pp1[j] = pp1[j-1] + pp1_[j];
	pp1_[j] = pp1[j-1];
    }

    PROTECT(x1 = allocVector(kind2type(clt[0]), pp1[n-1]));
    PROTECT(i1 = allocVector(INTSXP, pp1[n-1]));
    pi1 = INTEGER(i1);

#define CRSPARSE_SKEWPART_GENERAL(_DO_ASSIGN_,			\
				  _DO_ASSIGN_FROM_TRANSPOSE_,	\
				  _DO_INCR_FROM_TRANSPOSE_,	\
				  _DO_NEGATE_,			\
				  _DO_NEGATE_FROM_TRANSPOSE_)	\
    do {							\
	for (j = 0, k = 0, k_ = 0; j < n; ++j) {		\
	    kend = pp0[j];					\
	    kend_ = pp0_[j];					\
	    while (k < kend) {					\
		if (pi0[k] >= j) {				\
		    k = kend;					\
		    break;					\
		}						\
		while (k_ < kend_ && pi0_[k_] < pi0[k]) {	\
		    pi1[pp1_[j]] = pi0_[k_];			\
		    _DO_ASSIGN_FROM_TRANSPOSE_;			\
		    pi1[pp1_[pi0_[k_]]] = j;			\
		    _DO_NEGATE_FROM_TRANSPOSE_;			\
		    ++pp1_[j];					\
		    ++pp1_[pi0_[k_]];				\
		    ++k_;					\
		}						\
		pi1[pp1_[j]] = pi0[k];				\
		_DO_ASSIGN_;					\
		if (k_ < kend_ && pi0_[k_] == pi0[k]) {		\
		    _DO_INCR_FROM_TRANSPOSE_;			\
		    ++k_;					\
		}						\
		pi1[pp1_[pi0[k]]] = j;				\
		_DO_NEGATE_;					\
		++pp1_[j];					\
		++pp1_[pi0[k]];					\
		++k;						\
	    }							\
	    while (k_ < kend_) {				\
		if (pi0_[k_] >= j) {				\
		    k_ = kend_;					\
		    break;					\
		}						\
		pi1[pp1_[j]] = pi0_[k_];			\
		_DO_ASSIGN_FROM_TRANSPOSE_;			\
		pi1[pp1_[pi0_[k_]]] = j;			\
		_DO_NEGATE_FROM_TRANSPOSE_;			\
		++pp1_[j];					\
		++pp1_[pi0_[k_]];				\
		++k_;						\
	    }							\
	}							\
    } while (0)

    switch (clf[0]) {
    case 'n':
    {
	double *px1 = REAL(x1);
	CRSPARSE_SKEWPART_GENERAL(
	    px1[pp1_[j]]         =  0.5,
	    px1[pp1_[j]]         = -0.5,
	    px1[pp1_[j]]        -=  0.5,
	    px1[pp1_[pi0[k]]]    = -px1[pp1_[j]],
	    px1[pp1_[pi0_[k_]]]  = -px1[pp1_[j]]);
	break;
    }
    case 'l':
    {
	int *px0 = LOGICAL(x0), *px0_ = LOGICAL(x0_);
	double *px1 = REAL(x1);
	CRSPARSE_SKEWPART_GENERAL(
	    px1[pp1_[j]]         = (px0[k]   == NA_LOGICAL
				    ? NA_REAL : ((px0[k]   != 0) ?  0.5 : 0.0)),
	    px1[pp1_[j]]         = (px0_[k_] == NA_LOGICAL
				    ? NA_REAL : ((px0_[k_] != 0) ? -0.5 : 0.0)),
	    px1[pp1_[j]]        -= (px0_[k_] == NA_LOGICAL
				    ? NA_REAL : ((px0_[k_] != 0) ?  0.5 : 0.0)),
	    px1[pp1_[pi0[k]]]    = -px1[pp1_[j]],
	    px1[pp1_[pi0_[k_]]]  = -px1[pp1_[j]]);
	break;
    }
    case 'i':
    {
	int *px0 = INTEGER(x0), *px0_ = INTEGER(x0_);
	double *px1 = REAL(x1);
	CRSPARSE_SKEWPART_GENERAL(
	    px1[pp1_[j]]         = (px0[k]   == NA_INTEGER
				    ? NA_REAL :  0.5 * (double) px0[k]),
	    px1[pp1_[j]]         = (px0_[k_] == NA_INTEGER
				    ? NA_REAL : -0.5 * (double) px0_[k_]),
	    px1[pp1_[j]]        -= (px0_[k_] == NA_INTEGER
				    ? NA_REAL :  0.5 * (double) px0_[k_]),
	    px1[pp1_[pi0[k]]]    = -px1[pp1_[j]],
	    px1[pp1_[pi0_[k_]]]  = -px1[pp1_[j]]);
	break;
    }
    case 'd':
    {
	double *px0 = REAL(x0), *px0_ = REAL(x0_),
	    *px1 = REAL(x1);
	CRSPARSE_SKEWPART_GENERAL(
	    px1[pp1_[j]]         =  0.5 * px0[k],
	    px1[pp1_[j]]         = -0.5 * px0_[k_],
	    px1[pp1_[j]]        -=  0.5 * px0_[k_],
	    px1[pp1_[pi0[k]]]    = -px1[pp1_[j]],
	    px1[pp1_[pi0_[k_]]]  = -px1[pp1_[j]]);
	break;
    }
    case 'z':
    {
	Rcomplex *px0 = COMPLEX(x0), *px0_ = COMPLEX(x0_),
	    *px1 = COMPLEX(x1);
	CRSPARSE_SKEWPART_GENERAL(
	    do {
		px1[pp1_[j]].r          =  0.5 * px0[k].r;
		px1[pp1_[j]].i          =  0.5 * px0[k].i;
	    } while (0),
	    do {
		px1[pp1_[j]].r          = -0.5 * px0_[k_].r;
		px1[pp1_[j]].i          = -0.5 * px0_[k_].i;
	    } while (0),
	    do {
		px1[pp1_[j]].r         -=  0.5 * px0_[k_].r;
		px1[pp1_[j]].i         -=  0.5 * px0_[k_].i;
	    } while (0),
	    do {
		px1[pp1_[pi0[k]]].r     = -px1[pp1_[j]].r;
		px1[pp1_[pi0[k]]].i     = -px1[pp1_[j]].i;
	    } while (0),
	    do {
		px1[pp1_[pi0_[k_]]].r   = -px1[pp1_[j]].r;
		px1[pp1_[pi0_[k_]]].i   = -px1[pp1_[j]].i;
	    } while (0));
	break;
    }
    default:
	break;
    }

#undef CRSPARSE_SKEWPART_GENERAL
    
    SET_SLOT(to, Matrix_pSym, p1);
    SET_SLOT(to, iSym, i1);
    SET_SLOT(to, Matrix_xSym, x1);
    Free_FROM(pp1_, n);
    UNPROTECT(5);
    return to;
}

/* as(<[CR]sparseMatrix>, "TsparseMatrix") */
SEXP CRsparse_as_Tsparse(SEXP from)
{
    static const char *valid[] = { VALID_CRSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "CRsparse_as_Tsparse");
    const char *clf = valid[ivalid];

    char clt[] = "..TMatrix";
    clt[0] = clf[0];
    clt[1] = clf[1];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (clf[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (clf[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    else
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    if (clf[0] != 'n')
	SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
    
    SEXP iSym, jSym;
    int n, nnz, *pp = INTEGER(GET_SLOT(from, Matrix_pSym));
    if (clf[2] == 'C') {
	iSym = Matrix_iSym;
	jSym = Matrix_jSym;
	n = INTEGER(dim)[1];
    } else {
	iSym = Matrix_jSym;
	jSym = Matrix_iSym;
	n = INTEGER(dim)[0];
    }
    nnz = pp[n];
    
    SEXP j1 = PROTECT(allocVector(INTSXP, nnz));
    int j, k, kend, *pj1 = INTEGER(j1);
    SET_SLOT(to, iSym, GET_SLOT(from, iSym));
    SET_SLOT(to, jSym, j1);
    
    for (j = 0, k = 0; j < n; ++j) {
	kend = *(++pp);
	while (k < kend) {
	    *(pj1++) = j;
	    ++k;
	}
    }
    
    UNPROTECT(2);
    return to;
}

SEXP Tsparse_as_CRsparse(SEXP from, SEXP Csparse)
{
    static const char *valid[] = { VALID_TSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "Tsparse_as_CRsparse");
    const char *clf = valid[ivalid];
    int doC = (asLogical(Csparse) != 0);

    char clt[] = "...Matrix";
    clt[0] = clf[0];
    clt[1] = clf[1];
    clt[2] = (doC) ? 'C' : 'R';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	iSym = (doC) ? Matrix_iSym : Matrix_jSym,
	jSym = (doC) ? Matrix_jSym : Matrix_iSym,
	i0 = GET_SLOT(from, iSym),
	j0 = GET_SLOT(from, jSym),
	x0 = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue;
    int *pdim = INTEGER(dim), *pi0 = INTEGER(i0), *pj0 = INTEGER(j0),
	m = pdim[0],
	n = pdim[1],
	m_ = (doC) ? m : n,
	n_ = (doC) ? n : m,
	r_ = (m_ < n_) ? n_ : m_;
    R_xlen_t nnz0 = XLENGTH(i0), nnz1 = 0;

    /* FIXME? we would ideally only throw an error if the number
       of _unique_ (i,j) pairs exceeds INT_MAX ... 
    */
    if (nnz0 > INT_MAX)
	error(_("unable to coerce from TsparseMatrix to [CR]sparseMatrix"
		"when length of 'i' slot exceeds 2^31-1"));
    
    SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n_ + 1)),
	i1 = R_NilValue, x1 = R_NilValue;
    R_xlen_t k, kstart, kend, kend_, w = (R_xlen_t) m_ + r_ + m_;
    int *pp1 = INTEGER(p1), *pi1, *pj_, *workA, *workB, *workC, i, j;
    *(pp1++) = 0;
    Calloc_or_Alloca_TO(pj_, nnz0, int);
    Calloc_or_Alloca_TO(workA, w, int);
    workB = workA + m_;
    workC = workB + r_;
    
    /* 1. Tabulate column indices in workA[i]
       
          workA[.]: unused
          workB[.]: unused
          workC[.]: unused			  
    */
    
#define T_AS_CR_1				\
    do {					\
	Memzero(workA, m_);			\
	for (k = 0; k < nnz0; ++k)		\
	    ++workA[pi0[k]];			\
    } while (0)
    
    /* 2. Compute cumulative sum in workA[i], copying to workB[i]
       
          workA[i]: number of column indices listed for row i,
	  incl. duplicates
    */
    
#define T_AS_CR_2					\
    do {						\
	if (r_ > 0) {					\
	    workB[0] = 0;				\
	    for (i = 1; i < m_; ++i)			\
		workA[i] += (workB[i] = workA[i-1]);	\
	}						\
    } while (0)
    
    /* 3. Group column indices and data by row in pj_[k], px_[k]

          workA[i]: number of column indices listed for row <= i,
                    incl. duplicates
          workB[i]: number of column indices listed for row <  i,
                    incl. duplicates
    */
    
#define T_AS_CR_3(_XASSIGN_)					\
    do {							\
	for (k = 0; k < nnz0; ++k) {				\
	    pj_[workB[pi0[k]]] = pj0[k];			\
	    _XASSIGN_; /* px_[workB[pi0[k]]] = px0[k]; */	\
	    ++workB[pi0[k]];					\
	}							\
    } while (0)
    
    /* 4. Gather _unique_ column indices at the front of each group,
          aggregating data accordingly; record in workC[i] where the
	  unique column indices stop and the duplicates begin

	  workB[.]: unused
	    pj_[k]: column indices grouped by row, incl. duplicates, unsorted
	    px_[k]: corresponding data
    */
    
#define T_AS_CR_4(_XASSIGN_, _XINCR_)					\
    do {								\
	k = 0;								\
	for (j = 0; j < n_; ++j)					\
	    workB[j] = -1;						\
	for (i = 0; i < m_; ++i) {					\
	    kstart = k;							\
	    kend_ = k;							\
	    kend = workA[i];						\
	    while (k < kend) {						\
		if (workB[pj_[k]] < kstart) {				\
		    /* Have not yet seen this column index */		\
		    workB[pj_[k]] = kend_;				\
		    pj_[kend_] = pj_[k];				\
		    _XASSIGN_; /* px_[kend_] = px_[k]; */		\
		    ++kend_;						\
		} else {						\
		    /* Have already seen this column index */		\
		    _XINCR_; /* px_[workB[pj_[k]]] += px_[k]; */	\
		}							\
		++k;							\
	    }								\
	    workC[i] = kend_;						\
	    nnz1 += kend_ - kstart;					\
	}								\
    } while (0)
    
    /* 5. Tabulate _unique_ column indices in workB[j]

          workC[i]: pointer to first non-unique column index in row i  
            pi_[k]: column indices grouped by row, with unique indices in front
                    i.e., in positions workA[i-1] <= k < workC[i]
	    px_[k]: corresponding data, "cumulated" appropriately 
    */

#define T_AS_CR_5				\
    do {					\
	k = 0;					\
	Memzero(workB, n_);			\
	for (i = 0; i < m_; ++i) {		\
	    kend_ = workC[i];			\
	    while (k < kend_) {			\
		++workB[pj_[k]];		\
		++k;				\
	    }					\
	    k = workA[i];			\
	}					\
    } while (0)
    
    /* 6. Compute cumulative sum in pp1[j], copying to workB[j]
    
          workB[j]: number of nonzero elements in column j
    */

#define T_AS_CR_6				\
    do {					\
	for (j = 0; j < n_; ++j) {		\
	    pp1[j] = pp1[j-1] + workB[j];	\
	    workB[j] = pp1[j-1];		\
	}					\
    } while (0)

    /* 7. Pop unique (i,j) pairs from the unsorted stacks 0 <= i < m
          onto new stacks 0 <= j < n, which will be sorted 
       
	  workB[j]: number of nonzero elements in columns <  j
	    pp1[j]: number of nonzero elements in columns <= j
    */

#define T_AS_CR_7(_XASSIGN_)					\
    do {							\
	PROTECT(i1 = allocVector(INTSXP, nnz1));		\
	pi1 = INTEGER(i1);					\
	k = 0;							\
	for (i = 0; i < m_; ++i) {				\
	    kend_ = workC[i];					\
	    while (k < kend_) {					\
		pi1[workB[pj_[k]]] = i;				\
		_XASSIGN_; /* px1[workB[pj_[k]]] = px_[k]; */	\
		++workB[pj_[k]];				\
		++k;						\
	    }							\
	    k = workA[i];					\
	}							\
    } while (0)

#define T_AS_CR_N				\
    do {					\
	T_AS_CR_1;				\
	T_AS_CR_2;				\
	T_AS_CR_3();				\
	T_AS_CR_4(, );				\
	T_AS_CR_5;				\
	T_AS_CR_6;				\
	T_AS_CR_7();				\
    } while (0)
    
#define T_AS_CR_X(_CTYPE_, _PTR_, _SEXPTYPE_, _XINCR_)	\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1, *px_;		\
	Calloc_or_Alloca_TO(px_, nnz0, _CTYPE_);	\
	T_AS_CR_1;					\
	T_AS_CR_2;					\
	T_AS_CR_3(px_[workB[pi0[k]]] = px0[k]);		\
	T_AS_CR_4(px_[kend_] = px_[k], _XINCR_);	\
	T_AS_CR_5;					\
	T_AS_CR_6;					\
	PROTECT(x1 = allocVector(_SEXPTYPE_, nnz1));	\
	px1 = _PTR_(x1);				\
	T_AS_CR_7(px1[workB[pj_[k]]] = px_[k]);		\
	Free_FROM(px_, nnz0);				\
    } while (0)

#define T_AS_CR_CASES(_KIND_, _DO_N_, _DO_X_)				\
    do {								\
	switch (_KIND_) {						\
	case 'n':							\
	    _DO_N_;							\
	    break;							\
	case 'l':							\
	    _DO_X_(int, LOGICAL, LGLSXP,				\
		   do {							\
		       if (px_[k] != 0) {				\
			   if (px_[k] != NA_LOGICAL)			\
			       px_[workB[pj_[k]]] = 1;			\
			   else if (px_[workB[pj_[k]]] == 0)		\
			       px_[workB[pj_[k]]] = NA_LOGICAL;		\
		       }						\
		   } while (0));					\
	    break;							\
	case 'i':							\
	    _DO_X_(int, INTEGER, INTSXP,				\
		   /* FIXME: not detecting integer overflow here */	\
		   px_[workB[pj_[k]]] += px_[k]);			\
	    break;							\
	case 'd':							\
	    _DO_X_(double, REAL, REALSXP,				\
		   px_[workB[pj_[k]]] += px_[k]);			\
	    break;							\
	case 'z':							\
	    _DO_X_(Rcomplex, COMPLEX, CPLXSXP,				\
		   do {							\
		       px_[workB[pj_[k]]].r += px_[k].r;		\
		       px_[workB[pj_[k]]].i += px_[k].i;		\
		   } while (0));					\
	    break;							\
	default:							\
	    break;							\
	}								\
    } while (0)

    T_AS_CR_CASES(clf[0], T_AS_CR_N, T_AS_CR_X);
    Free_FROM(workA, w);
    Free_FROM(pj_, nnz0);
    
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    if (clf[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (clf[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    else
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    SET_SLOT(to, Matrix_pSym, p1);
    SET_SLOT(to, iSym, i1);
    if (clf[0] != 'n')
	SET_SLOT(to, Matrix_xSym, x1);
    
    UNPROTECT((clf[0] == 'n') ? 3 : 4);
    return to;
}

SEXP Tsparse_aggregate(SEXP from)
{
    static const char *valid[] = { VALID_TSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "Tsparse_aggregate");
    const char *cl = valid[ivalid];

    /* Need to behave as Tsparse_as_CRsparse(from, FALSE) 
       in order to get aggregated triplets sorted by column */
    SEXP to,
	dim = GET_SLOT(from, Matrix_DimSym),
	i0 = GET_SLOT(from, Matrix_jSym),
	j0 = GET_SLOT(from, Matrix_iSym),
	x0 = (cl[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue,
	i1 = R_NilValue,
	j1 = R_NilValue,
	x1 = R_NilValue;
    int *pdim = INTEGER(dim),
	*pi0 = INTEGER(i0), *pj0 = INTEGER(j0), *pi1, *pj1,
	m_ = pdim[1], n_ = pdim[0], r_ = (m_ < n_) ? n_ : m_;
    R_xlen_t nnz0 = XLENGTH(j0), nnz1 = 0;

    if (nnz0 > INT_MAX)
	error(_("unable to aggregate TsparseMatrix with 'i' slot "
		"of length exceeding 2^31-1"));
    
    R_xlen_t k, kstart, kend, kend_, w = (R_xlen_t) m_ + r_ + m_;
    int *pj_, *workA, *workB, *workC, i, j;
    Calloc_or_Alloca_TO(pj_, nnz0, int);
    Calloc_or_Alloca_TO(workA, w, int);
    workB = workA + m_;
    workC = workB + r_;
    
#define ALLOC_TRIPLET(_XASSIGN_)			\
    do {						\
	PROTECT(i1 = allocVector(INTSXP, nnz1));	\
	PROTECT(j1 = allocVector(INTSXP, nnz1));	\
	pi1 = INTEGER(i1);				\
	pj1 = INTEGER(j1);				\
							\
	k = 0;						\
	for (i = 0; i < m_; ++i) {			\
	    kend_ = workC[i];				\
	    while (k < kend_) {				\
		*(pi1++) = i;				\
		*(pj1++) = pj_[k];			\
		_XASSIGN_; /* *(px1++) = px_[k]; */	\
		++k;					\
	    }						\
	    k = workA[i];				\
	}						\
    } while (0)

#define T_AGGR_N				\
    do {					\
	T_AS_CR_1;				\
	T_AS_CR_2;				\
	T_AS_CR_3();				\
	T_AS_CR_4(, );				\
	if (nnz1 == nnz0)			\
	    return from;			\
	ALLOC_TRIPLET();			\
    } while (0)

#define T_AGGR_X(_CTYPE_, _PTR_, _SEXPTYPE_, _XINCR_)	\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1, *px_;		\
	Calloc_or_Alloca_TO(px_, nnz0, _CTYPE_);	\
	T_AS_CR_1;					\
	T_AS_CR_2;					\
	T_AS_CR_3(px_[workB[pi0[k]]] = px0[k]);		\
	T_AS_CR_4(px_[kend_] = px_[k], _XINCR_);	\
	if (nnz1 == nnz0)				\
	    return from;				\
	PROTECT(x1 = allocVector(_SEXPTYPE_, nnz1));	\
	px1 = _PTR_(x1);				\
	ALLOC_TRIPLET(*(px1++) = px_[k]);		\
	Free_FROM(px_, nnz0);				\
    } while (0)

    T_AS_CR_CASES(cl[0], T_AGGR_N, T_AGGR_X);
    Free_FROM(workA, w);
    Free_FROM(pj_, nnz0);

    PROTECT(to = NEW_OBJECT_OF_CLASS(cl));
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    if (cl[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (cl[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    else
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    SET_SLOT(to, Matrix_iSym, j1);
    SET_SLOT(to, Matrix_jSym, i1);
    if (cl[0] != 'n')
	SET_SLOT(to, Matrix_xSym, x1);

    UNPROTECT((cl[0] == 'n') ? 3 : 4);
    return to;
}

/* as(t(<[CR]sparseMatrix>), "[RC]sparseMatrix") */
SEXP tCRsparse_as_RCsparse(SEXP from)
{
    static const char *valid[] = { VALID_CRSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "tCRsparse_as_RCsparse");
    const char *clf = valid[ivalid];
    
    char clt[] = "...Matrix";
    clt[0] = clf[0];
    clt[1] = clf[1];
    clt[2] = (clf[2] == 'C') ? 'R' : 'C';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    
    if (m == n)	{
	SET_SLOT(to, Matrix_DimSym, dim);
    } else {
	pdim = INTEGER(GET_SLOT(to, Matrix_DimSym));
	pdim[0] = n;
	pdim[1] = m;
    }
    if (clf[1] == 's')
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    else
	set_reversed_DimNames(to, dimnames);
    SET_SLOT(to, Matrix_pSym, GET_SLOT(from, Matrix_pSym));
    if (clf[2] == 'R')
	SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_jSym));
    else
	SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_iSym));
    if (clf[0] != 'n')
	SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
    if (clf[1] != 'g')
	SET_SLOT(to, Matrix_uploSym,
		 mkString((*uplo_P(from) == 'U') ? "L" : "U"));
    if (clf[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    if (clf[1] == 's')
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));

    UNPROTECT(1);
    return to;
}

/* isDiagonal(<[CR]sparseMatrix>) */
#define CR_IS_DIAGONAL(_C_, _I_)					\
SEXP _C_ ## sparse_is_diagonal(SEXP obj)				\
{									\
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];	\
    if (pdim[1] != n)							\
	return ScalarLogical(0);					\
    int *pp = INTEGER(GET_SLOT(obj, Matrix_pSym));			\
    if (pp[n] > n)							\
	return ScalarLogical(0);					\
    int d, j,								\
	*pi = INTEGER(GET_SLOT(obj, Matrix_ ## _I_ ## Sym));		\
    for (j = 0; j < n; ++j)						\
	if ((d = pp[j+1] - pp[j]) > 1 || (d == 1 && *(pi++) != j))	\
	    return ScalarLogical(0);					\
    return ScalarLogical(1);						\
}

/* Csparse_is_diagonal() */
CR_IS_DIAGONAL(C, i)
/* Rsparse_is_diagonal() */
CR_IS_DIAGONAL(R, j)

/* isDiagonal(<TsparseMatrix>) */
SEXP Tsparse_is_diagonal(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    R_xlen_t nnz = XLENGTH(i);
    if (nnz > n)
	return ScalarLogical(0);
    int *pi = INTEGER(i),
	*pj = INTEGER(GET_SLOT(obj, Matrix_jSym));
    R_xlen_t k;
    for (k = 0; k < nnz; ++k)
	if (*(pi++) != *(pj++))
	    return ScalarLogical(0);
    return ScalarLogical(1);
}

#define RETURN_TRUE_OF_KIND(_KIND_)					\
    do {								\
	SEXP ans = PROTECT(allocVector(LGLSXP, 1));			\
	LOGICAL(ans)[0] = 1;						\
	setAttrib(ans, install("kind"), _KIND_);			\
	UNPROTECT(1);							\
	return ans;							\
    } while (0)

/* isTriangular(<.g[CR]Matrix>, upper) */
#define CR_IS_TRIANGULAR(_C_, _I_, _UPPER_, _LOWER_)			\
SEXP _C_ ## sparse_is_triangular(SEXP obj, SEXP upper)			\
{									\
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];	\
    if (pdim[1] != n)							\
	return ScalarLogical(0);					\
    int j, k, kend,							\
	*pp = INTEGER(GET_SLOT(obj, Matrix_pSym)),			\
	*pi = INTEGER(GET_SLOT(obj, Matrix_ ## _I_ ## Sym)),		\
	need_upper = asLogical(upper);					\
    ++pp;								\
    if (need_upper == NA_LOGICAL) {					\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if (_LOWER_)						\
		    goto opposite;					\
		++k;							\
	    }								\
	}								\
	RETURN_TRUE_OF_KIND(mkString("U"));				\
    opposite:								\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if (_UPPER_)						\
		    return ScalarLogical(0);				\
		++k;							\
	    }								\
	}								\
	RETURN_TRUE_OF_KIND(mkString("L"));				\
    } else if (need_upper != 0) {					\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if (_LOWER_)						\
		    return ScalarLogical(0);				\
		++k;							\
	    }								\
	}								\
    } else {								\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if (_UPPER_)						\
		    return ScalarLogical(0);				\
		++k;							\
	    }								\
	}								\
    }									\
    return ScalarLogical(1);						\
}

/* Csparse_is_triangular() */
CR_IS_TRIANGULAR(C, i, pi[k] < j, pi[k] > j)
/* Rsparse_is_triangular() */
CR_IS_TRIANGULAR(R, j, pi[k] > j, pi[k] < j)

/* isTriangular(<.gTMatrix>, upper) */
SEXP Tsparse_is_triangular(SEXP obj, SEXP upper)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    int *pi = INTEGER(i),
	*pj = INTEGER(GET_SLOT(obj, Matrix_jSym)),
	need_upper = asLogical(upper);
    R_xlen_t k, nnz = XLENGTH(i);
    if (need_upper == NA_LOGICAL) {
	for (k = 0; k < nnz; ++k)
	    if (pi[k] > pj[k])
		goto opposite;
	RETURN_TRUE_OF_KIND(mkString("U"));
    opposite:
	for (k = 0; k < nnz; ++k)
	    if (pi[k] < pj[k])
		return ScalarLogical(0);
	RETURN_TRUE_OF_KIND(mkString("L"));
    } else if (need_upper != 0) {
	for (k = 0; k < nnz; ++k)
	    if (pi[k] > pj[k])
		return ScalarLogical(0);
    } else {
	for (k = 0; k < nnz; ++k)
	    if (pi[k] < pj[k])
		return ScalarLogical(0);
    }
    return ScalarLogical(1);
}

#define CR_IS_SYMMETRIC_LOOP(_XCOND_)					\
    do {								\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if ((i = pi[k]) >= j) {					\
		    if (i == j)						\
			++pp_[j];					\
		    k = kend;						\
		    break;						\
		}							\
		if (pp_[i] == pp[i] || pi[pp_[i]] != j || (_XCOND_))	\
		    return ScalarLogical(0);				\
		++pp_[i];						\
		++pp_[j];						\
		++k;							\
	    }								\
	}								\
    } while (0)

/* isSymmetric(<.g[CR]Matrix>, tol = 0, checkDN) */
#define CR_IS_SYMMETRIC(_C_, _I_)					\
SEXP _C_ ## sparse_is_symmetric(SEXP obj, SEXP checkDN)			\
{									\
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];	\
    if (pdim[1] != n)							\
	return ScalarLogical(0);					\
    if (asLogical(checkDN) != 0 &&					\
	!DimNames_is_symmetric(GET_SLOT(obj, Matrix_DimNamesSym)))	\
	return ScalarLogical(0);					\
    int i, j, k, kend, *pp_,						\
	*pp = INTEGER(GET_SLOT(obj, Matrix_pSym)),			\
	*pi = INTEGER(GET_SLOT(obj, Matrix_ ## _I_ ## Sym));		\
    Calloc_or_Alloca_TO(pp_, n, int);					\
    Memcpy(pp_, pp, n);							\
    ++pp;								\
									\
    /* For all X[i,j] in "leading" triangle, */				\
    /* need that X[j,i] exists and X[j,i] == X[i,j] */			\
    if (R_has_slot(obj, Matrix_xSym)) {					\
	SEXP x = GET_SLOT(obj, Matrix_xSym);				\
	SEXPTYPE tx = TYPEOF(x);					\
	switch (tx) {							\
	case REALSXP:							\
	{								\
	    double *px = REAL(x);					\
	    CR_IS_SYMMETRIC_LOOP(					\
		ISNAN(px[pp_[i]])					\
		? !ISNAN(px[k])						\
		: (ISNAN(px[k]) || px[pp_[i]] != px[k]));		\
	    break;							\
	}								\
	case LGLSXP:							\
	{								\
	    int *px = LOGICAL(x);					\
	    CR_IS_SYMMETRIC_LOOP(					\
		px[pp_[i]] == NA_LOGICAL				\
		? (px[k] != NA_LOGICAL)					\
		: (px[k] == NA_LOGICAL || px[pp_[i]] != px[k]));	\
	    break;							\
	}								\
	case INTSXP:							\
	{								\
	    int *px = INTEGER(x);					\
	    CR_IS_SYMMETRIC_LOOP(					\
		px[pp_[i]] == NA_INTEGER				\
		? (px[k] != NA_INTEGER)					\
		: (px[k] == NA_INTEGER || px[pp_[i]] != px[k]));	\
	    break;							\
	}								\
	case CPLXSXP:							\
	{								\
	    Rcomplex *px = COMPLEX(x);					\
	    CR_IS_SYMMETRIC_LOOP(					\
		ISNAN(px[pp_[i]].r) || ISNAN(px[pp_[i]].i)		\
		? !(ISNAN(px[k].r) || ISNAN(px[k].i))			\
		: (ISNAN(px[k].r) || ISNAN(px[k].i) ||			\
		   px[pp_[i]].r != px[k].r || px[pp_[i]].i != px[k].i)); \
	    break;							\
	}								\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", tx, "[CR]sparse_is_symmetric"); \
	    break;							\
	}								\
    } else {								\
	CR_IS_SYMMETRIC_LOOP(0);					\
    }									\
									\
    /* Need upper, lower triangles to have same number of nonzero elements */ \
    for (j = 0; j < n; ++j)						\
	if (pp_[j] != pp[j])						\
	    return ScalarLogical(0);					\
									\
    Free_FROM(pp_, n);							\
    return ScalarLogical(1);						\
}

/* Csparse_is_symmetric() */
/* FIXME: not checking for real diagonal in complex case */
CR_IS_SYMMETRIC(C, i)
/* Rsparse_is_symmetric() */
/* FIXME: not checking for real diagonal in complex case */
CR_IS_SYMMETRIC(R, j)
