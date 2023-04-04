#include "dsyMatrix.h"

SEXP dsyMatrix_trf_(SEXP obj, int warn)
{
    SEXP val;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(val = get_factor(obj, "BunchKaufman"), &pid);
    if (!isNull(val)) {
	UNPROTECT(1);
	return val;
    }
    REPROTECT(val = NEW_OBJECT_OF_CLASS("BunchKaufman"), pid);

    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    SET_SLOT(val, Matrix_uploSym, uplo);
    
    if (n > 0) {
	R_xlen_t nn;
	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	    perm = PROTECT(allocVector(INTSXP, n)),
	    x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	    y = PROTECT(allocVector(REALSXP, nn = XLENGTH(x)));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	int *pperm = INTEGER(perm), lwork = -1, info;
	double *px = REAL(x), *py = REAL(y), tmp, *work;
	
#define DSYTRF_FINISH(_UL_)						\
	do {								\
	    Matrix_memset(py, 0, nn, sizeof(double));			\
	    F77_CALL(dlacpy)(&_UL_, pdim, pdim, px, pdim, py, pdim FCONE); \
	    F77_CALL(dsytrf)(&_UL_, pdim, py, pdim, pperm, &tmp, &lwork, \
			     &info FCONE);				\
	    lwork = (int) tmp;						\
	    Calloc_or_Alloca_TO(work, lwork, double);			\
	    F77_CALL(dsytrf)(&_UL_, pdim, py, pdim, pperm, work, &lwork, \
			     &info FCONE);				\
	    Free_FROM(work, lwork);					\
	    								\
	    if (info < 0)						\
		error(_("LAPACK '%s' gave error code %d"),		\
		      "dsytrf", info);					\
	    else if (info > 0 && warn > 0) {				\
		/* MJ: 'dsytrf' does not distinguish between singular, */ \
		/*     finite matrices and matrices containing NaN ... */ \
		/*     hence this message can mislead                  */ \
		if (warn > 1)						\
		    error  (_("LAPACK '%s': matrix is exactly singular, " \
			      "D[i,i]=0, i=%d"),			\
			    "dsytrf", info);				\
		else							\
		    warning(_("LAPACK '%s': matrix is exactly singular, " \
			      "D[i,i]=0, i=%d"),			\
			    "dsytrf", info);				\
	    }								\
	    								\
	    SET_SLOT(val, Matrix_DimSym, dim);				\
	    if (!isNull(dimnames))					\
		set_symmetrized_DimNames(val, dimnames, -1);		\
	    SET_SLOT(val, Matrix_permSym, perm);			\
	    SET_SLOT(val, Matrix_xSym, y);				\
	} while (0)

	DSYTRF_FINISH(ul);
	UNPROTECT(4);
    }
    
    set_factor(obj, "BunchKaufman", val);
    UNPROTECT(3);
    return val;
}

SEXP dsyMatrix_trf(SEXP obj, SEXP warn)
{
    return dsyMatrix_trf_(obj, asInteger(warn));
}

SEXP matrix_trf_(SEXP obj, int warn, char uplo)
{
    SEXP dim = PROTECT(getAttrib(obj, R_DimSymbol));
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("'matrix_trf()' requires a square matrix"));
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("BunchKaufman")),
	ul = PROTECT(mkString((uplo == 'U') ? "U" : "L"));
    SET_SLOT(val, Matrix_DimSym, ul);

    if (n > 0) {
	R_xlen_t nn = XLENGTH(obj);
	SEXP dimnames = PROTECT(getAttrib(obj, R_DimNamesSymbol)),
	    perm = PROTECT(allocVector(INTSXP, n)),
	    y = PROTECT(allocVector(REALSXP, nn));
	int *pperm = INTEGER(perm), lwork = -1, info;
	double *px = REAL(obj), *py = REAL(y), tmp, *work;
    
	DSYTRF_FINISH(uplo);

#undef DSYTRF_FINISH

	UNPROTECT(3);
    }
    
    UNPROTECT(3);
    return val;
}

SEXP matrix_trf(SEXP obj, SEXP warn, SEXP uplo)
{
    if (TYPEOF(obj) != REALSXP)
	ERROR_INVALID_TYPE("matrix", TYPEOF(obj), "matrix_trf");
    if (!isMatrix(obj))
	ERROR_INVALID_CLASS(obj, "matrix_trf");
    
    char ul = 'U';
    if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
	(uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
	((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
	error(_("invalid 'uplo' to 'matrix_trf()'; must be \"U\" or \"L\""));
    
    return matrix_trf_(obj, asInteger(warn), ul);
}

double get_norm_dsy(SEXP obj, const char *typstr)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *pdim = INTEGER(dim);
    double *px = REAL(x), norm, *work = NULL;
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    if (typstr[0] == 'I' || typstr[0] == 'O')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlansy)(typstr, ul, pdim, px, pdim, work FCONE FCONE);

    UNPROTECT(3);
    return norm;
}

SEXP dsyMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dsy(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dsyMatrix_rcond(SEXP obj)
{
    SEXP trf = PROTECT(dsyMatrix_trf_(obj, 2)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	perm = PROTECT(GET_SLOT(trf, Matrix_permSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));
    
    int *pdim = INTEGER(dim), *pperm = INTEGER(perm), info;
    double *px = REAL(x), norm = get_norm_dsy(obj, "O"), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    F77_CALL(dsycon)(ul, pdim, px, pdim, pperm, &norm, &rcond,
		     (double *) R_alloc((size_t) 2 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE);
    
    UNPROTECT(5);
    return ScalarReal(rcond);
}

SEXP dsyMatrix_determinant(SEXP obj, SEXP logarithm)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    SEXP res;
    if (n == 0) {
	int givelog = asLogical(logarithm), sign = 1;
	double modulus = (givelog) ? 0.0 : 1.0;
	res = as_det_obj(modulus, givelog, sign);
    } else {
	SEXP trf = PROTECT(dsyMatrix_trf_(obj, 0));
	res = BunchKaufman_determinant(trf, logarithm);
	UNPROTECT(1); /* trf */
    }
    return res;
}

SEXP dsyMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dsyMatrix")),
	trf = PROTECT(dsyMatrix_trf_(a, 2)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(trf, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	perm = PROTECT(GET_SLOT(trf, Matrix_permSym)),
	x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(trf, Matrix_xSym), &pid);
    REPROTECT(x = duplicate(x), pid);
    
    SET_SLOT(val, Matrix_DimSym, dim);
    SET_SLOT(val, Matrix_DimNamesSym, dimnames);
    SET_SLOT(val, Matrix_uploSym, uplo);
    SET_SLOT(val, Matrix_xSym, x);
    
    int *pdim = INTEGER(dim), *pperm = INTEGER(perm), info;
    double *px = REAL(x);
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    F77_CALL(dsytri)(ul, pdim, px, pdim, pperm,
		     (double *) R_alloc((size_t) pdim[0], sizeof(double)),
		     &info FCONE);
    
    UNPROTECT(7);
    return val;
}

SEXP dsyMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0)),
	adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	bdim = PROTECT(GET_SLOT(val, Matrix_DimSym));
    int *padim = INTEGER(adim), *pbdim = INTEGER(bdim);
    
    if (padim[0] != pbdim[0] || padim[0] < 1 || pbdim[1] < 1)
	error(_("dimensions of system to be solved are inconsistent"));
    
    SEXP trf = PROTECT(dsyMatrix_trf_(a, 2)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	perm = PROTECT(GET_SLOT(trf, Matrix_permSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym)),
	y = PROTECT(GET_SLOT(val, Matrix_xSym));
    
    int *pperm = INTEGER(perm), info;
    double *px = REAL(x), *py = REAL(y);
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    F77_CALL(dsytrs)(ul, pbdim, pbdim + 1, px, pbdim, pperm, py, pbdim,
		     &info FCONE);

    UNPROTECT(8);
    return val;
}

SEXP dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rtP)
{
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0));// incl. dimnames
    int rt = asLogical(rtP); /* if(rt), compute b %*% a,  else  a %*% b */
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym)),
	m = bdims[0], n = bdims[1];

    if ((rt && n != adims[0]) || (!rt && m != adims[0]))
	error(_("Matrices are not conformable for multiplication"));

    double one = 1., zero = 0.;
    R_xlen_t mn = m * (R_xlen_t)n;
    double *bcp, *vx = REAL(GET_SLOT(val, Matrix_xSym));
    Calloc_or_Alloca_TO(bcp, mn, double);
    Memcpy(bcp, vx, mn);

    if (m >=1 && n >= 1)
	F77_CALL(dsymm)(rt ? "R" :"L", uplo_P(a), &m, &n, &one,
			REAL(GET_SLOT(a, Matrix_xSym)), adims, bcp,
			&m, &zero, vx, &m FCONE FCONE);
    // add dimnames:
    int nd = rt ?
	1 : // v <- b %*% a : rownames(v) == rownames(b)  are already there
	0;  // v <- a %*% b : colnames(v) == colnames(b)  are already there
    SEXP nms = PROTECT(VECTOR_ELT(get_symmetrized_DimNames(a, -1), nd));
    SET_VECTOR_ELT(GET_SLOT(val, Matrix_DimNamesSym), nd, nms);
    Free_FROM(bcp, mn);
    UNPROTECT(2);
    return val;
}

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0

// this is very close to lsyMatrix_as_lsp*() in ./ldense.c  -- keep synced !
SEXP dsyMatrix_as_dspMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dspMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_pack(
	REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, (n*(n+1))/2)),
	REAL( GET_SLOT(from, Matrix_xSym)),
	n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	NUN);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    SET_SLOT(val, Matrix_factorSym,
	     duplicate(GET_SLOT(from, Matrix_factorSym)));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0

SEXP dsyMatrix_as_matrix(SEXP from, SEXP keep_dimnames)
{
    int n = INTEGER(GET_SLOT(from, Matrix_DimSym))[0];
    SEXP val = PROTECT(allocMatrix(REALSXP, n, n));
    R_xlen_t nsqr = n; nsqr *= n;

    ddense_unpacked_make_symmetric(Memcpy(REAL(val),
					  REAL(GET_SLOT(from, Matrix_xSym)),
					  nsqr),
				   from);
    if(asLogical(keep_dimnames))
	setAttrib(val, R_DimNamesSymbol, get_symmetrized_DimNames(from, -1));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */


