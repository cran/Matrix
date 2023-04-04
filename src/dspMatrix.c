#include "dspMatrix.h"

SEXP dspMatrix_trf_(SEXP obj, int warn)
{
    SEXP val;
    PROTECT_INDEX pidA;
    PROTECT_WITH_INDEX(val = get_factor(obj, "pBunchKaufman"), &pidA);
    if (!isNull(val)) {
	UNPROTECT(1);
	return val;
    }
    REPROTECT(val = NEW_OBJECT_OF_CLASS("pBunchKaufman"), pidA);

    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    SET_SLOT(val, Matrix_uploSym, uplo);
    
    if (n > 0) {
	PROTECT_INDEX pidB;
	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	    perm = PROTECT(allocVector(INTSXP, n)), x;
	PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pidB);
	REPROTECT(x = duplicate(x), pidB);
	char ul = *CHAR(STRING_ELT(uplo, 0));
	int *pperm = INTEGER(perm), info;
	double *px = REAL(x);
    
	F77_CALL(dsptrf)(&ul, pdim, px, pperm, &info FCONE);
    
	if (info < 0)
	    error(_("LAPACK '%s' gave error code %d"),
		  "dsptrf", info);
	else if (info > 0 && warn > 0) {
	    /* MJ: 'dsptrf' does not distinguish between singular, */
	    /*     finite matrices and matrices containing NaN ... */
	    /*     hence this message can mislead                  */
	    if (warn > 1)
		error  (_("LAPACK '%s': matrix is exactly singular, "
			  "D[i,i]=0, i=%d"),
			"dsptrf", info);
	    else
		warning(_("LAPACK '%s': matrix is exactly singular, "
			  "D[i,i]=0, i=%d"),
			"dsptrf", info);
	}

	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_permSym, perm);
	SET_SLOT(val, Matrix_xSym, x);
	UNPROTECT(3);
    }
    
    set_factor(obj, "pBunchKaufman", val);
    UNPROTECT(3);
    return val;
}

SEXP dspMatrix_trf(SEXP obj, SEXP warn)
{
    return dspMatrix_trf_(obj, asInteger(warn));
}

double get_norm_dsp(SEXP obj, const char *typstr)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *pdim = INTEGER(dim);
    double *px = REAL(x), norm, *work = NULL;
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    if (typstr[0] == 'I' || typstr[0] == 'O')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlansp)(typstr, ul, pdim, px, work FCONE FCONE);

    UNPROTECT(3);
    return norm;
}

SEXP dspMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dsp(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dspMatrix_rcond(SEXP obj)
{
    SEXP trf = PROTECT(dspMatrix_trf_(obj, 2)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	perm = PROTECT(GET_SLOT(trf, Matrix_permSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));
    
    int *pdim = INTEGER(dim), *pperm = INTEGER(perm), info;
    double *px = REAL(x), norm = get_norm_dsp(obj, "O"), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0));

    F77_CALL(dspcon)(ul, pdim, px, pperm, &norm, &rcond, 
		     (double *) R_alloc((size_t) 2 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE);
    
    UNPROTECT(5);
    return ScalarReal(rcond);
}

SEXP dspMatrix_determinant(SEXP obj, SEXP logarithm)
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
	SEXP trf = PROTECT(dspMatrix_trf_(obj, 0));
	res = BunchKaufman_determinant(trf, logarithm);
	UNPROTECT(1); /* trf */
    }
    return res;
}

SEXP dspMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dspMatrix")),
	trf = PROTECT(dspMatrix_trf_(a, 2)),
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

    F77_CALL(dsptri)(ul, pdim, px, pperm, 
		     (double *) R_alloc((size_t) pdim[0], sizeof(double)),
		     &info FCONE);
    
    UNPROTECT(7);
    return val;
}

SEXP dspMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0)),
	adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	bdim = PROTECT(GET_SLOT(val, Matrix_DimSym));
    int *padim = INTEGER(adim), *pbdim = INTEGER(bdim);
    
    if (padim[0] != pbdim[0] || padim[0] < 1 || pbdim[1] < 1)
	error(_("dimensions of system to be solved are inconsistent"));
    
    SEXP trf = PROTECT(dspMatrix_trf_(a, 2)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	perm = PROTECT(GET_SLOT(trf, Matrix_permSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym)),
	y = PROTECT(GET_SLOT(val, Matrix_xSym));
    
    int *pperm = INTEGER(perm), info;
    double *px = REAL(x), *py = REAL(y);
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    F77_CALL(dsptrs)(ul, pbdim, pbdim + 1, px, pperm, py, pbdim,
		     &info FCONE);
    
    UNPROTECT(8);
    return val;
}

SEXP dspMatrix_matrix_mm(SEXP a, SEXP b)
{
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0));
    int *bdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int i, ione = 1, n = bdims[0], nrhs = bdims[1];
    R_xlen_t nn = n * (R_xlen_t) nrhs;
    const char *uplo = uplo_P(a);
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)), one = 1., zero = 0.,
	*vx = REAL(GET_SLOT(val, Matrix_xSym)), *bx;

    Calloc_or_Alloca_TO(bx, nn, double);
    Memcpy(bx, vx, nn);
    if (bdims[0] != n)
	error(_("Matrices are not conformable for multiplication"));
    if (nrhs >= 1 && n >= 1) {
	R_xlen_t in;
	for (i = 0, in = 0; i < nrhs; i++, in += n) { // in := i * n (w/o overflow!)
	    F77_CALL(dspmv)(uplo, &n, &one, ax, bx + in, &ione,
			    &zero, vx + in, &ione FCONE);
	}
	Free_FROM(bx, nn);
    }
    UNPROTECT(1);
    return val;
}

/* MJ: no longer needed ... prefer more general packedMatrix_diag_[gs]et() */
#if 0

SEXP dspMatrix_getDiag(SEXP x)

{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(REALSXP, n));

    d_packed_getDiag(REAL(val), x, n);
    UNPROTECT(1);
    return val;
}

SEXP lspMatrix_getDiag(SEXP x)
{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(LGLSXP, n));

    l_packed_getDiag(LOGICAL(val), x, n);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_setDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return d_packed_setDiag(REAL(d), LENGTH(d), x, n);
}

SEXP lspMatrix_setDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return l_packed_setDiag(INTEGER(d), LENGTH(d), x, n);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general packedMatrix_unpack() */
#if 0

SEXP dspMatrix_as_dsyMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dsyMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_unpack(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n*n)),
		  REAL(GET_SLOT(from, Matrix_xSym)),
		  n,
		  *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
		  NUN);
    UNPROTECT(1);
    return val;
}

#endif /* MJ */
