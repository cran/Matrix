#include "dppMatrix.h"

SEXP dppMatrix_trf_(SEXP obj, int warn)
{
    SEXP val;
    PROTECT_INDEX pidA;
    PROTECT_WITH_INDEX(val = get_factor(obj, "pCholesky"), &pidA);
    if (!isNull(val)) {
	UNPROTECT(1);
	return val;
    }
    REPROTECT(val = NEW_OBJECT_OF_CLASS("pCholesky"), pidA);

    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    SET_SLOT(val, Matrix_uploSym, uplo);

    if (n > 0) {
	PROTECT_INDEX pidB;
	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)), x;
	PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pidB);
	REPROTECT(x = duplicate(x), pidB);
	char ul = *CHAR(STRING_ELT(uplo, 0));
	int info;
	double *px = REAL(x);
	
	F77_CALL(dpptrf)(&ul, pdim, px, &info FCONE);

	if (info < 0)
	    error(_("LAPACK '%s' gave error code %d"),
		  "dpptrf", info);
	else if (info > 0) {
	    if (warn > 1)
		error  (_("LAPACK '%s': leading minor of order %d is not "
			  "positive definite"),
			"dpptrf", info);
	    else if (warn > 0)
		warning(_("LAPACK '%s': leading minor of order %d is not "
			  "positive definite"),
			"dpptrf", info);
	    UNPROTECT(5);
	    return ScalarInteger(info);
	}
	
	SET_SLOT(val, Matrix_DimSym, dim);
	set_symmetrized_DimNames(val, dimnames, -1);
	SET_SLOT(val, Matrix_xSym, x);
	UNPROTECT(2);
    }

    set_factor(obj, "pCholesky", val);
    UNPROTECT(3);
    return val;
}

SEXP dppMatrix_trf(SEXP obj, SEXP warn)
{
    return dppMatrix_trf_(obj, asInteger(warn));
}

SEXP dppMatrix_rcond(SEXP obj)
{
    SEXP trf = PROTECT(dppMatrix_trf_(obj, 2)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));

    int *pdim = INTEGER(dim), info;
    double *px = REAL(x), norm = get_norm_dsp(obj, "O"), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0));

    F77_CALL(dppcon)(ul, pdim, px, &norm, &rcond,
		     (double *) R_alloc(3 * pdim[0], sizeof(double)),
		     (int *) R_alloc(pdim[0], sizeof(int)),
		     &info FCONE);

    UNPROTECT(4);
    return ScalarReal(rcond);
}

SEXP dppMatrix_solve(SEXP a)
{
    PROTECT_INDEX pid;
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dppMatrix")),
	trf = PROTECT(dppMatrix_trf_(a, 2)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(trf, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	x;
    PROTECT_WITH_INDEX(x = GET_SLOT(trf, Matrix_xSym), &pid);
    REPROTECT(x = duplicate(x), pid);
    
    SET_SLOT(val, Matrix_DimSym, dim);
    SET_SLOT(val, Matrix_DimNamesSym, dimnames);
    SET_SLOT(val, Matrix_xSym, x);
    SET_SLOT(val, Matrix_uploSym, uplo);

    int *pdim = INTEGER(dim), info;
    double *px = REAL(x);
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    F77_CALL(dpptri)(ul, pdim, px, &info FCONE);

    UNPROTECT(6);
    return val;
}

SEXP dppMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0)),
	adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	bdim = PROTECT(GET_SLOT(val, Matrix_DimSym));
    int *padim = INTEGER(adim), *pbdim = INTEGER(bdim);
    
    if (padim[0] != pbdim[0] || padim[0] < 1 || pbdim[1] < 1)
	error(_("dimensions of system to be solved are inconsistent"));
    
    SEXP trf = PROTECT(dppMatrix_trf_(a, 2)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym)),
	y = PROTECT(GET_SLOT(val, Matrix_xSym));
    
    int info;
    double *px = REAL(x), *py = REAL(y);
    const char *ul = CHAR(STRING_ELT(uplo, 0));

    F77_CALL(dpptrs)(ul, pbdim, pbdim + 1, px, py, pbdim, &info FCONE);
    
    UNPROTECT(7);
    return val;
}
