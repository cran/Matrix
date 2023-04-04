/* double (precision) Triangular Packed Matrices
 * Note: this means *square* {n x n} matrices
*/

#include "dtpMatrix.h"

double get_norm_dtp(SEXP obj, const char *typstr)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(obj, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *pdim = INTEGER(dim);
    double *px = REAL(x), norm, *work = NULL;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
    if (typstr[0] == 'I')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlantp)(typstr, ul, di, pdim, px,
			    work FCONE FCONE FCONE);

    UNPROTECT(4);
    return norm;
}

SEXP dtpMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dtp(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dtpMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(obj, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_rcond_type(CHAR(type));
    
    int *pdim = INTEGER(dim), info;
    double *px = REAL(x), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));

    F77_CALL(dtpcon)(typstr, ul, di, pdim, px, &rcond,
		     (double *) R_alloc((size_t) 3 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE FCONE FCONE);

    UNPROTECT(5);
    return ScalarReal(rcond);
}

SEXP dtpMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dtpMatrix")),
	dim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(a, Matrix_diagSym)),
	x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(a, Matrix_xSym), &pid);
    REPROTECT(x = duplicate(x), pid);
    
    SET_SLOT(val, Matrix_DimSym, dim);
    set_reversed_DimNames(val, dimnames);
    SET_SLOT(val, Matrix_uploSym, uplo);
    SET_SLOT(val, Matrix_diagSym, diag);
    SET_SLOT(val, Matrix_xSym, x);
    
    int *pdim = INTEGER(dim), info;
    double *px = REAL(x);
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
    F77_CALL(dtptri)(ul, di, pdim, px, &info FCONE FCONE);
    
    UNPROTECT(6);
    return val;
}

SEXP dtpMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0)),
	adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	bdim = PROTECT(GET_SLOT(val, Matrix_DimSym));
    int *padim = INTEGER(adim), *pbdim = INTEGER(bdim);
    
    if (padim[0] != pbdim[0] || padim[0] < 1 || pbdim[1] < 1)
	error(_("dimensions of system to be solved are inconsistent"));

    SEXP uplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(a, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(a, Matrix_xSym)),
	y = PROTECT(GET_SLOT(val, Matrix_xSym));
    
    int one = 1;
    double *px = REAL(x), *py = REAL(y);
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
#ifdef PRE_2013_08_30
    /* a^{-1} %*% b[, j] via BLAS 2 DTPSV */
    int j, n = pbdim[1];
    for (j = 0; j < n; ++j)
	F77_CALL(dtpsv)(ul, "N", di, pbdim, px, py + (R_xlen_t) j * pbdim[0],
			&one FCONE FCONE);
#else
    F77_CALL(dtptrs)(ul, "N", di, pbdim, pbdim + 1, px, py, pbdim,
		     &one FCONE FCONE);
#endif
    
    UNPROTECT(7);
    return val;
}

SEXP dtpMatrix_matrix_mm(SEXP x, SEXP y, SEXP right, SEXP trans)
{
    SEXP val = PROTECT(dense_as_general(y, 'd', 2, 0));
    int rt = asLogical(right); // if(rt), compute b %*% op(a), else op(a) %*% b
    int tr = asLogical(trans); // if(tr), op(a) = t(a), else op(a) = a
    /* Since 'x' is square (n x n ),   dim(x %*% y) = dim(y) */
    int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDim = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int m = yDim[0], n = yDim[1];
    int ione = 1;
    const char *uplo = uplo_P(x), *diag = diag_P(x);
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)),
	*vx = REAL(GET_SLOT(val, Matrix_xSym));

    if (yDim[0] != xDim[1])
    if ((rt && xDim[0] != n) || (!rt && xDim[1] != m))
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      xDim[0], xDim[1], yDim[0], yDim[1]);
    if (m < 1 || n < 1) {
/* 	error(_("Matrices with zero extents cannot be multiplied")); */
    } else /* BLAS */
	// go via BLAS 2  dtpmv(.); there is no dtpmm in Lapack!
	if(rt) {
	    error(_("right=TRUE is not yet implemented __ FIXME"));
	} else {
	    for (int j = 0; j < n; j++) // X %*% y[,j]
		F77_CALL(dtpmv)(uplo, /*trans = */ tr ? "T" : "N",
				diag, yDim, xx,
				vx + j * (size_t) m, &ione FCONE FCONE FCONE);
	}
    UNPROTECT(1);
    return val;
}

/* FIXME: This function should be removed and a rt argument added to
 * dtpMatrix_matrix_mm -- also to be more parallel to ./dtrMatrix.c code */
SEXP dgeMatrix_dtpMatrix_mm(SEXP x, SEXP y)
{
    SEXP val = PROTECT(duplicate(x));
    /* Since 'y' is square (n x n ),   dim(x %*% y) = dim(x) */
    int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDim = INTEGER(GET_SLOT(y, Matrix_DimSym));
    const char *uplo = uplo_P(y), *diag = diag_P(y);
    double *yx = REAL(GET_SLOT(y, Matrix_xSym)),
 	*vx = REAL(GET_SLOT(val, Matrix_xSym));

    if (yDim[0] != xDim[1])
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      xDim[0], xDim[1], yDim[0], yDim[1]);
    for (int i = 0; i < xDim[0]; i++)/* val[i,] := Y' %*% x[i,]  */
	F77_CALL(dtpmv)(uplo, "T", diag, yDim, yx,
			vx + i, /* incr = */ xDim FCONE FCONE FCONE);
    UNPROTECT(1);
    return val;
}

/* MJ: no longer needed ... prefer more general packedMatrix_diag_[gs]et() */
#if 0

// also applicable to dspMatrix , dppMatrix :
SEXP dtpMatrix_getDiag(SEXP x)
{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(REALSXP, n));

    tr_d_packed_getDiag(REAL(val), x, n);
    UNPROTECT(1);
    return val;
}

// also applicable to lspMatrix :
SEXP ltpMatrix_getDiag(SEXP x)
{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(LGLSXP, n));

    tr_l_packed_getDiag(LOGICAL(val), x, n);
    UNPROTECT(1);
    return val;
}

SEXP dtpMatrix_setDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return tr_d_packed_setDiag(REAL(d), LENGTH(d), x, n);
}

SEXP ltpMatrix_setDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return tr_l_packed_setDiag(INTEGER(d), LENGTH(d), x, n);
}

/* was unused, not replaced: */
SEXP dtpMatrix_addDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return tr_d_packed_addDiag(REAL(d), LENGTH(d), x, n);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general packedMatrix_unpack() */
#if 0

SEXP dtpMatrix_as_dtrMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_unpack(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n*n)),
		  REAL(GET_SLOT(from, Matrix_xSym)),
		  n,
		  *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
		  *CHAR(STRING_ELT(diag, 0)) == 'N' ? NUN : UNT);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */

