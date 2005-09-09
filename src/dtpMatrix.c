/* double (precision) Triangular Packed Matrices
 * Note: this means *square* {n x n} matrices
*/

#include "dtpMatrix.h"

SEXP dtpMatrix_validate(SEXP obj)
{
    SEXP val = triangularMatrix_validate(obj);
    if(isString(val))
	return(val);
    else {
	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
	if (dims[0] != packed_ncol(length(GET_SLOT(obj, Matrix_xSym))))
	    return(mkString(_("Incorrect length of 'x' slot")));
	return ScalarLogical(1);
    }
}

static
double get_norm(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    double *work = (double *) NULL;

    typnm[0] = norm_type(typstr);
    if (*typnm == 'I') {
	work = (double *) R_alloc(dims[0], sizeof(double));
    }
    return F77_CALL(dlantp)(typnm,
			    CHAR(asChar(GET_SLOT(obj, Matrix_uploSym))),
			    CHAR(asChar(GET_SLOT(obj, Matrix_diagSym))),
			    dims, REAL(GET_SLOT(obj, Matrix_xSym)), work);
}

SEXP dtpMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm(obj, CHAR(asChar(type))));
}

static
double set_rcond(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    SEXP rcv = GET_SLOT(obj, Matrix_rcondSym);
    double rcond = get_double_by_name(rcv, typnm);

    typnm[0] = rcond_type(typstr);
    if (R_IsNA(rcond)) {
	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), info;
	F77_CALL(dtpcon)(typnm,
			 CHAR(asChar(GET_SLOT(obj, Matrix_uploSym))),
			 CHAR(asChar(GET_SLOT(obj, Matrix_diagSym))),
			 dims, REAL(GET_SLOT(obj, Matrix_xSym)),
			 &rcond,
			 (double *) R_alloc(3*dims[0], sizeof(double)),
			 (int *) R_alloc(dims[0], sizeof(int)), &info);
	SET_SLOT(obj, Matrix_rcondSym,
		 set_double_by_name(rcv, rcond, typnm));
    }
    return rcond;
}

SEXP dtpMatrix_rcond(SEXP obj, SEXP type)
{
    return ScalarReal(set_rcond(obj, CHAR(asChar(type))));
}

SEXP dtpMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(duplicate(a));
    int info, *Dim = INTEGER(GET_SLOT(val, Matrix_DimSym));
    F77_CALL(dtptri)(CHAR(asChar(GET_SLOT(val, Matrix_uploSym))),
		     CHAR(asChar(GET_SLOT(val, Matrix_diagSym))),
		     Dim, REAL(GET_SLOT(val, Matrix_xSym)), &info);
    UNPROTECT(1);
    return val;
}

SEXP dtpMatrix_getDiag(SEXP x)
{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(REALSXP, n));

    if (*CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0)) == 'U') {
	int j;
	for (j = 0; j < n; j++) REAL(val)[j] = 1.;
    } else {
	packed_getDiag(REAL(val), x);
    }
    UNPROTECT(1);
    return val;
}

SEXP dtpMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP val = PROTECT(duplicate(b));
    int *Dim = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bDim = INTEGER(getAttrib(val, R_DimSymbol));
    char *uplo = CHAR(STRING_ELT(GET_SLOT(a, Matrix_uploSym), 0)),
	*diag = CHAR(STRING_ELT(GET_SLOT(a, Matrix_diagSym), 0));
    double *ax = REAL(GET_SLOT(a, Matrix_xSym));
    int ione = 1, j;

    if (bDim[0] != Dim[1])
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      Dim[0], Dim[1], bDim[0], bDim[1]);
    for (j = 0; j < bDim[1]; j++)
	F77_CALL(dtpsv)(uplo, "N", diag, bDim, ax,
			REAL(val) + j * bDim[0], &ione);
    UNPROTECT(1);
    return val;
}

SEXP dtpMatrix_dgeMatrix_mm(SEXP x, SEXP y)
{
    SEXP val = PROTECT(duplicate(y));
    /* Since 'x' is square (n x n ),   dim(x %*% y) = dim(y) */
    int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDim = INTEGER(GET_SLOT(y, Matrix_DimSym));
    int ione = 1, j;
    char *uplo = CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0)),
	*diag = CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)),
	*vx = REAL(GET_SLOT(val, Matrix_xSym));

    if (yDim[0] != xDim[1])
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      xDim[0], xDim[1], yDim[0], yDim[1]);
    for (j = 0; j < yDim[1]; j++) /* X %*% y[,j]  via BLAS 2 DTPMV(.) */
	F77_CALL(dtpmv)(uplo, "N", diag, yDim, xx,
			vx + j * yDim[0], &ione);
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_dtpMatrix_mm(SEXP x, SEXP y)
{
    SEXP val = PROTECT(duplicate(x));
    /* Since 'y' is square (n x n ),   dim(x %*% y) = dim(x) */
    int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDim = INTEGER(GET_SLOT(y, Matrix_DimSym));
    int i;
    char *uplo = CHAR(STRING_ELT(GET_SLOT(y, Matrix_uploSym), 0)),
	 *diag = CHAR(STRING_ELT(GET_SLOT(y, Matrix_diagSym), 0));
    double *yx = REAL(GET_SLOT(y, Matrix_xSym)),
 	*vx = REAL(GET_SLOT(val, Matrix_xSym));

    if (yDim[0] != xDim[1])
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      xDim[0], xDim[1], yDim[0], yDim[1]);
    for (i = 0; i < xDim[0]; i++)/* val[i,] := Y' %*% x[i,]  */
	F77_CALL(dtpmv)(uplo, "T", diag, yDim, yx,
			vx + i, /* incr = */ xDim);
    UNPROTECT(1);
    return val;
}

SEXP dtpMatrix_matrix_mm(SEXP x, SEXP y)
{
    SEXP val = PROTECT(duplicate(y));
    /* Since 'x' is square (n x n ),   dim(x %*% y) = dim(y) */
    int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDim = INTEGER(getAttrib(y, R_DimSymbol));
    int ione = 1, j;
    char *uplo = CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0)),
	*diag = CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

    if (yDim[0] != xDim[1])
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      xDim[0], xDim[1], yDim[0], yDim[1]);
    for (j = 0; j < yDim[1]; j++)
	F77_CALL(dtpmv)(uplo, "N", diag, yDim, xx,
			REAL(val) + j * yDim[0], &ione);
    UNPROTECT(1);
    return val;
}

SEXP dtpMatrix_as_dtrMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_rcondSym,
	     duplicate(GET_SLOT(from, Matrix_rcondSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    packed_to_full(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n*n)),
		   REAL(GET_SLOT(from, Matrix_xSym)), n,
		   *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW);
    UNPROTECT(1);
    return val;
}

