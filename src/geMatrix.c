#include "geMatrix.h"

SEXP geMatrix_validate(SEXP obj)
{
    SEXP x = GET_SLOT(obj, Matrix_xSym),
	Dim = GET_SLOT(obj, Matrix_DimSym);
    int m, n;

    if (length(Dim) != 2)
	return ScalarString(mkChar("Dim slot must have length 2"));
    m = INTEGER(Dim)[0]; n = INTEGER(Dim)[1];
    if (m < 0 || n < 0)
	return
	    ScalarString(mkChar("Negative value(s) in Dim"));
    if (length(x) != m * n)
    	return
	    ScalarString(mkChar("x slot must have length that is prod(Dim)"));
    return ScalarLogical(1);
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
    return F77_CALL(dlange)(typstr, dims, dims+1,
			    REAL(GET_SLOT(obj, Matrix_xSym)),
			    dims, work);
}

SEXP geMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm(obj, CHAR(asChar(type))));
}

static
double set_rcond(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    SEXP rcv = GET_SLOT(obj, install("rcond"));
    double rcond = get_double_by_name(rcv, typnm);

    typnm[0] = rcond_type(typstr);
    if (R_IsNA(rcond)) {
        SEXP LU = geMatrix_LU(obj);
	int *dims = INTEGER(GET_SLOT(LU, Matrix_DimSym)), info;
	double anorm = get_norm(obj, typstr);

	if (dims[0] != dims[1] || dims[0] < 1)
            error("rcond requires a square, non-empty matrix");
	F77_CALL(dgecon)(typnm,
			 dims, REAL(GET_SLOT(LU, Matrix_xSym)),
			 dims, &anorm, &rcond,
			 (double *) R_alloc(4*dims[0], sizeof(double)),
			 (int *) R_alloc(dims[0], sizeof(int)), &info);
	SET_SLOT(obj, install("rcond"),
		 set_double_by_name(rcv, rcond, typnm));
    }
    return rcond;
}

SEXP geMatrix_rcond(SEXP obj, SEXP type)
{
  return ScalarReal(set_rcond(obj, CHAR(asChar(type))));
}

SEXP geMatrix_crossprod(SEXP x)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("poMatrix")));
    int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int n = Dims[1];
    double one = 1.0, zero = 0.0;

    if ((*Dims) > 0 && n > 0) {
	vDims[0] = vDims[1] = n;
	SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, n * n));
	F77_CALL(dsyrk)(CHAR(asChar(GET_SLOT(val, Matrix_uploSym))),
		       "T", Dims + 1, Dims, &one,
		       REAL(GET_SLOT(x, Matrix_xSym)),
		       Dims, &zero,
		       REAL(GET_SLOT(val, Matrix_xSym)), &n);
    }
    UNPROTECT(1);
    return val;
}

SEXP geMatrix_geMatrix_crossprod(SEXP x, SEXP y)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix")));
    int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDims = INTEGER(GET_SLOT(y, Matrix_DimSym)),
	*vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int m = xDims[1], n = yDims[1];
    double one = 1.0, zero = 0.0;

    if ((*xDims) > 0 && (*yDims) > 0 && n > 0 && m > 0) {
	if (*xDims != *yDims)
	    error("Dimensions of x and y are not compatible for crossprod");
	vDims[0] = m; vDims[1] = n;
	SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
	F77_CALL(dgemm)("T", "N", xDims + 1, yDims + 1, xDims, &one,
			REAL(GET_SLOT(x, Matrix_xSym)), xDims,
			REAL(GET_SLOT(y, Matrix_xSym)), yDims,
			&zero, REAL(GET_SLOT(val, Matrix_xSym)),
			xDims + 1);
    }
    UNPROTECT(1);
    return val;
}

SEXP geMatrix_matrix_crossprod(SEXP x, SEXP y)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix")));
    int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDims = INTEGER(getAttrib(y, R_DimSymbol)),
	*vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int m = xDims[1], n = yDims[1];
    double one = 1.0, zero = 0.0;

    if (!(isMatrix(y) && isReal(y)))
	error("Argument y must be a numeric (real) matrix");
    if ((*xDims) > 0 && (*yDims) > 0 && n > 0 && m > 0) {
	if (*xDims != *yDims)
	    error("Dimensions of x and y are not compatible for crossprod");
	vDims[0] = m; vDims[1] = n;
	SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
	F77_CALL(dgemm)("T", "N", xDims + 1, yDims + 1, xDims, &one,
			REAL(GET_SLOT(x, Matrix_xSym)), xDims,
			REAL(y), yDims,
			&zero, REAL(GET_SLOT(val, Matrix_xSym)),
			xDims + 1);
    }
    UNPROTECT(1);
    return val;
}

SEXP geMatrix_getDiag(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int i, m = dims[0], nret = (m < dims[1]) ? m : dims[1];
    SEXP ret = PROTECT(allocVector(REALSXP, nret)),
	xv = GET_SLOT(x, Matrix_xSym);

    for (i = 0; i < nret; i++) {
	REAL(ret)[i] = REAL(xv)[i * (m + 1)];
    }
    UNPROTECT(1);
    return ret;
}

SEXP geMatrix_LU(SEXP x)
{
    SEXP val = get_factorization(x, "LU");
    int *dims, npiv, info;
    
    if (val != R_NilValue) return val;
    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    if (dims[0] < 1 || dims[1] < 1)
	error("Cannot factor a matrix with zero extents");
    npiv = (dims[0] <dims[1]) ? dims[0] : dims[1];
    val = PROTECT(NEW_OBJECT(MAKE_CLASS("LU")));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(x, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    SET_SLOT(val, install("pivot"), allocVector(INTSXP, npiv));
    F77_CALL(dgetrf)(dims, dims + 1, REAL(GET_SLOT(val, Matrix_xSym)),
		     dims, INTEGER(GET_SLOT(val, install("pivot"))),
		     &info);
    if (info) error("Lapack routine dgetrf returned error code %d", info);
    UNPROTECT(1);
    return set_factorization(x, val, "LU");
}

SEXP geMatrix_determinant(SEXP x, SEXP logarithm)
{
    int lg = asLogical(logarithm);
    SEXP lu = geMatrix_LU(x);
    int *dims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
	*jpvt = INTEGER(GET_SLOT(lu, install("pivot"))),
	i, n = dims[0], sign = 1;
    double *luvals = REAL(GET_SLOT(lu, Matrix_xSym)), modulus;

    if (n != dims[1])
	error("Determinant requires a square matrix");
    for (i = 0; i < n; i++) if (jpvt[i] != (i + 1)) sign = -sign;
    if (lg) {
	modulus = 0.0;
	for (i = 0; i < n; i++) {
	    double dii = luvals[i*(n + 1)]; /* ith diagonal element */
	    modulus += log(dii < 0 ? -dii : dii);
	    if (dii < 0) sign = -sign;
	}
    } else {
	modulus = 1.0;
	for (i = 0; i < n; i++)
	    modulus *= luvals[i*(n + 1)];
	if (modulus < 0) {
	    modulus = -modulus;
	    sign = -sign;
	}
    }
    return as_det_obj(modulus, lg, sign);
}

SEXP geMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix"))),
	lu = geMatrix_LU(a);
    int *dims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
	*pivot = INTEGER(GET_SLOT(lu, install("pivot")));
    double *x, tmp;
    int	info, lwork = -1;


    if (dims[0] != dims[1]) error("Solve requires a square matrix");
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(lu, Matrix_xSym)));
    x = REAL(GET_SLOT(val, Matrix_xSym));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(lu, Matrix_DimSym)));
    F77_CALL(dgetri)(dims, x, dims, pivot, &tmp, &lwork, &info);
    lwork = (int) tmp;
    F77_CALL(dgetri)(dims, x, dims, pivot,
		     (double *) R_alloc((size_t) lwork, sizeof(double)),
		     &lwork, &info);
    UNPROTECT(1);
    return val;
}

SEXP geMatrix_geMatrix_mm(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	*cdims,
	m = adims[0], n = bdims[1], k = adims[1];
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix")));
    char *trans = "N";
    double one = 1., zero = 0.;
    
    if (bdims[0] != k)
	error("Matrices are not conformable for multiplication");
    if (m < 1 || n < 1 || k < 1)
	error("Matrices with zero extents cannot be multiplied");
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    cdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    cdims[0] = m; cdims[1] = n;
    F77_CALL(dgemm)(trans, trans, adims, bdims+1, bdims, &one,
		    REAL(GET_SLOT(a, Matrix_xSym)), adims,
		    REAL(GET_SLOT(b, Matrix_xSym)), bdims,
		    &zero, REAL(GET_SLOT(val, Matrix_xSym)), adims);
    UNPROTECT(1);
    return val;
}

SEXP geMatrix_svd(SEXP x, SEXP nnu, SEXP nnv)
{
    int nu = asInteger(nnu),
	nv = asInteger(nnv),
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));
    SEXP val = PROTECT(allocVector(VECSXP, 3));

    if (dims[0] && dims[1]) {
	int m = dims[0], n = dims[1], mm = (m < n)?m:n,
	    lwork = -1, info;
	int *iwork = Calloc(8 * mm, int);
	double tmp, *work;
/* 	int bdspac = 3*m*m + 4*m, */
/* 	    wrkbl, maxwrk, minwrk, itmp, */
/* 	    ione = 1, iminus1 = -1; */
/* 	int i1, i2, i3; */

	SET_VECTOR_ELT(val, 0, allocVector(REALSXP, mm));
	SET_VECTOR_ELT(val, 1, allocMatrix(REALSXP, m, mm));
	SET_VECTOR_ELT(val, 2, allocMatrix(REALSXP, mm, n));
	F77_CALL(dgesdd)("S", &m, &n, xx, &m,
			 REAL(VECTOR_ELT(val, 0)),
			 REAL(VECTOR_ELT(val, 1)), &m,
			 REAL(VECTOR_ELT(val, 2)), &mm,
			 &tmp, &lwork, iwork, &info);
	lwork = (int) tmp;
/* 	F77_CALL(foo)(&i1, &i2, &i3); */
/* 	wrkbl = 3*m+(m+n)*i1; */
/* 	if (wrkbl < (itmp = 3*m + m*i2)) wrkbl = itmp; */
/* 	if (wrkbl < (itmp = 3*m + m*i3)) wrkbl = itmp; */
/* 	itmp = bdspac+3*m; */
/* 	maxwrk = (wrkbl > itmp) ? wrkbl : itmp; */
/* 	minwrk = 3*m + ((bdspac > n) ?  bdspac : n); */
	work = Calloc(lwork, double);
	F77_CALL(dgesdd)("S", &m, &n, xx, &m,
			 REAL(VECTOR_ELT(val, 0)),
			 REAL(VECTOR_ELT(val, 1)), &m,
			 REAL(VECTOR_ELT(val, 2)), &mm,
			 work, &lwork, iwork, &info);
	Free(iwork); Free(work);
    }
    UNPROTECT(1);
    return val;
}
