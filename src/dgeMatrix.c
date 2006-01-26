#include "dgeMatrix.h"

SEXP dMatrix_validate(SEXP obj)
{
    SEXP x = GET_SLOT(obj, Matrix_xSym),
	Dim = GET_SLOT(obj, Matrix_DimSym);
    int m, n;

    if (length(Dim) != 2)
	return mkString(_("Dim slot must have length 2"));
    m = INTEGER(Dim)[0]; n = INTEGER(Dim)[1];
    if (m < 0 || n < 0)
	return mkString(_("Negative value(s) in Dim"));
    if (!isReal(x))
	return mkString(_("x slot must be numeric \"double\""));
    return ScalarLogical(1);
}

SEXP dgeMatrix_validate(SEXP obj)
{
    SEXP val,
	fact = GET_SLOT(obj, Matrix_factorSym);

    if (isString(val = dense_nonpacked_validate(obj)))
	return(val);

    if (length(fact) > 0 && getAttrib(fact, R_NamesSymbol) == R_NilValue)
	return mkString(_("factors slot must be named list"));
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

SEXP dgeMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm(obj, CHAR(asChar(type))));
}

SEXP dgeMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP LU = dgeMatrix_LU(obj);
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(LU, Matrix_DimSym)), info;
    double anorm, rcond;

    if (dims[0] != dims[1] || dims[0] < 1)
	error(_("rcond requires a square, non-empty matrix"));
    typnm[0] = rcond_type(CHAR(asChar(type)));
    anorm = get_norm(obj, typnm);
    F77_CALL(dgecon)(typnm,
		     dims, REAL(GET_SLOT(LU, Matrix_xSym)),
		     dims, &anorm, &rcond,
		     (double *) R_alloc(4*dims[0], sizeof(double)),
		     (int *) R_alloc(dims[0], sizeof(int)), &info);
    return ScalarReal(rcond);
}

SEXP dgeMatrix_crossprod(SEXP x, SEXP trans)
{
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x) */
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix"))),
	nms = VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym), tr ? 0 : 1),
	vDnms = ALLOC_SLOT(val, Matrix_DimNamesSym, VECSXP, 2);
    int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*vDims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));
    int k = tr ? Dims[1] : Dims[0], n = tr ? Dims[0] : Dims[1];
    double *vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n * n)),
	one = 1.0, zero = 0.0;

    AZERO(vx, n * n);
    SET_SLOT(val, Matrix_uploSym, mkString("U"));
    ALLOC_SLOT(val, Matrix_factorSym, VECSXP, 0);
    vDims[0] = vDims[1] = n;
    SET_VECTOR_ELT(vDnms, 0, duplicate(nms));
    SET_VECTOR_ELT(vDnms, 1, duplicate(nms));
    F77_CALL(dsyrk)("U", tr ? "N" : "T", &n, &k,
		    &one, REAL(GET_SLOT(x, Matrix_xSym)), Dims,
		    &zero, vx, &n);

    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x,y) */
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDims = INTEGER(GET_SLOT(y, Matrix_DimSym)),
	*vDims;
    int m  = xDims[!tr],  n = yDims[!tr];/* -> result dim */
    int xd = xDims[ tr], yd = yDims[ tr];/* the conformable dims */
    double one = 1.0, zero = 0.0;

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    if (xd > 0 && yd > 0 && n > 0 && m > 0) {
	if (xd != yd)
	    error(_("Dimensions of x and y are not compatible for %s"),
		  tr ? "tcrossprod" : "crossprod");
	vDims[0] = m; vDims[1] = n;
	SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
	F77_CALL(dgemm)(tr ? "N" : "T", tr ? "T" : "N", &m, &n, &xd, &one,
			REAL(GET_SLOT(x, Matrix_xSym)), xDims,
			REAL(GET_SLOT(y, Matrix_xSym)), yDims,
			&zero, REAL(GET_SLOT(val, Matrix_xSym)), &m);
    }
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x,y) */
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDims = INTEGER(getAttrib(y, R_DimSymbol)),
	*vDims;
    int m  = xDims[!tr],  n = yDims[!tr];/* -> result dim */
    int xd = xDims[ tr], yd = yDims[ tr];/* the conformable dims */
    double one = 1.0, zero = 0.0;

    if (!(isMatrix(y) && isReal(y)))
	error(_("Argument y must be a numeric (real) matrix"));
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    if (xd > 0 && yd > 0 && n > 0 && m > 0) {
	if (xd != yd)
	    error(_("Dimensions of x and y are not compatible for %s"),
		  tr ? "tcrossprod" : "crossprod");
	vDims[0] = m; vDims[1] = n;
	SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
	F77_CALL(dgemm)(tr ? "N" : "T", tr ? "T" : "N", &m, &n, &xd, &one,
			REAL(GET_SLOT(x, Matrix_xSym)), xDims,
			REAL(y), yDims,
			&zero, REAL(GET_SLOT(val, Matrix_xSym)), &m);
    }
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_getDiag(SEXP x)
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

SEXP dgeMatrix_LU(SEXP x)
{
    SEXP val = get_factors(x, "LU");
    int *dims, npiv, info;

    if (val != R_NilValue) return val;
    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    if (dims[0] < 1 || dims[1] < 1)
	error(_("Cannot factor a matrix with zero extents"));
    npiv = (dims[0] <dims[1]) ? dims[0] : dims[1];
    val = PROTECT(NEW_OBJECT(MAKE_CLASS("LU")));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(x, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    F77_CALL(dgetrf)(dims, dims + 1, REAL(GET_SLOT(val, Matrix_xSym)),
		     dims,
		     INTEGER(ALLOC_SLOT(val, Matrix_permSym, INTSXP, npiv)),
		     &info);
    if (info < 0)
	error(_("Lapack routine %s returned error code %d"), "dgetrf", info);
    else if (info > 0)
	warning(_("Exact singularity detected during LU decomposition."));
    UNPROTECT(1);
    return set_factors(x, val, "LU");
}

SEXP dgeMatrix_determinant(SEXP x, SEXP logarithm)
{
    int lg = asLogical(logarithm);
    SEXP lu = dgeMatrix_LU(x);
    int *dims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
	*jpvt = INTEGER(GET_SLOT(lu, Matrix_permSym)),
	i, n = dims[0], sign = 1;
    double *luvals = REAL(GET_SLOT(lu, Matrix_xSym)), modulus;

    if (n != dims[1])
	error(_("Determinant requires a square matrix"));
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

SEXP dgeMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	lu = dgeMatrix_LU(a);
    int *dims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
	*pivot = INTEGER(GET_SLOT(lu, Matrix_permSym));
    double *x, tmp;
    int	info, lwork = -1;


    if (dims[0] != dims[1]) error(_("Solve requires a square matrix"));
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

SEXP dgeMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed)
{
    int cl = asLogical(classed);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	lu = dgeMatrix_LU(a);
    int *adims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol));
    int info, n = bdims[0], nrhs = bdims[1];
    int sz = n * nrhs;

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1 || *adims != adims[1])
	error(_("Dimensions of system to be solved are inconsistent"));
    Memcpy(INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2)), bdims, 2);
    F77_CALL(dgetrs)("N", &n, &nrhs, REAL(GET_SLOT(lu, Matrix_xSym)), &n,
		     INTEGER(GET_SLOT(lu, Matrix_permSym)),
		     Memcpy(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, sz)),
			    REAL(cl ? GET_SLOT(b, Matrix_xSym):b), sz),
		     &n, &info);
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_matrix_mm(SEXP a, SEXP b, SEXP classed, SEXP right)
{
    int cl = asLogical(classed);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol)),
	*cdims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));
    double one = 1., zero = 0.;

    if (asLogical(right)) {
	int m = bdims[0], n = adims[1], k = bdims[1];
	if (adims[0] != k)
	    error(_("Matrices are not conformable for multiplication"));
	if (m < 1 || n < 1 || k < 1)
	    error(_("Matrices with zero extents cannot be multiplied"));
	cdims[0] = m; cdims[1] = n;
	F77_CALL(dgemm) ("N", "N", &m, &n, &k, &one,
			 REAL(cl ? GET_SLOT(b, Matrix_xSym) : b), &m,
			 REAL(GET_SLOT(a, Matrix_xSym)), &k, &zero,
			 REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n)),
			 &m);
    } else {
	int m = adims[0], n = bdims[1], k = adims[1];

	if (bdims[0] != k)
	    error(_("Matrices are not conformable for multiplication"));
	if (m < 1 || n < 1 || k < 1)
	    error(_("Matrices with zero extents cannot be multiplied"));
	cdims[0] = m; cdims[1] = n;
	F77_CALL(dgemm)
	    ("N", "N", &m, &n, &k, &one, REAL(GET_SLOT(a, Matrix_xSym)),
	     &m, REAL(cl ? GET_SLOT(b, Matrix_xSym) : b), &k, &zero,
	     REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n)), &m);
    }
    UNPROTECT(1);
    return val;
}


SEXP dgeMatrix_svd(SEXP x, SEXP nnu, SEXP nnv)
{
    int /* nu = asInteger(nnu),
	   nv = asInteger(nnv), */
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

static double padec [] =   /*  Constants for matrix exponential calculation. */
{
  5.0000000000000000e-1,
  1.1666666666666667e-1,
  1.6666666666666667e-2,
  1.6025641025641026e-3,
  1.0683760683760684e-4,
  4.8562548562548563e-6,
  1.3875013875013875e-7,
  1.9270852604185938e-9,
};

/**
 * Matrix exponential - based on the code for Octave's expm function.
 *
 * @param x real square matrix to exponentiate
 *
 * @return matrix exponential of x
 */
SEXP dgeMatrix_exp(SEXP x)
{
    SEXP val = PROTECT(duplicate(x));
    int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int i, ilo, ilos, ihi, ihis, j, nc = Dims[1], sqpow;
    int ncp1 = Dims[1] + 1, ncsqr = nc * nc;
    int *pivot = Calloc(nc, int);
    int *iperm = Calloc(nc, int);
    double *dpp = Calloc(ncsqr, double), /* denominator power Pade' */
	*npp = Calloc(ncsqr, double), /* numerator power Pade' */
	*perm = Calloc(nc, double),
	*scale = Calloc(nc, double),
	*v = REAL(GET_SLOT(val, Matrix_xSym)),
	*work = Calloc(ncsqr, double), inf_norm, m1_j, /* (-1)^j */
	one = 1., trshift, zero = 0.;

    if (nc < 1 || Dims[0] != nc)
	error(_("Matrix exponential requires square, non-null matrix"));

    /* FIXME: Add special treatment for nc == 1 */

    /* Preconditioning 1.  Shift diagonal by average diagonal if positive. */
    trshift = 0;		/* determine average diagonal element */
    for (i = 0; i < nc; i++) trshift += v[i * ncp1];
    trshift /= nc;
    if (trshift > 0.) {		/* shift diagonal by -trshift */
	for (i = 0; i < nc; i++) v[i * ncp1] -= trshift;
    }

    /* Preconditioning 2. Balancing with dgebal. */
    F77_CALL(dgebal)("P", &nc, v, &nc, &ilo, &ihi, perm, &j);
    if (j) error(_("dgeMatrix_exp: LAPACK routine dgebal returned %d"), j);
    F77_CALL(dgebal)("S", &nc, v, &nc, &ilos, &ihis, scale, &j);
    if (j) error(_("dgeMatrix_exp: LAPACK routine dgebal returned %d"), j);

    /* Preconditioning 3. Scaling according to infinity norm */
    inf_norm = F77_CALL(dlange)("I", &nc, &nc, v, &nc, work);
    sqpow = (inf_norm > 0) ? (int) (1 + log(inf_norm)/log(2.)) : 0;
    if (sqpow < 0) sqpow = 0;
    if (sqpow > 0) {
	double scale_factor = 1.0;
	for (i = 0; i < sqpow; i++) scale_factor *= 2.;
	for (i = 0; i < ncsqr; i++) v[i] /= scale_factor;
    }

    /* Pade' approximation. Powers v^8, v^7, ..., v^1 */
    AZERO(npp, ncsqr);
    AZERO(dpp, ncsqr);
    m1_j = -1;
    for (j = 7; j >=0; j--) {
	double mult = padec[j];
	/* npp = m * npp + padec[j] *m */
	F77_CALL(dgemm)("N", "N", &nc, &nc, &nc, &one, v, &nc, npp, &nc,
			&zero, work, &nc);
	for (i = 0; i < ncsqr; i++) npp[i] = work[i] + mult * v[i];
	/* dpp = m * dpp * (m1_j * padec[j]) * m */
	mult *= m1_j;
	F77_CALL(dgemm)("N", "N", &nc, &nc, &nc, &one, v, &nc, dpp, &nc,
			&zero, work, &nc);
	for (i = 0; i < ncsqr; i++) dpp[i] = work[i] + mult * v[i];
	m1_j *= -1;
    }
    /* Zero power */
    for (i = 0; i < ncsqr; i++) dpp[i] *= -1.;
    for (j = 0; j < nc; j++) {
	npp[j * ncp1] += 1.;
	dpp[j * ncp1] += 1.;
    }

    /* Pade' approximation is solve(dpp, npp) */
    F77_CALL(dgetrf)(&nc, &nc, dpp, &nc, pivot, &j);
    if (j) error(_("dgeMatrix_exp: dgetrf returned error code %d"), j);
    F77_CALL(dgetrs)("N", &nc, &nc, dpp, &nc, pivot, npp, &nc, &j);
    if (j) error(_("dgeMatrix_exp: dgetrs returned error code %d"), j);
    Memcpy(v, npp, ncsqr);

    /* Now undo all of the preconditioning */
    /* Preconditioning 3: square the result for every power of 2 */
    while (sqpow--) {
	F77_CALL(dgemm)("N", "N", &nc, &nc, &nc, &one, v, &nc, v, &nc,
			&zero, work, &nc);
	Memcpy(v, work, ncsqr);
    }
    /* Preconditioning 2: apply inverse scaling */
    for (j = 0; j < nc; j++)
	for (i = 0; i < nc; i++)
	    v[i + j * nc] *= scale[i]/scale[j];
    /* Construct balancing permutation vector */
    for (i = 0; i < nc; i++) iperm[i] = i; /* identity permutation */
    /* Leading permutations applied in forward order */
    for (i = 0; i < (ilo - 1); i++) {
	int swapidx = (int) (perm[i]) - 1;
	int tmp = iperm[i];
	iperm[i] = iperm[swapidx];
	iperm[swapidx] = tmp;
    }
    /* Trailing permutations applied in reverse order */
    for (i = nc - 1; i >= ihi; i--) {
	int swapidx = (int) (perm[i]) - 1;
	int tmp = iperm[i];
	iperm[i] = iperm[swapidx];
	iperm[swapidx] = tmp;
    }
    /* Construct inverse balancing permutation vector */
    Memcpy(pivot, iperm, nc);
    for (i = 0; i < nc; i++) iperm[pivot[i]] = i;
    /* Apply inverse permutation */
    Memcpy(work, v, ncsqr);
    for (j = 0; j < nc; j++)
	for (i = 0; i < nc; i++)
	    v[i + j * nc] = work[iperm[i] + iperm[j] * nc];

    /* Preconditioning 1: Trace normalization */
    if (trshift > 0.) {
	double mult = exp(trshift);
	for (i = 0; i < ncsqr; i++) v[i] *= mult;
    }

    /* Clean up */
    Free(dpp); Free(npp); Free(perm); Free(iperm); Free(pivot); Free(scale); Free(work);
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_Schur(SEXP x, SEXP vectors)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int vecs = asLogical(vectors), info, izero = 0, lwork = -1, n = dims[0];
    double *work, tmp;
    char *nms[] = {"WR", "WI", "T", "Z", ""};
    SEXP val = PROTECT(Matrix_make_named(VECSXP, nms));

    if (n != dims[1] || n < 1)
	error(_("dgeMatrix_Schur: argument x must be a non-null square matrix"));
    SET_VECTOR_ELT(val, 0, allocVector(REALSXP, n));
    SET_VECTOR_ELT(val, 1, allocVector(REALSXP, n));
    SET_VECTOR_ELT(val, 2, allocMatrix(REALSXP, n, n));
    Memcpy(REAL(VECTOR_ELT(val, 2)), REAL(GET_SLOT(x, Matrix_xSym)), n * n);
    SET_VECTOR_ELT(val, 3, allocMatrix(REALSXP, vecs ? n : 0, vecs ? n : 0));
    F77_CALL(dgees)(vecs ? "V" : "N", "N", NULL, dims, (double *) NULL, dims, &izero,
		    (double *) NULL, (double *) NULL, (double *) NULL, dims,
		    &tmp, &lwork, (int *) NULL, &info);
    if (info) error(_("dgeMatrix_Schur: first call to dgees failed"));
    lwork = (int) tmp;
    work = Calloc(lwork, double);
    F77_CALL(dgees)(vecs ? "V" : "N", "N", NULL, dims, REAL(VECTOR_ELT(val, 2)), dims,
		    &izero, REAL(VECTOR_ELT(val, 0)), REAL(VECTOR_ELT(val, 1)),
		    REAL(VECTOR_ELT(val, 3)), dims, work, &lwork,
		    (int *) NULL, &info);
    if (info) error(_("dgeMatrix_Schur: dgees returned code %d"), info);
    Free(work);
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_colsums(SEXP x, SEXP naRmP, SEXP cols, SEXP mean)
{
    int keepNA = !asLogical(naRmP);
    int doMean = asLogical(mean);
    int useCols = asLogical(cols);
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int cnt = 0, i, j, n = dims[0], p = dims[1];
    SEXP ans = PROTECT(allocVector(REALSXP, (useCols) ? p : n));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *rx, sum;

    if (useCols) {
	cnt = n;
	for (j = 0; j < p; j++) {
	    rx = xx + n*j;
	    if (keepNA)
		for (sum = 0., i = 0; i < n; i++) sum += *rx++;
	    else {
		for (cnt = 0, sum = 0., i = 0; i < n; i++, rx++)
		    if (!ISNAN(*rx)) {cnt++; sum += *rx;}
	    }
	    if (doMean) {
		if (cnt > 0) sum /= cnt; else sum = NA_REAL;
	    }
	    REAL(ans)[j] = sum;
	}
    } else {
	double *rans = REAL(ans), *ra = rans, *rx = xx, *Cnt = NULL, *c;
	cnt = p;
	if (!keepNA && doMean) Cnt = Calloc(n, double);
	for (ra = rans, i = 0; i < n; i++) *ra++ = 0.0;
	for (j = 0; j < p; j++) {
	    ra = rans;
	    if (keepNA)
		for (i = 0; i < n; i++) *ra++ += *rx++;
	    else
		for (c = Cnt, i = 0; i < n; i++, ra++, rx++, c++)
		    if (!ISNAN(*rx)) {
			*ra += *rx;
			if (doMean) (*c)++;
		    }
	}
	if (doMean) {
	    if (keepNA)
		for (ra = rans, i = 0; i < n; i++)
		    *ra++ /= p;
	    else {
		for (ra = rans, c = Cnt, i = 0; i < n; i++, c++)
		    if (*c > 0) *ra++ /= *c; else *ra++ = NA_REAL;
		Free(Cnt);
	    }
	}
    }

    UNPROTECT(1);
    return ans;
}
