#include "dgCMatrix.h"

SEXP dgCMatrix_validate(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym),
	xslot = GET_SLOT(x, Matrix_xSym);
    int j,
	ncol = length(pslot) - 1,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow,
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    nrow = dims[0];
    if (length(islot) != length(xslot))
	return mkString(_("lengths of slots i and x must match"));
    if (length(pslot) <= 0)
	return mkString(_("slot p must have length > 0"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return mkString(_("last element of slot p must match length of slots i and x"));
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
    }
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    if (csc_unsorted_columns(ncol, xp, xi)) {
	csc_sort_columns(ncol, xp, xi, REAL(xslot));
    }
    return ScalarLogical(1);
}

SEXP csc_crossprod(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix"))), tmp;
    int *xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

    int j, *iVal, ncol = length(pslot) - 1, maxnz, nnz = 0, *pVal;
    double *xVal;

    SET_SLOT(ans, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    maxnz = (ncol * (ncol + 1))/2;
    iVal = Calloc(maxnz, int); xVal = Calloc(maxnz, double);
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, ncol + 1));
    tmp = GET_SLOT(ans, Matrix_pSym);
    pVal = INTEGER(tmp);
    for (j = 0; j < ncol; j++) {
	pVal[j] = nnz;
	if (xp[j] < xp[j+1]) {	/* column j contains some non-zeros */
	    int ind, jj;
	    double accum = 0.;
				/* diagonal elements */
	    for (ind = xp[j]; ind < xp[j+1]; ind++)
		accum += xx[ind] * xx[ind];
	    iVal[nnz] = j;
	    xVal[nnz] = accum;
	    nnz++;
				/* off-diagonals (lower triangle only) */
	    for (jj = j+1; jj < ncol; jj++) {
		int ind2;

		ind = xp[j];
		ind2 = xp[jj];
		accum = 0.;
		while (ind < xp[j+1] && ind2 < xp[jj+1]) {
		    if (xi[ind] < xi[ind2]) ind++;
		    else {
			if (xi[ind] > xi[ind2]) ind2++;
			else {
			    accum += xx[ind] * xx[ind2];
			    ind++; ind2++;
			}
		    }
		}
		if (accum != 0.) {
		    iVal[nnz] = jj;
		    xVal[nnz] = accum;
		    nnz++;
		}
	    }
	}
    }
    pVal[ncol] = nnz;

    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    Memcpy(INTEGER(GET_SLOT(ans, Matrix_iSym)), iVal, nnz);
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    Memcpy(REAL(GET_SLOT(ans, Matrix_xSym)), xVal, nnz);
    Free(iVal); Free(xVal); UNPROTECT(1);
    return dgCMatrix_set_Dim(ans, ncol);
}

SEXP csc_tcrossprod(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix")));
    int *xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

    int j, ntrip, *iVal, nrow = dims[0], ncol = dims[1], *jVal, nnz, pos;
    int *itmp, *ansp;
    double *xVal, *xtmp;

    SET_SLOT(ans, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
    ntrip = nrow;  		/* number of triplets */
    for (j = 0; j < ncol; j++) {
	int nzj = xp[j+1] - xp[j];
	ntrip += (nzj * (nzj - 1))/2;
    }
    iVal = Calloc(ntrip, int); jVal = Calloc(ntrip, int);
    xVal = Calloc(ntrip, double);
    for (j = 0; j < nrow; j++) {
	iVal[j] = jVal[j] = j;
	xVal[j] = 0.;
    }
    pos = nrow;
    for (j = 0; j < ncol; j++) {
	int k, kk, k2 = xp[j+1];
	for (k = xp[j]; k < k2; k++) {
	    int r1 = xi[k];
	    double x1 = xx[k];
	    xVal[r1] += x1 * x1;
	    for (kk = k + 1; kk < k2; kk++) {
		int r2 = xi[kk];
		double x2 = xx[kk];
		jVal[pos] = r1;
		iVal[pos] = r2;
		xVal[pos] = x1 * x2;
		pos++;
	    }
	}
    }
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, nrow + 1));
    ansp = INTEGER(GET_SLOT(ans, Matrix_pSym));
    itmp = Calloc(ntrip, int); xtmp = Calloc(ntrip, double);
    triplet_to_col(nrow, nrow, ntrip, iVal, jVal, xVal,
		   ansp, itmp, xtmp);
    nnz = ansp[nrow];
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    Memcpy(INTEGER(GET_SLOT(ans, Matrix_iSym)), itmp, nnz);
    Memcpy(REAL(GET_SLOT(ans, Matrix_xSym)), xtmp, nnz);
    dims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    dims[0] = dims[1] = nrow;
    Free(itmp); Free(xtmp); Free(iVal); Free(jVal); Free(xVal);
    UNPROTECT(1);
    return ans;
}

SEXP csc_matrix_crossprod(SEXP x, SEXP y)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym), ans;
    int j,
	*xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym)),
	xncol = length(pslot) - 1,
	xnrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	*ydims;
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

    if (!(isMatrix(y) && isReal(y))) error(_("y must be a numeric matrix"));
    ydims = INTEGER(getAttrib(y, R_DimSymbol));
    if (xnrow != ydims[0]) error(_("x and y must have the same number of rows"));
    ans = PROTECT(allocMatrix(REALSXP, xncol, ydims[1]));
    for (j = 0; j < ydims[1]; j++) {
	int i; double *ypt = REAL(y) + j * ydims[0];
	for(i = 0; i < xncol; i++) {
	    int ii; double accum = 0.;
	    for (ii = xp[i]; ii < xp[i+1]; ii++) {
		accum += xx[ii] * ypt[xi[ii]];
	    }
	    REAL(ans)[i + j * xncol] = accum;
	}
    }
    UNPROTECT(1);
    return ans;
}

SEXP compressed_to_dgTMatrix(SEXP x, SEXP colP)
{
    int col = asLogical(colP);
    SEXP indSym = col ? Matrix_iSym : Matrix_jSym;
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgTMatrix"))),
	indP = GET_SLOT(x, indSym),
	pP = GET_SLOT(x, Matrix_pSym);
    int npt = length(pP) - 1;

    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    SET_SLOT(ans, Matrix_xSym, duplicate(GET_SLOT(x, Matrix_xSym)));
    SET_SLOT(ans, indSym, duplicate(indP));
    expand_cmprPt(npt, INTEGER(pP),
		  INTEGER(ALLOC_SLOT(ans, col ? Matrix_jSym : Matrix_iSym,
				     INTSXP, length(indP))));
    UNPROTECT(1);
    return ans;
}

SEXP csc_to_matrix(SEXP x)
{
    SEXP ans, pslot = GET_SLOT(x, Matrix_pSym);
    int j, ncol = length(pslot) - 1,
	nrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	*xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *ax;

    ax = REAL(ans = PROTECT(allocMatrix(REALSXP, nrow, ncol)));
    for (j = 0; j < (nrow * ncol); j++) ax[j] = 0.;
    for (j = 0; j < ncol; j++) {
	int ind;
	for (ind = xp[j]; ind < xp[j+1]; ind++) {
	    ax[j * nrow + xi[ind]] = xx[ind];
	}
    }
    UNPROTECT(1);
    return ans;
}

SEXP csc_to_dgeMatrix(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	Dimslot = GET_SLOT(x, Matrix_DimSym);
    int *dims = INTEGER(Dimslot),
	*xp = INTEGER(GET_SLOT(x, Matrix_pSym)),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *ax;
    int j, nrow = dims[0], ncol = dims[1];

    SET_SLOT(ans, Matrix_DimSym, duplicate(Dimslot));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nrow*ncol));
    SET_SLOT(ans, Matrix_rcondSym, allocVector(REALSXP, 0));
    SET_SLOT(ans, Matrix_factorSym, allocVector(VECSXP, 0));
    ax = REAL(GET_SLOT(ans, Matrix_xSym));
    for (j = 0; j < (nrow * ncol); j++) ax[j] = 0.;
    for (j = 0; j < ncol; j++) {
	int ind;
	for (ind = xp[j]; ind < xp[j+1]; ind++) {
	    ax[j * nrow + xi[ind]] = xx[ind];
	}
    }
    UNPROTECT(1);
    return ans;
}

SEXP matrix_to_csc(SEXP A)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    int *adims = INTEGER(getAttrib(A, R_DimSymbol)), j,
	maxnz, nrow, ncol, nnz, *vp, *vi;

    double *vx;

    if (!(isMatrix(A) && isReal(A)))
	error(_("A must be a numeric matrix"));
    nrow = adims[0]; ncol = adims[1];
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    SET_SLOT(val, Matrix_pSym, allocVector(INTSXP, ncol + 1));
    vp = INTEGER(GET_SLOT(val, Matrix_pSym));
    maxnz = nrow * ncol;
    vi = Calloc(maxnz, int); vx = Calloc(maxnz, double);
    nnz = 0;
    for (j = 0; j < ncol; j++) {
	int i;
	vp[j] = nnz;
	for (i = 0; i < nrow; i++) {
	    double val = REAL(A)[i + j * nrow];
	    if (val != 0.) {
		vi[nnz] = i;
		vx[nnz] = val;
		nnz++;
	    }
	}
    }
    vp[ncol] = nnz;
    SET_SLOT(val, Matrix_iSym, allocVector(INTSXP, nnz));
    Memcpy(INTEGER(GET_SLOT(val, Matrix_iSym)), vi, nnz);
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, nnz));
    Memcpy(REAL(GET_SLOT(val, Matrix_xSym)), vx, nnz);
    Free(vi); Free(vx);
    UNPROTECT(1);
    return dgCMatrix_set_Dim(val, nrow);
}


SEXP dgTMatrix_to_csc(SEXP dgTMatrix)
{
    SEXP Tisl = GET_SLOT(dgTMatrix, Matrix_iSym);
    int *Ti = INTEGER(Tisl),
	*Tj = INTEGER(GET_SLOT(dgTMatrix, Matrix_jSym)),
	i, nrow, ncol,
	nz = length(Tisl);

    nrow = ncol = -1;
    for(i = 0; i < nz; i++) {
	if (Ti[i] > nrow) nrow = Ti[i];
	if (Tj[i] > ncol) ncol = Tj[i];
    }
    return triple_as_SEXP(nrow + 1, ncol + 1, nz, Ti, Tj,
			  REAL(GET_SLOT(dgTMatrix, Matrix_xSym)),
			  "dgCMatrix");
}

SEXP csc_getDiag(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym), ans;
    int *xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym)),
	j,
	ncol = length(pslot) - 1,
	nrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	ndiag;
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *diag;

    ndiag = (nrow < ncol) ? nrow : ncol;
    ans = PROTECT(allocVector(REALSXP, ndiag));
    diag = REAL(ans);
    for (j = 0; j < ndiag; j++) {
	int ind;
	diag[j] = 0.;
	for (ind = xp[j]; ind < xp[j+1]; ind++) {
	    if (xi[ind] == j) diag[j] = xx[ind];
	}
    }
    UNPROTECT(1);
    return ans;
}

SEXP csc_transpose(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))),
	islot = GET_SLOT(x, Matrix_iSym);
    int *adims,	*xdims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nnz = length(islot);

    adims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    adims[0] = xdims[1]; adims[1] = xdims[0];
    csc_compTr(xdims[0], xdims[1], nnz,
	       INTEGER(GET_SLOT(x, Matrix_pSym)), INTEGER(islot),
	       REAL(GET_SLOT(x, Matrix_xSym)),
	       INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, xdims[0] + 1)),
	       INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nnz)),
	       REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nnz)));
    UNPROTECT(1);
    return ans;
}

SEXP csc_matrix_mm(SEXP a, SEXP b)
{
    int *adim = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*ai = INTEGER(GET_SLOT(a, Matrix_iSym)),
	*ap = INTEGER(GET_SLOT(a, Matrix_pSym)),
	*bdim = INTEGER(getAttrib(b, R_DimSymbol));
    int j, k, m = adim[0], n = bdim[1], r = adim[1];
    double *ax = REAL(GET_SLOT(a, Matrix_xSym));
    SEXP val;

    if (bdim[0] != r)
	error(_("Matrices of sizes (%d,%d) and (%d,%d) cannot be multiplied"),
	      m, r, bdim[0], n);
    val = PROTECT(allocMatrix(REALSXP, m, n));
    for (j = 0; j < n; j++) {	/* across columns of b */
	double *ccol = REAL(val) + j * m,
	    *bcol = REAL(b) + j * r;

	for (k = 0; k < m; k++) ccol[k] = 0.; /* zero the accumulators */
	for (k = 0; k < r; k++) { /* across columns of a */
	    int kk, k2 = ap[k + 1];
	    for (kk = ap[k]; kk < k2; kk++) {
		ccol[ai[kk]] += ax[kk] * bcol[k];
	    }
	}
    }
    UNPROTECT(1);
    return val;
}

SEXP csc_col_permute(SEXP x, SEXP perm)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))), tmp;
    int *iperm, *prm, *vi, *vp, *xi, *xp, j, k, ncol, pos;
    double *vx, *xx;

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    tmp = GET_SLOT(x, Matrix_DimSym);
    SET_SLOT(val, Matrix_DimSym, duplicate(tmp));
    ncol = INTEGER(tmp)[1];
    if (!(isInteger(perm) && length(perm) == ncol))
	error(_("perm must be an integer vector of length %d"),
	      ncol);
    prm = INTEGER(perm);
    if (!R_ldl_valid_perm(ncol, prm))
	error(_("perm is not a valid 0-based permutation"));
    iperm = Calloc(ncol, int);
    for (j = 0; j < ncol; j++) iperm[prm[j]] = j;
    tmp = GET_SLOT(x, Matrix_pSym);
    xp = INTEGER(tmp);
    SET_SLOT(val, Matrix_pSym, duplicate(tmp));
    vp = INTEGER(GET_SLOT(val, Matrix_pSym));
    tmp = GET_SLOT(x, Matrix_iSym);
    xi = INTEGER(tmp);
    SET_SLOT(val, Matrix_iSym, duplicate(tmp));
    vi = INTEGER(GET_SLOT(val, Matrix_iSym));
    tmp = GET_SLOT(x, Matrix_xSym);
    xx = REAL(tmp);
    SET_SLOT(val, Matrix_xSym, duplicate(tmp));
    vx = REAL(GET_SLOT(val, Matrix_xSym));

    pos = vp[0] = 0;
    for (j = 0; j < ncol; j++) {
	int jj = iperm[j];
	int j1 = xp[jj], j2 = xp[jj+1];
	vp[j + 1] = vp[j] + (j2 - j1);
	for (k = j1; k < j2; k++) {
	    vi[pos] = xi[k];
	    vx[pos] = xx[k];
	    pos++;
	}
    }
    Free(iperm);
    UNPROTECT(1);
    return val;
}



