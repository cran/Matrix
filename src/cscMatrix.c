#include "cscMatrix.h"

SEXP csc_validate(SEXP x)
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
	return ScalarString(mkChar("lengths of slots i and x must match"));
    if (length(pslot) <= 0)
	return ScalarString(mkChar("slot p must have length > 0"));
    if (xp[0] != 0)
	return ScalarString(mkChar("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return ScalarString(
	    mkChar(
		"last element of slot p must match length of slots i and x"));
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return ScalarString(mkChar("slot p must be non-decreasing"));
    }
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return ScalarString(
		mkChar("all row indices must be between 0 and nrow-1"));
    }
    if (csc_unsorted_columns(ncol, xp, xi)) {
	csc_sort_columns(ncol, xp, xi, REAL(xslot));
    }
    return ScalarLogical(1);
}
    
SEXP csc_crossprod(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("sscMatrix"))), tmp;
    int *xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

    int j, *iVal, ncol = length(pslot) - 1, maxnz, nnz = 0, *pVal;
    double *xVal;
    
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
    return cscMatrix_set_Dim(ans, ncol);
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

    if (!(isMatrix(y) && isReal(y))) error("y must be a numeric matrix");
    ydims = INTEGER(getAttrib(y, R_DimSymbol));
    if (xnrow != ydims[0]) error("x and y must have the same number of rows");
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

SEXP csc_to_triplet(SEXP x)
{
    SEXP
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("tripletMatrix"))),
	dimslot = GET_SLOT(x, Matrix_DimSym),
	islot = GET_SLOT(x, Matrix_iSym), 
	pslot = GET_SLOT(x, Matrix_pSym);
    int *dims = INTEGER(dimslot), j, jj,
	*xp = INTEGER(pslot), *yj;
    
    SET_SLOT(ans, Matrix_iSym, duplicate(islot));
    SET_SLOT(ans, Matrix_DimSym, duplicate(dimslot));
    SET_SLOT(ans, Matrix_xSym, duplicate(GET_SLOT(x, Matrix_xSym)));
    SET_SLOT(ans, Matrix_jSym, allocVector(INTSXP, length(islot)));
    yj = INTEGER(GET_SLOT(ans, Matrix_jSym));
    jj = 0;
    for (j = 0; j < dims[1]; j++) {
	while (jj < xp[j + 1]) {
	    yj[jj] = j;
	    jj++;
	}
    }
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

SEXP csc_to_geMatrix(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix"))),
	Dimslot = GET_SLOT(x, Matrix_DimSym);
    int *dims = INTEGER(Dimslot),
	*xp = INTEGER(GET_SLOT(x, Matrix_pSym)),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *ax;
    int j, nrow = dims[0], ncol = dims[1];
		      
    SET_SLOT(ans, Matrix_DimSym, duplicate(Dimslot));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nrow*ncol));
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

SEXP csc_to_imagemat(SEXP x)
{
    SEXP ans, pslot = GET_SLOT(x, Matrix_pSym);
    int j, ncol = length(pslot) - 1,
	nrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	*xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*ax;
				/* ans is the transpose of the indicator
				   of non-zero */
    ax = INTEGER(ans = PROTECT(allocMatrix(INTSXP, ncol, nrow)));
    for (j = 0; j < (ncol * nrow); j++) ax[j] = 0;
    for (j = 0; j < ncol; j++) {
	int ind;

	for (ind = xp[j]; ind < xp[j+1]; ind++) {
				/* reverse rows of transpose */
	    ax[j + ncol*(nrow - 1 - xi[ind])] = 1;
	}
    }
    UNPROTECT(1);
    return ans;
}

SEXP matrix_to_csc(SEXP A)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("cscMatrix")));
    int *adims = INTEGER(getAttrib(A, R_DimSymbol)), j,
	maxnz, nrow, ncol, nnz, *vp, *vi;
    
    double *vx;

    if (!(isMatrix(A) && isReal(A)))
	error("A must be a numeric matrix");
    nrow = adims[0]; ncol = adims[1];
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
    return cscMatrix_set_Dim(val, nrow);
}

    
SEXP triplet_to_csc(SEXP triplet)
{
    SEXP Tisl = GET_SLOT(triplet, Matrix_iSym);
    int *Ti = INTEGER(Tisl),
	*Tj = INTEGER(GET_SLOT(triplet, Matrix_jSym)),
	i, nrow, ncol,
	nz = length(Tisl);

    nrow = ncol = -1;
    for(i = 0; i < nz; i++) {
	if (Ti[i] > nrow) nrow = Ti[i];
	if (Tj[i] > ncol) ncol = Tj[i];
    }
    return triple_as_SEXP(nrow + 1, ncol + 1, nz, Ti, Tj,
			  REAL(GET_SLOT(triplet, Matrix_xSym)),
			  "cscMatrix");
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
    SEXP
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("cscMatrix"))),
	islot = GET_SLOT(x, Matrix_iSym);
    int nnz = length(islot),
	*adims = INTEGER(GET_SLOT(ans, Matrix_DimSym)),
	*xdims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    adims[0] = xdims[1]; adims[1] = xdims[0];
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, xdims[0] + 1));
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    csc_components_transpose(xdims[0], xdims[1], nnz,
			     INTEGER(GET_SLOT(x, Matrix_pSym)),
			     INTEGER(islot),
			     REAL(GET_SLOT(x, Matrix_xSym)),
			     INTEGER(GET_SLOT(ans, Matrix_pSym)),
			     INTEGER(GET_SLOT(ans, Matrix_iSym)),
			     REAL(GET_SLOT(ans, Matrix_xSym)));
    UNPROTECT(1);
    return ans;
}
