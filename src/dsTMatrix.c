				/* Sparse symmetric matrices in triplet format */
#include "dsTMatrix.h"

SEXP dsTMatrix_validate(SEXP x) /* == lsTMatrix_validate */
{
    SEXP xxP = symmetricMatrix_validate(x);
    if(isString(xxP))
	return(xxP);
    else {
	SEXP xiP = GET_SLOT(x, Matrix_iSym),
	    xjP = GET_SLOT(x, Matrix_jSym);
	xxP = GET_SLOT(x, Matrix_xSym);
	if (length(xiP) != length(xjP) || length(xjP) != length(xxP))
	    return mkString(_("slots i, j and x must have the same length"));
	return ScalarLogical(1);
    }
}

SEXP dsTMatrix_as_dsyMatrix(SEXP x)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dsyMatrix"))),
	DimP = GET_SLOT(x, Matrix_DimSym),
	xiP = GET_SLOT(x, Matrix_iSym);
    int k, n = INTEGER(DimP)[1], nnz = length(xiP);
    int *xi = INTEGER(xiP), *xj = INTEGER(GET_SLOT(x, Matrix_jSym)),
	sz = n * n;
    double *tx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, sz)),
	*xx = REAL(GET_SLOT(x, Matrix_xSym));

    SET_SLOT(val, Matrix_DimSym, duplicate(DimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(x, Matrix_uploSym)));
    AZERO(tx, sz);
    for (k = 0; k < nnz; k++) tx[xi[k] + xj[k] * n] = xx[k];
    UNPROTECT(1);
    return val;
}

/* this corresponds to changing 'stype' of a cholmod_triplet; seems not available there */
SEXP dsTMatrix_as_dgTMatrix(SEXP x)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgTMatrix"))),
	dimP = GET_SLOT(x, Matrix_DimSym),
	xiP = GET_SLOT(x, Matrix_iSym);
    /* , uplo = GET_SLOT(x, Matrix_uploSym); */
    int i, nnz = length(xiP), n0d, nv,
	*xi = INTEGER(xiP),
	*xj = INTEGER(GET_SLOT(x, Matrix_jSym)),
	*vi, *vj;
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *vx;

    /* Find *length* of result slots:  = 2 * nnz - n0d; n0d := #{non-0 diagonals} :*/
    for(i = 0, n0d = 0; i < nnz; i++)
	if(xi[i] == xj[i]) n0d++ ;
    nv = 2 * nnz - n0d;

    vi = INTEGER(ALLOC_SLOT(val, Matrix_iSym, INTSXP, nv));
    vj = INTEGER(ALLOC_SLOT(val, Matrix_jSym, INTSXP, nv));
    vx =    REAL(ALLOC_SLOT(val, Matrix_xSym,REALSXP, nv));

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));

    /* copy the upper/lower triangle (including the diagonal) "at end" ([nv]): */
    nv = nnz - n0d;
    Memcpy(&vi[nv], xi, nnz);
    Memcpy(&vj[nv], xj, nnz);
    Memcpy(&vx[nv], xx, nnz);

    for(i = 0, nv = 0; i < nnz; i++) { /* copy the other triangle */
	if(xi[i] != xj[i]) { /* but not the diagonal */
	    vi[nv] = xj[i];
	    vj[nv] = xi[i];
	    vx[nv] = xx[i];
	    nv++;
	}
    }

    UNPROTECT(1);
    return val;
}
