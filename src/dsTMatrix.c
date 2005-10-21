				/* Sparse symmetric matrices in triplet format */
#include "dsTMatrix.h"

SEXP dsTMatrix_validate(SEXP x)
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

/* FIXME: no longer needed, with Tsparse_to_Csparse() ! */
SEXP dsTMatrix_as_dsCMatrix(SEXP x)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix"))),
	dimP = GET_SLOT(x, Matrix_DimSym),
	xiP = GET_SLOT(x, Matrix_iSym);
    int n = INTEGER(dimP)[0],
	nnz = length(xiP);
    int *ti = Calloc(nnz, int),
        *vp = INTEGER(ALLOC_SLOT(val, Matrix_pSym, INTSXP, n + 1));
    double *tx = Calloc(nnz, double);

    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(x, Matrix_uploSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    triplet_to_col(n, n, nnz, INTEGER(xiP),
		   INTEGER(GET_SLOT(x, Matrix_jSym)),
		   REAL(GET_SLOT(x, Matrix_xSym)),
		   vp, ti, tx);
    nnz = vp[n];
    Memcpy(INTEGER(ALLOC_SLOT(val, Matrix_iSym, INTSXP, nnz)), ti, nnz);
    Memcpy(   REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, nnz)), tx, nnz);
    Free(ti); Free(tx);
    UNPROTECT(1);
    return val;
}

/* this corresponds to changing 'stype' of a cholmod_triplet; seems not available there */
SEXP dsTMatrix_as_dgTMatrix(SEXP x)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgTMatrix"))),
	dimP = GET_SLOT(x, Matrix_DimSym),
	xiP = GET_SLOT(x, Matrix_iSym),
	uplo = GET_SLOT(x, Matrix_uploSym);
    int i, nnz = length(xiP);
    int *vi = INTEGER(ALLOC_SLOT(val, Matrix_iSym, INTSXP, 2 * nnz)),
	*vj = INTEGER(ALLOC_SLOT(val, Matrix_jSym, INTSXP, 2 * nnz)),
	*vx = INTEGER(ALLOC_SLOT(val, Matrix_xSym,REALSXP, 2 * nnz));

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));

    if(*CHAR(STRING_ELT(uplo, 0)) == 'U') { /* x stored in upper triangle */

	Memcpy(&vi[nnz], INTEGER(GET_SLOT(x, Matrix_iSym)), nnz);
	Memcpy(&vj[nnz], INTEGER(GET_SLOT(x, Matrix_jSym)), nnz);
	Memcpy(&vx[nnz],    REAL(GET_SLOT(x, Matrix_xSym)), nnz);

	for(i = 0; i < nnz; i++) { /* copy the other triangle */
	    vi[i] = INTEGER(GET_SLOT(x, Matrix_jSym))[i];
	    vj[i] = INTEGER(GET_SLOT(x, Matrix_iSym))[i];
	    vx[i] = INTEGER(GET_SLOT(x, Matrix_xSym))[i];
	}

    }
    else { /* x is stored in lower triangle */

	error("dsTMatrix_as_dgTMatrix(.) not yet for  uplo != 'U'");

    }

    UNPROTECT(1);
    return val;
}
