			/* Sparse triangular matrices in triplet format */
#include "dtTMatrix.h"

SEXP dtTMatrix_validate(SEXP x)
{
    return triangularMatrix_validate(x);
    /* see ./dtpMatrix.c as example to do more testing here */
}

SEXP dtTMatrix_as_dtrMatrix(SEXP x)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))),
	dimP = GET_SLOT(x, Matrix_DimSym),
	xiP = GET_SLOT(x, Matrix_iSym);
    int k, m = INTEGER(dimP)[0], n = INTEGER(dimP)[1], nnz = length(xiP);
    int *xi = INTEGER(xiP), *xj = INTEGER(GET_SLOT(x, Matrix_jSym)),
	sz = m * n;
    double *tx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, sz)),
	*xx = REAL(GET_SLOT(x, Matrix_xSym));

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(x, Matrix_uploSym)));
    SET_SLOT(val, Matrix_diagSym, duplicate(GET_SLOT(x, Matrix_diagSym)));
    AZERO(tx, sz);
    for (k = 0; k < nnz; k++) tx[xi[k] + xj[k] * m] = xx[k];
    UNPROTECT(1);
    return val;
}

SEXP dtTMatrix_as_dgCMatrix(SEXP x)
{
    cholmod_triplet *tx = as_cholmod_triplet(x);
    cholmod_sparse *cx = cholmod_triplet_to_sparse(tx, tx->nzmax, &c);
				/* chm_sparse_to_SEXP cholmod_frees cx */
    SEXP val = PROTECT(chm_sparse_to_SEXP(cx, 1, 0, "", R_NilValue));


    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(GET_SLOT(x, Matrix_DimNamesSym)));
    Free(tx);
    UNPROTECT(1);
    return val;

#if 0
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix"))),
	dimP = GET_SLOT(x, Matrix_DimSym),
	xiP = GET_SLOT(x, Matrix_iSym);
    int n = INTEGER(dimP)[0],
	nnz = length(xiP);
    int *ti = Calloc(nnz, int),
        *vp = INTEGER(ALLOC_SLOT(val, Matrix_pSym, INTSXP, n + 1));
    double *tx = Calloc(nnz, double);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(x, Matrix_uploSym)));
    SET_SLOT(val, Matrix_diagSym, duplicate(GET_SLOT(x, Matrix_diagSym)));

    triplet_to_col(n, n, nnz, INTEGER(xiP),
		   INTEGER(GET_SLOT(x, Matrix_jSym)),
		   REAL(GET_SLOT(x, Matrix_xSym)),
		   vp, ti, tx);
    nnz = vp[n];
    Memcpy(INTEGER(ALLOC_SLOT(val, Matrix_iSym,  INTSXP, nnz)), ti, nnz);
    Memcpy(   REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, nnz)), tx, nnz);
    Free(ti); Free(tx);
    UNPROTECT(1);
    return val;
#endif
}

