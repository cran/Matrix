			/* Sparse triangular matrices in triplet format */
#include "dtTMatrix.h"

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
/* Should generalize this, also for ltT -> lgC ... .. or drop completely
 * FIXME
 * int Rkind = (tx->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;
 */
    Free(tx);
				/* chm_sparse_to_SEXP cholmod_frees cx */
    return chm_sparse_to_SEXP(cx, 1, 0, /*Rkind*/ 0, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

