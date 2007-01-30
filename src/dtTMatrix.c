			/* Sparse triangular matrices in triplet format */
#include "dtTMatrix.h"
#include "dgTMatrix.h"

/* This should be use for *BOTH* triangular and symmetric Tsparse: */
SEXP tTMatrix_validate(SEXP x)
{
    SEXP val = xTMatrix_validate(x);/* checks x slot */
    if(isString(val))
	return(val);
    else {
	SEXP
	    islot = GET_SLOT(x, Matrix_iSym),
	    jslot = GET_SLOT(x, Matrix_jSym);
	int uploT = (*uplo_P(x) == 'U'),
	    k, nnz = length(islot),
	    *xj = INTEGER(jslot),
	    *xi = INTEGER(islot);

	/* Maybe FIXME: ">" should be ">="	for diag = 'U' (uplo = 'U') */
	if(uploT) {
	    for (k = 0; k < nnz; k++)
		if(xi[k] > xj[k])
		    return mkString(_("uplo='U' must not have sparse entries in lower diagonal"));
	}
	else {
	    for (k = 0; k < nnz; k++)
		if(xi[k] < xj[k])
		    return mkString(_("uplo='L' must not have sparse entries in upper diagonal"));
	}

	return ScalarLogical(1);
    }
}

/* SEXP dtTMatrix_as_dtrMatrix(SEXP x) ---> now in ./TMatrix_as.c */

/* Should generalize this, also for ltT -> lgC --
 * along the lines in ./TMatrix_as.c  ..... or drop completely : */
SEXP dtTMatrix_as_dgCMatrix(SEXP x)
{
    cholmod_triplet *tx = as_cholmod_triplet(x);
    cholmod_sparse *cx = cholmod_triplet_to_sparse(tx, tx->nzmax, &c);

 /* FIXME
 * int Rkind = (tx->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;
 */
    Free(tx);
				/* chm_sparse_to_SEXP cholmod_frees cx */
    return chm_sparse_to_SEXP(cx, 1, 0, /*Rkind*/ 0, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

