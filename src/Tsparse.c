				/* Sparse matrices in triplet form */
#include "Tsparse.h"
#ifdef USE_CHOLMOD
#include "chm_common.h"
#endif	/* USE_CHOLMOD */

SEXP Tsparse_validate(SEXP x)
{
    SEXP
	islot = GET_SLOT(x, Matrix_iSym),
	jslot = GET_SLOT(x, Matrix_jSym),
	dimslot = GET_SLOT(x, Matrix_DimSym);
    int j,
	*dims = INTEGER(dimslot),
	ncol, nrow, nnz = length(islot),
	*xj = INTEGER(jslot),
	*xi = INTEGER(islot);

    if (length(jslot) != nnz)
	return mkString(_("lengths of slots i and j must match"));
    if (length(dimslot) != 2)
	return mkString(_("slot Dim must have length 2"));
    nrow = dims[0]; ncol = dims[1];
    for (j = 0; j < nnz; j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
	if (xj[j] < 0 || xj[j] >= ncol)
	    return mkString(_("all column indices must be between 0 and ncol-1"));
    }
    return ScalarLogical(1);
}

SEXP Tsparse_to_Csparse(SEXP x)
{
#ifdef USE_CHOLMOD
    cholmod_triplet *chxt = as_cholmod_triplet(x);
    cholmod_sparse *chxs = cholmod_triplet_to_sparse(chxt, chxt->nnz, &c);

    Free(chxt);
    return chm_sparse_to_SEXP(chxs, 1);
#else
    error("General conversion requires CHOLMOD");
#endif	/* USE_CHOLMOD */
}

