/* Sparse matrices in triplet form */
#include "triplet.h"

SEXP triplet_validate(SEXP x)
{
    SEXP 
	islot = GET_SLOT(x, Matrix_iSym),
	jslot = GET_SLOT(x, Matrix_jSym),	
	xslot = GET_SLOT(x, Matrix_xSym),
	dimslot = GET_SLOT(x, Matrix_DimSym);
    int j,
	*dims = INTEGER(dimslot),
	ncol, nrow, nnz = length(islot),
	*xj = INTEGER(jslot),
	*xi = INTEGER(islot);
    
    if (length(xslot) != nnz || length(jslot) != nnz)
	return ScalarString(mkChar("lengths of slots i, j, and x must match"));
    if (length(dimslot) != 2)
	return ScalarString(mkChar("slot Dim must have length 2"));
    nrow = dims[0]; ncol = dims[1];
    for (j = 0; j < nnz; j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return ScalarString(
		mkChar("all row indices must be between 0 and nrow-1"));
	if (xj[j] < 0 || xj[j] >= ncol)
	    return ScalarString(
		mkChar("all column indices must be between 0 and ncol-1"));
    }
    return ScalarLogical(1);
}
