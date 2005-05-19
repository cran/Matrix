				/* Logical, sparse matrices in triplet form */
#include "lgTMatrix.h"

SEXP lgTMatrix_validate(SEXP x)
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

SEXP lgTMatrix_as_lgCMatrix(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lgCMatrix"))),
	xDim = GET_SLOT(x, Matrix_DimSym),
        xiP = GET_SLOT(x, Matrix_iSym);
    int m = INTEGER(xDim)[0], n = INTEGER(xDim)[1], nz = length(xiP);

    SET_SLOT(ans, Matrix_DimSym, duplicate(xDim));
    SET_SLOT(ans, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(x, Matrix_DimNamesSym)));
    triplet_to_col(m, n, nz, INTEGER(xiP),
		   INTEGER(GET_SLOT(x, Matrix_jSym)), (double *) NULL,
		   INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1)),
		   INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)),
		   (double *) NULL);
    UNPROTECT(1);
    return ans;
}
