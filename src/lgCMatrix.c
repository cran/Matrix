#include "lgCMatrix.h"

SEXP lgCMatrix_validate(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym);
    int j,
	ncol = length(pslot) - 1,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow,
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    nrow = dims[0];
    if (length(pslot) <= 0)
	return mkString(_("slot p must have length > 0"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return mkString(_("last element of slot p must match length of slot i"));
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
    }
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    if (csc_unsorted_columns(ncol, xp, xi))
	csc_sort_columns(ncol, xp, xi, (double *) NULL);

    return ScalarLogical(1);
}

/* very parallel to csc_to_matrix() in ./dgCMatrix.c */
SEXP lcsc_to_matrix(SEXP x)
{
    SEXP ans, pslot = GET_SLOT(x, Matrix_pSym);
    int j, ncol = length(pslot) - 1,
	nrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	*xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    int *ax;

    ax = LOGICAL(ans = PROTECT(allocMatrix(LGLSXP, nrow, ncol)));
    for (j = 0; j < (nrow * ncol); j++) ax[j] = 0;
    for (j = 0; j < ncol; j++) {
	int ind;
	for (ind = xp[j]; ind < xp[j+1]; ind++)
	    ax[j * nrow + xi[ind]] = 1;
    }
    UNPROTECT(1);
    return ans;
}

#ifdef _NEED_logical_to_csc_FIRST_
/* very parallel to matrix_to_csc() in ./dgCMatrix.c */
SEXP matrix_to_lcsc(SEXP A)
{
    if (!(isMatrix(A) && isLogical(A)))
	error(_("A must be a logical matrix"));
    return logical_to_csc(LOGICAL(A),
			  INTEGER(getAttrib(A, R_DimSymbol)));
}
#endif
