#include "lgCMatrix.h"

SEXP lgCMatrix_validate(SEXP x)
{
    /* Almost all is now done in Csparse_validate
     * *but* the checking of the 'x' slot */
    SEXP islot = GET_SLOT(x, Matrix_iSym),
	xslot = GET_SLOT(x, Matrix_xSym);

    if (length(islot) != length(xslot))
	return mkString(_("lengths of slots 'i' and 'x' must match"));

    return ScalarLogical(1);
}

SEXP lcsc_to_matrix(SEXP x)
{
    SEXP ans, pslot = GET_SLOT(x, Matrix_pSym);
    int j, ncol = length(pslot) - 1,
	nrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	*xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    int *xx = LOGICAL(GET_SLOT(x, Matrix_xSym)), *ax;

    ax = LOGICAL(ans = PROTECT(allocMatrix(LGLSXP, nrow, ncol)));
    for (j = 0; j < (nrow * ncol); j++) ax[j] = 0;
    for (j = 0; j < ncol; j++) {
	int ind;
	for (ind = xp[j]; ind < xp[j+1]; ind++)
	    ax[j * nrow + xi[ind]] = xx[ind];
    }
    UNPROTECT(1);
    return ans;
}

/* as above,  '1' instead of 'x' slot: */
SEXP ncsc_to_matrix(SEXP x)
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
