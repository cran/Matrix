				/* Sparse triangular matrices in triplet format */
#include "dtTMatrix.h"

SEXP dtTMatrix_validate(SEXP x)
{
    SEXP val;

    if (isString(val = check_scalar_string(GET_SLOT(x, Matrix_uploSym),
					   "LU", "uplo"))) return val;
    if (isString(val = check_scalar_string(GET_SLOT(x, Matrix_diagSym),
					   "NU", "diag"))) return val;
    return ScalarLogical(1);
}

SEXP dtTMatrix_as_dtrMatrix(SEXP x)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))),
	DimP = GET_SLOT(x, Matrix_DimSym),
	xiP = GET_SLOT(x, Matrix_iSym);
    int k, m = INTEGER(DimP)[0], n = INTEGER(DimP)[1], nnz = length(xiP);
    int *xi = INTEGER(xiP), *xj = INTEGER(GET_SLOT(x, Matrix_jSym)),
	sz = m * n;
    double *tx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, sz)),
	*xx = REAL(GET_SLOT(x, Matrix_xSym));
    
    SET_SLOT(val, Matrix_DimSym, duplicate(DimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(x, Matrix_uploSym)));
    SET_SLOT(val, Matrix_diagSym, duplicate(GET_SLOT(x, Matrix_diagSym)));
    AZERO(tx, sz);
    for (k = 0; k < nnz; k++) tx[xi[k] + xj[k] * m] = xx[k];
    UNPROTECT(1);
    return val;
}
