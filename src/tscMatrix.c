				/* Sparse triangular matrices */
#include "tscMatrix.h"

SEXP tsc_validate(SEXP x)
{
    return ScalarLogical(1);
}

SEXP tsc_transpose(SEXP x)
{
    SEXP
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("tscMatrix"))),
	islot = GET_SLOT(x, Matrix_iSym);
    int nnz = length(islot),
	*adims = INTEGER(GET_SLOT(ans, Matrix_DimSym)),
	*xdims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    adims[0] = xdims[1]; adims[1] = xdims[0];
    if (toupper(CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0]) == 'U')
	SET_SLOT(ans, Matrix_uploSym, ScalarString(mkChar("L")));
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, xdims[0] + 1));
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    csc_components_transpose(xdims[0], xdims[1], nnz,
			     INTEGER(GET_SLOT(x, Matrix_pSym)),
			     INTEGER(islot),
			     REAL(GET_SLOT(x, Matrix_xSym)),
			     INTEGER(GET_SLOT(ans, Matrix_pSym)),
			     INTEGER(GET_SLOT(ans, Matrix_iSym)),
			     REAL(GET_SLOT(ans, Matrix_xSym)));
    UNPROTECT(1);
    return ans;
}
