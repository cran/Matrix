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

SEXP tsc_to_triplet(SEXP x)
{
    SEXP ans;
    if (toupper(CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0))[0]) != 'U')
	ans = csc_to_triplet(x);
    else {			/* unit triangular matrix */
	SEXP islot = GET_SLOT(x, Matrix_iSym), 
	    pslot = GET_SLOT(x, Matrix_pSym);
	int *ai, *aj, j,
	    n = length(pslot) - 1,
	    nod = length(islot),
	    nout = n + nod,
	    *p = INTEGER(pslot);
	double *ax;
    
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("tripletMatrix")));
	SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
	SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nout));
	ai = INTEGER(GET_SLOT(ans, Matrix_iSym));
	Memcpy(ai, INTEGER(islot), nod);
	SET_SLOT(ans, Matrix_jSym, allocVector(INTSXP, nout));
	aj = INTEGER(GET_SLOT(ans, Matrix_jSym));
	SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nout));
	ax = REAL(GET_SLOT(ans, Matrix_xSym));
	Memcpy(ax, REAL(GET_SLOT(x, Matrix_xSym)), nod);
	for (j = 0; j < n; j++) {
	    int jj, npj = nod + j,  p2 = p[j+1];
	    aj[npj] = ai[npj] = j;
	    ax[npj] = 1.;
	    for (jj = p[j]; jj < p2; jj++) aj[jj] = j;
	}
	UNPROTECT(1);
    }
    return ans;
}
