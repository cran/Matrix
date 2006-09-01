#include "lsCMatrix.h"

/**
 * Check the validity of the slots of an lsCMatrix object
 *
 * @param x Pointer to an lsCMatrix object
 *
 * @return an SEXP that is either TRUE or a character string
 * describing the way in which the object failed the validity check
 */
SEXP lsCMatrix_validate(SEXP x)
{
    SEXP val = symmetricMatrix_validate(x);
    if(isString(val))
	return(val);
    else {
	/* FIXME needed? ltC* inherits from lgC* which does this in validate*/
	SEXP pslot = GET_SLOT(x, Matrix_pSym),
	    islot = GET_SLOT(x, Matrix_iSym);
	int
	    ncol = length(pslot) - 1,
	    *xp = INTEGER(pslot),
	    *xi = INTEGER(islot);

	if (csc_unsorted_columns(ncol, xp, xi))
	    csc_sort_columns(ncol, xp, xi, (double *) NULL);

	return ScalarLogical(1);
    }
}

#if 0				/* no longer used */
/**
 * Transpose an lsCMatrix
 *
 * @param x Pointer to an lsCMatrix object
 *
 * @return the transpose of x.  It represents the same matrix but is
 * stored in the opposite triangle.
 */
SEXP lsCMatrix_trans(SEXP x)
{
    SEXP Xi = GET_SLOT(x, Matrix_iSym),
	xDim = GET_SLOT(x, Matrix_DimSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lsCMatrix")));
    int n = INTEGER(xDim)[0], nz = length(Xi);
    int *xj = expand_cmprPt(n, INTEGER(GET_SLOT(x, Matrix_pSym)),
			    Calloc(nz, int)),
	*xi = Memcpy(Calloc(nz, int), Xi, nz);
    int up = uplo_P(x)[0] == 'U';

    SET_SLOT(ans, Matrix_DimSym, duplicate(xDim));
    SET_SLOT(ans, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(x, Matrix_DimNamesSym)));
    SET_SLOT(ans, Matrix_uploSym, mkString(up ? "L" : "U"));
    make_upper_triangular(up ? xj : xi, up ? xi : xj, nz);
    triplet_to_col(n, n, nz, xi, xj, (double *) NULL,
		   INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP,  n + 1)),
		   INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP,  nz)),
		   (double *) NULL);
    Free(xj); Free(xi);
    UNPROTECT(1);
    return ans;
}
#endif
