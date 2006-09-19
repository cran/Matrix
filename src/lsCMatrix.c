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
	/* Almost all is now done in Csparse_validate
	 * *but* the checking of the 'x' slot */
	SEXP islot = GET_SLOT(x, Matrix_iSym),
	    xslot = GET_SLOT(x, Matrix_xSym);

	if (length(islot) != length(xslot))
	    return mkString(_("lengths of slots 'i' and 'x' must match"));

	return ScalarLogical(1);
    }
}
