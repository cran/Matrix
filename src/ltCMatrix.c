				/* Sparse triangular logical matrices */
#include "ltCMatrix.h"

/**
 * Check the validity of the slots of an ltCMatrix object
 *
 * @param x Pointer to an ltCMatrix object
 *
 * @return an SEXP that is either TRUE or a character string
 * describing the way in which the object failed the validity check
 */

SEXP ltCMatrix_validate(SEXP x)
{
/* ltCMatrix extends both CsparseMatrix and triangularMatrix
 * which have their own _validate() each  */

    /* Almost all is now done in Csparse_validate
     * *but* the checking of the 'x' slot */
    SEXP islot = GET_SLOT(x, Matrix_iSym),
	xslot = GET_SLOT(x, Matrix_xSym);

    if (length(islot) != length(xslot))
	return mkString(_("lengths of slots 'i' and 'x' must match"));

    return ScalarLogical(1);
}
