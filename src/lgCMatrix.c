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
    if (csc_unsorted_columns(ncol, xp, xi)) {
      csc_sort_columns(ncol, xp, xi, (double *) NULL);
    }
    return ScalarLogical(1);
}

/** 
 * C := op(A) %*% op(B) + beta ^ C for logical sparse column-oriented matrices
 * 
 * @param tra nonzero if A is to be transposed
 * @param trb nonzero if B is to be transposed
 * @param m number of rows in C
 * @param n number of columns in C
 * @param k number of columns in A if tra == 0, otherwise number of
 *          rows in A
 * @param ai vector of row indices of TRUE elements in A
 * @param ap column pointers for A
 * @param bi vector of row indices of TRUE elements in B
 * @param bp column pointers for B
 * @param beta if non-zero existing TRUE elements in C are retained
 * @param ciP pointer to the column indices of TRUE elements in C.
 *        This may be reallocated if the number of TRUE elements in C changes
 * @param cp column pointers for C
 */
void Matrix_lgClgCmm(int tra, int trb, int m, int n, int k,
		     const int ai[], const int ap[],
		     const int bi[], const int bp[],
		     int beta, SEXP *CIP, int cp[])
{
    int annz = ap[tra ? m : k], bnnz = bp[trb ? k : n], cnnz = cp[n], extra = 0;
    int *ci = INTEGER(*CIP);
    
    if (beta) {
	if (tra) {
	    if (trb) {		/* t(A) %*% t(B) */
	    } else {		/* t(A) %*% B */
		
	    }
	} else {
	    if (trb) {		/* A %*% t(B) */
	    } else {		/* A %*% B */
		
	    }
	}
    }
}

SEXP lgCMatrix_lgCMatrix_mm(SEXP a, SEXP b, SEXP transa, SEXP transb)
{
    return R_NilValue;
}
