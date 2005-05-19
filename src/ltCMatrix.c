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
    SEXP val = check_scalar_string(GET_SLOT(x, Matrix_uploSym),
				   "LU", "uplo");
    int *Dim = INTEGER(GET_SLOT(x, Matrix_DimSym));

    if (isString(val)) return val;
    if (isString(val = check_scalar_string(GET_SLOT(x, Matrix_diagSym),
					   "NU", "diag"))) return val;
    if (Dim[0] != Dim[1])
	return mkString(_("Symmetric matrix must be square"));
    csc_check_column_sorting(x);
    return ScalarLogical(1);
}

/** 
 * Transpose an ltCMatrix
 * 
 * @param x Pointer to an ltCMatrix object
 * 
 * @return the transpose of x.  It represents the same matrix but is
 * stored in the opposite triangle.
 */
SEXP ltCMatrix_trans(SEXP x)
{
    SEXP Xi = GET_SLOT(x, Matrix_iSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("ltCMatrix"))),
	xdn = GET_SLOT(x, Matrix_DimNamesSym);
    SEXP adn = ALLOC_SLOT(ans, Matrix_DimNamesSym, VECSXP, 2);
    int *adims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	*xdims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	up = CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0] == 'U';
    int m = xdims[0], n = xdims[1], nz = length(Xi);
    int *xj = expand_cmprPt(n, INTEGER(GET_SLOT(x, Matrix_pSym)),
			    Calloc(nz, int));

    adims[0] = n; adims[1] = m;
    SET_VECTOR_ELT(adn, 0, VECTOR_ELT(xdn, 1));
    SET_VECTOR_ELT(adn, 1, VECTOR_ELT(xdn, 0));
    SET_SLOT(ans, Matrix_uploSym, mkString(up ? "L" : "U"));
    SET_SLOT(ans, Matrix_diagSym, duplicate(GET_SLOT(x, Matrix_diagSym)));
    triplet_to_col(n, m, nz, xj, INTEGER(Xi), (double *) NULL,
		   INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP,  m + 1)),
		   INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP,  nz)),
		   (double *) NULL);
    Free(xj);
    UNPROTECT(1);
    return ans;
}

/** 
 * Solve  one  of the matrix equations  op(A)*C = B, or
 * C*op(A) = B where A is a square ltCMatrix and B and C are lgCMatrix
 * objects.
 * 
 * @param side LFT or RGT
 * @param transa TRN or NTR
 * @param A pointer to an ltCMatrix object
 * @param B pointer to an lgCMatrix object
 * @param C pointer to an lgCMatrix object
 */
void
ltClgCsm(enum CBLAS_SIDE side, enum CBLAS_TRANSPOSE transa,
	 SEXP A, SEXP B, SEXP C)
{
    error(_("code not yet written"));
}
