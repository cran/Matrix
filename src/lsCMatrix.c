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
    SEXP val = check_scalar_string(GET_SLOT(x, Matrix_uploSym),
				   "LU", "uplo");
    int *Dim = INTEGER(GET_SLOT(x, Matrix_DimSym));

    if (isString(val)) return val;
    if (Dim[0] != Dim[1])
	return mkString(_("Symmetric matrix must be square"));
    csc_check_column_sorting(x);
    return ScalarLogical(1);
}

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
    int up = CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0] == 'U';

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

SEXP lsCMatrix_chol(SEXP x, SEXP pivot)
{
    int piv = asLogical(pivot);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lCholCMatrix")));
    int j, n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    int *Xi = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*Xp = INTEGER(GET_SLOT(x, Matrix_pSym)),
	*P, *Pinv, *Parent, *Lp;
    double *D = Calloc(n, double), *Tx, *Xx = Calloc(Xp[n], double);

    if (CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0] != 'U')
	error(_("Must have uplo == 'U' in x argument to lsCMatrix_chol"));
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_diagSym, mkString("U"));
    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    SET_SLOT(ans, Matrix_DimNamesSym, duplicate(GET_SLOT(x, Matrix_DimNamesSym)));
    P = INTEGER(ALLOC_SLOT(ans, Matrix_permSym, INTSXP, n));
    if (piv) {
	Pinv = Calloc(n, int);
	ssc_metis_order(n, Xp, Xi, P, Pinv);
    } else {
	int i;
	for (i = 0; i < n; i++) P[i] = i;
    }
    Lp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1));
    Parent = INTEGER(ALLOC_SLOT(ans, Matrix_ParentSym, INTSXP, n));
    R_ldl_symbolic(n, Xp, Xi, Lp, Parent,
		   (piv) ? P : (int *) NULL, (piv) ? Pinv : (int *) NULL);
    				/* Decompose the identity to get Li */
    for (j = 0; j < n; j++) {	/* Create an identity from Xp, Xi and Xx */
	int ii, ii2 = Xp[j + 1];
	for (ii = Xp[j]; ii < ii2; ii++)
	    Xx[ii] = (Xi[ii] == j) ? 1. : 0.;
    }
    Tx = Calloc(Lp[n], double);
    R_ldl_numeric(n, Xp, Xi, Xx, Lp, Parent,
		  INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, Lp[n])),
		  Tx, D, (piv) ? P : (int *) NULL,
		  (piv) ? Pinv : (int *) NULL);
    if (piv) Free(Pinv);
    Free(Xx); Free(Tx); Free(D);
    UNPROTECT(1);
    return ans;
}
