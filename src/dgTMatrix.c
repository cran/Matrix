/* Sparse matrices in triplet form */
#include "dgTMatrix.h"

SEXP dgTMatrix_validate(SEXP x)
{
    SEXP
	islot = GET_SLOT(x, Matrix_iSym),
	jslot = GET_SLOT(x, Matrix_jSym),
	xslot = GET_SLOT(x, Matrix_xSym),
	dimslot = GET_SLOT(x, Matrix_DimSym);
    int j,
	*dims = INTEGER(dimslot),
	ncol, nrow, nnz = length(islot),
	*xj = INTEGER(jslot),
	*xi = INTEGER(islot);

    if (length(xslot) != nnz || length(jslot) != nnz)
	return mkString(_("lengths of slots i, j, and x must match"));
    if (length(dimslot) != 2)
	return mkString(_("slot Dim must have length 2"));
    nrow = dims[0]; ncol = dims[1];
    for (j = 0; j < nnz; j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
	if (xj[j] < 0 || xj[j] >= ncol)
	    return mkString(_("all column indices must be between 0 and ncol-1"));
    }
    return ScalarLogical(1);
}

SEXP dgTMatrix_to_dgCMatrix(SEXP x)
{
    SEXP dd = GET_SLOT(x, Matrix_DimSym),
	iP = GET_SLOT(x, Matrix_iSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    int *dims = INTEGER(dd), nnz = length(iP);
    int *p, *ti = Calloc(nnz, int), m = dims[0], n = dims[1];
    double *tx = Calloc(nnz, double);

    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, n + 1));
    SET_SLOT(ans, Matrix_DimSym, duplicate(dd));
    p = INTEGER(GET_SLOT(ans, Matrix_pSym));
    triplet_to_col(m, n, nnz, INTEGER(iP),
		   INTEGER(GET_SLOT(x, Matrix_jSym)),
		   REAL(GET_SLOT(x, Matrix_xSym)),
		   p, ti, tx);
    nnz = p[n];
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    Memcpy(INTEGER(GET_SLOT(ans, Matrix_iSym)), ti, nnz);
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    Memcpy(REAL(GET_SLOT(ans, Matrix_xSym)), tx, nnz);

    Free(ti); Free(tx);
    UNPROTECT(1);
    return ans;
}

static void
insert_triplets_in_array(int m, int n, int nnz,
			 const int xi[], const int xj[], const double xx[],
			 double vx[])
{
    int i;
    memset(vx, 0, sizeof(double) * m * n);
    for (i = 0; i < nnz; i++) {
	vx[xi[i] + xj[i] * m] += xx[i];	/* allow redundant entries in x */
    }
}

SEXP dgTMatrix_to_dgeMatrix(SEXP x)
{
    SEXP dd = GET_SLOT(x, Matrix_DimSym),
	islot = GET_SLOT(x, Matrix_iSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));

    int *dims = INTEGER(dd),
	m = dims[0],
	n = dims[1];

    SET_SLOT(ans, Matrix_rcondSym, allocVector(REALSXP, 0));
    SET_SLOT(ans, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(ans, Matrix_DimSym, duplicate(dd));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, m * n));
    insert_triplets_in_array(m, n, length(islot),
			     INTEGER(islot), INTEGER(GET_SLOT(x, Matrix_jSym)),
			     REAL(GET_SLOT(x, Matrix_xSym)),
			     REAL(GET_SLOT(ans, Matrix_xSym)));
    UNPROTECT(1);
    return ans;
}

SEXP dgTMatrix_to_matrix(SEXP x)
{
    SEXP dd = GET_SLOT(x, Matrix_DimSym),
	islot = GET_SLOT(x, Matrix_iSym);
    int m = INTEGER(dd)[0],
	n = INTEGER(dd)[1];
    SEXP ans = PROTECT(allocMatrix(REALSXP, m, n));

    insert_triplets_in_array(m, n, length(islot),
			     INTEGER(islot), INTEGER(GET_SLOT(x, Matrix_jSym)),
			     REAL(GET_SLOT(x, Matrix_xSym)),
			     REAL(ans));
    UNPROTECT(1);
    return ans;
}
