#include "dgTMatrix.h"

#include "chm_common.h"
#include "Tsparse.h"

SEXP dgTMatrix_validate(SEXP x)
{
    SEXP
	islot = GET_SLOT(x, Matrix_iSym),
	xslot = GET_SLOT(x, Matrix_xSym);

    if (LENGTH(xslot) != LENGTH(islot))
	return mkString(_("lengths of slots i and x must match"));
    return ScalarLogical(1);
}

SEXP dgTMatrix_to_dgCMatrix(SEXP x)
{
    return Tsparse_to_Csparse(x);
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

SEXP graphNEL_as_dgTMatrix(SEXP x, SEXP symmetric)
{
    int sym = asLogical(symmetric);
    SEXP nodes = GET_SLOT(x, install("nodes")),
	edgeL = GET_SLOT(x, install("edgeL")),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS(sym
					    ? "dsTMatrix"
					    : "dgTMatrix")));
    int *ii, *jj, *dims, i, j, nnd = LENGTH(nodes), pos, totl;
    double *xx;

    totl = 0;
    for (i = 0; i < nnd; i++)
	totl += LENGTH(Matrix_getElement(VECTOR_ELT(edgeL, i), "edges"));
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = dims[1] = nnd;
    if (isString(nodes)) {
	SEXP dnms = ALLOC_SLOT(ans, Matrix_DimNamesSym, VECSXP, 2);
	SET_VECTOR_ELT(dnms, 0, duplicate(nodes));
	SET_VECTOR_ELT(dnms, 1, duplicate(nodes));
    }
    ii = Calloc(totl, int);
    jj = Calloc(totl, int);
    xx = Calloc(totl, double);
    pos = 0;
    for (i = 0; i < nnd; i++) {
	SEXP edg = VECTOR_ELT(edgeL, i);
	SEXP edges = Matrix_getElement(edg, "edges"),
	    weights = Matrix_getElement(edg, "weights");
	int *edgs = INTEGER(PROTECT(coerceVector(edges, INTSXP))),
	    nedg = LENGTH(edges);
	double *wts = REAL(weights);

	for (j = 0; j < nedg; j++) {
	    int j1 = edgs[j] - 1;
			/* symmetric case stores upper triangle only */
	    if ((!sym) || i <= j1) {
		ii[pos] = i;
		jj[pos] = j1;
		xx[pos] = wts[j];
		pos++;
	    }
	}
	UNPROTECT(1);
    }
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, pos)), ii, pos);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_jSym, INTSXP, pos)), jj, pos);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, pos)), xx, pos);

    Free(ii); Free(jj); Free(xx);
    UNPROTECT(1);
    return ans;
}
