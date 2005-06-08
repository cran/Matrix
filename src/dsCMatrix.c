#include "dsCMatrix.h"

SEXP dsCMatrix_validate(SEXP obj)
{
    SEXP val = check_scalar_string(GET_SLOT(obj, Matrix_uploSym),
				   "LU", "uplo");
    int *Dim = INTEGER(GET_SLOT(obj, Matrix_DimSym));

    if (isString(val)) return val;
    if (Dim[0] != Dim[1])
	return mkString(_("Symmetric matrix must be square"));
    csc_check_column_sorting(obj);
    return ScalarLogical(1);
}

SEXP dsCMatrix_chol(SEXP x, SEXP pivot)
{
    SEXP pSlot = GET_SLOT(x, Matrix_pSym), xorig = x;
    int *Ai = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*Ap = INTEGER(pSlot),
	*Lp, *Parent, info,
	lo = CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0] == 'L',
	n = length(pSlot)-1,
	nnz, piv = asLogical(pivot);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dCholCMatrix")));
    int *P, *Pinv = (int *) NULL;
    double *Ax;

    /* FIXME: Check if there is a Cholesky factorization.  If yes,
       check if the permutation status matches that of the call.  If
       so, return it. */

    if (lo) {
	x = PROTECT(ssc_transpose(x));
	Ai = INTEGER(GET_SLOT(x, Matrix_iSym));
	Ap = INTEGER(GET_SLOT(x, Matrix_pSym));
    }
    SET_SLOT(val, Matrix_uploSym, mkString("L"));
    SET_SLOT(val, Matrix_diagSym, mkString("U"));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    Parent = INTEGER(ALLOC_SLOT(val, Matrix_ParentSym, INTSXP, n));
    Lp = INTEGER(ALLOC_SLOT(val, Matrix_pSym, INTSXP, n + 1));
    P = INTEGER(ALLOC_SLOT(val, Matrix_permSym, INTSXP, n));
    if (piv) {
	SEXP trip = PROTECT(dsCMatrix_to_dgTMatrix(x));
	SEXP Ti = GET_SLOT(trip, Matrix_iSym);

	/* determine the permutation with Metis */
	Pinv = Calloc(n, int);
	ssc_metis_order(n, Ap, Ai, P, Pinv);
	/* create a symmetrized form of x */
	nnz = length(Ti);
	Ai = Calloc(nnz, int);
	Ax = Calloc(nnz, double);
	Ap = Calloc(n + 1, int);
	triplet_to_col(n, n, nnz, INTEGER(Ti),
		       INTEGER(GET_SLOT(trip, Matrix_jSym)),
		       REAL(GET_SLOT(trip, Matrix_xSym)),
		       Ap, Ai, Ax);
	UNPROTECT(1);
    } else {
	int i;
	for (i = 0; i < n; i++) P[i] = i;
	Ax = REAL(GET_SLOT(x, Matrix_xSym));

    }
    R_ldl_symbolic(n, Ap, Ai, Lp, Parent, (piv) ? P : (int *)NULL,
		   (piv) ? Pinv : (int *)NULL);
    nnz = Lp[n];
    SET_SLOT(val, Matrix_iSym, allocVector(INTSXP, nnz));
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, nnz));
    SET_SLOT(val, Matrix_DSym, allocVector(REALSXP, n));
    info = R_ldl_numeric(n, Ap, Ai, Ax, Lp, Parent,
			 INTEGER(GET_SLOT(val, Matrix_iSym)),
			 REAL(GET_SLOT(val, Matrix_xSym)),
			 REAL(GET_SLOT(val, Matrix_DSym)),
			 (piv) ? P : (int *)NULL,
			 (piv) ? Pinv : (int *)NULL);
    if (info != n)
	error(_("Leading minor of size %d (possibly after permutation) is indefinite"),
	      info + 1);
    if (piv) {
	Free(Pinv); Free(Ax); Free(Ai); Free(Ap);
    }
    UNPROTECT(lo ? 2 : 1);
    return set_factors(xorig, val, "Cholesky");
}

SEXP dsCMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed)
{
    int cl = asLogical(classed);
    SEXP Chol = get_factors(a, "Cholesky"), perm,
	bdP = cl ? GET_SLOT(b, Matrix_DimSym) : getAttrib(b, R_DimSymbol),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(bdP),
	*Li, *Lp, j, piv;
    int n = adims[1], ncol = bdims[1];
    double *Lx, *D, *in = REAL(cl ? GET_SLOT(b, Matrix_xSym) : b),
	*out = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n * ncol)),
	*tmp = (double *) NULL;

    if (!cl && !(isReal(b) && isMatrix(b)))
	error(_("Argument b must be a numeric matrix"));
    if (*adims != *bdims || ncol < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    if (Chol == R_NilValue) Chol = dsCMatrix_chol(a, ScalarLogical(1));
    SET_SLOT(val, Matrix_DimSym, duplicate(bdP));
    perm = GET_SLOT(Chol, Matrix_permSym);
    piv = length(perm);
    if (piv) tmp = Calloc(n, double);
    Li = INTEGER(GET_SLOT(Chol, Matrix_iSym));
    Lp = INTEGER(GET_SLOT(Chol, Matrix_pSym));
    Lx = REAL(GET_SLOT(Chol, Matrix_xSym));
    D = REAL(GET_SLOT(Chol, Matrix_DSym));
    for (j = 0; j < ncol; j++, in += n, out += n) {
	if (piv) R_ldl_perm(n, out, in, INTEGER(perm));
	else Memcpy(out, in, n);
	R_ldl_lsolve(n, out, Lp, Li, Lx);
	R_ldl_dsolve(n, out, D);
	R_ldl_ltsolve(n, out, Lp, Li, Lx);
	if (piv) R_ldl_permt(n, out, Memcpy(tmp, out, n), INTEGER(perm));
    }
    if (piv) Free(tmp);
    UNPROTECT(1);
    return val;
}

SEXP dsCMatrix_inverse_factor(SEXP A)
{
    return R_NilValue;		/* FIXME: Write this function. */
}

SEXP ssc_transpose(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix"))),
	islot = GET_SLOT(x, Matrix_iSym);
    int nnz = length(islot), *adims,
	*xdims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    adims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    adims[0] = xdims[1]; adims[1] = xdims[0];
    if (CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0] == 'U')
	SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    else
	SET_SLOT(ans, Matrix_uploSym, mkString("U"));
    csc_compTr(xdims[0], xdims[1], nnz,
	       INTEGER(GET_SLOT(x, Matrix_pSym)), INTEGER(islot),
	       REAL(GET_SLOT(x, Matrix_xSym)),
	       INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, xdims[0] + 1)),
	       INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nnz)),
	       REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nnz)));
    UNPROTECT(1);
    return ans;
}

SEXP dsCMatrix_to_dgTMatrix(SEXP x)
{
    SEXP
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgTMatrix"))),
	islot = GET_SLOT(x, Matrix_iSym),
	pslot = GET_SLOT(x, Matrix_pSym);
    int *ai, *aj, *iv = INTEGER(islot),
	j, jj, nnz = length(islot), nout,
	n = length(pslot) - 1,
	*p = INTEGER(pslot), pos;
    double *ax, *xv = REAL(GET_SLOT(x, Matrix_xSym));

    /* increment output count by number of off-diagonals */
    nout = nnz;
    for (j = 0; j < n; j++) {
	int p2 = p[j+1];
	for (jj = p[j]; jj < p2; jj++) {
	    if (iv[jj] != j) nout++;
	}
    }
    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nout));
    ai = INTEGER(GET_SLOT(ans, Matrix_iSym));
    SET_SLOT(ans, Matrix_jSym, allocVector(INTSXP, nout));
    aj = INTEGER(GET_SLOT(ans, Matrix_jSym));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nout));
    ax = REAL(GET_SLOT(ans, Matrix_xSym));
    pos = 0;
    for (j = 0; j < n; j++) {
	int p2 = p[j+1];
	for (jj = p[j]; jj < p2; jj++) {
	    int ii = iv[jj];
	    double xx = xv[jj];

	    ai[pos] = ii; aj[pos] = j; ax[pos] = xx; pos++;
	    if (ii != j) {
		aj[pos] = ii; ai[pos] = j; ax[pos] = xx; pos++;
	    }
	}
    }
    UNPROTECT(1);
    return ans;
}

SEXP dsCMatrix_ldl_symbolic(SEXP x, SEXP doPerm)
{
    SEXP Ax, Dims = GET_SLOT(x, Matrix_DimSym),
	ans = PROTECT(allocVector(VECSXP, 3)), tsc;
    int i, n = INTEGER(Dims)[0], nz, nza,
	*Ap, *Ai, *Lp, *Li, *Parent,
	doperm = asLogical(doPerm),
	*P = (int *) NULL, *Pinv = (int *) NULL;


    if (CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0] == 'L') {
	x = PROTECT(ssc_transpose(x));
    } else {
	x = PROTECT(duplicate(x));
    }
    Ax = GET_SLOT(x, Matrix_xSym);
    nza = length(Ax);
    Ap = INTEGER(GET_SLOT(x, Matrix_pSym));
    Ai = INTEGER(GET_SLOT(x, Matrix_iSym));
    if (doperm) {
	int *perm, *iperm = Calloc(n, int);

	SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, n));
	perm = INTEGER(VECTOR_ELT(ans, 2));
	ssc_metis_order(n, Ap, Ai, perm, iperm);
	ssc_symbolic_permute(n, 1, iperm, Ap, Ai);
	Free(iperm);
    }
    SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, n));
    Parent = INTEGER(VECTOR_ELT(ans, 0));
    SET_VECTOR_ELT(ans, 1, NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    tsc = VECTOR_ELT(ans, 1);
    SET_SLOT(tsc, Matrix_uploSym, mkString("L"));
    SET_SLOT(tsc, Matrix_diagSym, mkString("U"));
    SET_SLOT(tsc, Matrix_DimSym, Dims);
    SET_SLOT(tsc, Matrix_pSym, allocVector(INTSXP, n + 1));
    Lp = INTEGER(GET_SLOT(tsc, Matrix_pSym));
    R_ldl_symbolic(n, Ap, Ai, Lp, Parent, P, Pinv);
    nz = Lp[n];
    SET_SLOT(tsc, Matrix_iSym, allocVector(INTSXP, nz));
    Li = INTEGER(GET_SLOT(tsc, Matrix_iSym));
    SET_SLOT(tsc, Matrix_xSym, allocVector(REALSXP, nz));
    for (i = 0; i < nza; i++) REAL(Ax)[i] = 0.00001;
    for (i = 0; i < n; i++) REAL(Ax)[Ap[i+1]-1] = 10000.;
    i = R_ldl_numeric(n, Ap, Ai, REAL(Ax), Lp, Parent, Li,
		    REAL(GET_SLOT(tsc, Matrix_xSym)),
		    (double *) R_alloc(n, sizeof(double)), /* D */
		    P, Pinv);
    UNPROTECT(2);
    return ans;
}

SEXP dsCMatrix_metis_perm(SEXP x)
{
    SEXP pSlot = GET_SLOT(x, Matrix_pSym),
	ans = PROTECT(allocVector(VECSXP, 2));
    int n = length(pSlot) - 1;

    SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, n));
    SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, n));
    ssc_metis_order(n,
		    INTEGER(pSlot),
		    INTEGER(GET_SLOT(x, Matrix_iSym)),
		    INTEGER(VECTOR_ELT(ans, 0)),
		    INTEGER(VECTOR_ELT(ans, 1)));
    UNPROTECT(1);
    return ans;
}

