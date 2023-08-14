#include "dgCMatrix.h"
#include "cs_utils.h"
#include "chm_common.h"

/** Return a 2 column matrix  '' cbind(i, j) ''  of 0-origin index vectors (i,j)
 *  which entirely correspond to the (i,j) slots of
 *  as(x, "TsparseMatrix") :
 */
SEXP compressed_non_0_ij(SEXP x, SEXP colP)
{
    int col = asLogical(colP); /* 1 if "C"olumn compressed;  0 if "R"ow */
    SEXP ans, indSym = col ? Matrix_iSym : Matrix_jSym;
    SEXP indP = PROTECT(GET_SLOT(x, indSym)),
	 pP   = PROTECT(GET_SLOT(x, Matrix_pSym));
    int i, *ij;
    int nouter = INTEGER(GET_SLOT(x, Matrix_DimSym))[col ? 1 : 0],
	n_el   = INTEGER(pP)[nouter]; /* is only == length(indP), if the
				     inner slot is not over-allocated */

    ij = INTEGER(ans = PROTECT(allocMatrix(INTSXP, n_el, 2)));
    /* expand the compressed margin to 'i' or 'j' : */
    expand_cmprPt(nouter, INTEGER(pP), &ij[col ? n_el : 0]);
    /* and copy the other one: */
    if (col)
	for(i = 0; i < n_el; i++)
	    ij[i] = INTEGER(indP)[i];
    else /* row compressed */
	for(i = 0; i < n_el; i++)
	    ij[i + n_el] = INTEGER(indP)[i];

    UNPROTECT(3);
    return ans;
}

SEXP dgCMatrix_lusol(SEXP x, SEXP y)
{
    SEXP ycp = PROTECT((TYPEOF(y) == REALSXP) ?
		       duplicate(y) : coerceVector(y, REALSXP));
    CSP xc = AS_CSP__(x);
    R_CheckStack();

    if (xc->m != xc->n || xc->m <= 0)
	error(_("dgCMatrix_lusol requires a square, non-empty matrix"));
    if (LENGTH(ycp) != xc->m)
	error(_("Dimensions of system to be solved are inconsistent"));
    if (!cs_lusol(/*order*/ 1, xc, REAL(ycp), /*tol*/ 1e-7))
	error(_("cs_lusol failed"));

    UNPROTECT(1);
    return ycp;
}

// called from package MatrixModels's R code
SEXP dgCMatrix_qrsol(SEXP x, SEXP y, SEXP ord)
{
    /* FIXME: extend this to work in multivariate case, i.e. y a matrix with > 1 column ! */
    SEXP ycp = PROTECT((TYPEOF(y) == REALSXP) ?
		       duplicate(y) : coerceVector(y, REALSXP));
    CSP xc = AS_CSP(x); /* <--> x  may be  dgC* or dtC* */
    int order = asInteger(ord);
#ifdef _not_yet_do_FIXME__
    const char *nms[] = {"L", "coef", "Xty", "resid", ""};
    SEXP ans = PROTECT(Rf_mkNamed(VECSXP, nms));
#endif
    R_CheckStack();

    if (order < 0 || order > 3)
	error(_("dgCMatrix_qrsol(., order) needs order in {0,..,3}"));
    /* --> cs_amd()  ---  order 0: natural, 1: Chol, 2: LU, 3: QR */
    if (LENGTH(ycp) != xc->m)
	error(_("Dimensions of system to be solved are inconsistent"));
    /* FIXME?  Note that qr_sol() would allow *under-determined systems;
     *		In general, we'd need  LENGTH(ycp) = max(n,m)
     * FIXME also: multivariate y (see above)
     */
    if (xc->m < xc->n || xc->n <= 0)
	error(_("dgCMatrix_qrsol(<%d x %d>-matrix) requires a 'tall' rectangular matrix"),
		xc->m, xc->n);

    /* cs_qrsol(): Tim Davis (2006) .. "8.2 Using a QR factorization", p.136f , calling
     * -------      cs_sqr(order, ..), see  p.76 */
    /* MM: FIXME: write our *OWN* version of - the first case (m >= n) - of cs_qrsol()
     * ---------  which will  (1) work with a *multivariate* y
     *                        (2) compute coefficients properly, not overwriting RHS
     */
    if (!cs_qrsol(order, xc, REAL(ycp)))
	/* return value really is 0 or 1 - no more info there */
	error(_("cs_qrsol() failed inside dgCMatrix_qrsol()"));

    /* Solution is only in the first part of ycp -- cut its length back to n : */
    ycp = lengthgets(ycp, (R_xlen_t) xc->n);

    UNPROTECT(1);
    return ycp;
}

// called from package MatrixModels's R code:
SEXP dgCMatrix_cholsol(SEXP x, SEXP y)
{
    /* Solve Sparse Least Squares X %*% beta ~= y  with dense RHS y,
     * where X = t(x) i.e. we pass  x = t(X)  as argument,
     * via  "Cholesky(X'X)" .. well not really:
     * cholmod_factorize("x", ..) finds L in  X'X = L'L directly */
    CHM_SP cx = AS_CHM_SP(x);
    /* FIXME: extend this to work in multivariate case, i.e. y a matrix with > 1 column ! */
    SEXP y_ = PROTECT(coerceVector(y, REALSXP));
    CHM_DN cy = AS_CHM_DN(y_), rhs, cAns, resid;
    /* const -- but do not fit when used in calls: */
    double one[] = {1,0}, zero[] = {0,0}, neg1[] = {-1,0};
    const char *nms[] = {"L", "coef", "Xty", "resid", ""};
    SEXP ans = PROTECT(Rf_mkNamed(VECSXP, nms));
    R_CheckStack();

    size_t n = cx->ncol;/* #{obs.} {x = t(X) !} */
    if (n < cx->nrow || n <= 0)
	error(_("dgCMatrix_cholsol requires a 'short, wide' rectangular matrix"));
    if (cy->nrow != n)
	error(_("Dimensions of system to be solved are inconsistent"));
    rhs = cholmod_allocate_dense(cx->nrow, 1, cx->nrow, CHOLMOD_REAL, &c);
    /* cholmod_sdmult(A, transp, alpha, beta, X, Y, &c):
     *		Y := alpha*(A*X) + beta*Y or alpha*(A'*X) + beta*Y ;
     * here: rhs := 1 * x %*% y + 0 =  x %*% y =  X'y  */
    if (!(cholmod_sdmult(cx, 0 /* trans */, one, zero, cy, rhs, &c)))
	error(_("cholmod_sdmult error (rhs)"));
    CHM_FR L = cholmod_analyze(cx, &c);
    if (!cholmod_factorize(cx, L, &c))
	error(_("cholmod_factorize failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
/* FIXME: Do this in stages so an "effects" vector can be calculated */
    if (!(cAns = cholmod_solve(CHOLMOD_A, L, rhs, &c)))
	error(_("cholmod_solve (CHOLMOD_A) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    /* L : */
    SET_VECTOR_ELT(ans, 0, chm_factor_to_SEXP(L, 0));
    /* coef : */
            SET_VECTOR_ELT(ans, 1, allocVector(REALSXP,  cx->nrow));
    Memcpy(REAL(VECTOR_ELT(ans, 1)), (double*)(cAns->x), cx->nrow);
    /* X'y : */
/* FIXME: Change this when the "effects" vector is available */
            SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, cx->nrow));
    Memcpy(REAL(VECTOR_ELT(ans, 2)), (double*)(rhs->x), cx->nrow);
    /* resid := y */
    resid = cholmod_copy_dense(cy, &c);
    /* cholmod_sdmult(A, transp, alp, bet, X, Y, &c):
     *		Y := alp*(A*X) + bet*Y or alp*(A'*X) + beta*Y ;
     * here: resid := -1 * x' %*% coef + 1 * y = y - X %*% coef  */
    if (!(cholmod_sdmult(cx, 1/* trans */, neg1, one, cAns, resid, &c)))
	error(_("cholmod_sdmult error (resid)"));
    /* FIXME: for multivariate case, i.e. resid  *matrix* with > 1 column ! */
            SET_VECTOR_ELT(ans, 3,   allocVector(REALSXP, n));
    Memcpy(REAL(VECTOR_ELT(ans, 3)), (double*)(resid->x), n);

    cholmod_free_factor(&L, &c);
    cholmod_free_dense(&resid, &c);
    cholmod_free_dense(&rhs, &c);
    cholmod_free_dense(&cAns, &c);
    UNPROTECT(2);
    return ans;
}

/* MJ: unused */
#if 0

/**
 * Return a SuiteSparse QR factorization of the sparse matrix A
 *
 * @param Ap (pointer to) a [m x n] dgCMatrix
 * @param ordering integer SEXP specifying the ordering strategy to be used
 *	see SPQR/Include/SuiteSparseQR_definitions.h
 * @param econ integer SEXP ("economy"): number of rows of R and columns of Q
 *      to return. The default is m. Using n gives the standard economy form.
 *      A value less than the estimated rank r is set to r, so econ=0 gives the
 *      "rank-sized" factorization, where nrow(R)==nnz(diag(R))==r.
 * @param tol double SEXP: if tol <= -2 use SPQR's default,
 *                         if -2 < tol < 0, then no tol is used; otherwise,
 *      tol > 0, use as tolerance: columns with 2-norm <= tol treated as 0
 *
 *
 * @return SEXP  "SPQR" object with slots (Q, R, p, rank, Dim):
 *	Q: dgCMatrix; R: dgCMatrix  [subject to change to dtCMatrix FIXME ?]
 *	p: integer: 0-based permutation (or length 0 <=> identity);
 *	rank: integer, the "revealed" rank   Dim: integer, original matrix dim.
 */
SEXP dgCMatrix_SPQR(SEXP Ap, SEXP ordering, SEXP econ, SEXP tol)
{
/* SEXP ans = PROTECT(allocVector(VECSXP, 4)); */
    SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS("SPQR"));

    CHM_SP A = AS_CHM_SP(Ap), Q, R;
    SuiteSparse_long *E, rank;/* not always = int   FIXME  (Windows_64 ?) */

    if ((rank = SuiteSparseQR_C_QR(asInteger(ordering),
				   asReal(tol),/* originally had SPQR_DEFAULT_TOL */
				   (SuiteSparse_long)asInteger(econ),/* originally had 0 */
				   A, &Q, &R, &E, &cl)) == -1)
	error(_("SuiteSparseQR_C_QR returned an error code"));

    slot_dup(ans, Ap, Matrix_DimSym);
/*     SET_VECTOR_ELT(ans, 0, */
/* 		   chm_sparse_to_SEXP(Q, 0, 0, 0, "", R_NilValue)); */
    SET_SLOT(ans, Matrix_QSym,
	     chm_sparse_to_SEXP(Q, 0, 0, 0, "", R_NilValue));

    /* Also gives a dgCMatrix (not a dtC* *triangular*) :
     * may make sense if to be used in the "spqr_solve" routines .. ?? */
/*     SET_VECTOR_ELT(ans, 1, */
/* 		   chm_sparse_to_SEXP(R, 0, 0, 0, "", R_NilValue)); */
    SET_SLOT(ans, Matrix_RSym,
	     chm_sparse_to_SEXP(R, 0, 0, 0, "", R_NilValue));
    cholmod_free_sparse(&Al, &cl);
    cholmod_free_sparse(&R, &cl);
    cholmod_free_sparse(&Q, &cl);
    if (E) {
	SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, A->ncol));
	int *Er = INTEGER(VECTOR_ELT(ans, 2));
	for (int i = 0; i < A->ncol; i++) Er[i] = (int) E[i];
	R_Free(E);
    } else SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, 0));
    SET_VECTOR_ELT(ans, 3, ScalarInteger((int)rank));
    UNPROTECT(1);
    return ans;
}

#endif /* MJ */
