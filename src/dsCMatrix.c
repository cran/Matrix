#include "dsCMatrix.h"

SEXP dsCMatrix_chol(SEXP x, SEXP pivot)
{
    cholmod_factor
	*N = as_cholmod_factor(dsCMatrix_Cholesky(x, pivot,
						  ScalarLogical(FALSE),
						  ScalarLogical(FALSE)));
    /* Must use a copy; cholmod_factor_to_sparse modifies first arg. */
    cholmod_factor *Ncp = cholmod_copy_factor(N, &c);
    cholmod_sparse *L, *R;
    SEXP ans;

    L = cholmod_factor_to_sparse(Ncp, &c); cholmod_free_factor(&Ncp, &c);
    R = cholmod_transpose(L, /*values*/ 1, &c); cholmod_free_sparse(&L, &c);
    ans = PROTECT(chm_sparse_to_SEXP(R, /*cholmod_free*/ 1,
				     /*uploT*/ 1, /*Rkind*/ 0, /*diag*/ "N",
				     GET_SLOT(x, Matrix_DimNamesSym)));
    if (asLogical(pivot)) {
	SEXP piv = PROTECT(allocVector(INTSXP, N->n));
	int *dest = INTEGER(piv), *src = (int*)N->Perm, i;

	for (i = 0; i < N->n; i++) dest[i] = src[i] + 1;
	setAttrib(ans, install("pivot"), piv);
	/* FIXME: Because of the cholmod_factor -> S4 obj ->
	 * cholmod_factor conversions, the value of N->minor will
	 * always be N->n.  Change as_cholmod_factor and
	 * chm_factor_to_SEXP to keep track of Minor.
	 */
	setAttrib(ans, install("rank"), ScalarInteger((size_t) N->minor));
	UNPROTECT(1);
    }
    Free(N);
    UNPROTECT(1);
    return ans;
}

SEXP dsCMatrix_Cholesky(SEXP Ap, SEXP permP, SEXP LDLp, SEXP superP)
{
    char fname[12] = "spdCholesky"; /* template for factorization name */
    /* S|s : super or not
     * P|p : permuted or not
     * D|d :  LDL' or not (= LL')
     */
    int perm = asLogical(permP),
	LDL = asLogical(LDLp),
	super = asLogical(superP);
    SEXP Chol;
    cholmod_sparse *A;
    cholmod_factor *L;
    int sup, ll;

    if (super) fname[0] = 'S';
    if (perm) fname[1] = 'P';
    if (LDL) fname[2] = 'D';
    Chol = get_factors(Ap, fname);
    /* If Ap has already cached the factor, we return it immediately */
    if (Chol != R_NilValue) return Chol;
    A = as_cholmod_sparse(Ap);
    if (!A->stype)
	error("Non-symmetric matrix passed to dsCMatrix_chol");

    sup = c.supernodal;
    ll = c.final_ll;

    c.final_ll = !LDL;	/* leave as LL' or form LDL' */
    c.supernodal = super ? CHOLMOD_SUPERNODAL : CHOLMOD_SIMPLICIAL;

    if (perm) {
	L = cholmod_analyze(A, &c); /* get fill-reducing permutation */
    } else {			    /* require identity permutation */
	int nmethods = c.nmethods, ord0 = c.method[0].ordering,
	    postorder = c.postorder;
	c.nmethods = 1;
	c.method[0].ordering = CHOLMOD_NATURAL;
	c.postorder = FALSE;
	L = cholmod_analyze(A, &c);
	c.nmethods = nmethods; c.method[0].ordering = ord0;
	c.postorder = postorder;
    }
    if (!cholmod_factorize(A, L, &c))
	error(_("Cholesky factorization failed"));
    c.supernodal = sup;	/* restore previous setting */
    c.final_ll = ll;
    Free(A);
    Chol = set_factors(Ap, chm_factor_to_SEXP(L, 1), fname);
    return Chol;
}

static
SEXP get_factor_pattern(SEXP obj, char *pat, int offset)
{
    SEXP facs = GET_SLOT(obj, Matrix_factorSym), nms;
    int i;

    if (!LENGTH(facs)) /* Fast return for empty factor slot: */
	return R_NilValue;
    nms = getAttrib(facs, R_NamesSymbol);
    for (i = 0; i < LENGTH(nms); i++) {
	char *nm = CHAR(STRING_ELT(nms, i));
	if (strlen(nm) > offset && !strcmp(pat + offset, nm + offset))
	    return VECTOR_ELT(facs, i);
    }
    return R_NilValue;
}

SEXP dsCMatrix_Csparse_solve(SEXP a, SEXP b)
{
    SEXP Chol = get_factor_pattern(a, "...Cholesky", 3);
    cholmod_factor *L;
    cholmod_sparse *cx, *cb = as_cholmod_sparse(b);

    if (Chol == R_NilValue) /* compute (and cache) "sPDCholesky" */
	Chol = dsCMatrix_Cholesky(a,
				  ScalarLogical(1),  /* permuted  : "P" */
				  ScalarLogical(1),  /* LDL'	  : "D" */
				  ScalarLogical(0)); /* simplicial: "s" */
    L = as_cholmod_factor(Chol);
    cx = cholmod_spsolve(CHOLMOD_A, L, cb, &c);
    Free(cb); Free(L);
    return chm_sparse_to_SEXP(cx, /*cholmod_free*/ 1, /*uploT*/ 0,
			      /*Rkind*/ 0, /*diag*/ "N",
			      /*dimnames = */ R_NilValue);
}

SEXP dsCMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP Chol = get_factor_pattern(a, "...Cholesky", 3);
    cholmod_factor *L;
    cholmod_dense  *cx,
	*cb = as_cholmod_dense(PROTECT(mMatrix_as_dgeMatrix(b)));

    if (Chol == R_NilValue) /* compute (and cache) "sPDCholesky" */
	Chol = dsCMatrix_Cholesky(a,
				  ScalarLogical(1),  /* permuted  : "P" */
				  ScalarLogical(1),  /* LDL'      : "D" */
				  ScalarLogical(0)); /* simplicial: "s" */
    L = as_cholmod_factor(Chol);
    cx = cholmod_solve(CHOLMOD_A, L, cb, &c);
    Free(cb); Free(L);
    UNPROTECT(1);
    return chm_dense_to_SEXP(cx, 1, 0, /*dimnames = */ R_NilValue);
}

/* Needed for printing dsCMatrix objects */
/* FIXME: Create a more general version of this operation: also for lsC, (dsR?),..
*         e.g. make  compressed_to_dgTMatrix() in ./dgCMatrix.c work for dsC */
SEXP dsCMatrix_to_dgTMatrix(SEXP x)
{
    cholmod_sparse *A = as_cholmod_sparse(x);
    cholmod_sparse *Afull = cholmod_copy(A, /*stype*/ 0, /*mode*/ 1, &c);
    cholmod_triplet *At = cholmod_sparse_to_triplet(Afull, &c);

    if (!A->stype)
	error("Non-symmetric matrix passed to dsCMatrix_to_dgTMatrix");
    Free(A); cholmod_free_sparse(&Afull, &c);
    return chm_triplet_to_SEXP(At, 1, /*uploT*/ 0, /*Rkind*/ 0, "",
			       GET_SLOT(x, Matrix_DimNamesSym));
}
