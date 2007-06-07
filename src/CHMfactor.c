				/* CHOLMOD factors */
#include "CHMfactor.h"

SEXP CHMfactor_to_sparse(SEXP x)
{
    cholmod_factor *L = as_cholmod_factor(x), *Lcp;
    cholmod_sparse *Lm;

    Lcp = cholmod_copy_factor(L, &c); Free(L); /* next call changes Lcp */
    if (!(Lcp->is_ll))
	if (!cholmod_change_factor(Lcp->xtype, 1, 0, 1, 1, Lcp, &c))
	    error(_("cholmod_change_factor failed with status %d"), c.status);
    Lm = cholmod_factor_to_sparse(Lcp, &c); cholmod_free_factor(&Lcp, &c);
    return chm_sparse_to_SEXP(Lm, 1, -1, /*Rkind*/ 0, "N", R_NilValue);
}

SEXP CHMfactor_solve(SEXP a, SEXP b, SEXP system)
{
    cholmod_factor *L = as_cholmod_factor(a);
    SEXP bb = PROTECT(dup_mMatrix_as_dgeMatrix(b));
    cholmod_dense *B = as_cholmod_dense(bb), *X;
    int sys = asInteger(system);

    if (!(sys--))		/* -- to align with CHOLMOD definitions */
	error(_("system argument is not valid"));
    X = cholmod_solve(sys, L, B, &c);
    UNPROTECT(1); Free(L); Free(B);
    return chm_dense_to_SEXP(X, 1, 0/*Rkind*/, GET_SLOT(bb, Matrix_DimNamesSym));
}

SEXP CHMfactor_spsolve(SEXP a, SEXP b, SEXP system)
{
    cholmod_factor *L = as_cholmod_factor(a);
    cholmod_sparse *B = as_cholmod_sparse(b), *X;
    int sys = asInteger(system);

    if (!(sys--))		/* -- to align with CHOLMOD definitions */
	error(_("system argument is not valid"));
    X = cholmod_spsolve(sys, L, B, &c);
    Free(L); Free(B);
    return chm_sparse_to_SEXP(X, 1, 0/*uploT*/, 0/*Rkind*/, "",
			      GET_SLOT(b, Matrix_DimNamesSym));
}

