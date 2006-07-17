				/* CHOLMOD factors */
#include "CHMfactor.h"

SEXP CHMfactor_to_sparse(SEXP x)
{
    cholmod_factor *L = as_cholmod_factor(x), *Lcp;
    cholmod_sparse *Lm;

    Lcp = cholmod_copy_factor(L, &c); Free(L); /* next call changes Lcp */
    Lm = cholmod_factor_to_sparse(Lcp, &c); cholmod_free_factor(&Lcp, &c);
    return chm_sparse_to_SEXP(Lm, -1);
}

