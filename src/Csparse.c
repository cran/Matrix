			/* Sparse matrices in compress column-oriented form */
#include "Csparse.h"
#ifdef USE_CHOLMOD
#include "chm_common.h"
#endif	/* USE_CHOLMOD */

SEXP Csparse_validate(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym);
    int j, ncol = length(pslot) - 1,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow, *xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    nrow = dims[0];
    if (length(pslot) <= 0)
	return mkString(_("slot p must have length > 0"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return mkString(_("last element of slot p must match length of slots i and x"));
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
    }
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    return ScalarLogical(1);
}

SEXP Csparse_to_Tsparse(SEXP x)
{
#ifdef USE_CHOLMOD
    cholmod_sparse *chxs = as_cholmod_sparse(x);
    cholmod_triplet *chxt = cholmod_sparse_to_triplet(chxs, &c);

    Free(chxs);
    return chm_triplet_to_SEXP(chxt, 1);
#else
    error("General conversion requires CHOLMOD");
    return R_NilValue;		/* -Wall */
#endif	/* USE_CHOLMOD */
}

SEXP Csparse_transpose(SEXP x)
{
#ifdef USE_CHOLMOD
    cholmod_sparse *chx = as_cholmod_sparse(x);
    cholmod_sparse *chxt = cholmod_transpose(chx, (int) chx->xtype, &c);

    Free(chx);
    return chm_sparse_to_SEXP(chxt, 1);
#else
    error("General conversion requires CHOLMOD");
    return R_NilValue;		/* -Wall */
#endif	/* USE_CHOLMOD */
}


SEXP Csparse_Csparse_prod(SEXP a, SEXP b)
{
#ifdef USE_CHOLMOD
    cholmod_sparse *cha = as_cholmod_sparse(a), *chb = as_cholmod_sparse(b);
    cholmod_sparse *chc = cholmod_ssmult(cha, chb, 0, (int) cha->xtype, 1, &c);

    Free(cha); Free(chb);
    return chm_sparse_to_SEXP(chc, 1);
#else
    error("General multiplication requires CHOLMOD");
    return R_NilValue;		/* -Wall */
#endif	/* USE_CHOLMOD */
}

SEXP Csparse_dense_prod(SEXP a, SEXP b)
{
#ifdef USE_CHOLMOD
    cholmod_sparse *cha = as_cholmod_sparse(a);
    cholmod_dense *chb = as_cholmod_dense(b);
    cholmod_dense *chc = cholmod_allocate_dense(cha->nrow, chb->ncol,
						cha->nrow, chb->xtype, &c);
    double alpha = 1, beta = 0;

    cholmod_sdmult(cha, 0, &alpha, &beta, chb, chc, &c);
    Free(cha); Free(chb);
    return chm_dense_to_SEXP(chc, 1);
#else
    error("General multiplication requires CHOLMOD");
    return R_NilValue;		/* -Wall */
#endif	/* USE_CHOLMOD */
}

SEXP Csparse_crossprod(SEXP x, SEXP trans, SEXP triplet)
{
#ifdef USE_CHOLMOD
    int trip = asLogical(triplet),
	tr   = asLogical(trans); /* gets reversed because _aat is tcrossprod */
    cholmod_triplet
	*cht = trip ? as_cholmod_triplet(x) : (cholmod_triplet*) NULL;
    cholmod_sparse *chcp, *chxt,
	*chx = trip ? cholmod_triplet_to_sparse(cht, cht->nnz, &c)
	: as_cholmod_sparse(x);

    if (!tr)
	chxt = cholmod_transpose(chx, (int) chx->xtype, &c);
    chcp = cholmod_aat((!tr) ? chxt : chx, (int *) NULL, 0, chx->xtype, &c);
    if(!chcp)
	error("Csparse_crossprod(): error return from cholmod_aat()");

    if (trip) {
	cholmod_free_sparse(&chx, &c);
	Free(cht);
    } else {
	Free(chx);
    }
    if (!tr) cholmod_free_sparse(&chxt, &c);
    return chm_sparse_to_SEXP(chcp, 1);
#else
    error("General crossproduct requires CHOLMOD");
    return R_NilValue;		/* -Wall */
#endif	/* USE_CHOLMOD */
}

