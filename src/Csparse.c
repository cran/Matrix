			/* Sparse matrices in compressed column-oriented form */
#include "Csparse.h"
#include "chm_common.h"

SEXP Csparse_validate(SEXP x)
{
    /* NB: we do *NOT* check a potential 'x' slot here, at all */
    cholmod_sparse *chx = as_cholmod_sparse(x);
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym);
    int j, k, ncol = length(pslot) - 1,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow, sorted, *xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    nrow = dims[0];
    if (length(pslot) <= 0)
	return mkString(_("slot p must have length > 0"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return
	    mkString(_("last element of slot p must match length of slots i and x"));
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    sorted = TRUE;
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
	for (k = xp[j] + 1; k < xp[j + 1]; k++)
	    if (xi[k] < xi[k - 1]) sorted = FALSE;
    }
    if (!sorted) cholmod_sort(chx, &c);
    Free(chx);
    return ScalarLogical(1);
}

SEXP Csparse_to_dense(SEXP x)
{
    cholmod_sparse *chxs = as_cholmod_sparse(x);
    cholmod_dense *chxd = cholmod_sparse_to_dense(chxs, &c);

    Free(chxs);
    return chm_dense_to_SEXP(chxd, 1, Real_kind(x));
}

SEXP Csparse_to_nz_pattern(SEXP x, SEXP tri)
{
    cholmod_sparse *chxs = as_cholmod_sparse(x);
    cholmod_sparse
	*chxcp = cholmod_copy(chxs, chxs->stype, CHOLMOD_PATTERN, &c);
    int uploT = 0; char *diag = "";

    Free(chxs);
    if (asLogical(tri)) {	/* triangular sparse matrices */
	uploT = (strcmp(CHAR(asChar(GET_SLOT(x, Matrix_uploSym))), "U")) ?
	    -1 : 1;
	diag = CHAR(asChar(GET_SLOT(x, Matrix_diagSym)));
    }
    return chm_sparse_to_SEXP(chxcp, 1, uploT, 0, diag,
			      GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_to_matrix(SEXP x)
{
    cholmod_sparse *chxs = as_cholmod_sparse(x);
    cholmod_dense *chxd = cholmod_sparse_to_dense(chxs, &c);

    Free(chxs);
    return chm_dense_to_matrix(chxd, 1,
			       GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_to_Tsparse(SEXP x, SEXP tri)
{
    cholmod_sparse *chxs = as_cholmod_sparse(x);
    cholmod_triplet *chxt = cholmod_sparse_to_triplet(chxs, &c);
    int uploT = 0;
    char *diag = "";
    int Rkind = (chxs->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;

    Free(chxs);
    if (asLogical(tri)) {	/* triangular sparse matrices */
	uploT = (*uplo_P(x) == 'U') ? -1 : 1;
	diag = diag_P(x);
    }
    return chm_triplet_to_SEXP(chxt, 1, uploT, Rkind, diag,
			       GET_SLOT(x, Matrix_DimNamesSym));
}

/* this used to be called  sCMatrix_to_gCMatrix(..)   [in ./dsCMatrix.c ]: */
SEXP Csparse_symmetric_to_general(SEXP x)
{
    cholmod_sparse *chx = as_cholmod_sparse(x), *chgx;
    int Rkind = (chx->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;

    if (!(chx->stype))
	error(_("Nonsymmetric matrix in Csparse_symmetric_to_general"));
    chgx = cholmod_copy(chx, /* stype: */ 0, chx->xtype, &c);
    /* xtype: pattern, "real", complex or .. */
    Free(chx);
    return chm_sparse_to_SEXP(chgx, 1, 0, Rkind, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

#ifdef _not_yet_FIXME_
/* MM: This would seem useful; e.g. lsC* can hardly be coerced to ! */
SEXP Csparse_general_to_symmetric(SEXP x,
				  int stype)/*-1 : "L", +1 : "U" */
{
    cholmod_sparse *chx = as_cholmod_sparse(x), *chgx;
    int Rkind = (chx->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;

    chgx = cholmod_copy(chx, /* stype: */ stype, chx->xtype, &c);
    /* xtype: pattern, "real", complex or .. */
    Free(chx);
    return chm_sparse_to_SEXP(chgx, 1, 0, Rkind, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

#endif

SEXP Csparse_transpose(SEXP x, SEXP tri)
{
    cholmod_sparse *chx = as_cholmod_sparse(x);
    int Rkind = (chx->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;
    cholmod_sparse *chxt = cholmod_transpose(chx, (int) chx->xtype, &c);
    SEXP dn = PROTECT(duplicate(GET_SLOT(x, Matrix_DimNamesSym))), tmp;
    int uploT = 0; char *diag = "";

    Free(chx);
    tmp = VECTOR_ELT(dn, 0);	/* swap the dimnames */
    SET_VECTOR_ELT(dn, 0, VECTOR_ELT(dn, 1));
    SET_VECTOR_ELT(dn, 1, tmp);
    UNPROTECT(1);
    if (asLogical(tri)) {	/* triangular sparse matrices */
	uploT = (*uplo_P(x) == 'U') ? -1 : 1;
	diag = diag_P(x);
    }
    return chm_sparse_to_SEXP(chxt, 1, uploT, Rkind, diag, dn);
}

SEXP Csparse_Csparse_prod(SEXP a, SEXP b)
{
    cholmod_sparse *cha = as_cholmod_sparse(a),
	*chb = as_cholmod_sparse(b);
    cholmod_sparse *chc = cholmod_ssmult(cha, chb, 0, cha->xtype, 1, &c);
    SEXP dn = allocVector(VECSXP, 2);

    Free(cha); Free(chb);
    SET_VECTOR_ELT(dn, 0,	/* establish dimnames */
		   duplicate(VECTOR_ELT(GET_SLOT(a, Matrix_DimNamesSym), 0)));
    SET_VECTOR_ELT(dn, 1,
		   duplicate(VECTOR_ELT(GET_SLOT(b, Matrix_DimNamesSym), 1)));
    return chm_sparse_to_SEXP(chc, 1, 0, 0, "", dn);
}

SEXP Csparse_dense_prod(SEXP a, SEXP b)
{
    cholmod_sparse *cha = as_cholmod_sparse(a);
    cholmod_dense *chb = as_cholmod_dense(PROTECT(mMatrix_as_dgeMatrix(b)));
    cholmod_dense *chc =
	cholmod_allocate_dense(cha->nrow, chb->ncol, cha->nrow, chb->xtype, &c);
    double alpha[] = {1,0}, beta[] = {0,0};

    cholmod_sdmult(cha, 0, alpha, beta, chb, chc, &c);
    Free(cha); Free(chb);
    UNPROTECT(1);
    return chm_dense_to_SEXP(chc, 1, 0);
}

SEXP Csparse_dense_crossprod(SEXP a, SEXP b)
{
    cholmod_sparse *cha = as_cholmod_sparse(a);
    cholmod_dense *chb = as_cholmod_dense(PROTECT(mMatrix_as_dgeMatrix(b)));
    cholmod_dense *chc =
	cholmod_allocate_dense(cha->ncol, chb->ncol, cha->ncol, chb->xtype, &c);
    double alpha[] = {1,0}, beta[] = {0,0};

    cholmod_sdmult(cha, 1, alpha, beta, chb, chc, &c);
    Free(cha); Free(chb);
    UNPROTECT(1);
    return chm_dense_to_SEXP(chc, 1, 0);
}

SEXP Csparse_crossprod(SEXP x, SEXP trans, SEXP triplet)
{
    int trip = asLogical(triplet),
	tr   = asLogical(trans); /* gets reversed because _aat is tcrossprod */
    cholmod_triplet
	*cht = trip ? as_cholmod_triplet(x) : (cholmod_triplet*) NULL;
    cholmod_sparse *chcp, *chxt,
	*chx = trip ? cholmod_triplet_to_sparse(cht, cht->nnz, &c)
	: as_cholmod_sparse(x);
    SEXP dn = PROTECT(allocVector(VECSXP, 2));

    if (!tr)
	chxt = cholmod_transpose(chx, chx->xtype, &c);
    chcp = cholmod_aat((!tr) ? chxt : chx, (int *) NULL, 0, chx->xtype, &c);
    if(!chcp)
	error("Csparse_crossprod(): error return from cholmod_aat()");
    cholmod_band_inplace(0, chcp->ncol, chcp->xtype, chcp, &c);
    chcp->stype = 1;
    if (trip) {
	cholmod_free_sparse(&chx, &c);
	Free(cht);
    } else {
	Free(chx);
    }
    if (!tr) cholmod_free_sparse(&chxt, &c);
				/* create dimnames */
    SET_VECTOR_ELT(dn, 0,
		   duplicate(VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
					(tr) ? 1 : 0)));
    SET_VECTOR_ELT(dn, 1, duplicate(VECTOR_ELT(dn, 0)));
    UNPROTECT(1);
    return chm_sparse_to_SEXP(chcp, 1, 0, 0, "", dn);
}

SEXP Csparse_horzcat(SEXP x, SEXP y)
{
    cholmod_sparse *chx = as_cholmod_sparse(x),
	*chy = as_cholmod_sparse(y), *ans;
    int Rkind = 0; /* only for "d" - FIXME */

    ans = cholmod_horzcat(chx, chy, 1, &c);
    Free(chx); Free(chy);
    /* FIXME: currently drops dimnames */
    return chm_sparse_to_SEXP(ans, 1, 0, Rkind, "", R_NilValue);
}

SEXP Csparse_vertcat(SEXP x, SEXP y)
{
    cholmod_sparse *chx = as_cholmod_sparse(x),
	*chy = as_cholmod_sparse(y), *ans;
    int Rkind = 0; /* only for "d" - FIXME */

    ans = cholmod_vertcat(chx, chy, 1, &c);
    Free(chx); Free(chy);
    /* FIXME: currently drops dimnames */
    return chm_sparse_to_SEXP(ans, 1, 0, Rkind, "", R_NilValue);
}

SEXP Csparse_band(SEXP x, SEXP k1, SEXP k2)
{
    cholmod_sparse *chx = as_cholmod_sparse(x), *ans;
    int Rkind = (chx->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;

    ans = cholmod_band(chx, asInteger(k1), asInteger(k2), chx->xtype, &c);
    Free(chx);
    return chm_sparse_to_SEXP(ans, 1, 0, Rkind, "", R_NilValue);
}

SEXP Csparse_diagU2N(SEXP x)
{
    cholmod_sparse *chx = as_cholmod_sparse(x);
    cholmod_sparse *eye = cholmod_speye(chx->nrow, chx->ncol, chx->xtype, &c);
    double one[] = {1, 0};
    cholmod_sparse *ans = cholmod_add(chx, eye, one, one, TRUE, TRUE, &c);
    int uploT = (strcmp(CHAR(asChar(GET_SLOT(x, Matrix_uploSym))), "U")) ?
	-1 : 1;
    int Rkind = (chx->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;

    Free(chx); cholmod_free_sparse(&eye, &c);
    return chm_sparse_to_SEXP(ans, 1, uploT, Rkind, "N",
			      duplicate(GET_SLOT(x, Matrix_DimNamesSym)));
}

SEXP Csparse_submatrix(SEXP x, SEXP i, SEXP j)
{
    cholmod_sparse *chx = as_cholmod_sparse(x);
    int rsize = (isNull(i)) ? -1 : LENGTH(i),
	csize = (isNull(j)) ? -1 : LENGTH(j);
    int Rkind = (chx->xtype == CHOLMOD_REAL) ? Real_kind(x) : 0;

    if (rsize >= 0 && !isInteger(i))
	error(_("Index i must be NULL or integer"));
    if (csize >= 0 && !isInteger(j))
	error(_("Index j must be NULL or integer"));
    return chm_sparse_to_SEXP(cholmod_submatrix(chx, INTEGER(i), rsize,
						INTEGER(j), csize,
						TRUE, TRUE, &c),
			      1, 0, Rkind, "", R_NilValue);
}
