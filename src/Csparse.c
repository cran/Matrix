			/* Sparse matrices in compressed column-oriented form */
#include "Csparse.h"
#include "chm_common.h"

SEXP Csparse_validate(SEXP x)
{
    /* NB: we do *NOT* check a potential 'x' slot here, at all */
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym);
    Rboolean sorted, strictly;
    int j, k,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow = dims[0],
	ncol = dims[1],
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    if (length(pslot) != dims[1] + 1)
	return mkString(_("slot p must have length = ncol(.) + 1"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) < xp[ncol]) /* allow larger slots from over-allocation!*/
	return
	    mkString(_("last element of slot p must match length of slots i and x"));
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    sorted = TRUE; strictly = TRUE;
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
	if(sorted)
	    for (k = xp[j] + 1; k < xp[j + 1]; k++) {
		if (xi[k] < xi[k - 1])
		    sorted = FALSE;
		else if (xi[k] == xi[k - 1])
		    strictly = FALSE;
	    }
    }
    if (!sorted) {
	CHM_SP chx = AS_CHM_SP(x);
	R_CheckStack();

	cholmod_sort(chx, &c);
	/* Now re-check that row indices are *strictly* increasing
	 * (and not just increasing) within each column : */
	for (j = 0; j < ncol; j++) {
	    for (k = xp[j] + 1; k < xp[j + 1]; k++)
		if (xi[k] == xi[k - 1])
		    return mkString(_("slot i is not *strictly* increasing inside a column (even after cholmod_sort)"));
	}

    } else if(!strictly) {  /* sorted, but not strictly */
	return mkString(_("slot i is not *strictly* increasing inside a column"));
    }
    return ScalarLogical(1);
}

/* Called from ../R/Csparse.R : */
/* Can only return [dln]geMatrix (no symm/triang);
 * FIXME: replace by non-CHOLMOD code ! */
SEXP Csparse_to_dense(SEXP x)
{
    CHM_SP chxs = AS_CHM_SP(x);
    /* This loses the symmetry property, since cholmod_dense has none,
     * BUT, much worse (FIXME!), it also transforms CHOLMOD_PATTERN ("n") matrices
     * to numeric (CHOLMOD_REAL) ones : */
    CHM_DN chxd = cholmod_sparse_to_dense(chxs, &c);
    int Rkind = (chxs->xtype == CHOLMOD_PATTERN)? -1 : Real_kind(x);
    R_CheckStack();

    return chm_dense_to_SEXP(chxd, 1, Rkind, GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_to_nz_pattern(SEXP x, SEXP tri)
{
    CHM_SP chxs = AS_CHM_SP(x);
    CHM_SP chxcp = cholmod_copy(chxs, chxs->stype, CHOLMOD_PATTERN, &c);
    int tr = asLogical(tri);
    R_CheckStack();

    return chm_sparse_to_SEXP(chxcp, 1/*do_free*/,
			      tr ? ((*uplo_P(x) == 'U') ? 1 : -1) : 0,
			      0, tr ? diag_P(x) : "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_to_matrix(SEXP x)
{
    return chm_dense_to_matrix(cholmod_sparse_to_dense(AS_CHM_SP(x), &c),
			       1 /*do_free*/, GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_to_Tsparse(SEXP x, SEXP tri)
{
    CHM_SP chxs = AS_CHM_SP(x);
    CHM_TR chxt = cholmod_sparse_to_triplet(chxs, &c);
    int tr = asLogical(tri);
    int Rkind = (chxs->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    return chm_triplet_to_SEXP(chxt, 1,
			       tr ? ((*uplo_P(x) == 'U') ? 1 : -1) : 0,
			       Rkind, tr ? diag_P(x) : "",
			       GET_SLOT(x, Matrix_DimNamesSym));
}

/* this used to be called  sCMatrix_to_gCMatrix(..)   [in ./dsCMatrix.c ]: */
SEXP Csparse_symmetric_to_general(SEXP x)
{
    CHM_SP chx = AS_CHM_SP(x), chgx;
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    if (!(chx->stype))
	error(_("Nonsymmetric matrix in Csparse_symmetric_to_general"));
    chgx = cholmod_copy(chx, /* stype: */ 0, chx->xtype, &c);
    /* xtype: pattern, "real", complex or .. */
    return chm_sparse_to_SEXP(chgx, 1, 0, Rkind, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_general_to_symmetric(SEXP x, SEXP uplo)
{
    CHM_SP chx = AS_CHM_SP(x), chgx;
    int uploT = (*CHAR(asChar(uplo)) == 'U') ? 1 : -1;
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    chgx = cholmod_copy(chx, /* stype: */ uploT, chx->xtype, &c);
    /* xtype: pattern, "real", complex or .. */
    return chm_sparse_to_SEXP(chgx, 1, 0, Rkind, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_transpose(SEXP x, SEXP tri)
{
    /* TODO: lgCMatrix & igC* currently go via double prec. cholmod -
     *       since cholmod (& cs) lacks sparse 'int' matrices */
    CHM_SP chx = AS_CHM_SP(x);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    CHM_SP chxt = cholmod_transpose(chx, chx->xtype, &c);
    SEXP dn = PROTECT(duplicate(GET_SLOT(x, Matrix_DimNamesSym))), tmp;
    int tr = asLogical(tri);
    R_CheckStack();

    tmp = VECTOR_ELT(dn, 0);	/* swap the dimnames */
    SET_VECTOR_ELT(dn, 0, VECTOR_ELT(dn, 1));
    SET_VECTOR_ELT(dn, 1, tmp);
    UNPROTECT(1);
    return chm_sparse_to_SEXP(chxt, 1, /* SWAP 'uplo' for triangular */
			      tr ? ((*uplo_P(x) == 'U') ? -1 : 1) : 0,
			      Rkind, tr ? diag_P(x) : "", dn);
}

SEXP Csparse_Csparse_prod(SEXP a, SEXP b)
{
    CHM_SP cha = AS_CHM_SP(a), chb = AS_CHM_SP(b);
    CHM_SP chc = cholmod_ssmult(cha, chb, 0, cha->xtype, 1, &c);
    SEXP dn = allocVector(VECSXP, 2);
    R_CheckStack();

    SET_VECTOR_ELT(dn, 0,	/* establish dimnames */
		   duplicate(VECTOR_ELT(GET_SLOT(a, Matrix_DimNamesSym), 0)));
    SET_VECTOR_ELT(dn, 1,
		   duplicate(VECTOR_ELT(GET_SLOT(b, Matrix_DimNamesSym), 1)));
    return chm_sparse_to_SEXP(chc, 1, 0, 0, "", dn);
}

SEXP Csparse_Csparse_crossprod(SEXP a, SEXP b, SEXP trans)
{
    int tr = asLogical(trans);
    CHM_SP cha = AS_CHM_SP(a), chb = AS_CHM_SP(b), chTr, chc;
    SEXP dn = allocVector(VECSXP, 2);
    R_CheckStack();

    chTr = cholmod_transpose((tr) ? chb : cha, chb->xtype, &c);
    chc = cholmod_ssmult((tr) ? cha : chTr, (tr) ? chTr : chb,
			 0, cha->xtype, 1, &c);
    cholmod_free_sparse(&chTr, &c);

    SET_VECTOR_ELT(dn, 0,	/* establish dimnames */
		   duplicate(VECTOR_ELT(GET_SLOT(a, Matrix_DimNamesSym), (tr) ? 0 : 1)));
    SET_VECTOR_ELT(dn, 1,
		   duplicate(VECTOR_ELT(GET_SLOT(b, Matrix_DimNamesSym), (tr) ? 0 : 1)));
    return chm_sparse_to_SEXP(chc, 1, 0, 0, "", dn);
}

SEXP Csparse_dense_prod(SEXP a, SEXP b)
{
    CHM_SP cha = AS_CHM_SP(a);
    SEXP b_M = PROTECT(mMatrix_as_dgeMatrix(b));
    CHM_DN chb = AS_CHM_DN(b_M);
    CHM_DN chc = cholmod_allocate_dense(cha->nrow, chb->ncol, cha->nrow,
					chb->xtype, &c);
    SEXP dn = PROTECT(allocVector(VECSXP, 2));
    double one[] = {1,0}, zero[] = {0,0};
    R_CheckStack();

    cholmod_sdmult(cha, 0, one, zero, chb, chc, &c);
    SET_VECTOR_ELT(dn, 0,	/* establish dimnames */
		   duplicate(VECTOR_ELT(GET_SLOT(a, Matrix_DimNamesSym), 0)));
    SET_VECTOR_ELT(dn, 1,
		   duplicate(VECTOR_ELT(GET_SLOT(b_M, Matrix_DimNamesSym), 1)));
    UNPROTECT(2);
    return chm_dense_to_SEXP(chc, 1, 0, dn);
}

SEXP Csparse_dense_crossprod(SEXP a, SEXP b)
{
    CHM_SP cha = AS_CHM_SP(a);
    SEXP b_M = PROTECT(mMatrix_as_dgeMatrix(b));
    CHM_DN chb = AS_CHM_DN(b_M);
    CHM_DN chc = cholmod_allocate_dense(cha->ncol, chb->ncol, cha->ncol,
					chb->xtype, &c);
    SEXP dn = PROTECT(allocVector(VECSXP, 2));
    double one[] = {1,0}, zero[] = {0,0};
    R_CheckStack();

    cholmod_sdmult(cha, 1, one, zero, chb, chc, &c);
    SET_VECTOR_ELT(dn, 0,	/* establish dimnames */
		   duplicate(VECTOR_ELT(GET_SLOT(a, Matrix_DimNamesSym), 1)));
    SET_VECTOR_ELT(dn, 1,
		   duplicate(VECTOR_ELT(GET_SLOT(b_M, Matrix_DimNamesSym), 1)));
    UNPROTECT(2);
    return chm_dense_to_SEXP(chc, 1, 0, dn);
}

/* Computes   x'x  or  x x'  -- see Csparse_Csparse_crossprod above for  x'y and x y' */
SEXP Csparse_crossprod(SEXP x, SEXP trans, SEXP triplet)
{
    int trip = asLogical(triplet),
	tr   = asLogical(trans); /* gets reversed because _aat is tcrossprod */
    CHM_TR cht = trip ? AS_CHM_TR(x) : (CHM_TR) NULL;
    CHM_SP chcp, chxt,
	chx = trip ? cholmod_triplet_to_sparse(cht, cht->nnz, &c) : AS_CHM_SP(x);
    SEXP dn = PROTECT(allocVector(VECSXP, 2));
    R_CheckStack();

    if (!tr) chxt = cholmod_transpose(chx, chx->xtype, &c);
    chcp = cholmod_aat((!tr) ? chxt : chx, (int *) NULL, 0, chx->xtype, &c);
    if(!chcp) error(_("Csparse_crossprod(): error return from cholmod_aat()"));
    cholmod_band_inplace(0, chcp->ncol, chcp->xtype, chcp, &c);
    chcp->stype = 1;
    if (trip) cholmod_free_sparse(&chx, &c);
    if (!tr) cholmod_free_sparse(&chxt, &c);
    SET_VECTOR_ELT(dn, 0,	/* establish dimnames */
		   duplicate(VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
					(tr) ? 0 : 1)));
    SET_VECTOR_ELT(dn, 1, duplicate(VECTOR_ELT(dn, 0)));
    UNPROTECT(1);
    return chm_sparse_to_SEXP(chcp, 1, 0, 0, "", dn);
}

SEXP Csparse_drop(SEXP x, SEXP tol)
{
    CHM_SP chx = AS_CHM_SP(x);
    CHM_SP ans = cholmod_copy(chx, chx->stype, chx->xtype, &c);
    double dtol = asReal(tol);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    if(!cholmod_drop(dtol, ans, &c))
	error(_("cholmod_drop() failed"));
    return chm_sparse_to_SEXP(ans, 1, 0, Rkind, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_horzcat(SEXP x, SEXP y)
{
    CHM_SP chx = AS_CHM_SP(x), chy = AS_CHM_SP(y);
    int Rkind = 0; /* only for "d" - FIXME */
    R_CheckStack();

    /* FIXME: currently drops dimnames */
    return chm_sparse_to_SEXP(cholmod_horzcat(chx, chy, 1, &c),
			      1, 0, Rkind, "", R_NilValue);
}

SEXP Csparse_vertcat(SEXP x, SEXP y)
{
    CHM_SP chx = AS_CHM_SP(x), chy = AS_CHM_SP(y);
    int Rkind = 0; /* only for "d" - FIXME */
    R_CheckStack();

    /* FIXME: currently drops dimnames */
    return chm_sparse_to_SEXP(cholmod_vertcat(chx, chy, 1, &c),
			      1, 0, Rkind, "", R_NilValue);
}

SEXP Csparse_band(SEXP x, SEXP k1, SEXP k2)
{
    CHM_SP chx = AS_CHM_SP(x);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    CHM_SP ans = cholmod_band(chx, asInteger(k1), asInteger(k2), chx->xtype, &c);
    R_CheckStack();

    return chm_sparse_to_SEXP(ans, 1, 0, Rkind, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_diagU2N(SEXP x)
{
    if (*diag_P(x) != 'U') {/* "trivially fast" when there's no 'diag' slot at all */
	return (x);
    }
    else {
	CHM_SP chx = AS_CHM_SP(x);
	CHM_SP eye = cholmod_speye(chx->nrow, chx->ncol, chx->xtype, &c);
	double one[] = {1, 0};
	CHM_SP ans = cholmod_add(chx, eye, one, one, TRUE, TRUE, &c);
	int uploT = (*uplo_P(x) == 'U') ? 1 : -1;
	int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;

	R_CheckStack();
	cholmod_free_sparse(&eye, &c);
	return chm_sparse_to_SEXP(ans, 1, uploT, Rkind, "N",
				  GET_SLOT(x, Matrix_DimNamesSym));
    }
}

SEXP Csparse_submatrix(SEXP x, SEXP i, SEXP j)
{
    CHM_SP chx = AS_CHM_SP(x);
    int rsize = (isNull(i)) ? -1 : LENGTH(i),
	csize = (isNull(j)) ? -1 : LENGTH(j);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    if (rsize >= 0 && !isInteger(i))
	error(_("Index i must be NULL or integer"));
    if (csize >= 0 && !isInteger(j))
	error(_("Index j must be NULL or integer"));

    return chm_sparse_to_SEXP(cholmod_submatrix(chx, INTEGER(i), rsize,
						INTEGER(j), csize,
						TRUE, TRUE, &c),
			      1, 0, Rkind, "",
			      /* FIXME: drops dimnames */ R_NilValue);
}
