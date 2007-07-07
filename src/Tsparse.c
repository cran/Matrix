				/* Sparse matrices in triplet form */
#include "Tsparse.h"
#include "chm_common.h"

SEXP Tsparse_validate(SEXP x)
{
    /* NB: we do *NOT* check a potential 'x' slot here, at all */
    SEXP
	islot = GET_SLOT(x, Matrix_iSym),
	jslot = GET_SLOT(x, Matrix_jSym),
	dimslot = GET_SLOT(x, Matrix_DimSym);
    int j,
	nrow = INTEGER(dimslot)[0],
	ncol = INTEGER(dimslot)[1],
	nnz = length(islot),
	*xj = INTEGER(jslot),
	*xi = INTEGER(islot);

    if (length(jslot) != nnz)
	return mkString(_("lengths of slots i and j must match"));
    /* FIXME: this is checked in super class -- no need to do here: */
    if (length(dimslot) != 2)
	return mkString(_("slot Dim must have length 2"));

    for (j = 0; j < nnz; j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
	if (xj[j] < 0 || xj[j] >= ncol)
	    return mkString(_("all column indices must be between 0 and ncol-1"));
    }
    return ScalarLogical(1);
}

SEXP Tsparse_to_Csparse(SEXP x, SEXP tri)
{
    CHM_TR chxt = AS_CHM_TR(x);
    CHM_SP chxs = cholmod_triplet_to_sparse(chxt, chxt->nnz, &c);
    int tr = asLogical(tri);
    int Rkind = (chxt->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    return chm_sparse_to_SEXP(chxs, 1,
			      tr ? ((*uplo_P(x) == 'U') ? 1 : -1) : 0,
			      Rkind, tr ? diag_P(x) : "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}
