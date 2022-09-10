#include "dsyMatrix.h"

double get_norm_sy(SEXP obj, const char *typstr)
{
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    double *work = (double *) NULL;

    typnm[0] = La_norm_type(typstr);
    if (*typnm == 'I' || *typnm == 'O') {
        work = (double *) R_alloc(dims[0], sizeof(double));
    }
    return F77_CALL(dlansy)(typnm, uplo_P(obj),
			    dims, REAL(GET_SLOT(obj, Matrix_xSym)),
			    dims, work FCONE FCONE);
}

SEXP dsyMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm_sy(obj, CHAR(asChar(type))));
}

SEXP dsyMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP trf = dsyMatrix_trf(obj);
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), info;
    double anorm = get_norm_sy(obj, "O");
    double rcond;

    F77_CALL(dsycon)(uplo_P(trf), dims,
		     REAL   (GET_SLOT(trf, Matrix_xSym)), dims,
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     &anorm, &rcond,
		     (double *) R_alloc(2*dims[0], sizeof(double)),
		     (int *) R_alloc(dims[0], sizeof(int)), &info FCONE);
    return ScalarReal(rcond);
}

SEXP dsyMatrix_solve(SEXP a)
{
    SEXP trf = dsyMatrix_trf(a);
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dsyMatrix"));
    int *dims = INTEGER(GET_SLOT(trf, Matrix_DimSym)), info;

    slot_dup(val, trf, Matrix_uploSym);
    slot_dup(val, trf, Matrix_xSym);
    slot_dup(val, trf, Matrix_DimSym);
    F77_CALL(dsytri)(uplo_P(val), dims,
		     REAL(GET_SLOT(val, Matrix_xSym)), dims,
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     (double *) R_alloc((long) dims[0], sizeof(double)),
		     &info FCONE);
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP trf = dsyMatrix_trf(a),
	val = PROTECT(dense_as_general(b, 'd', 2, 0));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym)),
	info;

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    F77_CALL(dsytrs)(uplo_P(trf), adims, bdims + 1,
		     REAL(GET_SLOT(trf, Matrix_xSym)), adims,
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     REAL(GET_SLOT(val, Matrix_xSym)),
		     bdims, &info FCONE);
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP rtP)
{
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0));// incl. dimnames
    int rt = asLogical(rtP); /* if(rt), compute b %*% a,  else  a %*% b */
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym)),
	m = bdims[0], n = bdims[1];

    if ((rt && n != adims[0]) || (!rt && m != adims[0]))
	error(_("Matrices are not conformable for multiplication"));

    double one = 1., zero = 0.;
    R_xlen_t mn = m * (R_xlen_t)n;
    double *bcp, *vx = REAL(GET_SLOT(val, Matrix_xSym));
    Calloc_or_Alloca_TO(bcp, mn, double);
    Memcpy(bcp, vx, mn);

    if (m >=1 && n >= 1)
	F77_CALL(dsymm)(rt ? "R" :"L", uplo_P(a), &m, &n, &one,
			REAL(GET_SLOT(a, Matrix_xSym)), adims, bcp,
			&m, &zero, vx, &m FCONE FCONE);
    // add dimnames:
    int nd = rt ?
	1 : // v <- b %*% a : rownames(v) == rownames(b)  are already there
	0;  // v <- a %*% b : colnames(v) == colnames(b)  are already there
    SEXP nms = PROTECT(VECTOR_ELT(get_symmetrized_DimNames(a, -1), nd));
    SET_VECTOR_ELT(GET_SLOT(val, Matrix_DimNamesSym), nd, nms);
    Free_FROM(bcp, mn);
    UNPROTECT(2);
    return val;
}

SEXP dsyMatrix_trf(SEXP x)
{
    SEXP val = get_factor(x, "BunchKaufman");
    if (val != R_NilValue) return val;

    SEXP dimP = GET_SLOT(x, Matrix_DimSym),
	uploP = GET_SLOT(x, Matrix_uploSym);
    int n = INTEGER(dimP)[0];
    R_xlen_t nsqr = n; nsqr *= n; // nsqr = n^2 (w/o overflow !)
    const char *uplo = CHAR(STRING_ELT(uploP, 0));

    val = PROTECT(NEW_OBJECT_OF_CLASS("BunchKaufman"));
    SET_SLOT(val, Matrix_uploSym, duplicate(uploP));
    SET_SLOT(val, Matrix_diagSym, mkString("N"));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    double *vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, nsqr));
    AZERO(vx, nsqr, 0.0);
    F77_CALL(dlacpy)(uplo, &n, &n, REAL(GET_SLOT(x, Matrix_xSym)), &n, vx, &n FCONE);
    int *perm = INTEGER(ALLOC_SLOT(val, Matrix_permSym, INTSXP, n)),
	info, lwork = -1;
    double tmp, *work;
    F77_CALL(dsytrf)(uplo, &n, vx, &n, perm, &tmp, &lwork, &info FCONE);
    lwork = (int) tmp;
    Calloc_or_Alloca_TO(work, lwork, double);

    F77_CALL(dsytrf)(uplo, &n, vx, &n, perm, work, &lwork, &info FCONE);

    Free_FROM(work, lwork);
    if (info) error(_("Lapack routine dsytrf returned error code %d"), info);
    set_factor(x, "BunchKaufman", val);
    UNPROTECT(1);
    return val;
}

/** BunchKaufmann(<simple matrix>)
 */
SEXP matrix_trf(SEXP x, SEXP uploP)
{
    if (!(isReal(x) & isMatrix(x)))
	error(_("x must be a \"double\" (numeric) matrix"));
    SEXP dimP = getAttrib(x, R_DimSymbol);
    if(TYPEOF(dimP) == INTSXP)
	dimP = duplicate(dimP);
    else
        dimP = coerceVector(dimP, INTSXP);
    PROTECT(dimP);
    int *dims = INTEGER(dimP),
	n = dims[0];
    R_xlen_t nsqr = n; nsqr *= n; // nsqr = n^2 (w/o overflow !)

    if(n != dims[1])
	error(_("matrix_trf(x, *): matrix is not square"));
    /* In principle, we "should" check that the matrix is symmetric,
       OTOH, we only use its lower or upper (depending on 'uploP') triangular part */
    if(uploP == R_NilValue) {
	uploP = mkString("U"); // Default: if not specified, use "U"
    } else {
	if(TYPEOF(uploP) != STRSXP)
	    error(_("matrix_trf(*, uplo): uplo must be string"));
	uploP = duplicate(uploP); // as we "add" it to the result
    }
    PROTECT(uploP);
    const char *uplo = CHAR(STRING_ELT(uploP, 0));
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("BunchKaufman"));
    SET_SLOT(val, Matrix_uploSym, uploP);
    SET_SLOT(val, Matrix_diagSym, mkString("N"));
    SET_SLOT(val, Matrix_DimSym, dimP);
    double *vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, nsqr)); // n x n result matrix
    AZERO(vx, nsqr, 0.0);
    F77_CALL(dlacpy)(uplo, &n, &n, REAL(x), &n, vx, &n FCONE);
    int *perm = INTEGER(ALLOC_SLOT(val, Matrix_permSym, INTSXP, n)),
         info, lwork = -1;
    double tmp, *work;
    F77_CALL(dsytrf)(uplo, &n, vx, &n, perm, &tmp, &lwork, &info FCONE);
    lwork = (int) tmp;
    Calloc_or_Alloca_TO(work, lwork, double);
    F77_CALL(dsytrf)(uplo, &n, vx, &n, perm, work, &lwork, &info FCONE);
    Free_FROM(work, lwork);
    if (info) error(_("Lapack routine dsytrf returned error code %d"), info);
    UNPROTECT(3);
    return val;
}

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0

// this is very close to lsyMatrix_as_lsp*() in ./ldense.c  -- keep synced !
SEXP dsyMatrix_as_dspMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dspMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_pack(
	REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, (n*(n+1))/2)),
	REAL( GET_SLOT(from, Matrix_xSym)),
	n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	NUN);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    SET_SLOT(val, Matrix_factorSym,
	     duplicate(GET_SLOT(from, Matrix_factorSym)));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0

SEXP dsyMatrix_as_matrix(SEXP from, SEXP keep_dimnames)
{
    int n = INTEGER(GET_SLOT(from, Matrix_DimSym))[0];
    SEXP val = PROTECT(allocMatrix(REALSXP, n, n));
    R_xlen_t nsqr = n; nsqr *= n;

    ddense_unpacked_make_symmetric(Memcpy(REAL(val),
					  REAL(GET_SLOT(from, Matrix_xSym)),
					  nsqr),
				   from);
    if(asLogical(keep_dimnames))
	setAttrib(val, R_DimNamesSymbol, get_symmetrized_DimNames(from, -1));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */


