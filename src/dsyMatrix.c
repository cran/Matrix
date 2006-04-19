#include "dsyMatrix.h"

SEXP symmetricMatrix_validate(SEXP obj)
{
    SEXP val = GET_SLOT(obj, Matrix_DimSym);
    if (LENGTH(val) < 2)
	return mkString(_("'Dim' slot has length less than two"));
    if (INTEGER(val)[0] != INTEGER(val)[1])
        return mkString(_("Matrix is not square"));
    if (isString(val = check_scalar_string(GET_SLOT(obj, Matrix_uploSym),
					   "LU", "uplo"))) return val;
    return ScalarLogical(1);
}

SEXP dsyMatrix_validate(SEXP obj)
{
    /* since "dsy" inherits from "symmetric", and "dMatrix", only need this:*/
    return dense_nonpacked_validate(obj);
}

double get_norm_sy(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    double *work = (double *) NULL;

    typnm[0] = norm_type(typstr);
    if (*typnm == 'I' || *typnm == 'O') {
        work = (double *) R_alloc(dims[0], sizeof(double));
    }
    return F77_CALL(dlansy)(typnm, uplo_P(obj),
			    dims, REAL(GET_SLOT(obj, Matrix_xSym)),
			    dims, work);
}

SEXP dsyMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm_sy(obj, CHAR(asChar(type))));
}


SEXP dsyMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP trf = dsyMatrix_trf(obj);
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), info;
    double anorm = get_norm_sy(obj, "O");
    double rcond;

    typnm[0] = rcond_type(CHAR(asChar(type)));
    F77_CALL(dsycon)(uplo_P(trf), dims,
		     REAL   (GET_SLOT(trf, Matrix_xSym)), dims,
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     &anorm, &rcond,
		     (double *) R_alloc(2*dims[0], sizeof(double)),
		     (int *) R_alloc(dims[0], sizeof(int)), &info);
    return ScalarReal(rcond);
}

SEXP dsyMatrix_solve(SEXP a)
{
    SEXP trf = dsyMatrix_trf(a);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dsyMatrix")));
    int *dims = INTEGER(GET_SLOT(trf, Matrix_DimSym)), info;

    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(trf, Matrix_uploSym)));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(trf, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(trf, Matrix_DimSym)));
    F77_CALL(dsytri)(uplo_P(val), dims,
		     REAL(GET_SLOT(val, Matrix_xSym)), dims,
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     (double *) R_alloc((long) dims[0], sizeof(double)),
		     &info);
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_dgeMatrix_solve(SEXP a, SEXP b)
{
    SEXP trf = dsyMatrix_trf(a),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	info;

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(b, Matrix_DimSym)));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(b, Matrix_xSym)));
    F77_CALL(dsytrs)(uplo_P(trf), adims, bdims + 1,
		     REAL(GET_SLOT(trf, Matrix_xSym)), adims,
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     REAL(GET_SLOT(val, Matrix_xSym)),
		     bdims, &info);
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP trf = dsyMatrix_trf(a),
	val = PROTECT(duplicate(b));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(getAttrib(b, R_DimSymbol)),
	info;

    if (!(isReal(b) && isMatrix(b)))
	error(_("Argument b must be a numeric matrix"));
    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    F77_CALL(dsytrs)(uplo_P(trf), adims, bdims + 1,
		     REAL(GET_SLOT(trf, Matrix_xSym)), adims,
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     REAL(val), bdims, &info);
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_as_dgeMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(from, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym,
	     duplicate(GET_SLOT(from, Matrix_DimSym)));
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    make_d_matrix_symmetric(REAL(GET_SLOT(val, Matrix_xSym)), from);
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_as_matrix(SEXP from)
{
    int n = INTEGER(GET_SLOT(from, Matrix_DimSym))[0];
    SEXP val = PROTECT(allocMatrix(REALSXP, n, n));

    make_d_matrix_symmetric(Memcpy(REAL(val),
				   REAL(GET_SLOT(from, Matrix_xSym)), n * n),
			    from);
    setAttrib(val, R_DimNamesSymbol, GET_SLOT(from, Matrix_DimNamesSym));
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_dgeMatrix_mm(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	*cdims,
	m = adims[0], n = bdims[1], k = adims[1];
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    double one = 1., zero = 0.;

    if (bdims[0] != k)
	error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1 || k < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    cdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    cdims[0] = m; cdims[1] = n;
    F77_CALL(dsymm)("L", uplo_P(a), adims, bdims+1, &one,
		    REAL(GET_SLOT(a, Matrix_xSym)), adims,
		    REAL(GET_SLOT(b, Matrix_xSym)), bdims,
		    &zero, REAL(GET_SLOT(val, Matrix_xSym)), adims);
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_dgeMatrix_mm_R(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	*cdims,
	m = adims[0], n = bdims[1], k = adims[1];
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    double one = 1., zero = 0.;

    if (bdims[0] != k)
	error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1 || k < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    cdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    cdims[0] = m; cdims[1] = n;
    F77_CALL(dsymm)("R", uplo_P(a), adims, bdims+1, &one,
		    REAL(GET_SLOT(a, Matrix_xSym)), adims,
		    REAL(GET_SLOT(b, Matrix_xSym)), bdims,
		    &zero, REAL(GET_SLOT(val, Matrix_xSym)), adims);
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_trf(SEXP x)
{
    SEXP val = get_factors(x, "BunchKaufman"),
	dimP = GET_SLOT(x, Matrix_DimSym),
	uploP = GET_SLOT(x, Matrix_uploSym);
    int *dims = INTEGER(dimP), *perm, info;
    int lwork = -1, n = dims[0];
    char *uplo = CHAR(STRING_ELT(uploP, 0));
    double tmp, *vx, *work;

    if (val != R_NilValue) return val;
    dims = INTEGER(dimP);
    val = PROTECT(NEW_OBJECT(MAKE_CLASS("BunchKaufman")));
    SET_SLOT(val, Matrix_uploSym, duplicate(uploP));
    SET_SLOT(val, Matrix_diagSym, mkString("N"));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n * n));
    AZERO(vx, n * n);
    F77_CALL(dlacpy)(uplo, &n, &n, REAL(GET_SLOT(x, Matrix_xSym)), &n, vx, &n);
    perm = INTEGER(ALLOC_SLOT(val, Matrix_permSym, INTSXP, n));
    F77_CALL(dsytrf)(uplo, &n, vx, &n, perm, &tmp, &lwork, &info);
    lwork = (int) tmp;
    work = Calloc(lwork, double);
    F77_CALL(dsytrf)(uplo, &n, vx, &n, perm, work, &lwork, &info);
    if (info) error(_("Lapack routine dsytrf returned error code %d"), info);
    UNPROTECT(1);
    Free(work);
    return set_factors(x, val, "BunchKaufman");
}

SEXP dsyMatrix_as_dspMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dspMatrix"))),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    full_to_packed_double(
	REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, (n*(n+1))/2)),
	REAL(GET_SLOT(from, Matrix_xSym)), n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW, NUN);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    UNPROTECT(1);
    return val;
}
