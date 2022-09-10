#include "dpoMatrix.h"

SEXP dpoMatrix_validate(SEXP obj)
{
    double *x = REAL(GET_SLOT(obj, Matrix_xSym));
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    R_xlen_t pos = 0, np1 = (R_xlen_t) n + 1;
    
    /* Non-negative diagonal elements are necessary but _not_ sufficient */
    for (int i = 0; i < n; ++i, pos += np1)
	if (x[pos] < 0)
	    return mkString(_("'dpoMatrix' is not positive semidefinite"));
    return ScalarLogical(1);
}

SEXP corMatrix_validate(SEXP obj)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    SEXP sd = GET_SLOT(obj, install("sd"));
    if (XLENGTH(sd) != n)
	return mkString(_("length of 'sd' slot is not equal to n=Dim[1]"));
    double *psd = REAL(sd);
    for (int i = 0; i < n; ++i) {
	if (!R_FINITE(psd[i]))
	    return mkString(_("'sd' slot has nonfinite elements"));
	if (psd[i] < 0)
	    return mkString(_("'sd' slot has negative elements"));
    }
    return ScalarLogical(1);
}

SEXP dpoMatrix_chol(SEXP x)
{
    SEXP val = get_factor(x, "Cholesky"),
	dimP = GET_SLOT(x, Matrix_DimSym),
	uploP = GET_SLOT(x, Matrix_uploSym);
    const char *uplo = CHAR(STRING_ELT(uploP, 0));
    int *dims = INTEGER(dimP), info;
    int n = dims[0];
    const R_xlen_t n2 = ((R_xlen_t)n) * n; // = n^2
    double *vx;

    if (val != R_NilValue) return val;// use x@factors$Cholesky if available
    dims = INTEGER(dimP);
    val = PROTECT(NEW_OBJECT_OF_CLASS("Cholesky"));
    SET_SLOT(val, Matrix_uploSym, duplicate(uploP));
    SET_SLOT(val, Matrix_diagSym, mkString("N"));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    set_symmetrized_DimNames(val, GET_SLOT(x, Matrix_DimNamesSym), -1);
    vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n2));
    AZERO(vx, n2, 0.0);
    F77_CALL(dlacpy)(uplo, &n, &n, REAL(GET_SLOT(x, Matrix_xSym)),
		     &n, vx, &n FCONE);
    if (n > 0) {
	F77_CALL(dpotrf)(uplo, &n, vx, &n, &info  FCONE);
	if (info) {
	    if(info > 0)
		error(_("the leading minor of order %d is not positive definite"),
		      info);
	    else /* should never happen! */
		error(_("Lapack routine %s returned error code %d"), "dpotrf", info);
	}
    }
    set_factor(x, "Cholesky", val);
    UNPROTECT(1);
    return val;
}

SEXP dpoMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP Chol = dpoMatrix_chol(obj);
    const char typnm[] = {'O', '\0'};	/* always use the one norm */
    int *dims = INTEGER(GET_SLOT(Chol, Matrix_DimSym)), info;
    double anorm = get_norm_sy(obj, typnm), rcond;

    F77_CALL(dpocon)(uplo_P(Chol),
		     dims, REAL(GET_SLOT(Chol, Matrix_xSym)),
		     dims, &anorm, &rcond,
		     (double *) R_alloc(3*dims[0], sizeof(double)),
		     (int *) R_alloc(dims[0], sizeof(int)), &info FCONE);
    return ScalarReal(rcond);
}

SEXP dpoMatrix_solve(SEXP x)
{
    SEXP Chol = dpoMatrix_chol(x);
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dpoMatrix"));
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)), info;

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    slot_dup(val, Chol, Matrix_uploSym);
    slot_dup(val, Chol, Matrix_xSym);
    slot_dup(val, Chol, Matrix_DimSym);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(x, Matrix_DimNamesSym)));
    F77_CALL(dpotri)(uplo_P(val), dims,
		     REAL(GET_SLOT(val, Matrix_xSym)), dims, &info FCONE);
    UNPROTECT(1);
    return val;
}

SEXP dpoMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP Chol = dpoMatrix_chol(a),
	val = PROTECT(dense_as_general(b, 'd', 2, 0));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym)),
	info;

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    F77_CALL(dpotrs)(uplo_P(Chol), adims, bdims + 1,
		     REAL(GET_SLOT(Chol, Matrix_xSym)), adims,
		     REAL(GET_SLOT(val, Matrix_xSym)),
		     bdims, &info FCONE);
    UNPROTECT(1);
    return val;
}
