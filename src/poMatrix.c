#include "poMatrix.h"

SEXP poMatrix_chol(SEXP x)
{
    SEXP val = get_factorization(x, "Cholesky");
    int *dims, info;
    
    if (val != R_NilValue) return val;
    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    val = PROTECT(NEW_OBJECT(MAKE_CLASS("Cholesky")));
    SET_SLOT(val, install("uplo"), duplicate(GET_SLOT(x, install("uplo"))));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(x, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    F77_CALL(dpotrf)(CHAR(asChar(GET_SLOT(val, Matrix_uploSym))), dims,
		     REAL(GET_SLOT(val, Matrix_xSym)), dims, &info);
    if (info) error("Lapack routine dpotrf returned error code %d", info);
    UNPROTECT(1);
    return set_factorization(x, val, "Cholesky");
}

static
double set_rcond(SEXP obj, char *typstr)
{
    char typnm[] = {'O', '\0'};	/* always use the one norm */
    SEXP rcv = GET_SLOT(obj, install("rcond"));
    double rcond = get_double_by_name(rcv, typnm);

    if (R_IsNA(rcond)) {
        SEXP Chol = poMatrix_chol(obj);
	int *dims = INTEGER(GET_SLOT(Chol, Matrix_DimSym)), info;
	double anorm = get_norm_sy(obj, typnm);

	F77_CALL(dpocon)(CHAR(asChar(GET_SLOT(Chol, Matrix_uploSym))),
			 dims, REAL(GET_SLOT(Chol, Matrix_xSym)),
			 dims, &anorm, &rcond,
			 (double *) R_alloc(3*dims[0], sizeof(double)),
			 (int *) R_alloc(dims[0], sizeof(int)), &info);
	SET_SLOT(obj, install("rcond"),
		 set_double_by_name(rcv, rcond, typnm));
    }
    return rcond;
}

SEXP poMatrix_rcond(SEXP obj, SEXP type)
{
  return ScalarReal(set_rcond(obj, CHAR(asChar(type))));
}

SEXP poMatrix_solve(SEXP x)
{
    SEXP Chol = poMatrix_chol(x);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("poMatrix")));
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)), info;

    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(x, Matrix_uploSym)));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(Chol, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    F77_CALL(dpotri)(CHAR(asChar(GET_SLOT(val, Matrix_uploSym))),
		     dims, REAL(GET_SLOT(val, Matrix_xSym)), dims, &info);
    UNPROTECT(1);
    return val;
}

SEXP poMatrix_geMatrix_solve(SEXP a, SEXP b)
{
    SEXP Chol = poMatrix_chol(a),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	info;

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error("Dimensions of system to be solved are inconsistent");
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(b, Matrix_DimSym)));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(b, Matrix_xSym)));
    F77_CALL(dpotrs)(CHAR(asChar(GET_SLOT(Chol, Matrix_uploSym))),
		     adims, bdims + 1,
		     REAL(GET_SLOT(Chol, Matrix_xSym)), adims,
		     REAL(GET_SLOT(val, Matrix_xSym)),
		     bdims, &info);
    UNPROTECT(1);
    return val;
}

SEXP poMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP Chol = poMatrix_chol(a),
	val = PROTECT(duplicate(b));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(getAttrib(b, R_DimSymbol)),
	info;

    if (!(isReal(b) && isMatrix(b)))
	error("Argument b must be a numeric matrix");
    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error("Dimensions of system to be solved are inconsistent");
    F77_CALL(dpotrs)(CHAR(asChar(GET_SLOT(Chol, Matrix_uploSym))),
		     adims, bdims + 1,
		     REAL(GET_SLOT(Chol, Matrix_xSym)), adims,
		     REAL(val), bdims, &info);
    UNPROTECT(1);
    return val;
}