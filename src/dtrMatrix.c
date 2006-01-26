/* double (precision) TRiangular Matrices */

#include "dtrMatrix.h"

SEXP triangularMatrix_validate(SEXP obj)
{
    SEXP val = GET_SLOT(obj, Matrix_DimSym);

    if (LENGTH(val) < 2)
	return mkString(_("'Dim' slot has length less than two"));
    if (INTEGER(val)[0] != INTEGER(val)[1])
        return mkString(_("Matrix is not square"));
    if (isString(val = check_scalar_string(GET_SLOT(obj, Matrix_uploSym),
					   "LU", "uplo"))) return val;
    if (isString(val = check_scalar_string(GET_SLOT(obj, Matrix_diagSym),
					   "NU", "diag"))) return val;
    return ScalarLogical(1);
}

SEXP dtrMatrix_validate(SEXP obj)
{
    /* since "dtr" inherits from "triangular", and "dMatrix", only need this:*/
    return dense_nonpacked_validate(obj);
}


static
double get_norm(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    double *work = (double *) NULL;

    typnm[0] = norm_type(typstr);
    if (*typnm == 'I') {
	work = (double *) R_alloc(dims[0], sizeof(double));
    }
    return F77_CALL(dlantr)(typnm, uplo_P(obj), diag_P(obj), dims, dims+1,
			    REAL(GET_SLOT(obj, Matrix_xSym)), dims, work);
}


SEXP dtrMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm(obj, CHAR(asChar(type))));
}

SEXP dtrMatrix_rcond(SEXP obj, SEXP type)
{
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), info;
    double rcond;

    typnm[0] = rcond_type(CHAR(asChar(type)));
    F77_CALL(dtrcon)(typnm, uplo_P(obj), diag_P(obj), dims,
		     REAL(GET_SLOT(obj, Matrix_xSym)), dims, &rcond,
		     (double *) R_alloc(3*dims[0], sizeof(double)),
		     (int *) R_alloc(dims[0], sizeof(int)), &info);
    return ScalarReal(rcond);
}

SEXP dtrMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(duplicate(a));
    int info, *Dim = INTEGER(GET_SLOT(val, Matrix_DimSym));
    F77_CALL(dtrtri)(uplo_P(val), diag_P(val), Dim,
		     REAL(GET_SLOT(val, Matrix_xSym)), Dim, &info);
    UNPROTECT(1);
    return val;
}

SEXP dtrMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed)
{
    int cl = asLogical(classed);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol));
    int n = bdims[0], nrhs = bdims[1];
    int sz = n * nrhs;
    double one = 1.0;

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1 || *adims != adims[1])
	error(_("Dimensions of system to be solved are inconsistent"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)), bdims, 2);
    F77_CALL(dtrsm)("L", uplo_P(a),
		    "N", diag_P(a),
		    &n, &nrhs, &one, REAL(GET_SLOT(a, Matrix_xSym)), &n,
		    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, sz)),
			   REAL(cl ? GET_SLOT(b, Matrix_xSym):b), sz), &n);
    UNPROTECT(1);
    return ans;
}

SEXP dtrMatrix_matrix_mm(SEXP a, SEXP b, SEXP classed, SEXP right)
{
    int cl = asLogical(classed), rt = asLogical(right);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol)),
	*cdims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));
    int m, n, sz;
    double one = 1.;

    if (!cl && !(isReal(b) && isMatrix(b)))
	error(_("Argument b must be a numeric matrix"));
    if (adims[0] != adims[1]) error(_("dtrMatrix in %*% must be square"));
    m = rt ? bdims[0] : adims[0];
    n = rt ? adims[1] : bdims[1];
    if ((rt && (adims[0] != m)) || (!rt && (bdims[0] != m)))
	    error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    cdims[0] = m; cdims[1] = n; sz = m * n;
    F77_CALL(dtrmm)(rt ? "R" : "L", uplo_P(a),
		    "N", diag_P(a), &m, &n,
		    &one, REAL(GET_SLOT(a, Matrix_xSym)), rt ? &n : &m,
		    Memcpy(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, sz)),
			   REAL(cl ? GET_SLOT(b, Matrix_xSym) : b), sz),
		    rt ? &m : &n);
    UNPROTECT(1);
    return val;
}

SEXP dtrMatrix_as_dgeMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));

    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(from, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(from, Matrix_DimSym)));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    make_d_matrix_triangular(REAL(GET_SLOT(val, Matrix_xSym)), from);
    UNPROTECT(1);
    return val;
}

SEXP dtrMatrix_as_matrix(SEXP from)
{
    int *Dim = INTEGER(GET_SLOT(from, Matrix_DimSym));
    int m = Dim[0], n = Dim[1];
    SEXP val = PROTECT(allocMatrix(REALSXP, m, n));
    make_d_matrix_triangular(Memcpy(REAL(val),
				    REAL(GET_SLOT(from, Matrix_xSym)), m * n),
			     from);
    setAttrib(val, R_DimNamesSymbol, GET_SLOT(from, Matrix_DimNamesSym));
    UNPROTECT(1);
    return val;
}

SEXP dtrMatrix_getDiag(SEXP x)
{
    int i, n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    SEXP ret = PROTECT(allocVector(REALSXP, n)),
	xv = GET_SLOT(x, Matrix_xSym);

    if ('U' == diag_P(x)[0]) {
	for (i = 0; i < n; i++) REAL(ret)[i] = 1.;
    } else {
	for (i = 0; i < n; i++) {
	    REAL(ret)[i] = REAL(xv)[i * (n + 1)];
	}
    }
    UNPROTECT(1);
    return ret;
}

SEXP dtrMatrix_dgeMatrix_mm_R(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	m = adims[0], n = bdims[1], k = adims[1];
    SEXP val = PROTECT(duplicate(b));
    double one = 1.;

    if (bdims[0] != k)
	error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1 || k < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    F77_CALL(dtrmm)("R", uplo_P(a), "N", diag_P(a), adims, bdims+1, &one,
		    REAL(GET_SLOT(a, Matrix_xSym)), adims,
		    REAL(GET_SLOT(val, Matrix_xSym)), bdims);
    UNPROTECT(1);
    return val;
}

SEXP dtrMatrix_as_dtpMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dtpMatrix"))),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    full_to_packed_double(
	REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, (n*(n+1))/2)),
	REAL(GET_SLOT(from, Matrix_xSym)), n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	*CHAR(STRING_ELT(diag, 0)) == 'U' ? UNT : NUN);
    UNPROTECT(1);
    return val;
}
