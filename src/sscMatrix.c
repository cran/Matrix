#include "sscMatrix.h"

SEXP sscMatrix_validate(SEXP obj)
{
    SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
    int *Dim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    char *val;
    
    if (length(uplo) != 1)
	return ScalarString(mkChar("uplo slot must have length 1"));
    val = CHAR(STRING_ELT(uplo, 0));
    if (strlen(val) != 1) 
    	return ScalarString(mkChar("uplo[1] must have string length 1"));
    if (toupper(*val) != 'U' && toupper(*val) != 'L')
    	return ScalarString(mkChar("uplo[1] must be \"U\" or \"L\""));
    if (Dim[0] != Dim[1])
	return ScalarString(mkChar("Symmetric matrix must be square"));
    csc_check_column_sorting(obj);
    return ScalarLogical(1);
}

SEXP sscMatrix_chol(SEXP x, SEXP pivot)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("sscChol")));
    taucs_ccs_matrix *tm =
	csc_taucs_ptr(x, TAUCS_DOUBLE | TAUCS_LOWER | TAUCS_SYMMETRIC);
    int nnz, piv = asLogical(pivot);
    int *dims;
    void *L;

    if (piv) {
	int *iperm, *perm;

	SET_SLOT(val, Matrix_permSym, allocVector(INTSXP, tm->n));
	SET_SLOT(val, Matrix_ipermSym, allocVector(INTSXP, tm->n));
	perm = INTEGER(GET_SLOT(val, Matrix_permSym));
	iperm = INTEGER(GET_SLOT(val, Matrix_ipermSym));
	ssc_metis_order(tm->n, (tm->colptr)[tm->n], tm->colptr,
			tm->rowind, perm, iperm);
	tm = taucs_dccs_permute_symmetrically(tm, perm, iperm);
    }
    if (!(L = taucs_ccs_factor_llt_mf(tm)))
	error("Matrix is not positive definite");
    if (piv) taucs_dccs_free(tm);
    tm = taucs_supernodal_factor_to_ccs(L);
    taucs_supernodal_factor_free(L);
    nnz = tm->colptr[tm->n];
    SET_SLOT(val, Matrix_pSym, allocVector(INTSXP, tm->n + 1));
    Memcpy(INTEGER(GET_SLOT(val, Matrix_pSym)), tm->colptr, tm->n + 1);
    SET_SLOT(val, Matrix_iSym, allocVector(INTSXP, nnz));
    Memcpy(INTEGER(GET_SLOT(val, Matrix_iSym)), tm->rowind, nnz);
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, nnz));
    Memcpy(REAL(GET_SLOT(val, Matrix_xSym)), tm->values.d, nnz);
    dims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    dims[0] = dims[1] = tm->n;
    taucs_dccs_free(tm);
    UNPROTECT(1);
    return set_factorization(x, val, "Cholesky");
}

SEXP sscMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP Chol = get_factorization(a, "Cholesky"),
	val = PROTECT(duplicate(b));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(getAttrib(b, R_DimSymbol)),
	j, n = adims[1];
    taucs_ccs_matrix* tm;
    double *in = REAL(b), *out = REAL(val) ;

    if (!(isReal(b) && isMatrix(b)))
	error("Argument b must be a numeric matrix");
    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error("Dimensions of system to be solved are inconsistent");
    if (Chol == R_NilValue) Chol = sscMatrix_chol(a, ScalarLogical(1.));
    tm = csc_taucs_ptr(Chol, TAUCS_DOUBLE | TAUCS_LOWER | TAUCS_TRIANGULAR);
    if (!length(GET_SLOT(Chol, Matrix_permSym))) {
	for (j = 0; j < bdims[1]; j++, in += n, out += n) {
	    int errcode = taucs_dccs_solve_llt(tm, out, in);
	    if (errcode)
		error("taucs_solve returned error code %d for column %d",
		      errcode, j + 1);
	}
    } else {
	int *iperm = INTEGER(GET_SLOT(Chol, Matrix_ipermSym));
	double *tmpIn = (double *) R_alloc(n, sizeof(double)),
	    *tmpOut = (double *) R_alloc(n, sizeof(double));

	for (j = 0; j < bdims[1]; j++, in += n, out += n) {
	    int errcode, i;
				/* permute y */
	    for (i = 0; i < n; i++) tmpIn[iperm[i]] = in[i];
				/* solve */
	    errcode = taucs_dccs_solve_llt(tm, tmpOut, tmpIn);
	    if (errcode)
		error("taucs_solve returned error code %d for column %d",
		      errcode, j + 1);
				/* inverse permute b */
	    for (i = 0; i < n; i++) out[i] = tmpOut[iperm[i]];
	}
    }
    UNPROTECT(1);
    return val;
}

SEXP sscMatrix_inverse_factor(SEXP A)
{
    return mat_from_taucs(taucs_ccs_factor_xxt(
	csc_taucs_ptr(A, TAUCS_DOUBLE | TAUCS_LOWER | TAUCS_SYMMETRIC)));
}

SEXP ssc_transpose(SEXP x)
{
    SEXP
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("sscMatrix"))),
	islot = GET_SLOT(x, Matrix_iSym);
    int nnz = length(islot),
	*adims = INTEGER(GET_SLOT(ans, Matrix_DimSym)),
	*xdims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    adims[0] = xdims[1]; adims[1] = xdims[0];
    if (toupper(CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0]) == 'U')
	SET_SLOT(ans, Matrix_uploSym, ScalarString(mkChar("L")));
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, xdims[0] + 1));
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    csc_components_transpose(xdims[0], xdims[1], nnz,
			     INTEGER(GET_SLOT(x, Matrix_pSym)),
			     INTEGER(islot),
			     REAL(GET_SLOT(x, Matrix_xSym)),
			     INTEGER(GET_SLOT(ans, Matrix_pSym)),
			     INTEGER(GET_SLOT(ans, Matrix_iSym)),
			     REAL(GET_SLOT(ans, Matrix_xSym)));
    UNPROTECT(1);
    return ans;
}
