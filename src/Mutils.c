#include "Mutils.h"
#include "triplet_to_col.h"
#include <R_ext/Lapack.h>

SEXP
    Matrix_DSym,
    Matrix_DIsqrtSym,
    Matrix_DimSym,
    Matrix_GpSym,
    Matrix_LIiSym,
    Matrix_LIpSym,
    Matrix_LIxSym,
    Matrix_LiSym,
    Matrix_LpSym,
    Matrix_LxSym,
    Matrix_OmegaSym,
    Matrix_ParentSym,
    Matrix_RXXSym,
    Matrix_RZXSym,
    Matrix_XtXSym,
    Matrix_ZtXSym,
    Matrix_bVarSym,
    Matrix_devianceSym,
    Matrix_devCompSym,
    Matrix_diagSym,
    Matrix_iSym,
    Matrix_ipermSym,
    Matrix_jSym,
    Matrix_matSym,
    Matrix_ncSym,
    Matrix_pSym,
    Matrix_permSym,
    Matrix_statusSym,
    Matrix_uploSym,
    Matrix_xSym,
    Matrix_zSym;

SEXP Matrix_init(void)
{
    Matrix_DSym = install("D");
    Matrix_DIsqrtSym = install("DIsqrt");
    Matrix_DimSym = install("Dim");
    Matrix_GpSym = install("Gp");
    Matrix_LIiSym = install("LIi");
    Matrix_LIpSym = install("LIp");
    Matrix_LIxSym = install("LIx");
    Matrix_LiSym = install("Li");
    Matrix_LpSym = install("Lp");
    Matrix_LxSym = install("Lx");
    Matrix_OmegaSym = install("Omega");
    Matrix_ParentSym = install("Parent");
    Matrix_RXXSym = install("RXX");
    Matrix_RZXSym = install("RZX");
    Matrix_XtXSym = install("XtX");
    Matrix_ZtXSym = install("ZtX");
    Matrix_bVarSym = install("bVar");
    Matrix_devianceSym = install("deviance");
    Matrix_devCompSym = install("devComp");
    Matrix_diagSym = install("diag");
    Matrix_iSym = install("i");
    Matrix_ipermSym = install("iperm");
    Matrix_jSym = install("j");
    Matrix_matSym = install("mat");
    Matrix_ncSym = install("nc");
    Matrix_pSym = install("p");
    Matrix_permSym = install("perm");
    Matrix_statusSym = install("status");
    Matrix_uploSym = install("uplo");
    Matrix_xSym = install("x");
    Matrix_zSym = install("z");
    return R_NilValue;
}

char norm_type(char *typstr)
{
    char typup;
    
    if (strlen(typstr) != 1)
	error("argument type[1]='%s' must be a character string of string length 1",
	      typstr);
    typup = toupper(*typstr);
    if (typup == '1') typup = 'O'; /* aliases */
    if (typup == 'E') typup = 'F';
    if (typup != 'M' && typup != 'O' && typup != 'I' && typup != 'F')
	error("argument type[1]='%s' must be one of 'M','1','O','I','F' or 'E'",
	      typstr);
    return typup;
}

char rcond_type(char *typstr)
{
    char typup;
    
    if (strlen(typstr) != 1)
	error("argument type[1]='%s' must be a character string of string length 1",
	      typstr);
    typup = toupper(*typstr);
    if (typup == '1') typup = 'O'; /* alias */
    if (typup != 'O' && typup != 'I')
	error("argument type[1]='%s' must be one of '1','O', or 'I'",
	      typstr);
    return typup;
}

double get_double_by_name(SEXP obj, char *nm)
{
    SEXP nms = getAttrib(obj, R_NamesSymbol);
    int i, len = length(obj);
    
    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error("object must be a named, numeric vector");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    return REAL(obj)[i];
	}
    }
    return R_NaReal;
}

SEXP
set_double_by_name(SEXP obj, double val, char *nm)
{
    SEXP nms = getAttrib(obj, R_NamesSymbol);
    int i, len = length(obj);

    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error("object must be a named, numeric vector");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    REAL(obj)[i] = val;
	    return obj;
	}
    }
    {SEXP nx = PROTECT(allocVector(REALSXP, len + 1)),
	 nnms = allocVector(STRSXP, len + 1);

    setAttrib(nx, R_NamesSymbol, nnms);
    for (i = 0; i < len; i++) {
	REAL(nx)[i] = REAL(obj)[i];
	SET_STRING_ELT(nnms, i, duplicate(STRING_ELT(nms, i)));
    }
    REAL(nx)[len] = val;
    SET_STRING_ELT(nnms, len, mkChar(nm));
    UNPROTECT(1);
    return nx;
    }
}

SEXP as_det_obj(double val, int log, int sign)
{
    SEXP det = PROTECT(allocVector(VECSXP, 2)),
	nms = allocVector(STRSXP, 2),
	vv = ScalarReal(val);

    setAttrib(det, R_NamesSymbol, nms);
    SET_STRING_ELT(nms, 0, mkChar("modulus"));
    SET_STRING_ELT(nms, 1, mkChar("sign"));
    setAttrib(vv, install("logarithm"), ScalarLogical(log));
    SET_VECTOR_ELT(det, 0, vv);
    SET_VECTOR_ELT(det, 1, ScalarInteger(sign));
    setAttrib(det, R_ClassSymbol, ScalarString(mkChar("det")));
    UNPROTECT(1);
    return det;
}

SEXP get_factorization(SEXP obj, char *nm)
{
    SEXP fac = GET_SLOT(obj, install("factorization")),
	nms = getAttrib(fac, R_NamesSymbol);
    int i, len = length(fac);
    
    if ((!isNewList(fac)) || (length(fac) > 0 && nms == R_NilValue))
	error("factorization slot must be a named list");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    return VECTOR_ELT(fac, i);
	}
    }
    return R_NilValue;
}

SEXP set_factorization(SEXP obj, SEXP val, char *nm)
{
    SEXP fac = GET_SLOT(obj, install("factorization")),
	nms = getAttrib(fac, R_NamesSymbol), nfac, nnms;
    int i, len = length(fac);

    if ((!isNewList(fac)) || (length(fac) > 0 && nms == R_NilValue))
	error("factorization slot must be a named list");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    SET_VECTOR_ELT(fac, i, val);
	    return val;
	}
    }
    nfac = allocVector(VECSXP, len + 1);
    nnms = allocVector(STRSXP, len + 1);
    setAttrib(nfac, R_NamesSymbol, nnms);
    for (i = 0; i < len; i++) {
	SET_VECTOR_ELT(nfac, i, VECTOR_ELT(fac, i));
	SET_STRING_ELT(nnms, i, duplicate(STRING_ELT(nms, i)));
    }
    SET_VECTOR_ELT(nfac, len, val);
    SET_STRING_ELT(nnms, len, mkChar(nm));
    SET_SLOT(obj, install("factorization"), nfac);
    return val;
}

SEXP cscMatrix_set_Dim(SEXP x, int nrow)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    dims[0] = nrow;
    dims[1] = length(GET_SLOT(x, Matrix_pSym)) - 1;
    return x;
}

int csc_unsorted_columns(int ncol, const int p[], const int i[])
{
    int j;
    for (j = 0; j < ncol; j++) {
	int ind, lst = p[j+1] - 1;
	for (ind = p[j]; ind < lst; ind++) {
	    if (i[ind] > i[ind+1]) return 1;
	}
    }
    return 0;
}

void csc_sort_columns(int ncol, const int p[], int i[], double x[])
{
    int j, maxdiff, *ord;
    double *dd;
	
    maxdiff = -1;
    for (j = 0; j < ncol; j++) {
	int diff = p[j+1] - p[j];
	if (diff > maxdiff) maxdiff = diff;
    }
    ord = Calloc(maxdiff, int);
    dd = Calloc(maxdiff, double);
    for (j = 0; j < ncol; j++) {
	int cLen = p[j+1] - p[j];
	if (cLen > 1) {
	    int k, offset = p[j];
	    for (k = 0; k < cLen; k++) ord[k] = k;
	    R_qsort_int_I(i + offset, ord, 1, cLen);
	    for (k = 0; k < cLen; k++) dd[k] = x[ord[k] + offset];
	    Memcpy(x + offset, dd, cLen);
	}
    }
    Free(ord); Free(dd);
}

SEXP csc_check_column_sorting(SEXP m)
{
    int *mp = INTEGER(GET_SLOT(m, Matrix_pSym)),
	*mi = INTEGER(GET_SLOT(m, Matrix_iSym)),
	ncol = INTEGER(GET_SLOT(m, Matrix_DimSym))[1];

    if (csc_unsorted_columns(ncol, mp, mi))
	csc_sort_columns(ncol, mp, mi, REAL(GET_SLOT(m, Matrix_xSym)));
    return m;
}

SEXP triple_as_SEXP(int nrow, int ncol, int nz,
		    const int Ti [], const int Tj [], const double Tx [],
		    char *Rclass)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS(Rclass)));
    int *Ai, *Ap;
    double *Ax;

    SET_SLOT(val, Matrix_pSym, allocVector(INTSXP, ncol + 1));
    Ap = INTEGER(GET_SLOT(val, Matrix_pSym));
    Ai = Calloc(nz, int); Ax = Calloc(nz, double);
    triplet_to_col(nrow, ncol, nz, Ti, Tj, Tx, Ap, Ai, Ax);
    nz = Ap[ncol];
    SET_SLOT(val, Matrix_iSym, allocVector(INTSXP, nz));
    Memcpy(INTEGER(GET_SLOT(val, Matrix_iSym)), Ai, nz); Free(Ai);
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, nz));
    Memcpy(REAL(GET_SLOT(val, Matrix_xSym)), Ax, nz); Free(Ax);
    UNPROTECT(1);
    return cscMatrix_set_Dim(val, nrow);
}    

void csc_components_transpose(int m, int n, int nnz,
			      const int xp[], const int xi[],
			      const double xx[],
			      int ap[], int ai[], double ax[])
{
    int k, kk,
	*ind = (int *) R_alloc(nnz, sizeof(int)),
	*aj = (int *) R_alloc(nnz, sizeof(int));

    Memcpy(aj, xi, nnz);	/* copy xi into aj and sort */
    for (k = 0; k < nnz; k++) ind[k] = k;
    R_qsort_int_I(aj, ind, 1, nnz);

    ap[0] = 0; kk = 0;		/* generate ap from aj */
    for (k = 1; k < m; k++) {
	while (aj[kk] < k) kk++;
	ap[k] = kk;
    }
    ap[m] = nnz;

    for (k = 0; k < n; k++) { /* overwrite aj with (implicit) xj */
	for (kk = xp[k]; kk < xp[k+1]; kk++) aj[kk] = k;
    }
    for (k = 0; k < nnz; k++) {	/* write ax and ai from xx and xj */
	kk = ind[k];
	ax[k] = xx[kk];
	ai[k] = aj[kk];
    }
    if (csc_unsorted_columns(m, ap, ai)) csc_sort_columns(m, ap, ai, ax);
}

void ssc_symbolic_permute(int n, int upper, const int perm[],
			  int Ap[], int Ai[])
{
    int
	j, k,
	nnz = Ap[n],
	*Aj = Calloc(nnz, int),
	*ord = Calloc(nnz, int),
	*ii = Calloc(nnz, int);
    
    for (j = 0; j < n; j++) {
	int pj = perm[j];
	for (k = Ap[j]; k < Ap[j+1]; k++) {
	    Aj[k] = pj;
	}
    }
    for (k = 0; k < nnz; k++) {
	Ai[k] = perm[Ai[k]];
	ord[k] = k;
	if ((upper && Ai[k] > Aj[k]) || (!upper && Ai[k] < Aj[k])) {
	    int tmp = Ai[k]; Ai[k] = Aj[k]; Aj[k] = tmp;
	}
    }
    R_qsort_int_I(Aj, ord, 1, nnz); /* sort Aj carrying along ind */

    k = nnz - 1;
    for (j = n - 1; j >= 0; j--) {	/* generate new Ap */
	for(; Aj[k] >= j && k >= 0; k--) Ap[j] = k;
    }
    for (k = 0; k < nnz; k++) ii[k] = Ai[ord[k]];
    Memcpy(Ai, ii, nnz);
    for (j = 0; j < n; j++) R_isort(Ai + Ap[j], Ap[j+1] - Ap[j]);
    Free(Aj); Free(ord); Free(ii);
}
    
    
/** 
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 * 
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
double *
nlme_symmetrize(double *a, const int nc)
{
    int i, j;

    for (i = 1; i < nc; i++) 
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

/** 
 * Check the error code returned by an Lapack routine and create an
 * appropriate error message.
 * 
 * @param info Error code as returned from the Lapack routine
 * @param laName Character string containing the name of the Lapack routine
 */
void
nlme_check_Lapack_error(int info, const char *laName)
{
    if (info != 0) {
        if (info > 0)
            error("error code %d from Lapack routine %s", info, laName);
        error("argument no. %d to Lapack routine %s is illegal",
              -info, laName);
    }
}

/** 
 * Calculate the inner product of vec(nlev*D^{-1} - A'A)/2 and the
 * pdgradient array regarded as a nc*nc by plen matrix.  This
 * calculation is used in several of the LMEgradient methods.
 * 
 * @param factor The nc by nc factor of the pdMat object
 * @param A The nc by nc matrix A from the LME decomposition.
 * @param nlev The number of groups associated with the random effect
 * @param nc The number of columns in the matrix
 * @param pdgradient A pdgradient object of dimension nc by nc by plen
 * @param value An array of length plen in which the gradient will be  returned
 * 
 * @return value, with the LME gradient
 */
double *
LMEgradient(const double* factor, const double* A, const int nlev,
	     const int nc, const double* pdgradient, const int plen,
	     double* value)
{
    int info, ncsq = nc*nc, one_i = 1;
    double nlev2_d = ((double) nlev)/2., mhalf_d = -0.5, one_d = 1.0,
	zero_d = 0.0;
    double *fact = Calloc(nc * nc, double);
    
    F77_CALL(dlacpy)("U", &nc, &nc, factor, &nc, fact, &nc);
    F77_CALL(dpotri)("U", &nc, fact, &nc, &info);
    nlme_check_Lapack_error(info, "dpotri");
    F77_CALL(dsyrk)("U", "T", &nc, &nc, &mhalf_d, A, &nc, &nlev2_d,
		    fact, &nc);
    nlme_symmetrize(fact, nc);
    F77_CALL(dgemv)("T", &ncsq, &plen, &one_d, pdgradient, &ncsq,
		    fact, &one_i, &zero_d, value, &one_i);
    Free(fact);
    return value;
}
