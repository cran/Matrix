#include "Mutils.h"
#include "triplet_to_col.h"
#include <R_ext/Lapack.h>

char norm_type(char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(
	    _("argument type[1]='%s' must be a character string of string length 1"),
	    typstr);
    typup = toupper(*typstr);
    if (typup == '1') typup = 'O'; /* aliases */
    if (typup == 'E') typup = 'F';
    if (typup != 'M' && typup != 'O' && typup != 'I' && typup != 'F')
	error(_("argument type[1]='%s' must be one of 'M','1','O','I','F' or 'E'"),
	      typstr);
    return typup;
}

char rcond_type(char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = toupper(*typstr);
    if (typup == '1') typup = 'O'; /* alias */
    if (typup != 'O' && typup != 'I')
	error(_("argument type[1]='%s' must be one of '1','O', or 'I'"),
	      typstr);
    return typup;
}

double get_double_by_name(SEXP obj, char *nm)
{
    SEXP nms = getAttrib(obj, R_NamesSymbol);
    int i, len = length(obj);

    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error(_("object must be a named, numeric vector"));
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
    {
	SEXP nx = PROTECT(allocVector(REALSXP, len + 1)),
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
    setAttrib(det, R_ClassSymbol, mkString("det"));
    UNPROTECT(1);
    return det;
}

SEXP get_factors(SEXP obj, char *nm)
{
    SEXP fac = GET_SLOT(obj, Matrix_factorSym),
	nms = getAttrib(fac, R_NamesSymbol);
    int i, len = length(fac);

    if ((!isNewList(fac)) || (length(fac) > 0 && nms == R_NilValue))
	error("factors slot must be a named list");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    return VECTOR_ELT(fac, i);
	}
    }
    return R_NilValue;
}

SEXP set_factors(SEXP obj, SEXP val, char *nm)
{
    SEXP fac = GET_SLOT(obj, Matrix_factorSym),
	nms = getAttrib(fac, R_NamesSymbol), nfac, nnms;
    int i, len = length(fac);

    if ((!isNewList(fac)) || (length(fac) > 0 && nms == R_NilValue))
	error("factors slot must be a named list");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    SET_VECTOR_ELT(fac, i, duplicate(val));
	    return val;
	}
    }
    nfac = PROTECT(allocVector(VECSXP, len + 1));
    nnms = PROTECT(allocVector(STRSXP, len + 1));
    setAttrib(nfac, R_NamesSymbol, nnms);
    for (i = 0; i < len; i++) {
	SET_VECTOR_ELT(nfac, i, VECTOR_ELT(fac, i));
	SET_STRING_ELT(nnms, i, duplicate(STRING_ELT(nms, i)));
    }
    SET_VECTOR_ELT(nfac, len, duplicate(val));
    SET_STRING_ELT(nnms, len, mkChar(nm));
    SET_SLOT(obj, Matrix_factorSym, nfac);
    UNPROTECT(2);
    return val;
}

SEXP dgCMatrix_set_Dim(SEXP x, int nrow)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    dims[0] = nrow;
    dims[1] = length(GET_SLOT(x, Matrix_pSym)) - 1;
    return x;
}

/**
 * Check for unsorted columns in the row indices
 *
 * @param ncol number of columns
 * @param p column pointers
 * @param i row indices
 *
 * @return 0 if all columns are sorted, otherwise 1
 */
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

/**
 * Sort the columns in a sparse column-oriented matrix so that each
 * column is in increasing order of row index.
 *
 * @param ncol number of columns
 * @param p column pointers
 * @param i row indices
 * @param x values of nonzero elements
 */
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
    if (x) dd = Calloc(maxdiff, double);
    for (j = 0; j < ncol; j++) {
	int cLen = p[j+1] - p[j];
	if (cLen > 1) {
	    int k, offset = p[j];
	    for (k = 0; k < cLen; k++) ord[k] = k;
	    R_qsort_int_I(i + offset, ord, 1, cLen);
	    if (x) {
	      for (k = 0; k < cLen; k++) dd[k] = x[ord[k] + offset];
	      Memcpy(x + offset, dd, cLen);
	    }
	}
    }
    Free(ord);
    if (x) Free(dd);
}

/**
 * Check for sorted columns in an object that inherits from the
 * dgCMatrix class.  Resort the columns if necessary.
 *
 * @param m pointer to an object that inherits from the dgCMatrix class
 *
 * @return m with the columns sorted by increasing row index
 */
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
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    UNPROTECT(1);
    return dgCMatrix_set_Dim(val, nrow);
}

/* Create the components of the transpose of a csc matrix from its components */

void csc_compTr(int m, int n, int nnz,
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
	for(; k >= 0 && Aj[k] >= j; k--) Ap[j] = k;
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
            error(_("error code %d from Lapack routine %s"), info, laName);
        error(_("argument no. %d to Lapack routine %s is illegal"),
              -info, laName);
    }
}

void make_array_triangular(double *to, SEXP from)
{
    int i, j, *dims = INTEGER(GET_SLOT(from, Matrix_DimSym));
    int n = dims[0], m = dims[1];

    if (*CHAR(asChar(GET_SLOT(from, Matrix_uploSym))) == 'U') {
	for (j = 0; j < n; j++) {
	    for (i = j+1; i < m; i++) {
		to[i + j*m] = 0.;
	    }
	}
    } else {
	for (j = 1; j < n; j++) {
	    for (i = 0; i < j && i < m; i++) {
		to[i + j*m] = 0.;
	    }
	}
    }
    if (*CHAR(asChar(GET_SLOT(from, Matrix_diagSym))) == 'U') {
	j = (n < m) ? n : m;
	for (i = 0; i < j; i++) {
	    to[i * (m + 1)] = 1.;
	}
    }
}

/**
 * Replace the value of a slot or subslot of an object in place.  This
 * routine purposely does not copy the value of obj.  Use with caution.
 *
 * @param obj object with slot to be replaced
 * @param names vector of names.  The last element is the name of the slot to replace.  The leading elements are the names of slots and subslots of obj.
 * @param value the replacement value for the slot
 *
 * @return obj, with the named slot modified in place.
 */
SEXP
nlme_replaceSlot(SEXP obj, SEXP names, SEXP value)
{
    int lnm1 = length(names) - 1;

    if (lnm1 >= 0) {
	SEXP comp = obj;
	int i;

	for (i = 0; i < lnm1; i++) {
	    comp = GET_SLOT(comp, install(CHAR(STRING_ELT(names, i))));
	}
	SET_SLOT(comp, install(CHAR(STRING_ELT(names, lnm1))), value);
    }
    return obj;
}

/**
 * Produce a weighted copy of the matrices in MLin in the storage
 * allocated to MLout
 *
 * @param MLin input matrix list
 * @param wts real vector of weights
 * @param adjst adjusted response
 * @param MLout On input a list of matrices of the same dimensions as MLin.
 *
 * @return MLout with its contents overwritten by a weighted copy of
 * MLin according to wts with adjst overwriting the response.
 */
SEXP nlme_weight_matrix_list(SEXP MLin, SEXP wts, SEXP adjst, SEXP MLout)
{
    int i, j, n, nf;
    SEXP lastM;

    if (!(isNewList(MLin) && isReal(wts) && isReal(adjst) && isNewList(MLout)))
	error(_("Incorrect argument type"));
    nf = length(MLin);
    if (length(MLout) != nf)
	error(_("Lengths of MLin (%d) and MLout (%d) must match"), nf,
	      length(MLout));
    n = length(wts);
    if (length(adjst) != n)
	error(_("Expected adjst to have length %d, got %d"), n, length(adjst));
    for (i = 0; i < nf; i++) {
	SEXP Min = VECTOR_ELT(MLin, i),
	    Mout = VECTOR_ELT(MLout, i);
	int *din, *dout, k, nc;

	if (!(isMatrix(Min) && isReal(Min)))
	    error(_("component %d of MLin is not a numeric matrix"), i + 1);
	din = INTEGER(getAttrib(Min, R_DimSymbol));
	nc = din[1];
	if (din[0] != n)
	    error(_("component %d of MLin has %d rows, expected %d"), i + 1,
		  din[0], n);
	if (!(isMatrix(Mout) && isReal(Mout)))
	    error(_("component %d of MLout is not a numeric matrix"), i + 1);
	dout = INTEGER(getAttrib(Mout, R_DimSymbol));
	if (dout[0] != n)
	    error(_("component %d of MLout has %d rows, expected %d"), i + 1,
		  dout[0], n);
	if (dout[1] != nc)
	    error(_("component %d of MLout has %d columns, expected %d"), i + 1,
		  dout[1], nc);
	for (k = 0; k < nc; k++) {
	    for (j = 0; j < n; j++) {
		REAL(Mout)[j + k * n] = REAL(Min)[j + k * n] * REAL(wts)[j];
	    }
	}
    }
    lastM = VECTOR_ELT(MLout, nf - 1);
    j = INTEGER(getAttrib(lastM, R_DimSymbol))[1] - 1;
    for (i = 0; i < n; i++)
	REAL(lastM)[j*n + i] = REAL(adjst)[i] * REAL(wts)[i];
    return MLout;
}

/**
 * Create a named vector of type TYP
 *
 * @param TYP a vector SEXP type (e.g. REALSXP)
 * @param names names of list elements with null string appended
 *
 * @return pointer to a named vector of type TYP
 */
SEXP
Matrix_make_named(int TYP, char **names)
{
    SEXP ans, nms;
    int i, n;

    for (n = 0; strlen(names[n]) > 0; n++) {}
    ans = PROTECT(allocVector(TYP, n));
    nms = PROTECT(allocVector(STRSXP, n));
    for (i = 0; i < n; i++) SET_STRING_ELT(nms, i, mkChar(names[i]));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/**
 * Allocate a 3-dimensional array
 *
 * @param mode The R mode (e.g. INTSXP)
 * @param nrow number of rows
 * @param ncol number of columns
 * @param nface number of faces
 *
 * @return A 3-dimensional array of the indicated dimensions and mode
 */
SEXP alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface)
{
    SEXP s, t;
    int n;

    if (nrow < 0 || ncol < 0 || nface < 0)
	error(_("negative extents to 3D array"));
    if ((double)nrow * (double)ncol * (double)nface > INT_MAX)
	error(_("alloc3Darray: too many elements specified"));
    n = nrow * ncol * nface;
    PROTECT(s = allocVector(mode, n));
    PROTECT(t = allocVector(INTSXP, 3));
    INTEGER(t)[0] = nrow;
    INTEGER(t)[1] = ncol;
    INTEGER(t)[2] = nface;
    setAttrib(s, R_DimSymbol, t);
    UNPROTECT(2);
    return s;
}

/**
 * Expand a column of a compressed, sparse, column-oriented matrix.
 *
 * @param dest array to hold the result
 * @param m number of rows in the matrix
 * @param j index (0-based) of column to expand
 * @param Ap array of column pointers
 * @param Ai array of row indices
 * @param Ax array of non-zero values
 *
 * @return dest
 */
double *expand_csc_column(double *dest, int m, int j,
			  const int Ap[], const int Ai[], const double Ax[])
{
    int k, k2 = Ap[j + 1];

    for (k = 0; k < m; k++) dest[k] = 0.;
    for (k = Ap[j]; k < k2; k++) dest[Ai[k]] = Ax[k];
    return dest;
}

#define Matrix_Error_Bufsiz    4096

SEXP check_scalar_string(SEXP sP, char *vals, char *nm)
{
    SEXP val = ScalarLogical(1);
    char *buf, *str;
    /* only allocate when needed: in good case, none is needed */
#define SPRINTF buf = Calloc(Matrix_Error_Bufsiz, char); sprintf

    if (length(sP) != 1) {
	SPRINTF(buf, _("'%s' slot must have length 1"), nm);
    } else {
	str = CHAR(STRING_ELT(sP, 0));
	if (strlen(str) != 1) {
	    SPRINTF(buf, _("'%s' must have string length 1"), nm);
	} else {
	    int i, len, match;
	    for (i = 0, len = strlen(vals), match = 0; i < len; i++) {
		if (str[0] == vals[i])
		    return R_NilValue;
	    }
	    SPRINTF(buf, _("'%s' must be in '%s'"), nm, vals);
	}
    }
    /* 'error' returns : */
    val = mkString(buf);
    Free(buf);
    return val;
#undef SPRINTF
}

double *packed_to_full(double *dest, const double *src, int n,
		       enum CBLAS_UPLO uplo)
{
    int i, j, pos = 0;

    AZERO(dest, n*n);
    for (j = 0; j < n; j++) {
	switch(uplo) {
	case UPP:
	    for (i = 0; i <= j; i++) dest[i + j * n] = src[pos++];
	    break;
	case LOW:
	    for (i = j; i < n; i++) dest[i + j * n] = src[pos++];
	    break;
	default:
	    error(_("'uplo' must be UPP or LOW"));
	}
    }
    return dest;
}

double *full_to_packed(double *dest, const double *src, int n,
		       enum CBLAS_UPLO uplo, enum CBLAS_DIAG diag)
{
    int i, j, pos = 0;

    for (j = 0; j < n; j++) {
	switch(uplo) {
	case UPP:
	    for (i = 0; i <= j; i++)
		dest[pos++] = (i == j && diag == UNT) ? 1. : src[i + j * n];
	    break;
	case LOW:
	    for (i = j; i < n; i++)
		dest[pos++] = (i == j && diag == UNT) ? 1. : src[i + j * n];
	    break;
	default:
	    error(_("'uplo' must be UPP or LOW"));
	}
    }
    return dest;
}

/**
 * Copy the diagonal elements of the packed array x to dest
 *
 * @param dest vector of length ncol(x)
 * @param x pointer to an object representing a packed array
 *
 * @return dest
 */
double *packed_getDiag(double *dest, SEXP x)
{
    int j, n = *INTEGER(GET_SLOT(x, Matrix_DimSym)), pos;
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

    if (*CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0)) == 'U') {
	for (pos = 0, j = 0; j < n; pos += ++j) dest[j] = xx[pos];
    } else {
	for (pos = 0, j = 0; j < n; pos += (n - j), j++) dest[j] = xx[pos];
    }
    return dest;
}

SEXP Matrix_expand_pointers(SEXP pP)
{
    int n = length(pP) - 1;
    int *p = INTEGER(pP);
    SEXP ans = PROTECT(allocVector(INTSXP, p[n]));

    expand_cmprPt(n, p, INTEGER(ans));
    UNPROTECT(1);
    return ans;
}

	
