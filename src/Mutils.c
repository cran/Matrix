#include "Mutils.h"
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

/*MM: this is useful for all the ..CMatrix classes
  (and ..R by [0] <-> [1]): */
SEXP dgCMatrix_set_Dim(SEXP x, int nrow)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    dims[0] = nrow;
    dims[1] = length(GET_SLOT(x, Matrix_pSym)) - 1;
    return x;
}



/**  The following two csc_ functions are identically usable for rcs__
 *
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
    double *dd = (double *) NULL;

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

/* Fill in the "trivial remainder" in  n*m  array ;
 *  typically the 'x' slot of a "dtrMatrix" :
 * But should be usable for double/logical/int/complex : */

#define MAKE_TRIANGULAR_BODY(_TO_, _FROM_, _ZERO_, _ONE_)	\
{								\
    int i, j, *dims = INTEGER(GET_SLOT(_FROM_, Matrix_DimSym));	\
    int n = dims[0], m = dims[1];				\
								\
    if (*uplo_P(_FROM_) == 'U') {				\
	for (j = 0; j < n; j++)					\
	    for (i = j+1; i < m; i++)				\
		_TO_[i + j*m] = _ZERO_;				\
    } else {							\
	for (j = 1; j < n; j++)					\
	    for (i = 0; i < j && i < m; i++)			\
		_TO_[i + j*m] = _ZERO_;				\
    }								\
    if (*diag_P(_FROM_) == 'U') {				\
	j = (n < m) ? n : m;					\
	for (i = 0; i < j; i++)					\
	    _TO_[i * (m + 1)] = _ONE_;				\
    }								\
}

void make_d_matrix_triangular(double *to, SEXP from)
    MAKE_TRIANGULAR_BODY(to, from, 0., 1.)
void make_i_matrix_triangular(int *to, SEXP from)
    MAKE_TRIANGULAR_BODY(to, from, 0, 1)


/* Should work for double/logical/int/complex : */
#define MAKE_SYMMETRIC_BODY(_TO_, _FROM_)			\
{								\
    int i, j, n = INTEGER(GET_SLOT(_FROM_, Matrix_DimSym))[0];	\
								\
    if (*uplo_P(_FROM_) == 'U') {				\
	for (j = 0; j < n; j++)					\
	    for (i = j+1; i < n; i++)				\
		_TO_[i + j*n] = _TO_[j + i*n];			\
    } else {							\
	for (j = 1; j < n; j++)					\
	    for (i = 0; i < j && i < n; i++)			\
		_TO_[i + j*n] = _TO_[j + i*n];			\
    }								\
}

void make_d_matrix_symmetric(double *to, SEXP from)
    MAKE_SYMMETRIC_BODY(to, from)

void make_i_matrix_symmetric(int *to, SEXP from)
    MAKE_SYMMETRIC_BODY(to, from)


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

SEXP dense_nonpacked_validate(SEXP obj)
{
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if ((dims[0] * dims[1]) != length(GET_SLOT(obj, Matrix_xSym)))
	return mkString(_("length of x slot != prod(Dim)"));
    return ScalarLogical(1);
}


#define PACKED_TO_FULL(TYPE)						\
TYPE *packed_to_full_ ## TYPE(TYPE *dest, const TYPE *src,		\
		        int n, enum CBLAS_UPLO uplo)			\
{									\
    int i, j, pos = 0;							\
									\
    AZERO(dest, n*n);							\
    for (j = 0; j < n; j++) {						\
	switch(uplo) {							\
	case UPP:							\
	    for (i = 0; i <= j; i++) dest[i + j * n] = src[pos++];	\
	    break;							\
	case LOW:							\
	    for (i = j; i < n; i++) dest[i + j * n] = src[pos++];	\
	    break;							\
	default:							\
	    error(_("'uplo' must be UPP or LOW"));			\
	}								\
    }									\
    return dest;							\
}

PACKED_TO_FULL(double)
PACKED_TO_FULL(int)

#define FULL_TO_PACKED(TYPE)						\
TYPE *full_to_packed_ ## TYPE(TYPE *dest, const TYPE *src, int n,	\
		      enum CBLAS_UPLO uplo, enum CBLAS_DIAG diag)	\
{									\
    int i, j, pos = 0;							\
									\
    for (j = 0; j < n; j++) {						\
	switch(uplo) {							\
	case UPP:							\
	    for (i = 0; i <= j; i++)					\
		dest[pos++] = (i == j && diag== UNT) ? 1 : src[i + j*n]; \
	    break;							\
	case LOW:							\
	    for (i = j; i < n; i++)					\
		dest[pos++] = (i == j && diag== UNT) ? 1 : src[i + j*n]; \
	    break;							\
	default:							\
	    error(_("'uplo' must be UPP or LOW"));			\
	}								\
    }									\
    return dest;							\
}

FULL_TO_PACKED(double)
FULL_TO_PACKED(int)



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

    if (*uplo_P(x) == 'U') {
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


/**
 * Return the element of a given name from a named list
 *
 * @param list
 * @param nm name of desired element
 *
 * @return element of list with name nm
 */
SEXP
Matrix_getElement(SEXP list, char *nm) {
    SEXP names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < LENGTH(list); i++)
	if (!strcmp(CHAR(STRING_ELT(names, i)), nm))
	    return(VECTOR_ELT(list, i));
    return R_NilValue;
}

/**
 * Allocate a real classed matrix
 *
 * @param class character string of the type of Matrix to allocate
 * @param nrow number of rows
 * @param ncol number of columns
 *
 * @return pointer to a classed real matrix
 */
SEXP alloc_real_classed_matrix(char *class, int nrow, int ncol)
{
    SEXP val = NEW_OBJECT(MAKE_CLASS(class));
    int *dims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));

    dims[0] = nrow; dims[1] = ncol;
    ALLOC_SLOT(val, Matrix_xSym, REALSXP, nrow * ncol);
    return val;
}

SEXP alloc_dgeMatrix(int m, int n, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = m; dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, m * n);
    dn = ALLOC_SLOT(ans, Matrix_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP alloc_dpoMatrix(int n, char *uplo, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * n);
    SET_SLOT(ans, Matrix_uploSym, mkString(uplo));
    dn = ALLOC_SLOT(ans, Matrix_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP alloc_dtrMatrix(int n, char *uplo, char *diag, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * n);
    SET_SLOT(ans, Matrix_uploSym, mkString(uplo));
    SET_SLOT(ans, Matrix_diagSym, mkString(diag));
    dn = ALLOC_SLOT(ans, Matrix_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP alloc_dsCMatrix(int n, int nz, char *uplo, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz);
    ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz);
    ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1);
    SET_SLOT(ans, Matrix_uploSym, mkString(uplo));
    dn = ALLOC_SLOT(ans, Matrix_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

/**
 * Zero a square matrix of size nc then copy a vector to the diagonal
 *
 * @param dest destination array of length nc * nc
 * @param src diagonal elements in an array of length nc
 * @param nc number of columns (and rows) in the matrix
 *
 * @return dest
 */

static double *
install_diagonal(double *dest, SEXP A)
{
    int nc = INTEGER(GET_SLOT(A, Matrix_DimSym))[0];
    int i, ncp1 = nc + 1, unit = *diag_P(A) == 'U';
    double *ax = REAL(GET_SLOT(A, Matrix_xSym));

    AZERO(dest, nc * nc);
    for (i = 0; i < nc; i++)
	dest[i * ncp1] = (unit) ? 1. : ax[i];
    return dest;
}

/**  Duplicate a ddenseMatrix or a numeric matrix or vector
 *  as a dgeMatrix.
 *  This is for the many *_matrix_{prod,crossprod,tcrossprod,etc.}
 *  functions that work with both classed and unclassed matrices.
 *
 *
 * @param A	  either a ddenseMatrix object or a matrix object
 */

SEXP dup_mMatrix_as_dgeMatrix(SEXP A)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	ad = R_NilValue , an = R_NilValue;	/* -Wall */
    char *cl = class_P(A),
	*valid[] = {"_NOT_A_CLASS_", "dgeMatrix", "dtrMatrix",
		    "dsyMatrix", "dpoMatrix", "ddiMatrix",
		    "dtpMatrix", "dspMatrix", "dppMatrix",
		    /* sub classes of those above:*/
		    /* dtr */ "Cholesky", "LDL", "BunchKaufman",
		    /* dtp */ "pCholesky", "pBunchKaufman",
		    /* dpo */ "corMatrix",
		    ""};
    int ctype = Matrix_check_class(cl, valid), nprot = 1, sz;
    double *ansx;

    if (ctype > 0) {		/* a ddenseMatrix object */
	ad = GET_SLOT(A, Matrix_DimSym);
	an = GET_SLOT(A, Matrix_DimNamesSym);
    }
    else if (ctype < 0) {	/* not a (recognized) classed matrix */
	if (isMatrix(A)) {	/* "matrix" */
	    ad = getAttrib(A, R_DimSymbol);
	    an = getAttrib(A, R_DimNamesSymbol);
	} else { /* maybe "numeric" (incl "integer","logical") -->  (n x 1) */
	    int* dd = INTEGER(ad = PROTECT(allocVector(INTSXP, 2)));
	    nprot++;
	    dd[0] = LENGTH(A); dd[1] = 1;
	    an = R_NilValue;
	}
	if (isInteger(A) || isLogical(A)) {
	    A = PROTECT(coerceVector(A, REALSXP));
	    nprot++;
	}
	if (!isReal(A))
	    error(_("invalid class `%s' to dup_mMatrix_as_dgeMatrix"), cl);
	ctype = 0;
    }

    SET_SLOT(ans, Matrix_DimSym, duplicate(ad));
    SET_SLOT(ans, Matrix_DimNamesSym, (LENGTH(an) == 2) ? duplicate(an) :
	     allocVector(VECSXP, 2));
    sz = INTEGER(ad)[0] * INTEGER(ad)[1];
    ansx = REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, sz));
    switch(ctype) {
    case 0:			/* unclassed real matrix */
	Memcpy(ansx, REAL(A), sz);
	break;
    case 1:			/* dgeMatrix */
	Memcpy(ansx, REAL(GET_SLOT(A, Matrix_xSym)), sz);
	break;
    case 2:			/* dtrMatrix   and subclasses */
    case 9: case 10: case 11:   /* ---	Cholesky, LDL, BunchKaufman */
	Memcpy(ansx, REAL(GET_SLOT(A, Matrix_xSym)), sz);
	make_d_matrix_triangular(ansx, A);
	break;
    case 3:			/* dsyMatrix */
    case 4:			/* dpoMatrix  + subclass */
    case 14:	 		/* ---	corMatrix */
	Memcpy(ansx, REAL(GET_SLOT(A, Matrix_xSym)), sz);
	make_d_matrix_symmetric(ansx, A);
	break;
    case 5:			/* ddiMatrix */
	install_diagonal(ansx, A);
	break;
    case 6:			/* dtpMatrix  + subclasses */
    case 12: case 13: 		/* ---	pCholesky, pBunchKaufman */
	packed_to_full_double(ansx, REAL(GET_SLOT(A, Matrix_xSym)),
			      INTEGER(ad)[0],
			      *uplo_P(A) == 'U' ? UPP : LOW);
	make_d_matrix_triangular(ansx, A);
	break;
    case 7:			/* dspMatrix */
    case 8:			/* dppMatrix */
	packed_to_full_double(ansx, REAL(GET_SLOT(A, Matrix_xSym)),
			      INTEGER(ad)[0],
			      *uplo_P(A) == 'U' ? UPP : LOW);
	make_d_matrix_symmetric(ansx, A);
	break;
    default:
	error(_("unexpected ctype = %d in dup_mMatrix_as_dgeMatrix"), ctype);
    }
    UNPROTECT(nprot);
    return ans;
}
