#include "syMatrix.h"

SEXP syMatrix_validate(SEXP obj)
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
    return ScalarLogical(1);
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
    return F77_CALL(dlansy)(typnm,
			    CHAR(asChar(GET_SLOT(obj, Matrix_uploSym))),
			    dims, REAL(GET_SLOT(obj, Matrix_xSym)),
			    dims, work);
}

SEXP syMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm_sy(obj, CHAR(asChar(type))));
}

static
double set_rcond_sy(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    SEXP rcv = GET_SLOT(obj, install("rcond"));
    double rcond;

    typnm[0] = rcond_type(typstr);
    rcond = get_double_by_name(rcv, typnm);

/* FIXME: Need a factorization here. */
    if (R_IsNA(rcond)) {
	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), info;
	double anorm = get_norm_sy(obj, "O");

	error("Code for set_rcond_sy not yet written");
	F77_CALL(dsycon)(CHAR(asChar(GET_SLOT(obj, Matrix_uploSym))),
			 dims, REAL(GET_SLOT(obj, Matrix_xSym)),
			 dims, INTEGER(GET_SLOT(obj, install("pivot"))),
			 &anorm, &rcond, 
			 (double *) R_alloc(2*dims[0], sizeof(double)),
			 (int *) R_alloc(dims[0], sizeof(int)), &info);
	SET_SLOT(obj, install("rcond"),
		 set_double_by_name(rcv, rcond, typnm));
    }
    return rcond;
}

SEXP syMatrix_rcond(SEXP obj, SEXP type)
{
/* FIXME: This is a stub */
/*     return ScalarReal(set_rcond_sy(obj, CHAR(asChar(type)))); */
    return ScalarReal(NA_REAL);
}

static 
void make_symmetric(double *to, SEXP from, int n)
{
    int i, j;
    if (toupper(*CHAR(asChar(GET_SLOT(from, Matrix_uploSym)))) == 'U') {
	for (j = 0; j < n; j++) {
	    for (i = j+1; i < n; i++) {
		to[i + j*n] = to[j + i*n];
	    }
	}
    } else {
	for (j = 1; j < n; j++) {
	    for (i = 0; i < j && i < n; i++) {
		to[i + j*n] = to[j + i*n];
	    }
	}
    }
}
    
SEXP syMatrix_solve(SEXP a)
{
/* FIXME: Write the code */
    error("code for syMatrix_solve not yet written");
    return R_NilValue;
}

SEXP syMatrix_matrix_solve(SEXP a, SEXP b)
{
/* FIXME: Write the code */
    error("code for syMatrix_matrix_solve not yet written");
    return R_NilValue;
}

SEXP syMatrix_as_geMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix"))),
	rcondSym = install("rcond");
    
    SET_SLOT(val, rcondSym, duplicate(GET_SLOT(from, rcondSym)));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(from, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym,
	     duplicate(GET_SLOT(from, Matrix_DimSym)));
    make_symmetric(REAL(GET_SLOT(val, Matrix_xSym)), from,
		   INTEGER(GET_SLOT(val, Matrix_DimSym))[0]);
    UNPROTECT(1);
    return val;
}

SEXP syMatrix_as_matrix(SEXP from)
{
    int n = INTEGER(GET_SLOT(from, Matrix_DimSym))[0];
    SEXP val = PROTECT(allocMatrix(REALSXP, n, n));
    
    make_symmetric(Memcpy(REAL(val),
			  REAL(GET_SLOT(from, Matrix_xSym)), n * n),
		   from, n);
    UNPROTECT(1);
    return val;
}

SEXP syMatrix_geMatrix_mm(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	*cdims,
	m = adims[0], n = bdims[1], k = adims[1];
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix")));
    double one = 1., zero = 0.;
    
    if (bdims[0] != k)
	error("Matrices are not conformable for multiplication");
    if (m < 1 || n < 1 || k < 1)
	error("Matrices with zero extents cannot be multiplied");
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    cdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    cdims[0] = m; cdims[1] = n;
    F77_CALL(dsymm)("L", CHAR(asChar(GET_SLOT(a, Matrix_uploSym))),
		    adims, bdims+1, &one,
		    REAL(GET_SLOT(a, Matrix_xSym)), adims,
		    REAL(GET_SLOT(b, Matrix_xSym)), bdims,
		    &zero, REAL(GET_SLOT(val, Matrix_xSym)), adims);
    UNPROTECT(1);
    return val;
}

SEXP syMatrix_geMatrix_mm_R(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	*cdims,
	m = adims[0], n = bdims[1], k = adims[1];
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix")));
    double one = 1., zero = 0.;
    
    if (bdims[0] != k)
	error("Matrices are not conformable for multiplication");
    if (m < 1 || n < 1 || k < 1)
	error("Matrices with zero extents cannot be multiplied");
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    cdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    cdims[0] = m; cdims[1] = n;
    F77_CALL(dsymm)("R", CHAR(asChar(GET_SLOT(a, Matrix_uploSym))),
		    adims, bdims+1, &one,
		    REAL(GET_SLOT(a, Matrix_xSym)), adims,
		    REAL(GET_SLOT(b, Matrix_xSym)), bdims,
		    &zero, REAL(GET_SLOT(val, Matrix_xSym)), adims);
    UNPROTECT(1);
    return val;
}
    
