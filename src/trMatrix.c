#include "trMatrix.h"

SEXP trMatrix_validate(SEXP obj)
{
    SEXP uplo = GET_SLOT(obj, install("uplo")),
	diag = GET_SLOT(obj, install("uplo"));
    char *val;
    
    if (length(uplo) != 1)
	return ScalarString(mkChar("uplo slot must have length 1"));
    if (length(diag) != 1)
	return ScalarString(mkChar("diag slot must have length 1"));
    val = CHAR(STRING_ELT(uplo, 0));
    if (strlen(val) != 1) 
    	return ScalarString(mkChar("uplo[1] must have string length 1"));
    if (toupper(*val) != 'U' && toupper(*val) != 'L')
    	return ScalarString(mkChar("uplo[1] must be \"U\" or \"L\""));
    val = CHAR(STRING_ELT(diag, 0));
    if (strlen(val) != 1) 
    	return ScalarString(mkChar("diag[1] must have string length 1"));
    if (toupper(*val) != 'U' && toupper(*val) != 'N')
    	return ScalarString(mkChar("diag[1] must be \"U\" or \"N\""));
    return ScalarLogical(1);
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
    return F77_CALL(dlantr)(typnm,
			    CHAR(asChar(GET_SLOT(obj, install("uplo")))),
			    CHAR(asChar(GET_SLOT(obj, install("diag")))),
			    dims, dims+1,
			    REAL(GET_SLOT(obj, install("x"))),
			    dims, work);
}


SEXP trMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm(obj, CHAR(asChar(type))));
}

static
double set_rcond(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    SEXP rcv = GET_SLOT(obj, install("rcond"));
    double rcond = get_double_by_name(rcv, typnm);

    typnm[0] = rcond_type(typstr);
    if (R_IsNA(rcond)) {
	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), info;
	F77_CALL(dtrcon)(typnm,
			 CHAR(asChar(GET_SLOT(obj, install("uplo")))),
			 CHAR(asChar(GET_SLOT(obj, install("diag")))),
			 dims, REAL(GET_SLOT(obj, install("x"))),
			 dims, &rcond,
			 (double *) R_alloc(3*dims[0], sizeof(double)),
			 (int *) R_alloc(dims[0], sizeof(int)), &info);
	SET_SLOT(obj, install("rcond"),
		 set_double_by_name(rcv, rcond, typnm));
    }
    return rcond;
}

SEXP trMatrix_rcond(SEXP obj, SEXP type)
{
    return ScalarReal(set_rcond(obj, CHAR(asChar(type))));
}

SEXP trMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(duplicate(a));
    int info, *Dim = INTEGER(GET_SLOT(val, Matrix_DimSym));
    F77_CALL(dtrtri)(CHAR(asChar(GET_SLOT(val, install("uplo")))),
		     CHAR(asChar(GET_SLOT(val, install("diag")))),
		     Dim, REAL(GET_SLOT(val, install("x"))), Dim, &info);
    UNPROTECT(1);
    return val;
}

void make_array_triangular(double *to, SEXP from)
{
    int i, j, *dims = INTEGER(GET_SLOT(from, Matrix_DimSym));
    int n = dims[0], m = dims[1];

    if (toupper(*CHAR(asChar(GET_SLOT(from, install("uplo"))))) == 'U') {
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
    if (toupper(*CHAR(asChar(GET_SLOT(from, install("diag"))))) == 'U') {
	j = (n < m) ? n : m;
	for (i = 0; i < j; i++) {
	    to[i * (m + 1)] = 1.;
	}
    }
}
    
SEXP trMatrix_as_geMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("geMatrix")));
    
    SET_SLOT(val, install("rcond"),
	     duplicate(GET_SLOT(from, install("rcond"))));
    SET_SLOT(val, install("x"), duplicate(GET_SLOT(from, install("x"))));
    SET_SLOT(val, Matrix_DimSym,
	     duplicate(GET_SLOT(from, Matrix_DimSym)));
    make_array_triangular(REAL(GET_SLOT(val, install("x"))), from);
    UNPROTECT(1);
    return val;
}

SEXP trMatrix_as_matrix(SEXP from)
{
    int *Dim = INTEGER(GET_SLOT(from, Matrix_DimSym));
    int m = Dim[0], n = Dim[1];
    SEXP val = PROTECT(allocMatrix(REALSXP, m, n));
    
    make_array_triangular(Memcpy(REAL(val),
				 REAL(GET_SLOT(from, install("x"))), m * n),
			  from);
    UNPROTECT(1);
    return val;
}

SEXP trMatrix_getDiag(SEXP x)
{
    int i, n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    SEXP ret = PROTECT(allocVector(REALSXP, n)),
	xv = GET_SLOT(x, install("x"));

    if ('U' == toupper(CHAR(STRING_ELT(GET_SLOT(x, install("diag")), 0))[0])) {
	for (i = 0; i < n; i++) REAL(ret)[i] = 1.;
    } else {
	for (i = 0; i < n; i++) {
	    REAL(ret)[i] = REAL(xv)[i * (n + 1)];
	}
    }
    UNPROTECT(1);
    return ret;
}
    
SEXP trMatrix_geMatrix_mm(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	m = adims[0], n = bdims[1], k = adims[1];
    SEXP val = PROTECT(duplicate(b));
    double one = 1.;
    
    if (bdims[0] != k)
	error("Matrices are not conformable for multiplication");
    if (m < 1 || n < 1 || k < 1)
	error("Matrices with zero extents cannot be multiplied");
    F77_CALL(dtrmm)("L", CHAR(asChar(GET_SLOT(a, Matrix_uploSym))), "N",
		    CHAR(asChar(GET_SLOT(a, Matrix_diagSym))),
		    adims, bdims+1, &one,
		    REAL(GET_SLOT(a, Matrix_xSym)), adims,
		    REAL(GET_SLOT(val, Matrix_xSym)), bdims);
    UNPROTECT(1);
    return val;
}

SEXP trMatrix_geMatrix_mm_R(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	m = adims[0], n = bdims[1], k = adims[1];
    SEXP val = PROTECT(duplicate(b));
    double one = 1.;
    
    if (bdims[0] != k)
	error("Matrices are not conformable for multiplication");
    if (m < 1 || n < 1 || k < 1)
	error("Matrices with zero extents cannot be multiplied");
    F77_CALL(dtrmm)("R", CHAR(asChar(GET_SLOT(a, Matrix_uploSym))), "N",
		    CHAR(asChar(GET_SLOT(a, Matrix_diagSym))),
		    adims, bdims+1, &one,
		    REAL(GET_SLOT(a, Matrix_xSym)), adims,
		    REAL(GET_SLOT(val, Matrix_xSym)), bdims);
    UNPROTECT(1);
    return val;
}
