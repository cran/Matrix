#include "geMutils.h"

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

SEXP Matrix_init(void)
{
    Matrix_DimSym = install("Dim");
    Matrix_diagSym = install("diag");
    Matrix_iSym = install("i");
    Matrix_pSym = install("p");
    Matrix_uploSym = install("uplo");
    Matrix_xSym = install("x");
    Matrix_zSym = install("z");
    return R_NilValue;
}

SEXP cscMatrix_set_Dim(SEXP x, int nrow)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    dims[0] = nrow;
    dims[1] = length(GET_SLOT(x, Matrix_pSym)) - 1;
    return x;
}
