#include "LU.h"

SEXP LU_expand(SEXP x)
{
    SEXP L = PROTECT(NEW_OBJECT(MAKE_CLASS("trMatrix"))),
	U = PROTECT(NEW_OBJECT(MAKE_CLASS("trMatrix"))),
	val = PROTECT(allocVector(VECSXP, 2)),
	nms = PROTECT(allocVector(STRSXP, 2)),
	lux = GET_SLOT(x, Matrix_xSym),
	dd = GET_SLOT(x, Matrix_DimSym);

    SET_STRING_ELT(nms, 0, mkChar("L"));
    SET_STRING_ELT(nms, 1, mkChar("U"));
    setAttrib(val, R_NamesSymbol, nms);
    SET_VECTOR_ELT(val, 0, L);
    SET_VECTOR_ELT(val, 1, U);
    SET_SLOT(L, Matrix_xSym, duplicate(lux));
    SET_SLOT(L, Matrix_DimSym, dd);
    SET_SLOT(L, Matrix_uploSym, ScalarString(mkChar("L")));
    SET_SLOT(L, Matrix_diagSym, ScalarString(mkChar("U")));
    make_array_triangular(REAL(GET_SLOT(L, Matrix_xSym)), L);
    SET_SLOT(U, Matrix_xSym, duplicate(lux));
    SET_SLOT(U, Matrix_DimSym, dd);
    make_array_triangular(REAL(GET_SLOT(U, Matrix_xSym)), U);
    UNPROTECT(4);
    return val;
}

    
	     
