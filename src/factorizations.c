#include "factorizations.h"

SEXP LU_validate(SEXP obj)
{
    return ScalarLogical(1);
}

SEXP Cholesky_validate(SEXP obj)
{
    return ScalarLogical(1);
}

SEXP SVD_validate(SEXP obj)
{
    return ScalarLogical(1);
}

SEXP LU_expand(SEXP x)
{
    char *nms[] = {"L", "U", ""};
    SEXP L, U, val = PROTECT(Matrix_make_named(VECSXP, nms)),
	lux = GET_SLOT(x, Matrix_xSym),
	dd = GET_SLOT(x, Matrix_DimSym);

    SET_VECTOR_ELT(val, 0, NEW_OBJECT(MAKE_CLASS("dtrMatrix")));
    L = VECTOR_ELT(val, 0);
    SET_VECTOR_ELT(val, 1, NEW_OBJECT(MAKE_CLASS("dtrMatrix")));
    U = VECTOR_ELT(val, 1);
    SET_SLOT(L, Matrix_xSym, duplicate(lux));
    SET_SLOT(L, Matrix_DimSym, dd);
    SET_SLOT(L, Matrix_uploSym, mkString("L"));
    SET_SLOT(L, Matrix_diagSym, mkString("U"));
    make_array_triangular(REAL(GET_SLOT(L, Matrix_xSym)), L);
    SET_SLOT(U, Matrix_xSym, duplicate(lux));
    SET_SLOT(U, Matrix_DimSym, dd);
    SET_SLOT(U, Matrix_uploSym, mkString("U"));
    SET_SLOT(U, Matrix_diagSym, mkString("N"));
    make_array_triangular(REAL(GET_SLOT(U, Matrix_xSym)), U);
    UNPROTECT(1);
    return val;
}
