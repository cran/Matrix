#include "flame.h"
#include "R_ext/Lapack.h"
#include "FLAME/FLAME.h"

FLA_Obj *R_to_FLA_copy(SEXP Ain)
{
    FLA_Obj *val = (FLA_Obj *) R_alloc((long) 1, sizeof(FLA_Obj));
    int *Adims, m, n;

    if (!(isReal(Ain) && isMatrix(Ain)))
	error("A must be a numeric matrix");
    Adims = INTEGER(coerceVector(getAttrib(Ain, R_DimSymbol), INTSXP));
    m = Adims[0]; n = Adims[1];
    val->datatype = FLA_DOUBLE; val->m = m; val->n = n; val->ldim = m;
    val->buffer = (void *) R_alloc(m * n, sizeof(double));
    Memcpy((double *)val->buffer, REAL(Ain), m * n);
    return val;
}

FLA_Obj *R_to_FLA_inPlace(SEXP Ain)
{
    FLA_Obj *val = (FLA_Obj *) R_alloc((long) 1, sizeof(FLA_Obj));
    int *Adims, m, n;

    if (!(isReal(Ain) && isMatrix(Ain)))
	error("A must be a numeric matrix");
    Adims = INTEGER(coerceVector(getAttrib(Ain, R_DimSymbol), INTSXP));
    m = Adims[0]; n = Adims[1];
    val->datatype = FLA_DOUBLE;
    val->m = m; val->n = n; val->ldim = m;
    val->buffer = (void *) REAL(Ain);
    return val;
}

SEXP R_FLA_Init()
{
    FLA_Init();
    return ScalarLogical(1);
}

SEXP R_FLA_Finalize()
{
    FLA_Finalize();
    return ScalarLogical(1);
}

int FLA_Abort( char *msg, int line, char * fname)
{
    error(msg, line, fname);
    return 0;
}


SEXP lsq_Chol_flame(SEXP Xin, SEXP yin)
{
    FLA_Obj *X = R_to_FLA_inPlace(Xin), *y = R_to_FLA_inPlace(yin),
	XpX, Xpy;
    int n = FLA_Obj_length(*X), p = FLA_Obj_width(*X), k = FLA_Obj_width(*y);
    SEXP ans = PROTECT(allocMatrix(REALSXP, p, k));
    
    FLA_Obj_create(FLA_DOUBLE, p, p, &XpX);
    FLA_Syrk(FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, ONE, *X, ZERO, XpX);
    FLA_Obj_create_without_buffer(FLA_DOUBLE, p, k, &Xpy);
    FLA_Obj_attach_buffer(REAL(ans), p, &Xpy);
    FLA_Gemm(FLA_TRANSPOSE, FLA_NO_TRANSPOSE, ONE, *X, *y, ZERO, Xpy);
    FLA_Chol_eager(FLA_LOWER_TRIANGULAR, XpX, RFLAME_CHOL_NB);
    FLA_Trsv(FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, 
	     XpX, Xpy);
    FLA_Trsv(FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, 
	     XpX, Xpy);
    FLA_Obj_free(&XpX);
    UNPROTECT(1);
    return ans;
}

SEXP lsq_QR_flame(SEXP Xin, SEXP yin)
{
    FLA_Obj *X = R_to_FLA_copy(Xin), *y = R_to_FLA_copy(yin), S;
    int n = FLA_Obj_length(*X), p = FLA_Obj_width(*X), k = FLA_Obj_width(*y);
    SEXP ans = PROTECT(allocMatrix(REALSXP, p, k));
    
/*     FLA_Obj_create(X->datatype, p, RFLAME_QR_NB, &S); */
/*     FLA_QR_right_rec(*X, S, RFLAME_QR_NB); */
    FLA_Obj_create(FLA_DOUBLE, p, 1, &S);
    FLA_QR_right_unb(*X, S);
    FLA_Obj_free(&S); UNPROTECT(1);
    return ans;
}

