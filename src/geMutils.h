#include <Rinternals.h>
#include <Rdefines.h>

SEXP Matrix_DimSym, Matrix_xSym, Matrix_uploSym, Matrix_diagSym,
    Matrix_pSym, Matrix_iSym, Matrix_zSym;
char norm_type(char *typstr);
char rcond_type(char *typstr);
double get_double_by_name(SEXP obj, char *nm);
SEXP set_double_by_name(SEXP obj, double val, char *nm);
SEXP as_det_obj(double val, int log, int sign);
SEXP get_factorization(SEXP obj, char *nm);
SEXP set_factorization(SEXP obj, SEXP val, char *nm);
SEXP cscMatrix_set_Dim(SEXP x, int nrow);
