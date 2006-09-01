#ifndef CS_UTILS_H
#define CS_UTILS_H

#include "cs.h"
#include "Mutils.h"

cs *Matrix_as_cs(SEXP x);
css *Matrix_as_css(SEXP x);
csn *Matrix_as_csn(SEXP x);
SEXP Matrix_cs_to_SEXP(cs *A, char *cl, int dofree);
SEXP Matrix_css_to_SEXP(css *S, char *cl, int dofree, int m, int n);
SEXP Matrix_csn_to_SEXP(csn *N, char *cl, int dofree);

#endif
