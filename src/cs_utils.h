#ifndef CS_UTILS_H
#define CS_UTILS_H

#include "cs.h"
#include "Mutils.h"

cs *Matrix_as_cs(SEXP x);
SEXP Matrix_cs_to_SEXP(cs *A, const char *cl, int dofree);

#endif
