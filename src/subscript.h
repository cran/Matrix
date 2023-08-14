#ifndef MATRIX_SUBSCRIPT_H
#define MATRIX_SUBSCRIPT_H

#include "Mutils.h"

SEXP R_subscript_1ary    (SEXP x, SEXP i);
SEXP R_subscript_1ary_mat(SEXP x, SEXP i);
SEXP R_subscript_2ary    (SEXP x, SEXP i, SEXP j);

#endif
