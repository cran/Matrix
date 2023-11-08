#ifndef MATRIX_OBJECTS_H
#define MATRIX_OBJECTS_H

#include <Rinternals.h>

SEXP R_Matrix_nonvirtual(SEXP, SEXP);
SEXP R_Matrix_kind (SEXP);
SEXP R_Matrix_shape(SEXP);
SEXP R_Matrix_repr (SEXP);

#endif /* MATRIX_OBJECTS_H */
