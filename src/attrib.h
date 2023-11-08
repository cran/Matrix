#ifndef MATRIX_ATTRIB_H
#define MATRIX_ATTRIB_H

#include <Rinternals.h>

SEXP R_DimNames_is_symmetric(SEXP);
SEXP R_symDN(SEXP);
SEXP R_revDN(SEXP);
SEXP R_set_factor(SEXP, SEXP, SEXP, SEXP);

#endif /* MATRIX_ATTRIB_H */
