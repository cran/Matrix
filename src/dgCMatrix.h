#ifndef MATRIX_CSC_H
#define MATRIX_CSC_H

#include <Rdefines.h>
#include "Mutils.h"
#include "R_ldl.h"
#include "triplet_to_col.h"

SEXP csc_crossprod(SEXP x);
SEXP csc_tcrossprod(SEXP x);
SEXP csc_matrix_crossprod(SEXP x, SEXP y, SEXP classed);
SEXP dgCMatrix_validate(SEXP x);
SEXP compressed_to_dgTMatrix(SEXP x, SEXP colP);
SEXP csc_to_matrix(SEXP x);
SEXP csc_to_dgeMatrix(SEXP x);
SEXP matrix_to_csc(SEXP A);
SEXP dgTMatrix_to_csc(SEXP dgTMatrix);
SEXP csc_getDiag(SEXP x);
SEXP csc_transpose(SEXP x);
SEXP csc_matrix_mm(SEXP a, SEXP b);
SEXP csc_col_permute(SEXP x, SEXP perm);

#endif
