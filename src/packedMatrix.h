#ifndef MATRIX_PAMATRIX_H
#define MATRIX_PAMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

void fast_symmetric_DimNames(SEXP dn, SEXP *vec, SEXP *nm);
SEXP packedMatrix_t(SEXP obj);
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms);
SEXP packedMatrix_diag_set(SEXP obj, SEXP val);
SEXP packedMatrix_sub1(SEXP obj, SEXP index);
SEXP packedMatrix_sub1_mat(SEXP obj, SEXP index);
SEXP packedMatrix_sub2(SEXP obj, SEXP index1, SEXP index2, SEXP drop);

#endif
