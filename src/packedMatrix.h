#ifndef MATRIX_PACKEDMATRIX_H
#define MATRIX_PACKEDMATRIX_H

#include "Mutils.h"

SEXP packedMatrix_unpack(SEXP from, SEXP strict);
SEXP packedMatrix_force_symmetric(SEXP from, SEXP uplo_to);

SEXP packedMatrix_is_triangular(SEXP obj, SEXP upper);
SEXP packedMatrix_is_symmetric(SEXP obj, SEXP checkDN);
SEXP packedMatrix_is_diagonal(SEXP obj);

SEXP packedMatrix_transpose(SEXP from);
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms);
SEXP packedMatrix_diag_set(SEXP obj, SEXP val);

SEXP packedMatrix_symmpart(SEXP from);
SEXP packedMatrix_skewpart(SEXP from);

#endif /* MATRIX_PACKEDMATRIX_H */
