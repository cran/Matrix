#ifndef MATRIX_UNPACKEDMATRIX_H
#define MATRIX_UNPACKEDMATRIX_H

#include "Mutils.h"

SEXP unpackedMatrix_force_symmetric(SEXP from, SEXP uplo_to);

SEXP unpackedMatrix_is_triangular(SEXP obj, SEXP upper);
SEXP matrix_is_triangular(SEXP obj, SEXP upper);

SEXP unpackedMatrix_is_symmetric(SEXP obj, SEXP checkDN);
SEXP matrix_is_symmetric(SEXP obj, SEXP checkDN);

SEXP unpackedMatrix_is_diagonal(SEXP obj);
SEXP matrix_is_diagonal(SEXP obj);

SEXP unpackedMatrix_transpose(SEXP from);
SEXP unpackedMatrix_diag_get(SEXP obj, SEXP nms);
SEXP unpackedMatrix_diag_set(SEXP obj, SEXP val);

SEXP unpackedMatrix_symmpart(SEXP from);
SEXP matrix_symmpart(SEXP from);

SEXP unpackedMatrix_skewpart(SEXP from);
SEXP matrix_skewpart(SEXP from);

#endif /* MATRIX_UNPACKEDMATRIX_H */
