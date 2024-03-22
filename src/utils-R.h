#ifndef MATRIX_UTILS_R_H
#define MATRIX_UTILS_R_H

#include <Rinternals.h>

SEXP R_Matrix_version(void);

SEXP R_index_triangle(SEXP, SEXP, SEXP, SEXP);
SEXP R_index_diagonal(SEXP, SEXP, SEXP);

SEXP R_nnz(SEXP, SEXP, SEXP);

SEXP R_all0(SEXP);
SEXP R_any0(SEXP);

SEXP Mmatrix(SEXP);

SEXP compressed_non_0_ij(SEXP, SEXP);
SEXP Matrix_expand_pointers(SEXP);
SEXP m_encodeInd (SEXP, SEXP, SEXP, SEXP);
SEXP m_encodeInd2(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP Matrix_rle_i(SEXP, SEXP);
SEXP Matrix_rle_d(SEXP, SEXP);

#endif /* MATRIX_UTILS_R_H */
