#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP R_dense_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP diag);
SEXP R_dense_as_kind(SEXP from, SEXP kind);
SEXP R_dense_as_matrix(SEXP from, SEXP ndense);
SEXP R_geMatrix_as_matrix(SEXP from, SEXP ndense);
SEXP R_dense_as_vector(SEXP from, SEXP ndense);
SEXP R_geMatrix_as_vector(SEXP from, SEXP ndense);
SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2);

SEXP lsq_dense_Chol(SEXP X, SEXP y);
SEXP lsq_dense_QR(SEXP X, SEXP y);
SEXP lapack_qr(SEXP Xin, SEXP tl);

/* MJ: no longer needed ... prefer (un)?packedMatrix_(symm|skew)part() */
#if 0
SEXP ddense_symmpart(SEXP x);
SEXP ddense_skewpart(SEXP x);
#endif
/* MJ */

/* MJ: no longer needed ... prefer R_dense_as_sparse() above */
#if 0
SEXP dense_to_Csparse(SEXP x);
#endif /* MJ */

/* MJ: no longer needed ... prefer (un)?packedMatrix_force_symmetric() */
#if 0
SEXP dense_to_symmetric(SEXP x, SEXP uplo, SEXP symm_test);
#endif /* MJ */

/* MJ: no longer needed ... prefer R_dense_band() above */
#if 0
SEXP dense_band(SEXP x, SEXP k1, SEXP k2);
#endif /* MJ */

/* TODO: compare with macros in ./Mutils.h */

#define VALID_DDENSE							\
"dgeMatrix", "dtrMatrix", "dsyMatrix", "dtpMatrix", "dspMatrix"

#define VALID_LDENSE							\
"lgeMatrix", "ltrMatrix", "lsyMatrix", "ltpMatrix", "lspMatrix"

#define VALID_NDENSE							\
"ngeMatrix", "ntrMatrix", "nsyMatrix", "ntpMatrix", "nspMatrix"

#endif
