#ifndef MATRIX_CHM_COMMON_H
#define MATRIX_CHM_COMMON_H

#include "cholmod-etc.h"

cholmod_factor  *   sexp_as_cholmod_factor (cholmod_factor  *, SEXP);
cholmod_sparse  *   sexp_as_cholmod_sparse (cholmod_sparse  *, SEXP,
                                           Rboolean, Rboolean);
cholmod_triplet *   sexp_as_cholmod_triplet(cholmod_triplet *, SEXP,
                                            Rboolean);
cholmod_dense   *   sexp_as_cholmod_dense  (cholmod_dense   *, SEXP);
cholmod_dense   *numeric_as_cholmod_dense  (cholmod_dense   *, double *,
                                            int, int);

SEXP cholmod_factor_as_sexp (cholmod_factor  *, int);
SEXP cholmod_sparse_as_sexp (cholmod_sparse  *, int,
                             int, int, const char *, SEXP);
SEXP cholmod_triplet_as_sexp(cholmod_triplet *, int,
                             int, int, const char *, SEXP);
SEXP cholmod_dense_as_sexp  (cholmod_dense   *, int);

double          cholmod_factor_ldetA (cholmod_factor *);
cholmod_factor *cholmod_factor_update(cholmod_factor *, cholmod_sparse *,
                                      double);

int  R_cholmod_start (cholmod_common *);
int  R_cholmod_finish(cholmod_common *);

SEXP R_cholmod_common_envini(SEXP);
void R_cholmod_common_envset(void);
void R_cholmod_common_envget(void);

#endif /* MATRIX_CHM_COMMON_H */
