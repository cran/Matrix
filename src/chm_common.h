#ifndef CHM_COMMON_H
#define CHM_COMMON_H

#include "CHOLMOD/Include/cholmod.h"
#include "Mutils.h"

cholmod_common c;

cholmod_sparse *as_cholmod_sparse(SEXP x);
cholmod_triplet *as_cholmod_triplet(SEXP x);
cholmod_dense *as_cholmod_dense(SEXP x);

SEXP chm_sparse_to_SEXP(cholmod_sparse *a, int free);
SEXP chm_triplet_to_SEXP(cholmod_triplet *a, int free);
SEXP chm_dense_to_SEXP(cholmod_dense *a, int free);

#endif
