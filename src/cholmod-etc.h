#ifndef MATRIX_CHOLMOD_ETC_H
#define MATRIX_CHOLMOD_ETC_H

#include <Rinternals.h>
#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"

extern cholmod_common c ;
extern cholmod_common cl;

cholmod_factor *M2CHF(SEXP, int);
cholmod_sparse *M2CHS(SEXP, int);
cholmod_dense  *M2CHD(SEXP, int);

SEXP CHF2M(cholmod_factor *, int);
SEXP CHS2M(cholmod_sparse *, int, char);
SEXP CHD2M(cholmod_dense  *, int, char);

#endif /* MATRIX_CHOLMOD_ETC_H */
