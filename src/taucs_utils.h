#ifndef MATRIX_TAUCS_UTILS_H
#define MATRIX_TAUCS_UTILS_H

#include "Mutils.h"
#include "taucs/taucs.h"

taucs_ccs_matrix* csc_taucs_ptr(SEXP A, int flags);
SEXP mat_from_taucs(taucs_ccs_matrix *tm);

#endif
