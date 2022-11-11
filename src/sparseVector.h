#ifndef MATRIX_SPVECTOR_H
#define MATRIX_SPVECTOR_H

#include "Mutils.h"

#define SPV_SUB(_KIND_, _CTYPE_)				\
_CTYPE_ _KIND_ ## sparseVector_sub(int64_t i, int nnz_v,	\
				   double* v_i, _CTYPE_ *v_x,	\
				   int64_t len_v)
SPV_SUB(n, int);
SPV_SUB(l, int);
SPV_SUB(i, int);
SPV_SUB(d, double);
SPV_SUB(z, Rcomplex);
#undef SPV_SUB

#endif /* MATRIX_SPVECTOR_H */
