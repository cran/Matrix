#ifndef MATRIX_H
#define MATRIX_H
#include <Rdefines.h>
#include <Rconfig.h>
#include "cholmod.h"

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

cholmod_sparse* M_as_cholmod_sparse(SEXP x);
cholmod_dense* M_as_cholmod_dense(SEXP x);
cholmod_dense* M_numeric_as_chm_dense(double *v, int n);
cholmod_factor* M_as_cholmod_factor(SEXP x);
SEXP M_chm_factor_to_SEXP(cholmod_factor *f, int dofree);
SEXP M_dpoMatrix_chol(SEXP x);

#endif /* MATRIX_H */
