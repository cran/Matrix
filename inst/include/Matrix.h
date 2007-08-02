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

CHM_SP M_as_cholmod_sparse(CHM_SP ans, SEXP x);
CHM_TR M_as_cholmod_triplet(CHM_TR ans, SEXP x);
CHM_DN M_as_cholmod_dense(CHM_DN ans, SEXP x);
CHM_DN M_numeric_as_chm_dense(CHM_DN ans, double *v, int nr, int nc);
CHM_FR M_as_cholmod_factor(CHM_FR ans, SEXP x);

#define AS_CHM_DN(x) M_as_cholmod_dense((CHM_DN)alloca(sizeof(cholmod_dense)), x )
#define AS_CHM_FR(x) M_as_cholmod_factor((CHM_FR)alloca(sizeof(cholmod_factor)), x )
#define AS_CHM_SP(x) M_as_cholmod_sparse((CHM_SP)alloca(sizeof(cholmod_sparse)), x )
#define AS_CHM_TR(x) M_as_cholmod_triplet((CHM_TR)alloca(sizeof(cholmod_triplet)), x )
#define N_AS_CHM_DN(x,nr,nc) M_numeric_as_chm_dense((CHM_DN)alloca(sizeof(cholmod_dense)), x , nr, nc )

SEXP M_chm_factor_to_SEXP(CHM_FR f, int dofree);
SEXP M_chm_sparse_to_SEXP(CHM_SP a, int dofree, int uploT, int Rkind,
			  char *diag, SEXP dn);
SEXP M_chm_triplet_to_SEXP(CHM_TR a, int dofree, int uploT, int Rkind,
			   const char* diag, SEXP dn);

SEXP M_dpoMatrix_chol(SEXP x);

#endif /* MATRIX_H */
