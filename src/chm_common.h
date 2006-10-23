#ifndef CHM_COMMON_H
#define CHM_COMMON_H

#include "UFconfig/UFconfig.h"
#include "CHOLMOD/Include/cholmod.h"
#include "Mutils.h"

extern cholmod_common c;

cholmod_sparse *as_cholmod_sparse(SEXP x);
cholmod_triplet *as_cholmod_triplet(SEXP x);
cholmod_dense *as_cholmod_dense(SEXP x);
cholmod_dense *numeric_as_chm_dense(double *v, int n);
cholmod_factor *as_cholmod_factor(SEXP x);

int R_cholmod_start(cholmod_common *Common);
void R_cholmod_error(int status, char *file, int line, char *message);

SEXP chm_factor_to_SEXP(cholmod_factor *f, int dofree);
SEXP chm_sparse_to_SEXP(cholmod_sparse *a, int dofree,
			int uploT, int Rkind, char *diag, SEXP dn);
SEXP chm_triplet_to_SEXP(cholmod_triplet *a, int dofree,
 			 int uploT, int Rkind, char* diag, SEXP dn);
SEXP chm_dense_to_SEXP(cholmod_dense *a, int dofree, int Rkind);
/* 		       int uploST, char *diag, SEXP dn); */
SEXP chm_dense_to_matrix(cholmod_dense *a, int dofree, SEXP dn);

SEXP CHMfactor_validate(SEXP obj);
SEXP CHMsimpl_validate(SEXP obj);
SEXP CHMsuper_validate(SEXP obj);

#endif
