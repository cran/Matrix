#ifndef CHM_COMMON_H
#define CHM_COMMON_H

#include "UFconfig/UFconfig.h"
#include "CHOLMOD/Include/cholmod.h"
#include "Mutils.h"

typedef struct cholmod_common_struct  *CHM_CM ;
typedef struct cholmod_dense_struct   *CHM_DN ;
typedef struct cholmod_factor_struct  *CHM_FR ;
typedef struct cholmod_sparse_struct  *CHM_SP ;
typedef struct cholmod_triplet_struct *CHM_TR ;

extern cholmod_common c;

CHM_SP as_cholmod_sparse (CHM_SP ans, SEXP x);
CHM_TR as_cholmod_triplet(CHM_TR ans, SEXP x);
CHM_DN as_cholmod_dense  (CHM_DN ans, SEXP x);
CHM_DN as_cholmod_x_dense(CHM_DN ans, SEXP x);
CHM_DN numeric_as_chm_dense(CHM_DN ans, double *v, int nr, int nc);
CHM_FR as_cholmod_factor (CHM_FR ans, SEXP x);

#define AS_CHM_DN(x) as_cholmod_dense  ((CHM_DN)alloca(sizeof(cholmod_dense)), x )
#define AS_CHM_FR(x) as_cholmod_factor ((CHM_FR)alloca(sizeof(cholmod_factor)), x )
#define AS_CHM_SP(x) as_cholmod_sparse ((CHM_SP)alloca(sizeof(cholmod_sparse)), x )
#define AS_CHM_TR(x) as_cholmod_triplet((CHM_TR)alloca(sizeof(cholmod_triplet)), x )
#define N_AS_CHM_DN(x,nr,nc) M_numeric_as_chm_dense((CHM_DN)alloca(sizeof(cholmod_dense)), x , nr, nc )

int R_cholmod_start(CHM_CM Common);
void R_cholmod_error(int status, char *file, int line, char *message);

SEXP chm_factor_to_SEXP(CHM_FR f, int dofree);
SEXP chm_sparse_to_SEXP(CHM_SP a, int dofree, int uploT, int Rkind,
			const char *diag, SEXP dn);
SEXP chm_triplet_to_SEXP(CHM_TR a, int dofree, int uploT, int Rkind,
			 const char* diag, SEXP dn);
SEXP chm_dense_to_SEXP(CHM_DN a, int dofree, int Rkind, SEXP dn);
/* 		       int uploST, char *diag, SEXP dn); */
SEXP chm_dense_to_matrix(CHM_DN a, int dofree, SEXP dn);

SEXP CHMfactor_validate(SEXP obj);
SEXP CHMsimpl_validate(SEXP obj);
SEXP CHMsuper_validate(SEXP obj);

#endif
