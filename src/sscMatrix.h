#ifndef MATRIX_SSC_H
#define MATRIX_SSC_H

#include <Rdefines.h>
#include "Mutils.h"
#include "R_ldl.h"
#include "triplet_to_col.h"

SEXP sscMatrix_validate(SEXP x);
SEXP sscMatrix_chol(SEXP x, SEXP pivot);
SEXP sscMatrix_inverse_factor(SEXP A);
SEXP sscMatrix_matrix_solve(SEXP a, SEXP b);
SEXP ssc_transpose(SEXP x);
SEXP sscMatrix_to_triplet(SEXP x);
SEXP sscMatrix_ldl_symbolic(SEXP x, SEXP doPerm);
extern void ssc_metis_order(int n, const int Tp [], const int Ti [],
			    int Perm[], int iPerm[]);


#endif
