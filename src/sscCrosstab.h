#ifndef MATRIX_SSCCROSSTAB_H
#define MATRIX_SSCCROSSTAB_H

#include "Mutils.h"
#include "ldl.h"

SEXP sscCrosstab(SEXP flist, SEXP upper);
extern void ssc_metis_order(int n, const int Tp [], const int Ti [],
			    int Perm[], int iPerm[]);
SEXP sscCrosstab_groupedPerm(SEXP ctab);
SEXP sscCrosstab_project2(SEXP ctab);

#endif
