#ifndef MATRIX_SSCCROSSTAB_H
#define MATRIX_SSCCROSSTAB_H

#include "Mutils.h"
#include "ldl.h"

SEXP sscCrosstab(SEXP flist, SEXP upper);
SEXP sscCrosstab_L_LI_sizes(SEXP ctab, SEXP permexp);
SEXP sscCrosstab_groupedPerm(SEXP ctab);

#endif
