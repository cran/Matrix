#ifndef MATRIX_CHMFACTOR_H
#define MATRIX_CHMFACTOR_H

#include "chm_common.h"
#include "Mutils.h"

double chm_factor_ldetL2(CHM_FR f);
CHM_FR chm_factor_update(CHM_FR f, CHM_SP A, double fac);

/* MJ: unused or no longer needed, with replacement
   in ./factorizations.c or in ../R/factorizations.R
*/
#if 0
SEXP CHMfactor_to_sparse(SEXP x);
SEXP CHMfactor_solve(SEXP a, SEXP b, SEXP type);
SEXP CHMfactor_spsolve(SEXP a, SEXP b, SEXP type);
SEXP CHMfactor_updown(SEXP update, SEXP C, SEXP L);
SEXP CHMfactor_ldetL2(SEXP x);
SEXP CHMfactor_update(SEXP object, SEXP parent, SEXP mult);
SEXP destructive_CHM_update(SEXP object, SEXP parent, SEXP mult);
SEXP CHMfactor_ldetL2up(SEXP x, SEXP parent, SEXP mult);
#endif /* MJ */

#endif
