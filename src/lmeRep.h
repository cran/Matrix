#ifndef MATRIX_LMEUTILS_H
#define MATRIX_LMEUTILS_H

#include "Mutils.h"
#include "triplet_to_col.h"
#include <R_ext/Lapack.h>
#include <R_ext/Constants.h>

SEXP lmeRep_validate(SEXP x);
SEXP lmeRep_create(SEXP facs, SEXP ncv);
SEXP lmeRep_update_mm(SEXP x, SEXP facs, SEXP mmats);
SEXP lmeRep_initial(SEXP x);
SEXP lmeRep_factor(SEXP x);
SEXP lmeRep_invert(SEXP x);
SEXP lmeRep_sigma(SEXP x, SEXP REML);
SEXP lmeRep_coef(SEXP x, SEXP Unc);
SEXP lmeRep_coefGets(SEXP x, SEXP coef, SEXP Unc);
SEXP lmeRep_fixef(SEXP x);
SEXP lmeRep_ranef(SEXP x);
SEXP lmeRep_ECMEsteps(SEXP x, SEXP nsteps, SEXP REMLp, SEXP Verbp);
SEXP lmeRep_gradient(SEXP x, SEXP REMLp, SEXP Uncp);


#endif
