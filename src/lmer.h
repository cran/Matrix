#ifndef MATRIX_LMEUTILS_H
#define MATRIX_LMEUTILS_H

#include "Mutils.h"
#include "triplet_to_col.h"
#include "dgBCMatrix.h"
#include "bCrosstab.h"
#include <R_ext/Lapack.h>
#include <R_ext/Constants.h>

SEXP lmer_validate(SEXP x);
SEXP lmer_update_mm(SEXP x, SEXP mmats);
SEXP lmer_create(SEXP flist, SEXP mmats);
SEXP lmer_inflate(SEXP x);
SEXP lmer_initial(SEXP x);
SEXP lmer_factor(SEXP x);
SEXP lmer_invert(SEXP x);
SEXP lmer_sigma(SEXP x, SEXP REML);
SEXP lmer_coef(SEXP x, SEXP Unc);
SEXP lmer_coefGets(SEXP x, SEXP coef, SEXP Unc);
SEXP lmer_fixef(SEXP x);
SEXP lmer_ranef(SEXP x);
SEXP lmer_ECMEsteps(SEXP x, SEXP nsteps, SEXP REMLp, SEXP Verbp);
SEXP lmer_fitted(SEXP x, SEXP mmats, SEXP useRf);
SEXP lmer_gradient(SEXP x, SEXP REMLp, SEXP Uncp);
SEXP lmer_variances(SEXP x);
SEXP lmer_Crosstab(SEXP flist);
SEXP lmer_firstDer(SEXP x, SEXP val);

#endif
