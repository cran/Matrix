#ifndef MATRIX_LMEUTILS_H
#define MATRIX_LMEUTILS_H

#include "Mutils.h"
#include "triplet_to_col.h"
#include "dgBCMatrix.h"
#include "dgCMatrix.h"
#include "Metis_utils.h"
#include "R_ldl.h"
#include "dsCMatrix.h"
#include "dtCMatrix.h"
#include "lgCMatrix.h"
#include "lCholCMatrix.h"
#include "Rmath.h"
#include <R_ext/Lapack.h>
#include <R_ext/Constants.h>

SEXP Matrix_rWishart(SEXP ns, SEXP df, SEXP scal);
SEXP lmer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp);
SEXP lmer_validate(SEXP x);
SEXP lmer_update_mm(SEXP x, SEXP mmats);
SEXP lmer_create(SEXP flist, SEXP mmats, SEXP method);
SEXP lmer_inflate(SEXP x);
SEXP lmer_initial(SEXP x);
SEXP lmer_factor(SEXP x);
SEXP lmer_invert(SEXP x);
SEXP lmer_sigma(SEXP x, SEXP REML);
SEXP lmer_coef(SEXP x, SEXP pType);
SEXP lmer_coefGets(SEXP x, SEXP coef, SEXP pType);
SEXP lmer_fixef(SEXP x);
SEXP lmer_ranef(SEXP x);
SEXP lmer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp);
SEXP lmer_fitted(SEXP x, SEXP mmats, SEXP useRf);
SEXP lmer_gradient(SEXP x, SEXP pType);
SEXP lmer_variances(SEXP x);
SEXP lmer_Crosstab(SEXP flist);
SEXP lmer_firstDer(SEXP x, SEXP val);
SEXP lmer_secondDer(SEXP x);
SEXP glmer_MCMCsamp(SEXP GSpt, SEXP b, SEXP fixedp, SEXP varcp,
		    SEXP savebp, SEXP nsampp);
SEXP glmer_PQL(SEXP GSp);
SEXP glmer_ranef_update(SEXP GSp, SEXP fixed, SEXP varc, SEXP b);
SEXP glmer_devAGQ(SEXP pars, SEXP GSp, SEXP nAGQp);
SEXP glmer_finalize(SEXP GSpt);
SEXP glmer_fixed_update(SEXP GSp, SEXP b, SEXP fixed);
SEXP glmer_init(SEXP rho);

#endif
