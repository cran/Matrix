#ifndef MATRIX_LMER_H
#define MATRIX_LMER_H

#include "Mutils.h"
#include "triplet_to_col.h"
#include "dgCMatrix.h"
#include "dpoMatrix.h"
#include "R_ldl.h"
#include "dsCMatrix.h"
#include "dtCMatrix.h"
#include "lgCMatrix.h"
#include "lCholCMatrix.h"
#include "Rmath.h"
#include "chm_common.h"
#include <R_ext/Lapack.h>
#include <R_ext/Constants.h>
#include <R_ext/Utils.h>

SEXP Matrix_rWishart(SEXP ns, SEXP df, SEXP scal);
/* SEXP glmer_MCMCsamp(SEXP GSpt, SEXP b, SEXP fixedp, SEXP varcp, */
/* 		    SEXP savebp, SEXP nsampp); */
SEXP glmer_PQL(SEXP GSp);
/* SEXP glmer_bhat(SEXP pars, SEXP GSp); */
/* SEXP glmer_devAGQ(SEXP pars, SEXP GSp, SEXP nAGQp); */
SEXP glmer_devLaplace(SEXP pars, SEXP GSp);
SEXP glmer_finalize(SEXP GSpt);
/* SEXP glmer_fixed_update(SEXP GSp, SEXP b, SEXP fixed); */
SEXP glmer_init(SEXP rho);
/* SEXP glmer_ranef_update(SEXP GSp, SEXP fixed, SEXP varc, SEXP b); */
SEXP mer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp);
SEXP mer_Hessian(SEXP x);
SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp);
SEXP mer_coef(SEXP x, SEXP pType);
SEXP mer_coefGets(SEXP x, SEXP coef, SEXP pType);
SEXP mer_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP method,
		 SEXP nc, SEXP cnames, SEXP useS, SEXP call,
		 SEXP family);
SEXP mer_dtCMatrix(SEXP x);
SEXP mer_dtCMatrix_inv(SEXP x);
SEXP mer_factor(SEXP x);
SEXP mer_fitted(SEXP x);
SEXP mer_fixef(SEXP x);
SEXP mer_gradComp(SEXP x);
SEXP mer_gradient(SEXP x, SEXP pType);
SEXP mer_initial(SEXP x);
SEXP mer_postVar(SEXP x);
SEXP mer_ranef(SEXP x);
SEXP mer_secondary(SEXP x);
SEXP mer_sigma(SEXP x, SEXP REML);
SEXP mer_simulate(SEXP x, SEXP nsimP);
SEXP mer_update_ZXy(SEXP x);
SEXP mer_update_y(SEXP x, SEXP ynew);
SEXP mer_validate(SEXP x);
SEXP Zt_create(SEXP fl, SEXP Ztl);

#endif
