#ifndef MATRIX_CHMFACTOR_H
#define MATRIX_CHMFACTOR_H

#include "chm_common.h"
#include "Mutils.h"

double chm_factor_ldetL2(CHM_FR f);
CHM_FR chm_factor_update(CHM_FR f, CHM_SP A, double fac);

#endif
