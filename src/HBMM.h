#ifndef HBMM_H
#define HBMM_H

#include <Rdefines.h>
#include "Mutils.h"

#if 0
SEXP Matrix_writeHarwellBoeing(SEXP obj, SEXP filename, SEXP type);
#endif

SEXP Matrix_writeMatrixMarket(SEXP obj, SEXP filename, SEXP type);

#endif
