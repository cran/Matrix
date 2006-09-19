#ifndef HBMM_H
#define HBMM_H

#include <Rdefines.h>
#include "Mutils.h"

SEXP Matrix_writeHarwellBoeing(SEXP obj, SEXP filename, SEXP type);
SEXP Matrix_writeMatrixMarket(SEXP obj, SEXP filename, SEXP type);

#endif
