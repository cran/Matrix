#ifndef MATRIX_METIS_UTILS_H
#define MATRIX_METIS_UTILS_H

#include <Rdefines.h>
#include "metis.h"

void ssc_metis_order(int n, const int Tp [], const int Ti [],
		     idxtype* perm, idxtype* iperm);

#endif

