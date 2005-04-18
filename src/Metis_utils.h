#ifndef MATRIX_METIS_UTILS_H
#define MATRIX_METIS_UTILS_H

#include "Mutils.h"
void ssc_metis_order(int n, const int Tp [], const int Ti [],
		     int perm[], int iperm[]);

#endif

