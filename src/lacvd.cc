#include "lafnames.h"
#include LA_COL_VECTOR_DOUBLE_H

LaColVectorDouble& LaColVectorDouble::ref(SEXP x)
{				// create a reference to the data
    int n = LENGTH(x);
    LaColVectorDouble tmp(REAL(coerceVector(x, REALSXP)), n);
    return ref(tmp);
}
