#include "lafnames.h"
#include LA_ROW_VECTOR_DOUBLE_H

LaRowVectorDouble& LaRowVectorDouble::ref(SEXP x)
{				// create a reference to the data
    int n = LENGTH(x);
    LaRowVectorDouble tmp(REAL(coerceVector(x, REALSXP)), n);
    return ref(tmp);
}
