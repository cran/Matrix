//
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.
//
// Modifications Copyright (C) 2000-2000 the R Development Core Team

#include "lafnames.h"
#include LA_SPD_MAT_DOUBLE_H

SEXP LaSpdMatDouble::asSEXP() const
{
    int n = size(0);
    SEXP val = allocMatrix(REALSXP, n, n);
    F77_CALL(dlacpy)('L', n, n, &(*this)(0,0), gdim(0),
		     REAL(val), n);
    for (int i = 1; i < n; i++) // symmetrize the result
	for (int j = 0; j < i; j++)
	    REAL(val)[i * n + j] = REAL(val)[j * n + i];
    SEXP classes = allocVector(STRSXP, 3);
    STRING(classes)[0] = mkChar("PositiveDefinite");
    STRING(classes)[1] = mkChar("Hermitian");
    STRING(classes)[2] = mkChar("Matrix");
    setAttrib(val, R_ClassSymbol, classes);
    return val;
}
