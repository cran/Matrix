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
    SEXP val = PROTECT(LaSymmMatDouble::asSEXP());
    SEXP classes = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(classes, 0, mkChar("PositiveDefinite"));
    SET_STRING_ELT(classes, 1, mkChar("Hermitian"));
    SET_STRING_ELT(classes, 2, mkChar("Matrix"));
    setAttrib(val, R_ClassSymbol, classes);
    UNPROTECT(2);
    return val;
}
