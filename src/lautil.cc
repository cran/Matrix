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
#include LA_GEN_MAT_DOUBLE_H
#include LA_VECTOR_INT_H
#include LA_SYMM_MAT_DOUBLE_H

#include LA_UTIL_H
#include "lapackd.h"

int LaEnvBlockSize(const char *fname, const LaGenMatDouble &A)
{
    char *opts = "U";
    
    int one = 1;
    int M = A.size(0);
    int N = A.size(1);
    int junk = -1;

    return F77_CALL(ilaenv)(one, fname, opts, M, N, junk, junk);
}

int LaEnvBlockSize(const char *fname, const LaSymmMatDouble &A)
{
    char *opts = "U";
    int one = 1;
    int M = A.size(0);
    int N = A.size(1);
    int junk = -1;

    return F77_CALL(ilaenv)(one, fname, opts, M, N, junk, junk);
}


double Mach_eps_double()
{
    return F77_CALL(dlamch)('e');
}

void LaSwap(LaGenMatDouble &A, LaVectorInt &ipiv)
{
    F77_CALL(dlaswp)(A.size(1), &A(0,0), A.gdim(0), ipiv.start(),
		     ipiv.end(), &ipiv(0), ipiv.inc());
}


