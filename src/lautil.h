// -*- c++ -*-
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

#ifndef _LA_UTIL_H_
#define _LA_UTIL_H_

// only callable from C-Lapack due to added ftnlen parameters by f2c; 
extern "C" {
    int F77_NAME(ilaenv)(const int& i, const char* name, const char* opts, 
			 const int& n1, const int& n2, const int& n3,
			 const int& n4);
}

#ifdef _LA_GEN_MAT_DOUBLE_H_
#ifdef _LA_VECTOR_LONG_INT_H_
    void LaSwap(LaGenMatDouble &A, LaVectorLongInt &ipiv);
#endif
    int LaEnvBlockSize(const char *fname, const LaGenMatDouble &A);
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
    int LaEnvBlockSize(const char *fname, const LaSymmMatDouble &A);
#endif
double Mach_eps_double();
#endif // _LA_UTIL_H_
