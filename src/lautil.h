//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

// General Utilities that should have been included in C++
//

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
