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

#ifndef _BLAS_PP_H_
#define _BLAS_PP_H_

#include "laexcp.h"
#include "blas1++.h"
#include "blas2++.h"
#include "blas3++.h"
#include <complex.h>		// to declare abs
//inline double abs(double x) {return (x < 0.0) ? -x : x;}

				// Vector/Vector operators
#ifdef _LA_VECTOR_DOUBLE_H_

inline LaVectorDouble operator*(const LaVectorDouble &x, double a)
{
    int N = x.size();
    LaVectorDouble t(N);
    
    for (int i=0; i<N; i++)
    {
        t(i) = a * x(i);
    }
    
    return t;
}

inline LaVectorDouble operator*(double a, const LaVectorDouble &x)
{
    return operator*(x,a);
}

inline double operator*(const LaVectorDouble &dx, 
			const LaVectorDouble &dy)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    return F77_CALL(ddot)(dx.size(), &dx(0), dx.inc(), &dy(0), dy.inc());
}

inline LaVectorDouble operator+(const LaVectorDouble &dx, 
                                const LaVectorDouble &dy)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    
    LaVectorDouble tmp((int) dx.size());
    tmp = dy;
    
    F77_CALL(daxpy)(dx.size(), 1.0, &dx(0), dx.inc(), &tmp(0), dy.inc());
    return tmp;
}

inline LaVectorDouble operator-(const LaVectorDouble &dx, 
                                const LaVectorDouble &dy)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
//    int incx = dx.inc(), incy = dy.inc(), n = dx.size();
    
    LaVectorDouble tmp(dx.size());
    tmp = dx;
    
    F77_CALL(daxpy)(dx.size(), -1.0, &dy(0), dx.inc(), &tmp(0), dy.inc());
    return tmp;
}
				// Matrix/Vector operators
#ifdef _LA_GEN_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaGenMatDouble &A, 
                                const LaVectorDouble &dx)
{
    LaVectorDouble dy(A.size(0));
    dy = 0.0;
    
    F77_CALL(dgemv)('N', A.size(0), A.size(1), 1.0, &A(0,0), A.gdim(0),
		    &dx(0), dx.inc(), 0.0, &dy(0), dy.inc()); 
    return dy; 
    
}
#endif

#ifdef _LA_BAND_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaBandMatDouble &A, 
                                const LaVectorDouble &dx)
{
//    char trans = 'N';
//    double alpha = 1.0, beta = 0.0;
//    int M = A.size(0), N = A.size(1), lda = A.gdim(0),
//        kl = A.subdiags(), ku = A.superdiags(); 
    
    LaVectorDouble dy(A.size(1));
//    int incx = dx.inc(), incy = dy.inc();
    
    F77_CALL(dgbmv)('N', A.size(0), A.size(1), A.subdiags(),
		    A.superdiags(), 1.0, &A(0,0), A.gdim(0),
		    &dx(0), dx.inc(), 0.0, &dy(0), dy.inc());
    return dy;
}
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaSymmMatDouble &A, 
                                const LaVectorDouble &dx)
{
//    char uplo = 'L';
//    double alpha = 1.0, beta = 0.0;
//    int N = A.size(1), lda = A.gdim(0);
    
    LaVectorDouble dy(A.size(1));
//    int incx = dx.inc(), incy = dy.inc();
    
    F77_CALL(dsymv)('L', A.size(1), 1.0, &A(0,0), A.gdim(0), &dx(0),
		    dx.inc(), 0.0, &dy(0), dy.inc());
    return dy;
}
#endif

#ifdef _LA_SYMM_BAND_MAT_DOUBLE_H_ 
inline LaVectorDouble operator*(const LaSymmBandMatDouble &A, 
				const LaVectorDouble &dx) 
{
//    char uplo = 'L';
//    double alpha = 1.0, beta = 0.0;
//    int N = A.size(1), lda = A.gdim(0), k = A.subdiags();
    
    LaVectorDouble dy(A.size(1));
//    int incx = dx.inc(), incy = dy.inc();
    
    F77_CALL(dsbmv)('L', A.size(1), A.subdiags(), 1.0, &A(0,0),
		    A.gdim(0), &dx(0), dx.inc(), 0.0, &dy(0), dy.inc());
    return dy;
}
#endif

#ifdef _LA_SPD_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaSpdMatDouble &AP, 
                                const LaVectorDouble &dx)
{
//    char uplo = 'L';
//    double alpha = 1.0, beta = 0.0;
//    int N = AP.size(1), incx = dx.inc(); 
//    int lda = AP.gdim(0);
    
    LaVectorDouble dy(AP.size(0));
//    int incy = dy.inc();
    
    F77_CALL(dsymv)('L', AP.size(1), 1.0, &AP(0,0), AP.gdim(0),
		    &dx(0), dx.inc(), 0.0, &dy(0), dy.inc());
    return dy;
}
#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaLowerTriangMatDouble &A, 
                                const LaVectorDouble &dx)
{
//    char uplo = 'L', trans = 'N', diag = 'N';
//    int N = A.size(1), lda = A.gdim(0),
//        incx = dx.inc();
    
    LaVectorDouble dy(dx);
    
    F77_CALL(dtrmv)('L', 'N', 'N', A.size(1), &A(0,0), A.gdim(0),
		    &dy(0), dy.inc());
    return dy;
}
#endif

#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaUpperTriangMatDouble &A, 
                                const LaVectorDouble &dx)
{
//    char uplo = 'U', trans = 'N', diag = 'N';
//    int N = A.size(1), lda = A.gdim(0),
//        incx = dx.inc();
    
    LaVectorDouble dy(dx);
    
    F77_CALL(dtrmv)('U', 'N', 'N', A.size(1), &A(0,0), A.gdim(0),
		    &dy(0), dy.inc());
    return dy;
}
#endif

#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaUnitLowerTriangMatDouble &A,
                                const LaVectorDouble &dx)
{
//    char uplo = 'L', trans = 'N', diag = 'U';
//    int N = A.size(1), lda = A.gdim(0),
//        incx = dx.inc();
    
    LaVectorDouble dy(dx);
    
    F77_CALL(dtrmv)('L', 'N', 'U', A.size(1), &A(0,0),
		    A.gdim(0), &dy(0), dy.inc());
    return dy;
}
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaUnitUpperTriangMatDouble &A,
                                const LaVectorDouble &dx)
{
//    char uplo = 'U', trans = 'N', diag = 'U';
//    int N = A.size(1), lda = A.gdim(0),
//        incx = dx.inc();

    LaVectorDouble dy(dx);
    
    F77_CALL(dtrmv)('U', 'N', 'U', A.size(1), &A(0,0), A.gdim(0),
		    &dy(0), dy.inc());
    return dy;
}
#endif
				// Matrix/Matrix operators

inline LaGenMatDouble operator*(const LaGenMatDouble &A, 
                                const LaGenMatDouble &B)
{
//    char t = 'N';
//    int m = A.size(0), k = A.size(1), n = B.size(1);
//    int lda = A.gdim(0), ldb = B.gdim(0);
    
    LaGenMatDouble C(A.size(0), B.size(1));

//    int ldc = A.gdim(0);
    
    C = 0.0;
    
    F77_CALL(dgemm)('N', 'N', A.size(0), B.size(1), A.size(1),
		    1.0, &A(0,0), A.gdim(0), &B(0,0), B.gdim(0),
		    1.0, &C(0,0), A.gdim(0));
    return C;
}

#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaUnitLowerTriangMatDouble &A,
                                const LaGenMatDouble &B)
{
//    char side = 'L', uplo = 'L', transa = 'N', diag = 'U';
//    double alpha = 1.0;
//    int m = B.size(0), n = B.size(1),
//	lda = A.gdim(0), ldb = B.gdim(0);
    
    LaGenMatDouble C(B);
    
    F77_CALL(dtrmm)('L', 'L', 'N', 'U', B.size(0), B.size(1), 1.0,
		    &A(0,0), A.gdim(0), &C(0,0), B.gdim(0));
    return C;
}
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaUnitUpperTriangMatDouble &A,
                                const LaGenMatDouble &B)
{
//    char side = 'L', uplo = 'U', transa = 'N', diag = 'U';
//    double alpha = 1.0;
//    int m = B.size(0), n = B.size(1),
//	lda = A.gdim(0), ldb = B.gdim(0);
    
    LaGenMatDouble C(B);
    
    F77_CALL(dtrmm)('L', 'U', 'N', 'U', B.size(0), B.size(1), 1.0,
		    &A(0,0), A.gdim(0), &C(0,0), B.gdim(0));
    
    return C;
}
#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaLowerTriangMatDouble &A,
                                const LaGenMatDouble &B)
{
//    char side = 'L', uplo = 'L', transa = 'N', diag = 'N';
//    double alpha = 1.0;
//    int m = B.size(0), n = B.size(1),
//	lda = A.gdim(0), ldb = B.gdim(0);
    
    LaGenMatDouble C(B);
    
    F77_CALL(dtrmm)('L', 'L', 'N', 'N', B.size(0), B.size(1), 1.0,
		    &A(0,0), A.gdim(0), &C(0,0), B.gdim(0));
    
    return C;
}
#endif

#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaUpperTriangMatDouble &A,
                                const LaGenMatDouble &B)
{
//    char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
//    double alpha = 1.0;
//    int m = B.size(0), n = B.size(1),
//	lda = A.gdim(0), ldb = B.gdim(0);
    
    LaGenMatDouble C(B);
    
    F77_CALL(dtrmm)('L', 'U', 'N', 'N', B.size(0), B.size(1), 1.0,
		    &A(0,0), A.gdim(0), &C(0,0), B.gdim(0));
    
    return C;
}
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaSymmMatDouble &A, 
                                const LaGenMatDouble &B)
{
//    char side = 'L', uplo = 'L';
//    double alpha = 1.0, beta = 1.0;
    LaGenMatDouble C(B.size(1),A.size(1));
//    int m = C.size(0), n = C.size(1), lda = A.gdim(0),
//	ldb = B.gdim(0), ldc = C.gdim(0);
    
    F77_CALL(dsymm)('L', 'L', C.size(0), C.size(1), 1.0,
		    &A(0,0), A.gdim(0), &B(0,0), B.gdim(0), 1.0,
		    &C(0,0), C.gdim(0));
    
    return C;
}
#endif

#ifdef _LA_SYMM_TRIDIAG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaSymmTridiagMatDouble& A, 
                                const LaVectorDouble& X)
{
    int M = A.size();
    int N = X.size();
    LaVectorDouble R(M);
    
    R(0) = ((A.diag(0)(0) * X(0)) + (A.diag(-1)(0) * X(1)));
    
    for (int i = 1; i < M-2; i++)
    {
        R(i) = ((A.diag(-1)(i-1) * X(i-1)) +
                (A.diag(0)(i) * X(i)) +
                (A.diag(-1)(i) * X(i+1)));
    }
    
    R(M-1) = ((A.diag(0)(M-1) * X(N-1)) + (A.diag(-1)(M-2) * X(N-2)));
    
    return R;
}
#endif

#ifdef  _LA_TRIDIAG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaTridiagMatDouble& A, 
                                const LaVectorDouble& X)
{
    int M = A.size();
    int N = X.size();
    LaVectorDouble R(M);
    
    R(0) = ((A.diag(0)(0) * X(0)) + (A.diag(-1)(0) * X(1)));
    
    for (int i = 1; i < M-2; i++)
    {
        R(i) = ((A.diag(-1)(i-1) * X(i-1)) +
                (A.diag(0)(i) * X(i)) +
                (A.diag(1)(i) * X(i+1)));
    }
    
    R(M-1) = ((A.diag(0)(M-1) * X(N-1)) + (A.diag(1)(M-2) * X(N-2)));
    
    return R;
}
#endif

inline LaGenMatDouble operator-(const LaGenMatDouble &A, 
                                const LaGenMatDouble &B)
{
#ifndef HPPA
    const char fname[] = "operator-(A,B)";
#else
    char *fname = NULL;
#endif

    int M = A.size(0);
    int N = A.size(1);
    
    if (M != B.size(0) || N != B.size(1))
    {
       throw(LaException(fname, "matrices non-conformant."));
    }
    
    LaGenMatDouble C(M,N);
    
    // slow mode
    // we'll hook the BLAS in later
    
    for (int i = 0;  i < M; i++)
        for(int j = 0; j < N; j++)
            C(i,j) = A(i,j) - B(i,j);
    return C;
}

inline LaGenMatDouble operator+(const LaGenMatDouble &A, 
                                const LaGenMatDouble &B)
{
#ifndef HPPA
    const char fname[] = "operator+(A,B)";
#else
    char *fname = NULL;
#endif

    int M = A.size(0);
    int N = A.size(1);
    
    if (M != B.size(0) || N != B.size(1))
    {
        throw(LaException(fname, "matrices non-conformant."));
    }
    
    LaGenMatDouble C(M,N);
    
    // slow mode
    // we'll hook the BLAS in later
    
    for (int i = 0;  i < M; i++)
        for(int j = 0; j < N; j++)
            C(i,j) = A(i,j) + B(i,j);
    
    return C;
}
				// Matrix/Vector Norms
inline double Norm_Inf(const LaVectorDouble &x)
{   
    int index = Blas_Index_Max(x);
    return abs(x(index));
}

inline double Norm_Inf(const LaGenMatDouble &A)
{
    int M = A.size(0);
    int index;
    
    // max row-sum
    
    LaVectorDouble R(M);
    
    for (int i = 0; i < M; i++)
        R(i) = Blas_Norm1( A( i, LaIndex() ));
    
    index = Blas_Index_Max(R);
    return R(index);
}

#ifdef _LA_BAND_MAT_DOUBLE_H_
inline double Norm_Inf(const LaBandMatDouble &A)
{
//    int kl = A.subdiags(), ku = A.superdiags(); 
    int N=A.size(1);
    int M=N;
    
    // slow version
    
    LaVectorDouble R(M);
    for (int i = 0; i < M; i++)
    {
        R(i) = 0.0;
        for (int j = 0; j < N; j++)
            R(i) += abs(A(i,j));
    }
    
    double max = R(0);
    
    // report back largest row sum
    for (int i = 1; i < M; i++)
        if (R(i) > max) max = R(i);
    return max;
}
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
inline double Norm_Inf(const LaSymmMatDouble &S)
{
    int N = S.size(0); // square matrix
    
    // slow version
    
    LaVectorDouble R(N);
    
    for (int i = 0; i < N; i++)
    {
        R(i) = 0.0;
        for (int j = 0; j < N; j++)
            R(i) += abs(S(i,j));
    }
    
    double max = R(0);
    
    // report back largest row sum
    for (int i = 1; i < N; i++)
        if (R(i) > max) max=R(i);
    
    return max;
}
#endif

#ifdef _LA_SPD_MAT_DOUBLE_H_
inline double Norm_Inf(const LaSpdMatDouble &S)
{
    int N = S.size(0);		//SPD matrices are square
    
				// slow version
    LaVectorDouble R(N);
    for (int i = 0; i < N; i++) {
        R(i) = 0.0;
        for (int j = 0; j < N; j++) { R(i) += abs(S(i,j)); }
    }
    
    double max = R(0);
				// report back largest row sum
    for (int i = 1; i < N; i++) { if (R(i) > max) max = R(i); }
    return max;
}
#endif

#ifdef _LA_SYMM_TRIDIAG_MAT_DOUBLE_H_
inline double Norm_Inf(const LaSymmTridiagMatDouble &S)
{
    int N = S.size();		// S is square
    LaVectorDouble R(N);

    R(0) = abs(S(0,0)) + abs(S(0,1));
    
    for (int i=1; i<N-1; i++)
    {
        R(i) = abs(S(i,i-1)) + abs(S(i,i)) + abs(S(i,i+1));
    }
    
    R(N-1) = abs(S(N-1,N-2)) + abs(S(N-1,N-1));
    
    return Norm_Inf(R);
}
#endif

#ifdef _LA_TRIDIAG_MAT_DOUBLE_H_
inline double Norm_Inf(const LaTridiagMatDouble &T)
{
    int N = T.size();		// T is square
    LaVectorDouble R(N);
    
    R(0) = abs(T(0,0)) + abs(T(0,1));
    
    for (int i=1; i<N-1; i++)
    {
        R(i) = abs(T(i,i-1)) + abs(T(i,i)) + abs(T(i,i+1));
    }

    R(N-1) = abs(T(N-1,N-2)) + abs(T(N-1,N-1));

    return Norm_Inf(R);
}
#endif

#endif // LA_VECTOR_DOUBLE_H
#endif // _BLAS_PP_H_
