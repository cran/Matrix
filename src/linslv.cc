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
// Modifications Copyright (C) 2000-2000, 2002 the R Development Core Team

#include <iostream>
#include "lapackd.h"
#include "lafnames.h"
#include LA_GEN_MAT_DOUBLE_H
#include LA_VECTOR_DOUBLE_H
#include LA_VECTOR_INT_H 
#include LA_SPD_MAT_DOUBLE_H
#include LA_SYMM_MAT_DOUBLE_H
#include LA_EXCEPTION_H

#include LA_SOLVE_DOUBLE_H
#include LA_UTIL_H

#ifdef length
#undef length
#endif

#ifdef append
#undef append
#endif

void LaLinearSolve(const LaGenMatDouble& A, LaGenMatDouble& X,
		   LaGenMatDouble& B)
{
    int M = A.size(0), N = A.size(1);

    if ( M == N ) 
        LaLULinearSolve(A,X,B);
    else 
        LaQRLinearSolve(A,X,B);
}   
    
void LaLinearSolve(const LaSpdMatDouble &A, LaGenMatDouble& X, 
		   LaGenMatDouble& B )
{
    LaCholLinearSolve(A, X, B );
}

#if 0
void LaLinearSolve(const LaSymmMatDouble &A, LaGenMatDouble& X, 
		   LaGenMatDouble& B )
{
    LaCholLinearSolve(A, X, B );
}
#endif

void LaLinearSolveIP(LaSpdMatDouble &A, LaGenMatDouble& X, LaGenMatDouble& B )
{
    LaCholLinearSolveIP(A, X, B );
}

void LaLinearSolveIP(LaGenMatDouble& A, LaGenMatDouble& X,
		     LaGenMatDouble& B)
{
    int M = A.size(0), N = A.size(1);

    if ( M == N ) 
        LaLULinearSolveIP(A,X,B);
    else 
        LaQRLinearSolveIP(A,X,B);
}   

void LaLULinearSolve(const  LaGenMatDouble& A, LaGenMatDouble& X, 
		     LaGenMatDouble& B )
{
    LaGenMatDouble A1(A);
    LaLULinearSolveIP(A1, X, B);
}

// General LU Solver
// 
//                        N x N               N x nrhs          N x nrhs
//
void LaLULinearSolveIP( LaGenMatDouble& A, LaGenMatDouble& X, LaGenMatDouble& B )
{
#ifndef HPPA
    const char fname[] = "LaLULinearSolveIP(LaGenMatDouble &A, &X, &B)";
#else
    char *fname = NULL;  // HP C++ does not support string initalization!
#endif

    // let's not worry about non-unit column strides for the moment
    if ( A.inc(0) != 1 || A.inc(1) != 1)
        throw(LaException(fname, "A is  non-contiguous."));
    
    if (!(X.size(0) == B.size(0) && X.size(1) == B.size(1)))
        throw(LaException(fname, "X and B are non-conformant."));
    
    X.inject(B);            // will throw exception if not conformant
    
    
    // in the future this can call the linear least square routines
    // to handle non-square matrices
    
    if (A.size(0) != A.size(1))
        throw(LaException(fname, "Square matrix expected.\n"));
    
    if (A.size(1) != X.size(0))
        throw(LaException(fname, "A and X are non-comformant."));
    
    int info;

    LaVectorInt ipiv(A.size(0));        
    F77_CALL(dgesv)(A.size(0), X.size(1), &A(0,0), A.gdim(0) * A.inc(0),
		    &ipiv(0), &X(0,0), X.inc(0) * X.gdim(0), info);
    
    if (info != 0)
        throw(LaException(fname, "Internal error in LAPACK: SGESV()"));
}

void LaQRLinearSolve(const LaGenMatDouble& A, LaGenMatDouble& X, 
		     LaGenMatDouble& B )
{
    LaGenMatDouble A1(A);
    LaQRLinearSolveIP(A1, X, B);
}


// General QR solver
//
//                          M x N              N x nrhs           M  x nrhs
//
void LaQRLinearSolveIP(LaGenMatDouble& A, LaGenMatDouble& X, LaGenMatDouble& B )
{
#ifndef HPPA
    const char fname[] = "LaQRLinearSolveIP(LaGenMatDouble &A, &X, &B)";
#else
    char *fname = NULL;  // HP C++ does not support string initalization!
#endif
    
    // let's not worry about non-unit column strides for the moment
    if ( A.inc(0) != 1 || A.inc(1) != 1)
        throw(LaException(fname, "A is  non-contiguous."));

    if (!(  A.size(0) == B.size(0) &&
            A.size(1) == X.size(0) &&
            X.size(1) == B.size(1) ))
        throw(LaException(fname, "input matrices are non-conformant."));

    int info;
    int M = A.size(0), N = A.size(1);

    int nrhs = X.size(1);

#define MAX3(A,B,C) (A>B ? (A>C ? A : C) : (B>C ? B : C))
#define MIN(A,B) (A < B ? A : B )

    int nb = LaEnvBlockSize("DGELSV", A);
    int W = MIN(M,N) + nb * MAX3(M, N, nrhs);
    
    LaVectorDouble work(W);         
    // typically is A non-square, so we need to create tmp X becuase is
    //  X is N x nrhs, while B is M x nrhs.  We need to make copies of
    //  these so that the routine won't corrupt data around X and B
    
    if (M != N) {
        int Mm =  MAX3(M,N,1);
        LaGenMatDouble Xtmp(Mm, nrhs);
	
#undef MIN
#undef MAX3
	
        Xtmp(LaIndex(0,M-1),LaIndex(0,nrhs-1)) = B(LaIndex(0,M-1), LaIndex(0,nrhs-1));
        F77_CALL(dgels)('N', A.size(0), A.size(1), X.size(1), &A(0,0),
			A.inc(0) * A.gdim(0), &Xtmp(0,0), 
			X.inc(0) * X.gdim(0), &work(0), W, info);
        X.copy(Xtmp(LaIndex(0,N-1), LaIndex(0,nrhs-1)));
    } else {
        F77_CALL(dgels) ('N', A.size(0), A.size(1), X.size(1), &A(0,0),
			A.inc(0) * A.gdim(0), &X(0,0), 
			X.inc(0) * X.gdim(0), &work(0), W, info);
    }
    // this shouldn't really happen.
    //
    if (info < 0)
        throw(LaException(fname, "Internal error in LAPACK: SGELS()"));

}

void  LaCholLinearSolve( const LaSpdMatDouble& A, LaGenMatDouble& X,
			 LaGenMatDouble& B )
{
    LaSpdMatDouble A1(A);
    LaCholLinearSolveIP(A1, X, B);
}

#if 0
void  LaCholLinearSolve( const LaSymmMatDouble& A, LaGenMatDouble& X,
			 LaGenMatDouble& B )
{
    LaSymmMatDouble A1(A);
    LaCholLinearSolveIP(A1, X, B);
}
#endif

// A is NxN, X is MxN and B is MxN
//
void LaCholLinearSolveIP( LaSpdMatDouble& A, LaGenMatDouble& X, 
			  LaGenMatDouble& B )
{
#ifndef HPPA
    const char fname[] = "LaCholLinearSolveIP(LaSpdMatDouble &A, &X, &B)";
#else
    char *fname = NULL;  // HP C++ does not support string initalization!
#endif
    
    // let's not worry about non-unit column strides for the moment
    if ( A.inc(0) != 1 || A.inc(1) != 1)
        throw(LaException(fname, "A is  non-contiguous."));
    
    if (!(X.size(0) == B.size(0) && X.size(1) == B.size(1)))
        throw(LaException(fname, "X and B are non-conformant."));

    X.inject(B);            // will throw exception if not conformant


    // in the future this can call the linear least square routines
    // to handle non-square matrices
    
    if (A.size(0) != A.size(1))
        throw(LaException(fname, "Square matrix expected.\n"));
    
    if (A.size(1) != X.size(0))
        throw(LaException(fname, "A and X are non-comformant."));
    
    int info;
    
    F77_CALL(dposv)('L', A.size(0), X.size(0), &A(0,0),
		    A.inc(0) * A.gdim(0), &X(0,0), X.inc(0) * X.gdim(0),
		    info);
    
    // this shouldn't really happen.
    //
    if (info < 0)
        throw(LaException(fname, "Internal error in LAPACK: dposv()"));
    if (info > 0)
        throw(LaException(fname, "A is not symmetric-positive-definite."));
}
