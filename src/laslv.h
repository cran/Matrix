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

#ifdef _LA_GEN_MAT_DOUBLE_H_
void LaLinearSolve( const LaGenMatDouble& A, LaGenMatDouble& X, 
    LaGenMatDouble& B );

void LaLinearSolveIP( LaGenMatDouble& A, LaGenMatDouble& X, LaGenMatDouble& B );


void LaLULinearSolve( const LaGenMatDouble& A, LaGenMatDouble& X, 
    LaGenMatDouble& B );
void LaLULinearSolveIP( LaGenMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );

void LaQRLinearSolve( const LaGenMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );
void LaQRLinearSolveIP( LaGenMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );



#ifdef _LA_SPD_MAT_DOUBLE_H_

void LaLinearSolve( const LaSpdMatDouble& A, LaGenMatDouble& X, 
    LaGenMatDouble& B );
void LaLinearSolveIP( LaSpdMatDouble& A, LaGenMatDouble& X, LaGenMatDouble& B );

void LaCholLinearSolve( const LaSpdMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );
void LaCholLinearSolveIP( LaSpdMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_

void LaLinearSolve( const LaSymmMatDouble& A, LaGenMatDouble& X, 
    LaGenMatDouble& B );
void LaLinearSolveIP( LaSymmMatDouble& A, LaGenMatDouble& X, LaGenMatDouble& B );

void LaCholLinearSolve( const LaSymmMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );
void LaCholLinearSolveIP( LaSymmMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );

// Eigenvalue problems 

void LaEigSolve(const LaSymmMatDouble &S, LaVectorDouble &eigvals);
void LaEigSolve(const LaSymmMatDouble &S, LaVectorDouble &eigvals, 
    LaGenMatDouble &eigvec);
void LaEigSolveIP(LaSymmMatDouble &S, LaVectorDouble &eigvals);
void LaEigSolveVecIP(LaSymmMatDouble &S, LaVectorDouble &eigvals);
#endif



#endif
