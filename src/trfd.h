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

#ifndef _LA_TRIDIAG_FACT_DOUBLE_H_
#define _LA_TRIDIAG_FACT_DOUBLE_H_

#include LA_VECTOR_INT_H
#include LA_TRIDIAG_MAT_DOUBLE_H

#include "lapack.h"

class LaTridiagFactDouble
{
    LaTridiagMatDouble T_;
    LaVectorInt pivot_;
    int size_;
public:

    // constructor

    LaTridiagFactDouble();
    LaTridiagFactDouble(int);
    LaTridiagFactDouble(LaTridiagFactDouble &);
    ~LaTridiagFactDouble();

    LaTridiagMatDouble& T() { return T_; }
    LaVectorInt& pivot() { return pivot_; }
    int size() { return size_; }
    LaVectorDouble diag(int);

    // operators

    LaTridiagFactDouble& ref(LaTridiagMatDouble &);
    LaTridiagFactDouble& ref(LaTridiagFactDouble &);
    LaTridiagFactDouble& copy(const LaTridiagMatDouble &);
    LaTridiagFactDouble& copy(const LaTridiagFactDouble &);

};



    // constructor/destructor functions

inline LaTridiagFactDouble::LaTridiagFactDouble():T_(),pivot_(),size_(0)
{}


inline LaTridiagFactDouble::LaTridiagFactDouble(int N):T_(N),pivot_(N),size_(N)
{}


inline LaTridiagFactDouble::LaTridiagFactDouble(LaTridiagFactDouble &F)
{
  T_.copy(F.T_);
  pivot_.copy(F.pivot_);
  size_ = F.size_;
}

inline LaTridiagFactDouble::~LaTridiagFactDouble()
{}

    // member functions

inline LaVectorDouble LaTridiagFactDouble::diag(int k)
{
    return T_.diag(k);
}
    
    // operators


inline LaTridiagFactDouble& LaTridiagFactDouble::ref(LaTridiagFactDouble& F)
{
    T_.ref(F.T_);
    pivot_.ref(F.pivot_);
    size_ = F.size_;
    
    return *this;
}

inline LaTridiagFactDouble& LaTridiagFactDouble::ref(LaTridiagMatDouble& A)
{
    T_.ref(A);

    return *this;
}

inline LaTridiagFactDouble& LaTridiagFactDouble::copy(const LaTridiagFactDouble& F)
{
    T_.copy(F.T_);
    pivot_.copy(F.pivot_);
    size_ = F.size_;
    
    return *this;
}

inline LaTridiagFactDouble& LaTridiagFactDouble::copy(const LaTridiagMatDouble& A)
{
    T_.copy(A);

    return *this;
}

inline void LaTridiagMatFactorize(LaTridiagMatDouble &A,
                                 LaTridiagFactDouble &AF)
{
    int info = 0;
    AF.copy(A);
    
//    double *DL = &AF.diag(-1)(0), *D = &AF.diag(0)(0),
//         *DU = &AF.diag(1)(0), *DU2 = &AF.diag(2)(0);
    F77_CALL(dgttrf)(A.size(), &AF.diag(-1)(0), &AF.diag(0)(0),
		     &AF.diag(1)(0), &AF.diag(2)(0),
		     &(AF.pivot()(0)), info);
// Shouldn't this do something with info?
}


inline void LaLinearSolve(LaTridiagFactDouble &AF, LaGenMatDouble &X,
                        LaGenMatDouble &B)
{
//    char trans = 'N';
//    int N = AF.size(), nrhs = X.size(1), ldb = B.size(0), info = 0;
//    double *DL = &AF.diag(-1)(0), *D = &AF.diag(0)(0),
//         *DU =  &AF.diag(1)(0), *DU2 = &AF.diag(2)(0);
    int info;
    X.inject(B);
    F77_CALL(dgttrs)('N', AF.size(), X.size(1),  &AF.diag(-1)(0),
		    &AF.diag(0)(0), &AF.diag(1)(0), &AF.diag(2)(0),
		    &(AF.pivot()(0)), &X(0,0), B.size(0), info);
// Shouldn't this do something with info?

}

#endif 
// _LA_TRIDIAG_FACT_DOUBLE_H_
