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

#ifndef _LA_SYMM_FACT_DOUBLE_H_
#define _LA_SYMM_FACT_DOUBLE_H_

#include "lafnames.h"
#include LA_VECTOR_INT_H
#include LA_SYMM_MAT_DOUBLE_H


#include "lapack.h"

class LaSymmFactDouble
{
    LaSymmMatDouble             S_;
    LaVectorInt             pivot_;
    int                      info_;
    char                     uplo_;
    int                     size_;
    int                     gdim_;

public:

    // constructor

    inline LaSymmFactDouble();
    inline LaSymmFactDouble(int,int);
    inline LaSymmFactDouble(const LaSymmFactDouble &);
    inline ~LaSymmFactDouble();

    // extraction functions for components

    inline LaSymmMatDouble& S() { return S_; }
    inline LaVectorInt& pivot() { return pivot_; }
    inline int info() { return info_; }
    inline char uplo(){ return uplo_; }
    inline int size() { return size_; }
    inline int gdim() { return gdim_; }

    // operators

    inline LaSymmFactDouble ref(LaSymmFactDouble &);
    inline LaSymmFactDouble ref(LaSymmMatDouble &);
    inline LaSymmFactDouble& copy(const LaSymmFactDouble &);
    inline LaSymmFactDouble& copy(const LaSymmMatDouble &);

};



    // constructor/destructor functions

inline LaSymmFactDouble::LaSymmFactDouble():S_(),pivot_(),info_(0),uplo_('L')
{}


inline LaSymmFactDouble::LaSymmFactDouble(int n, int m):S_(n,m),pivot_(n*m),
                    info_(0),uplo_('L')
{}


inline LaSymmFactDouble::LaSymmFactDouble(const LaSymmFactDouble &F)
{
    S_.copy(F.S_);
    pivot_.copy(F.pivot_);
    info_ = F.info_;
    uplo_ = F.uplo_;
    size_ = F.size_;
    gdim_ = F.gdim_;
}

inline LaSymmFactDouble::~LaSymmFactDouble()
{}

    // operators


inline LaSymmFactDouble LaSymmFactDouble::ref(LaSymmFactDouble& F)
{
    S_.ref(F.S_);
    pivot_.ref(F.pivot_);
    info_ = F.info_;
    uplo_ = F.uplo_;
    size_ = F.size_;
    gdim_ = F.gdim_;
    
    return *this;
}

inline LaSymmFactDouble& LaSymmFactDouble::copy(const LaSymmFactDouble& F)
{
    S_.copy(F.S_);
    pivot_.copy(F.pivot_);
    info_ = F.info_;
    uplo_ = F.uplo_;
    size_ = F.size_;
    gdim_ = F.gdim_;
    
    return *this;
}

inline LaSymmFactDouble LaSymmFactDouble::ref(LaSymmMatDouble &G)
{
    S_.ref(G);
    info_ = 0;
    uplo_ = 'L';
    size_ = G.size(0);
    gdim_ = G.gdim(0);

    return *this;
}

inline LaSymmFactDouble& LaSymmFactDouble::copy(const LaSymmMatDouble &G)
{
    S_.copy(G);
    info_ = 0;
    uplo_ = 'L';
    size_ = G.size(0);
    gdim_ = G.gdim(0);

    return *this;
}

#if 0
inline void LaSymmMatFactorize(LaSymmMatDouble &A, LaSymmFactDouble &AF)
{
    char UPLO = 'L';
    int N = A.size(0), LDA = A.gdim(0), info = 0;
//    int M = DSYTRF;
//    int NB = F77NAME(get_nb)(&N,&M);

    int lwork = N*NB;
    double *WORK = new double[lwork];
    LaVectorInt piv(N);
    AF.pivot().copy(piv); // make copies of A and pivot information
    AF.copy(A);

    F77_CALL(dsytrf)('L', N, &(AF.S()(0,0)), A.gdim(0), &(AF.pivot()(0)),
		     work, lwork, &info);

    delete[] WORK;
}
#endif
inline void LaLinearSolve(LaSymmFactDouble &AF, LaGenMatDouble &X,
                           LaGenMatDouble &B)
{
//    char uplo = 'L';
//    int N = AF.size(), nrhs = X.size(1), lda = AF.gdim(),
//            ldb = B.size(0), info = 0;
    int info;

    X.inject(B);
    F77_CALL(dsytrs)('L', AF.size(), X.size(1), &(AF.S()(0,0)), AF.gdim(),
            &(AF.pivot()(0)), &X(0,0), B.size(0), info);

}

#endif 
// _LA_SYMM_FACT_DOUBLE_H_
