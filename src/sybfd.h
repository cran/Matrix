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

#ifndef _LA_SYMM_BAND_FACT_DOUBLE_H_
#define _LA_SYMM_BAND_FACT_DOUBLE_H_

#include LA_VECTOR_INT_H
#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_SYMM_BAND_MAT_DOUBLE_H
#include "lapack.h"

class LaSymmBandFactDouble
{
    LaSymmBandMatDouble         S_;
    int                      info_;
    char                     uplo_;

public:

    // constructor

    inline LaSymmBandFactDouble();
    inline LaSymmBandFactDouble(int,int);
    inline LaSymmBandFactDouble(LaSymmBandFactDouble &);
    inline ~LaSymmBandFactDouble();

    // extraction functions for components

    inline LaSymmBandMatDouble& S();
    inline int& info();
    inline char& uplo();

    // operators

    inline LaSymmBandFactDouble ref(LaSymmBandFactDouble &);
    inline LaSymmBandFactDouble ref(LaSymmBandMatDouble &);
    inline LaSymmBandFactDouble& copy(LaSymmBandFactDouble &);
    inline LaSymmBandFactDouble& copy(LaSymmBandMatDouble &);

};



    // constructor/destructor functions

inline LaSymmBandFactDouble::LaSymmBandFactDouble():S_(),
            info_(0),uplo_('L')
{}


inline LaSymmBandFactDouble::LaSymmBandFactDouble(int n, int m):S_(n,m),
            info_(0),uplo_('L')
{}


inline LaSymmBandFactDouble::LaSymmBandFactDouble(LaSymmBandFactDouble &F)
{
  S_.ref(F.S_);
  info_ = F.info_;
  uplo_ = F.uplo_;
}

inline LaSymmBandFactDouble::~LaSymmBandFactDouble()
{}

    // member functions

inline LaSymmBandMatDouble& LaSymmBandFactDouble::S()
{
    return S_;
}

inline int& LaSymmBandFactDouble::info()
{
    return info_;
}

inline char& LaSymmBandFactDouble::uplo()
{
    return uplo_;
}

    
    // operators


inline LaSymmBandFactDouble LaSymmBandFactDouble::ref(LaSymmBandFactDouble& F)
{
    S_.ref(F.S_);
    info_ = F.info_;
    uplo_ = F.uplo_;
    
    return *this;
}

inline LaSymmBandFactDouble LaSymmBandFactDouble::ref(LaSymmBandMatDouble &G)
{
  S_.ref(G);
  info_ = 0;
  uplo_ = 'L';

  return *this;
}

inline LaSymmBandFactDouble& LaSymmBandFactDouble::copy(LaSymmBandFactDouble& F)
{
    S_.copy(F.S_);
    info_ = F.info_;
    uplo_ = F.uplo_;
    
    return *this;
}

inline LaSymmBandFactDouble& LaSymmBandFactDouble::copy(LaSymmBandMatDouble &G)
{
  S_.copy(G);
  info_ = 0;
  uplo_ = 'L';

  return *this;
}

inline void LaSymmBandMatFactorize(const LaSymmBandMatDouble &A,
                                    LaSymmBandFactDouble& AF)
{
    char uplo = 'L';
    integer n = A.size(0), kd = A.subdiags(), lda = A.gdim(0),
            info = 0;

    AF.S().copy(A);
    F77NAME(dpbtrf)(&uplo, &n, &kd, &(AF.S()(0,0)), &lda, &info);
}


inline void LaLinearSolve(LaSymmBandFactDouble &AF, LaGenMatDouble &X,
                            LaGenMatDouble &B)
{
    char uplo = 'L';
    integer n = AF.S().size(1), kd = AF.S().subdiags(), ldaf = AF.S().gdim(0),
            info = 0, nrhs = X.size(1), ldb = B.size(0);

    X.inject(B);
    F77NAME(dpbtrs)(&uplo, &n, &kd, &nrhs, &(AF.S()(0,0)), &ldaf, 
                    &X(0,0), &ldb, &info);
}

#endif 
// _LA_SYMM_BAND_FACT_DOUBLE_H_
