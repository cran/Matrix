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


#ifndef _LA_GEN_FACT_DOUBLE_H
#define _LA_GEN_FACT_DOUBLE_H

#include "lafnames.h"
#include LA_VECTOR_INT_H
#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_UPPER_TRIANG_MAT_DOUBLE_H
#include LA_UTIL_H
#include "lapack.h"
#include "factor.h"

class LaLUFactorDouble : public Factor
{
    LaUnitLowerTriangMatDouble  L_;
    LaUpperTriangMatDouble      U_;
    LaVectorInt             pivot_;
    bool                 singular_;

public:
				// constructor
    LaLUFactorDouble()
	: L_(), U_(), pivot_() { singular_ = true; };
    inline LaLUFactorDouble(const LaGenMatDouble&);
    inline LaLUFactorDouble(const LaLUFactorDouble&);

    virtual ~LaLUFactorDouble() { };

				// extractor methods for components
    LaUnitLowerTriangMatDouble& L()
	{ return L_; };
    LaUpperTriangMatDouble& U()
	{ return U_; };
    LaVectorInt& pivot()
	{ return pivot_; };
    bool isSingular()
	{ return singular_; };

				// linear equation solvers
    inline LaGenMatDouble& solve() const;// inverse
    inline LaMatDouble& solve(LaMatDouble& B) const; // in-place solution
    inline LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const;
				// operators
    inline LaLUFactorDouble& ref(const LaLUFactorDouble&);
    inline LaLUFactorDouble& ref(const LaGenMatDouble&);
};



// constructor/destructor functions

inline LaLUFactorDouble::LaLUFactorDouble(const LaGenMatDouble& A) :
    pivot_(min(A.size(0), A.size(1)))
{
    LaGenMatDouble A1;
    A1.copy(A);
    ref(A1);
}

inline LaLUFactorDouble::LaLUFactorDouble(const LaLUFactorDouble& F)
{
  L_.ref(F.L_);
  U_.ref(F.U_);
  pivot_.ref(F.pivot_);
  singular_ = F.singular_;
}

// operators
inline LaLUFactorDouble& LaLUFactorDouble::ref(const LaGenMatDouble& A)
{
    if (A.size(0) != A.size(1))
	throw(LaException("LaLUFactorDouble::ref(const LaGenMatDouble&)",
			  "non-square input matrix"));
    if(A.inc(0) != 1 || A.inc(1) != 1)	
	throw(LaException("LaLUFactorDouble::ref(const LaGenMatDouble&)",
			  "input matrix has non unit increment"));
    pivot_.resize(A.size(0));
    L_.ref(A);
    U_.ref(A);
    int info;
    F77_CALL(dgetrf)(A.size(0), A.size(0), &A(0, 0), A.gdim(0),
		     &pivot_(0), info);
    if (info < 0)
	throw(LaException("LaLUFactorDouble::ref(const LaGenMatDouble&)",
			  "illegal input"));
    singular_ = info > 0;
    return *this;
}

inline LaLUFactorDouble& LaLUFactorDouble::ref(const LaLUFactorDouble& F)
{

    L_.ref(F.L_);
    U_.ref(F.U_);
    pivot_.ref(F.pivot_);
    singular_ = F.singular_;
    
    return *this;
}

inline LaGenMatDouble& LaLUFactorDouble::solve() const
{
    if (singular_)
	throw(LaException("singular matrix"));
    int info;
    LaGenMatDouble& ans = *(new LaGenMatDouble());
    ans.copy(U_);
    int lwork = U_.size(0) *
	F77_NAME(ilaenv)(1, "DGETRI", "", U_.size(0), -1, -1, -1);
    VectorDouble work(lwork);
    F77_CALL(dgetri)(U_.size(0), &ans(0,0), U_.gdim(0), &pivot_(0),
		     &work(0), lwork, info);
    if (info < 0)
	throw(LaException("LaLUFactorDouble::ref(const LaGenMatDouble&)",
			  "illegal input"));
    return ans;
}

inline LaMatDouble& LaLUFactorDouble::solve(LaMatDouble& B) const
{
    if (singular_)
	throw(LaException("singular matrix"));

    dynamic_cast<LaGenMatDouble&>(B);

    int info;
    F77_CALL(dgetrs)('N', L_.size(0), B.size(1), &U_(0,0),
		     L_.gdim(0), &pivot_(0), &B(0,0), B.size(0), info);
    return B;
}

inline LaMatDouble& LaLUFactorDouble::solve(LaMatDouble& X, const LaMatDouble& B ) const
{
//    dynamic_cast<LaGenMatDouble&>(B);
    dynamic_cast<LaGenMatDouble&>(X);
    X.inject(B);
    return solve(X);
}

#if 0
inline void LaGenMatFactorize(LaGenMatDouble &GM, LaLUFactorDouble &GF)
{
    integer m = GM.size(0), n = GM.size(1), lda = GM.gdim(0);
    integer info=0;

    F77NAME(dgetrf)(&m, &n, &GM(0,0), &lda, &(GF.pivot()(0)), &info);
}

inline void LaGenMatFactorizeUnblocked(LaGenMatDouble &A, LaLUFactorDouble &F)
{
    integer m = A.size(0), n=A.size(1), lda = A.gdim(0);
    integer info=0;

    F77NAME(dgetf2)(&m, &n, &A(0,0), &lda, &(F.pivot()(0)), &info);
}
#endif

#endif
