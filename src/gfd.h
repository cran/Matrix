// -*- c++ -*-
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//  Modifications Copyright (C) 2000-2000 the R Development Core Team
//


#ifndef _LA_GEN_FACT_DOUBLE_H
#define _LA_GEN_FACT_DOUBLE_H

#include "lafnames.h"
#include LA_VECTOR_INT_H
#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_UPPER_TRIANG_MAT_DOUBLE_H

#include "lapack.h"
#include "solvable.h"

class LaLUFactorDouble : public solvable
{
    LaUnitLowerTriangMatDouble  L_;
    LaUpperTriangMatDouble      U_;
    LaVectorInt             pivot_;
    int                      info_;
    int                 transpose_;

public:
				// constructor
    LaLUFactorDouble()
	: L_(), U_(), pivot_() { info_ = 0; transpose_ = 0; };
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
    int& info()
	{ return info_; };
    int& transpose()
	{ return transpose_; };

				// linear equation solvers
    inline LaMatrix& solve() const;// inverse
    inline LaMatrix& solve(LaMatrix& B) const; // in-place solution
    inline LaMatrix& solve(LaMatrix& X, const LaMatrix& B) const;
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
  info_ = F.info_;
  transpose_ = F.transpose_;
}

// operators
inline LaLUFactorDouble& LaLUFactorDouble::ref(const LaGenMatDouble& A)
{
    assert(A.size(0) == A.size(1) && A.inc(0) == 1 && A.inc(1) == 1);
    LaVectorInt pivot(A.size(0));
    pivot_.ref(pivot);
    L_.ref(A);
    U_.ref(A);
    transpose_ = 0;
    F77_CALL(dgetrf)(A.size(0), A.size(0), &A(0, 0), A.gdim(0),
		     &pivot(0), info_);
    return *this;
}

inline LaLUFactorDouble& LaLUFactorDouble::ref(const LaLUFactorDouble& F)
{

    L_.ref(F.L_);
    U_.ref(F.U_);
    pivot_.ref(F.pivot_);
    info_ = F.info_;
    transpose_ = F.transpose_;
    
    return *this;
}

inline LaMatrix& LaLUFactorDouble::solve() const
{
    return *new LaGenMatDouble(0,0);
}

inline LaMatrix& LaLUFactorDouble::solve(LaMatrix& B) const
{
    int info;
    dynamic_cast<LaGenMatDouble&>(B);
    F77_CALL(dgetrs)('N', L_.size(0), B.size(1), &U_(0,0),
		     L_.gdim(0), &pivot_(0), &B(0,0), B.size(0), info);
    return B;
}

inline LaMatrix& LaLUFactorDouble::solve(LaMatrix& X, const LaMatrix& B ) const
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
