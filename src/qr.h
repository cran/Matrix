// -*- c++ -*-
//
//  Copyright (C) 2000-2000 the R Development Core Team
//



#ifndef _LA_QR_FACT_DOUBLE_H
#define _LA_QR_FACT_DOUBLE_H

#include "lafnames.h"
#include LA_VECTOR_INT_H
#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_UPPER_TRIANG_MAT_DOUBLE_H
#include LA_VECTOR_DOUBLE_H
#include LA_UTIL_H

#include "lapackd.h"
#include "solvable.h"

class LaQRFactorDouble : public solvable
{
    LaGenMatDouble             qr_;
    LaUpperTriangMatDouble      R_;  // always a reference to qr_
    LaVectorDouble          qraux_;
    LaVectorInt             pivot_;
    int                      rank_;

public:
				// constructor
    LaQRFactorDouble()
	: qr_(), R_(), qraux_(), pivot_() { rank_ = -1; R_.ref(qr_); }
    inline LaQRFactorDouble(const LaGenMatDouble&);
    inline LaQRFactorDouble(const LaQRFactorDouble&);

    virtual ~LaQRFactorDouble() { };

				// extractor methods for components
    LaUpperTriangMatDouble& R()
	{ return R_; };
    LaGenMatDouble& qr()
	{ return qr_; };
    LaVectorInt& pivot()
	{ return pivot_; };
    int rank(double);
    inline int rank() { return rank(1.0e5); }
    LaVectorDouble& qraux()
	{ return qraux_; }
				// linear equation solvers
    inline LaMatrix& solve() const;// inverse
    LaMatrix& solve(LaMatrix& B) const; // in-place solution
    LaMatrix& solve(LaMatrix& X, const LaMatrix& B) const;

				// operators
    LaQRFactorDouble& ref(const LaQRFactorDouble&);
    inline LaQRFactorDouble& ref(const LaGenMatDouble&);
};



// constructor/destructor functions

inline LaQRFactorDouble::LaQRFactorDouble(const LaGenMatDouble& A)
    : qr_(), R_(), qraux_(), pivot_()
{
    LaGenMatDouble A1;
    A1.copy(A);
    ref(A1);
}

inline LaQRFactorDouble::LaQRFactorDouble(const LaQRFactorDouble& F)
{
  qr_.ref(F.qr_);
  R_.ref(qr_);
  qraux_.ref(F.qraux_);
  pivot_.ref(F.pivot_);
  rank_ = F.rank_;
}

// operators
inline LaQRFactorDouble& LaQRFactorDouble::ref(const LaQRFactorDouble& F)
{

    qr_.ref(F.qr_);
    R_.ref(qr_);
    pivot_.ref(F.pivot_);
    qraux_ = F.qraux_;
    rank_ = F.rank_;
    
    return *this;
}

inline LaMatrix& LaQRFactorDouble::solve() const
{
    return *new LaGenMatDouble(0,0);
}

#endif
