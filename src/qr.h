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

#ifndef _LA_QR_FACT_DOUBLE_H
#define _LA_QR_FACT_DOUBLE_H

#include "lafnames.h"
#include LA_VECTOR_INT_H
#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_UPPER_TRIANG_MAT_DOUBLE_H
#include LA_VECTOR_DOUBLE_H
#include LA_UTIL_H

#include "lapackd.h"
#include "factor.h"

class LaQRFactorDouble : public Factor
{
    LaGenMatDouble             qr_;
    LaUpperTriangMatDouble      R_;  // always a reference to qr_
    LaVectorDouble          qraux_;
    LaVectorInt             pivot_;
    int                      rank_;

public:
				// constructor
    LaQRFactorDouble()
	: rank_(-1) { R_.ref(qr_); }
    inline explicit LaQRFactorDouble(const LaGenMatDouble&);
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
    int rank()
	{ return rank(1.0e6); };
    bool isSingular()
	{ return rank() < min(qr_.size(0), qr_.size(1)); };
    LaVectorDouble& qraux()
	{ return qraux_; }
				// linear equation solvers
    LaGenMatDouble* solve() const;// inverse
    LaMatDouble& solve(LaMatDouble& B) const; // in-place solution
    LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const;

    LaMatDouble& applyQ(LaMatDouble& y, bool left = true,
		     bool transpose = true) const;

				// operators
    LaQRFactorDouble& ref(const LaQRFactorDouble&);
    inline LaQRFactorDouble& ref(const LaGenMatDouble&);
};



// constructor/destructor functions

inline LaQRFactorDouble::LaQRFactorDouble(const LaGenMatDouble& A)
{
    LaGenMatDouble A1;
    A1.copy(A);
    ref(A1);
}

inline LaQRFactorDouble::LaQRFactorDouble(const LaQRFactorDouble& F)
{
    ref(F);
}

// operators
inline LaQRFactorDouble& LaQRFactorDouble::ref(const LaQRFactorDouble& F)
{

    qr_.ref(F.qr_);
    R_.ref(qr_);
    pivot_.ref(F.pivot_);
    qraux_.ref(F.qraux_);
    rank_ = F.rank_;
    
    return *this;
}

#endif
