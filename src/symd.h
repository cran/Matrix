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

#include "bunch_kaufman.h"

#ifndef _LA_SYMM_MAT_DOUBLE_H_
#define _LA_SYMM_MAT_DOUBLE_H_

#include "lafnames.h"
#include "tgmd.h"
#include "eigen.h"

class LaSymmMatDouble : public LaMatDouble
{
    LaTriangMatDouble data_;

protected:

    mutable LaSymmFactor* factor_;

    void setFactor(LaSymmFactor *factor) const
	{ clearDecomposition(); factor_ = factor; }
public:
				// constructors
    LaSymmMatDouble(char uplo = 'U')
	: data_(uplo)
	{  factor_ = 0; };
    LaSymmMatDouble(int i, int j, char uplo = 'U')
	: data_(i, j, uplo)
	{  factor_ = 0; };
    LaSymmMatDouble(double* d, int i, int j, char uplo = 'U')
	: data_(d, i, j, uplo)
	{  factor_ = 0; };
    LaSymmMatDouble(const LaSymmMatDouble& A)
	: data_(A.data_)
	{  factor_ = 0; };
    explicit LaSymmMatDouble(SEXP s, char uplo = 'U')
	: data_(s, uplo)
	{  factor_ = 0; };
				// destructor
    ~LaSymmMatDouble()
	{ delete factor_; };

    char uplo() const { return data_.uplo(); }
    const LaSymmFactor& factor() const
	{
	    if (factor_ == 0)
		throw(LaException("No factor present"));
	    return *factor_;
	}
    int size(int d) const	// submatrix size
	{ return data_.size(d); }	
    int gdim(int d) const	// global dimensions
	{ return data_.gdim(d); }
    LaIndex index(int d) const	// return indices of matrix.
	{ return data_.index(d); }
    int ref_count() const	// return ref_count of matrix.
	{ return data_.ref_count(); }
    double* addr() const	// return address of matrix.
	{ return data_.addr(); };

				// operators
    double& operator()(int i,int j) 
	{
	    if (uplo() == 'U') {
		if (i < j)
		    return data_(i, j);
		else return data_(j, i);
	    } else {
		if (i < j)
		    return data_(j, i);
		else return data_(i, j);
	    }
	}
    double& operator()(int i, int j) const
	{
	    if (uplo() == 'U') {
		if (i < j)
		    return data_(i, j);
		else return data_(j, i);
	    } else {
		if (i < j)
		    return data_(j, i);
		else return data_(i, j);
	    }
	}
    LaMatDouble& operator=(double s)
	{ clearDecomposition(); data_ = s; return *this; };

    LaSymmMatDouble& inject(const LaMatDouble& A)
	{ clearDecomposition(); data_.inject(A); return *this; }
    LaSymmMatDouble& resize(const LaMatDouble& A)
	{ return resize(A.size(0), A.size(1)); }
    LaSymmMatDouble& resize(int m, int n)
	{ clearDecomposition(); data_.resize(m, n); return *this; };
    LaSymmMatDouble& ref(const LaSymmMatDouble& A)
	{ clearDecomposition(); data_.ref(A.data_); return *this; };
    LaSymmMatDouble& ref(SEXP s)
	{ clearDecomposition(); data_.ref(s); return *this; };
    LaSymmMatDouble& copy(const LaMatDouble& A)
	{ clearDecomposition(); data_.copy(A); return *this; };
    inline LaSymmMatDouble* clone() const;
   
    ostream &printMatrix(ostream &) const;

    operator LaGenMatDouble();
    operator LaTriangMatDouble()
	{
	    LaTriangMatDouble ans;
	    ans.ref(data_);
	    return ans;
	}
//    operator LaLowerTriangMatDouble();
//    operator LaUpperTriangMatDouble();

    inline const LaSymmFactor& doDecomposition() const;
    inline void clearDecomposition() const;

    inline double rcond(char which) const;

				// linear equation solvers
    inline LaSymmMatDouble* solve() const;   // inverse
    inline LaMatDouble& solve(LaMatDouble& B) const;
    inline LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const;
				// eigenvalues/eigenvectors
    inline LaSymmEigenDouble* eigen(bool leftEV = true, bool rightEV = true)
	{ return new LaSymmEigenDouble(*this, uplo(), leftEV | rightEV); }
				// matrix norms, etc.
    double norm(char) const;
    SEXP asSEXP() const;
};

inline void LaSymmMatDouble::clearDecomposition() const
{
    if (factor_ != 0) {
	delete factor_;
	factor_ = 0;
    }
}


inline const LaSymmFactor& LaSymmMatDouble::doDecomposition() const
{
    clearDecomposition();
    setFactor(new LaBunchKaufmanFactorDouble());
    LaSymmMatDouble tmp;
    tmp.copy(*this);
    factor_->ref(tmp);
    return factor();
}

inline double LaSymmMatDouble::rcond(char which) const
{
    if (factor_ == 0)
	doDecomposition();

    return factor().rcond(norm('O'));
}

inline LaSymmMatDouble* LaSymmMatDouble::solve() const   // inverse
{
    if (factor_ == 0)
	doDecomposition();
    return dynamic_cast<LaSymmMatDouble*>(factor().solve());
}

inline LaMatDouble& LaSymmMatDouble::solve(LaMatDouble& B) const
{
    if (factor_ == 0)
	doDecomposition();
    return factor().solve(B);
}

inline LaMatDouble& LaSymmMatDouble::solve(LaMatDouble& X, const LaMatDouble& B) const
{
    if (factor_ == 0)
	doDecomposition();
    return factor().solve(X, B);
}

inline LaSymmMatDouble* LaSymmMatDouble::clone() const
{
    LaTriangMatDouble* tmp = data_.clone();
    LaSymmMatDouble* ans = new LaSymmMatDouble();
    ans->data_.ref(*tmp);
    delete tmp;
    return ans;
}

#endif
