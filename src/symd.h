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

#ifndef _LA_SYMM_MAT_DOUBLE_H_
#define _LA_SYMM_MAT_DOUBLE_H_

#include "lafnames.h"
#include LA_LOWER_TRIANG_MAT_DOUBLE_H

class LaSymmMatDouble : public LaMatDouble
{
protected:
    LaLowerTriangMatDouble lower_data_;
public:
				// constructors
    LaSymmMatDouble() : lower_data_() { }
    LaSymmMatDouble(int i, int j) : lower_data_(i,j) { }
    LaSymmMatDouble(double* d, int i, int j) : lower_data_(d,i,j) { }
    LaSymmMatDouble(const LaSymmMatDouble& A)
	{ lower_data_.copy(A.lower_data_); }
    LaSymmMatDouble(SEXP s) : lower_data_(s) { };
				// destructor
    ~LaSymmMatDouble() { };

    int size(int d) const	// submatrix size
	{ return lower_data_.size(d); }
    int gdim(int d) const	// global dimensions
	{ return lower_data_.gdim(d); }
    LaIndex index(int d) const	// return indices of matrix.
        { return lower_data_.index(d); }
    int ref_count() const	// return ref_count of matrix.
        { return lower_data_.ref_count(); }
    double* addr() const	// return address of matrix.
        { return lower_data_.addr(); }

				// operators
    double& operator()(int i,int j) 
	{ if (i < j) return lower_data_(j,i); else return lower_data_(i,j); }
    double& operator()(int i, int j) const
	{ if (i < j) return lower_data_(j,i); else return lower_data_(i,j); }
    LaMatDouble& operator=(double s)
	{ lower_data_ = s; return *this; }

    LaSymmMatDouble& inject(const LaMatDouble& A)
	{ return copy(A); }
    LaSymmMatDouble& resize(const LaMatDouble& A)
	{ return resize(A.size(0), A.size(1)); }
    inline LaSymmMatDouble& resize(int m, int n)
	{ lower_data_.resize(m, n); return *this; }
    LaSymmMatDouble& ref(const LaSymmMatDouble& A)
	{ lower_data_.ref(A.lower_data_); return *this; }
    LaSymmMatDouble& ref(SEXP s)
	{ lower_data_.ref(s); return *this; }
    LaSymmMatDouble& copy(const LaMatDouble &);
   
    ostream &printMatrix(ostream &) const;

    operator LaGenMatDouble();
    operator LaLowerTriangMatDouble();
				// linear equation solvers
    LaSymmMatDouble& solve() const;	// inverse
    LaMatDouble& solve(LaMatDouble& B) const; // in-place solution
    LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const;
    
				// matrix norms, etc.
    double norm(char) const;
    double rcond(char which) const;
    SEXP asSEXP() const;
};

#endif
