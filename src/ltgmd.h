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
// Modifications Copyright (C) 2000-2002 the R Development Core Team

#include "lafnames.h"
#ifndef _LA_GEN_MAT_DOUBLE_H_
#include LA_GEN_MAT_DOUBLE_H
#endif

#ifndef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
#define _LA_LOWER_TRIANG_MAT_DOUBLE_H_

class LaLowerTriangMatDouble : public LaMatDouble
{
    LaGenMatDouble data_;
    static double outofbounds_;

public:

    // constructors
    LaLowerTriangMatDouble()
        : data_() { shallow_ = 0; }
    LaLowerTriangMatDouble(int i, int j)
        : data_(i, j) { }
    LaLowerTriangMatDouble(double* d, int i, int j)
        : data_(d,i,j) { }
    explicit LaLowerTriangMatDouble(SEXP s)
        : data_(s) { }
    LaLowerTriangMatDouble(const LaLowerTriangMatDouble& A)
        : data_()
	{ data_.ref(A.data_); }

    // destructor
    ~LaLowerTriangMatDouble() { }

    int size(int d) const          // submatrix size
	{ return data_.size(d); }
    int gdim(int d) const          // global dimensions
	{ return data_.gdim(d); }
    LaIndex index(int d) const     // return indices of matrix.
        { return data_.index(d); }
    int ref_count() const          // return ref_count of matrix.
        { return data_.ref_count(); }
    const double* addr() const           // return address of matrix.
        { return data_.addr(); }
    double* addr()               // return address of matrix.
        { return data_.addr(); }

    // operators
    double& operator()(int i,int j)
	{
        if (i < j) {
            outofbounds_ = 0;
            return outofbounds_;
        } else return data_(i,j);
    }
    double operator()(int i, int j) const
	{
        if (i < j) {
            outofbounds_ = 0;
            return outofbounds_;
        } else return data_(i,j);
    }
    LaMatDouble& operator=(double); 
//      operator LaGenMatDouble()
//  	{ LaGenMatDouble G; G.ref((*this).data_); return G; };

    LaLowerTriangMatDouble& inject(const LaMatDouble& A)
	{ data_.inject(A); return *this; }
    LaLowerTriangMatDouble& resize(const LaMatDouble& A)
	{ return resize(A.size(0), A.size(1)); }
    LaLowerTriangMatDouble& resize(int m, int n)
	{ data_.resize(m, n); return *this; }
    LaLowerTriangMatDouble& ref(const LaLowerTriangMatDouble& A)
	{ data_.ref(A.data_); return *this; }
    LaLowerTriangMatDouble& ref(const LaGenMatDouble& A)
	{ data_.ref(A); return *this; }
    LaLowerTriangMatDouble& ref(SEXP s)
	{ data_.ref(s); return *this; }
    LaLowerTriangMatDouble& copy(const LaMatDouble &);
    inline LaLowerTriangMatDouble* clone() const;
				// linear equation solvers
    LaLowerTriangMatDouble* solve() const;	// inverse
    LaMatDouble& solve(LaMatDouble& B) const; // in-place solution
    LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const;
				// matrix norms
    double norm(char which) const;
    double rcond(char which) const;
    SEXP asSEXP() const;

    std::ostream& printMatrix(std::ostream&) const;
};

inline LaLowerTriangMatDouble* LaLowerTriangMatDouble::clone() const
{
    LaGenMatDouble* tmp = data_.clone();
    LaLowerTriangMatDouble* ans = new LaLowerTriangMatDouble();
    ans->data_.ref(*tmp);
    delete tmp;
    return ans;
}

#endif 
// _LA_LOWER_TRIANG_MAT_DOUBLE_H_
