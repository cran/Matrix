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

#include "lafnames.h"
#ifndef _LA_GEN_MAT_DOUBLE_H_
#include LA_GEN_MAT_DOUBLE_H
#endif

#ifndef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
#define _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_

class LaUnitLowerTriangMatDouble : public LaMatDouble
{
    LaGenMatDouble data_;
    static double outofbounds_;
    
public:
				// constructors
    LaUnitLowerTriangMatDouble()
	: data_() { *info_ = 0; };
    LaUnitLowerTriangMatDouble(int i, int j)
	: data_(i, j) { *info_ = 0; };
    LaUnitLowerTriangMatDouble(double* d, int i, int j)
	: data_(d,i,j) { *info_ = 0; };
    LaUnitLowerTriangMatDouble(SEXP s)
	: data_(s) { };
    LaUnitLowerTriangMatDouble(const LaUnitLowerTriangMatDouble& A)
	{ data_.copy(A.data_); };
				// destructor
    ~LaUnitLowerTriangMatDouble() { };

    int size(int d) const	// submatrix size
	{ return data_.size(d); }
    int gdim(int d) const	// global dimensions
	{ return data_.gdim(d); }
    LaIndex index(int d) const	// return indices of matrix.
        { return data_.index(d); }
    int ref_count() const	// return ref_count of matrix.
        { return data_.ref_count(); }
    double* addr() const	// return address of matrix.
        { return data_.addr(); }

				// operators
    double& operator()(int,int);
    double& operator()(int,int) const;
    LaMatDouble& operator=(double); 
    operator LaGenMatDouble()
	{ LaGenMatDouble G; G.ref((*this).data_); return G; };
    
    LaUnitLowerTriangMatDouble& inject(const LaMatDouble& A)
	{ data_.inject(A); return *this; }
    LaUnitLowerTriangMatDouble& resize(const LaMatDouble& A)
	{ return resize(A.size(0), A.size(1)); }
    LaUnitLowerTriangMatDouble& resize(int m, int n)
	{ data_.resize(m, n); return *this; }
    LaUnitLowerTriangMatDouble& ref(const LaUnitLowerTriangMatDouble& A)
	{ data_.ref(A.data_); return *this; }
    LaUnitLowerTriangMatDouble& ref(const LaGenMatDouble& A)
	{ data_.ref(A); return *this; }
    LaUnitLowerTriangMatDouble& ref(SEXP s)
	{ data_.ref(s); return *this; }
    LaUnitLowerTriangMatDouble& copy(const LaMatDouble &);
				// linear equation solvers
    LaUnitLowerTriangMatDouble& solve() const;	// inverse
    LaMatDouble& solve(LaMatDouble& B) const; // in-place solution
    LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const;

    double norm(char which) const;
    double rcond(char which) const;
    SEXP asSEXP() const;

    ostream &printMatrix(ostream &) const;
};

#endif 
// _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
