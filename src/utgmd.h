// -*- c++ -*-
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
// Modifications Copyright (C) 2000-2000 the R Development Core Team

#include "lafnames.h"
#ifndef _LA_GEN_MAT_DOUBLE_H_
#include LA_GEN_MAT_DOUBLE_H
#endif

#ifndef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
#define _LA_UPPER_TRIANG_MAT_DOUBLE_H_

class LaUpperTriangMatDouble : public LaMatDouble
{
    LaGenMatDouble data_;
    static double outofbounds_;

public:

    // constructors
    LaUpperTriangMatDouble()
	{ *info_ = 0; };
    LaUpperTriangMatDouble(int i, int j)
	: data_(i, j) { *info_ = 0; };
    LaUpperTriangMatDouble(double* d, int i, int j)
	: data_(d,i,j) { *info_ = 0; };
    LaUpperTriangMatDouble(SEXP s)
	: data_(s) { };
    LaUpperTriangMatDouble(const LaUpperTriangMatDouble& A)
	{ data_.copy(A.data_); };

    // destructor
    ~LaUpperTriangMatDouble() { };

    int size(int d) const          // submatrix size
	{ return data_.size(d); }
    int gdim(int d) const          // global dimensions
	{ return data_.gdim(d); }
    LaIndex index(int d) const     // return indices of matrix.
        { return data_.index(d); }
    int ref_count() const          // return ref_count of matrix.
        { return data_.ref_count(); }
    double* addr() const           // return address of matrix.
        { return data_.addr(); }

    // operators
    inline double& operator()(int i, int j)
	{ if (i > j) return outofbounds_; else return data_(i,j); }
    inline double& operator()(int i, int j) const
	{ if (i > j) return outofbounds_; else return data_(i,j); }
    LaMatDouble& operator=(double); 
    operator LaGenMatDouble()
	{ LaGenMatDouble G; G.ref((*this).data_); return G; };

    LaMatrix& inject(const LaMatrix& A)
	{ data_.inject(A); return *this; }
    LaMatrix& resize(const LaMatrix& A)
	{ return resize(A.size(0), A.size(1)); }
    LaMatrix& resize(int m, int n)
	{ data_.resize(m, n); return *this; }
    LaMatrix& ref(const LaUpperTriangMatDouble& A)
	{ data_.ref(A.data_); return *this; }
    LaMatrix& ref(const LaGenMatDouble& A)
	{ data_.ref(A); return *this; }
    LaMatrix& ref(SEXP s)
	{ data_.ref(s); return *this; }
    LaMatrix& copy(const LaMatrix &);
				// linear equation solvers
    LaMatrix& solve() const;	// inverse
    LaMatrix& solve(LaMatrix& B) const; // in-place solution
    LaMatrix& solve(LaMatrix& X, const LaMatrix& B) const;
				// matrix norms
    double norm(char which) const;

    ostream &printMatrix(ostream &) const;
};
#endif 
// _LA_UPPER_TRIANG_MAT_DOUBLE_H_
