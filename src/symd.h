//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
// Modifications Copyright (C) 2000-2000 the R Development Core Team

#ifndef _LA_SYMM_MAT_DOUBLE_H_
#define _LA_SYMM_MAT_DOUBLE_H_

#include "lafnames.h"
#include LA_LOWER_TRIANG_MAT_DOUBLE_H

class LaSymmMatDouble : public LaMatDouble
{
    LaLowerTriangMatDouble lower_data_;

public:

    // constructors
    LaSymmMatDouble()
	: lower_data_() { };
    LaSymmMatDouble(int i, int j)
	: lower_data_(i,j) { }
    LaSymmMatDouble(double* d, int i, int j)
	: lower_data_(d,i,j) { };
    LaSymmMatDouble(const LaSymmMatDouble& A)
	{ lower_data_.copy(A.lower_data_); };
    LaSymmMatDouble(SEXP s)
	: lower_data_(s) { };

    // destructor
    ~LaSymmMatDouble() { };

    int size(int d) const          // submatrix size
	{ return lower_data_.size(d); }
    int gdim(int d) const          // global dimensions
	{ return lower_data_.gdim(d); }
    LaIndex index(int d) const     // return indices of matrix.
        { return lower_data_.index(d); }
    int ref_count() const          // return ref_count of matrix.
        { return lower_data_.ref_count(); }
    double* addr() const           // return address of matrix.
        { return lower_data_.addr(); }

    // operators
    double& operator()(int i,int j) 
	{ if (i < j) return lower_data_(j,i); else return lower_data_(i,j); }
    double& operator()(int i, int j) const
	{ if (i < j) return lower_data_(j,i); else return lower_data_(i,j); }
    LaMatDouble& operator=(double s)
	{ lower_data_ = s; return *this; }

    LaMatrix& inject(const LaMatrix& A)
	{ return copy(A); }
    LaMatrix& resize(const LaMatrix& A)
	{ return resize(A.size(0), A.size(1)); }
    inline LaMatrix& resize(int m, int n)
	{ lower_data_.resize(m, n); return *this; }
    LaMatrix& ref(const LaSymmMatDouble& A)
	{ lower_data_.ref(A.lower_data_); return *this; }
    LaMatrix& ref(SEXP s)
	{ lower_data_.ref(s); return *this; }
    LaMatrix& copy(const LaMatrix &);
   
    ostream &printMatrix(ostream &) const;

    operator LaGenMatDouble();
    operator LaLowerTriangMatDouble();

				// linear equation solvers
    LaMatrix& solve() const;	// inverse
    LaMatrix& solve(LaMatrix& B) const; // in-place solution
    LaMatrix& solve(LaMatrix& X, const LaMatrix& B) const;
				// matrix norms
    double norm(char) const;
    
};

#endif
