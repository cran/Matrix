//  -*- c++ -*-
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//  Modifications Copyright (C) 2000-2000 the R Development Core Team


#ifndef _LA_SPD_MAT_DOUBLE_H_
#define _LA_SPD_MAT_DOUBLE_H_

#include LA_SYMM_MAT_DOUBLE_H

class LaSpdMatDouble : public LaSymmMatDouble
{
public:
    // constructors and destructors
    LaSpdMatDouble()
	: LaSymmMatDouble() { };
    LaSpdMatDouble(int i, int j)
	: LaSymmMatDouble(i, j) { };
    LaSpdMatDouble(double* d, int i, int j)
	: LaSymmMatDouble(d, i, j) { };
    LaSpdMatDouble(const LaSpdMatDouble& A)
	: LaSymmMatDouble(A) { };
    LaSpdMatDouble(SEXP s)
	: LaSymmMatDouble(s) { };
    ~LaSpdMatDouble() { }

    // operators
    LaMatDouble& operator=(double s)
	{ LaSymmMatDouble::operator=(s); return *this; }
    LaMatrix& ref(const LaSpdMatDouble& A)
	{ LaSymmMatDouble::ref(A); return *this; }
    LaMatrix& ref(SEXP s)
	{ LaSymmMatDouble::ref(s); return *this; }
    LaMatrix& inject(const LaMatrix& A)
	{ LaSymmMatDouble::inject(A); return *this; }
    LaMatrix& copy(const LaMatrix& A)
	{ LaSymmMatDouble::copy(A); return *this; }

				// linear equation solvers
    LaMatrix& solve() const {throw(LaException("",""));}	// inverse
    LaMatrix& solve(LaMatrix& B) const  {throw(LaException("",""));}; // in-place solution
    LaMatrix& solve(LaMatrix& X, const LaMatrix& B) const  {throw(LaException("",""));};
				// matrix norms
    double norm(char) const  {throw(LaException("",""));};
};

#endif 
// _LA_SPD_MAT_DOUBLE_H_
