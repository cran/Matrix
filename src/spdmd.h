//  -*- c++ -*-
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
// Modifications Copyright (C) 2000-2000, 2002 the R Development Core Team

#include "cholesky.h"

#ifndef _LA_SPD_MAT_DOUBLE_H_
#define _LA_SPD_MAT_DOUBLE_H_

#include LA_SYMM_MAT_DOUBLE_H
#include "vi.h"

class LaSpdMatDouble : public LaSymmMatDouble
{
public:
				// constructors and destructors
    LaSpdMatDouble(char uplo = 'U')
	: LaSymmMatDouble(uplo)
	{ };
    LaSpdMatDouble(int i, int j, char uplo = 'U')
	: LaSymmMatDouble(i, j, uplo)
	{ };
    LaSpdMatDouble(double* d, int i, int j, char uplo = 'U')
	: LaSymmMatDouble(d, i, j, uplo)
	{ };
    LaSpdMatDouble(const LaSymmMatDouble& A)
	: LaSymmMatDouble(A)
	{ };
    explicit LaSpdMatDouble(SEXP s, char uplo = 'U')
	: LaSymmMatDouble(s, uplo)
	{ };
    ~LaSpdMatDouble() { }
	
                                // operators
    LaSpdMatDouble& operator=(double s)
	{ LaSymmMatDouble::operator=(s); return *this; }
    LaSpdMatDouble& ref(const LaSpdMatDouble& A)
	{ LaSymmMatDouble::ref(A); return *this; }
    LaSpdMatDouble& ref(SEXP s)
	{ LaSymmMatDouble::ref(s); return *this; }
    LaSpdMatDouble& inject(const LaMatDouble& A)
	{ LaSymmMatDouble::inject(A); return *this; }
    LaSpdMatDouble& copy(const LaMatDouble& A)
	{ LaSymmMatDouble::copy(A); return *this; }
    inline LaSpdMatDouble* clone() const;

				// linear equation solvers
    inline LaSpdMatDouble* solve() const;   // inverse
    inline LaMatDouble& solve(LaMatDouble& B) const
        { return LaSymmMatDouble::solve(B); }
    inline LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const
        { return LaSymmMatDouble::solve(X, B); }
				// matrix norms
//    double norm(char) const;
//    inline double rcond(char which) const;
    SEXP asSEXP() const;

};

inline LaSpdMatDouble* LaSpdMatDouble::solve() const   // inverse
{
    if (factor_ == 0)
	doDecomposition();
    return dynamic_cast<LaSpdMatDouble*>(factor().solve());
}

inline LaSpdMatDouble* LaSpdMatDouble::clone() const
{
    LaSymmMatDouble* tmp = LaSymmMatDouble::clone();
    LaSpdMatDouble* ans = new LaSpdMatDouble(*tmp);
    delete tmp;
    return ans;
}

#endif 

// _LA_SPD_MAT_DOUBLE_H_
