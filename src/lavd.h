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
// Modifications Copyright (C) 2000-2001 the R Development Core Team

#ifndef _LA_VECTOR_DOUBLE_H_
#define _LA_VECTOR_DOUBLE_H_

#include "lafnames.h"

#ifndef _LA_GEN_MAT_DOUBLE_H_
#include LA_GEN_MAT_DOUBLE_H
#endif

#include "laexcp.h"

// a vector is simply an nx1 or 1xn, matrix, only that it can
// be constructed and accessed by a single dimension

class LaVectorDouble: public LaGenMatDouble
{
 public:
    
    inline LaVectorDouble();
    inline LaVectorDouble(int);
    inline LaVectorDouble(int, int);  
    inline LaVectorDouble(double*, int);
//    inline LaVectorDouble(double*, int, int);
    inline LaVectorDouble(const LaGenMatDouble&);
    inline explicit LaVectorDouble(SEXP);
    
    inline int size() const;
    inline int inc() const;
    inline LaIndex index() const;
    inline int start() const;
    inline int end() const;
    
    inline int size(int d) const
        { if (d != 0) throw LaException("Vectors are one dimensional"); return size();}
    inline int inc(int d) const
        { if (d != 0) throw LaException("Vectors are one dimensional"); return inc();}
    inline LaIndex index(int d) const
        { if (d != 0) throw LaException("Vectors are one dimensional"); return index();}
    inline int start(int d) const
        { if (d != 0) throw LaException("Vectors are one dimensional"); return start();}
    inline int end(int d) const
        { if (d != 0) throw LaException("Vectors are one dimensional"); return end();}
    
    inline LaVectorDouble& ref(const LaGenMatDouble &);
    LaVectorDouble& ref(SEXP);
    inline LaVectorDouble& inject(const LaMatDouble &);
    inline LaVectorDouble& copy(const LaMatDouble &);
    inline LaVectorDouble* clone() const;
    LaVectorDouble& resize(int m)
	{ LaGenMatDouble::resize(m, 1); return *this; };
    LaVectorDouble& resize(int m, int n)
        { if (m!= 1 || n != 1) throw LaException("Vectors are one dimensional"); resize(m*n); return *this; }
    LaVectorDouble& resize(const LaVectorDouble &A)
        { resize(A.size()); return *this; }
    LaVectorDouble& resize(const LaMatDouble &A)
        { resize(A.size(0)*A.size(1)); return *this; }

    const double& operator()(int i, int j) const
	{ return LaGenMatDouble::operator()(i,j); }
    double& operator()(int i, int j)
	{ return LaGenMatDouble::operator()(i,j); }

    inline double& operator()(int i);
    inline const double& operator()(int i) const ;
    inline LaVectorDouble operator()(const LaIndex&);
    inline LaVectorDouble operator()(const LaIndex&) const;
    
    inline LaVectorDouble& operator=(double);

    SEXP asSEXP() const;
};

// NOTE: we default to column vectors, since matrices are column
//  oriented.

inline LaVectorDouble::LaVectorDouble() : LaGenMatDouble(0,1) {}
inline LaVectorDouble::LaVectorDouble(int i) : LaGenMatDouble(i,1) {}

// NOTE: one shouldn't be using this method to initalize, but
// it is here so that the constructor can be overloaded with 
// a runtime test.
//
inline LaVectorDouble::LaVectorDouble(int m, int n) : LaGenMatDouble(m,n)
{
    if (!(n==1 || m==1)) throw(LaException("assert failed : n==1 || m==1"));
}

inline LaVectorDouble::LaVectorDouble(double *d, int m) : 
    LaGenMatDouble(d,m,1) {}

inline LaVectorDouble::LaVectorDouble(SEXP s)
    : LaGenMatDouble()
{
    LaVectorDouble tmp(REAL(coerceVector(s, REALSXP)), LENGTH(s));
    copy(tmp);
}

#if 0
inline LaVectorDouble::LaVectorDouble(double *d, int m, int n) : 
    LaGenMatDouble(d,m,n) {}
#endif

inline LaVectorDouble::LaVectorDouble(const LaGenMatDouble& G) : 
    LaGenMatDouble(G)
{
    if (!(G.size(0)==1 || G.size(1)==1)) throw(LaException("assert failed : G.size(0)==1 || G.size(1)==1"));

}

//note that vectors can be either stored columnwise, or row-wise
// this will handle the 0x0 case as well.

inline int LaVectorDouble::size() const 
{
    return LaGenMatDouble::size(0)*LaGenMatDouble::size(1);
}

inline double& LaVectorDouble::operator()(int i)
{
    if (LaGenMatDouble::size(0)==1 )
	return LaGenMatDouble::operator()(0,i);
    return LaGenMatDouble::operator()(i,0);
}

inline const double& LaVectorDouble::operator()(int i) const
{
    if (LaGenMatDouble::size(0)==1 )
	return LaGenMatDouble::operator()(0,i);
    return LaGenMatDouble::operator()(i,0);
}

inline LaVectorDouble LaVectorDouble::operator()(const LaIndex& I)
{
    if (LaGenMatDouble::size(0)==1)
        return LaGenMatDouble::operator()(LaIndex(0,0),I);
    return LaGenMatDouble::operator()(I,LaIndex(0,0));
}

inline LaVectorDouble LaVectorDouble::operator()(const LaIndex& I) const
{
    if (LaGenMatDouble::size(0)==1)
        return LaGenMatDouble::operator()(LaIndex(0,0),I);
    return LaGenMatDouble::operator()(I,LaIndex(0,0));
}

inline LaVectorDouble& LaVectorDouble::copy(const LaMatDouble &A)
{
				//make sure rhs is a vector
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));
    LaGenMatDouble::copy(A);
    return *this;
}

inline LaVectorDouble* LaVectorDouble::clone() const
{
    LaGenMatDouble* tmp = LaGenMatDouble::clone();
    LaVectorDouble* ans = new LaVectorDouble();
    ans->ref(*tmp);
    delete tmp;
    return ans;
}

inline LaVectorDouble& LaVectorDouble::ref(const LaGenMatDouble &A)
{
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));
    LaGenMatDouble::ref(A);
    return *this;
}

inline LaVectorDouble& LaVectorDouble::operator=(double d)
{
    LaGenMatDouble::operator=(d);
    return *this;
}

inline LaVectorDouble& LaVectorDouble::inject(const LaMatDouble& A)
{
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));
    LaGenMatDouble::inject(A);
    return *this;
}

inline int LaVectorDouble::inc() const
{
    if (LaGenMatDouble::size(1)==1 )
	return LaGenMatDouble::inc(0);
    return LaGenMatDouble::inc(1);
}

inline LaIndex LaVectorDouble::index() const
{
    if (LaGenMatDouble::size(1)==1 )
	return LaGenMatDouble::index(0);
    return LaGenMatDouble::index(1);
}

inline int LaVectorDouble::start() const
{
    if (LaGenMatDouble::size(1)==1 )
	return LaGenMatDouble::start(0);
    return LaGenMatDouble::start(1);
}

inline int LaVectorDouble::end() const
{
    if (LaGenMatDouble::size(1)==1 )
	return LaGenMatDouble::end(0);
    return LaGenMatDouble::end(1);
}
#endif // _LA_VECTOR_DOUBLE_H_

