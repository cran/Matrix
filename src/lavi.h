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

#ifndef _LA_VECTOR_INT_H_
#define _LA_VECTOR_INT_H_

#include "lafnames.h"

#ifndef _LA_GEN_MAT_INT_H_
#include LA_GEN_MAT_INT_H
#endif

#include "laexcp.h"

extern "C" {
#include <R.h>
#include <Rinternals.h>
}

// a vector is simply an nx1 or 1xn, matrix, only that it can
// be constructed and accessed by a single dimension

class LaVectorInt: public LaGenMatInt
{
public:

    inline LaVectorInt();
    inline LaVectorInt(int);
    inline LaVectorInt(int, int);  
    inline LaVectorInt(int*, int);
    inline LaVectorInt(int*, int, int);
    inline LaVectorInt(const LaGenMatInt&);
    ~LaVectorInt() { };

    inline int size() const;
    inline int inc() const;
    inline LaIndex index() const;
    inline int start() const;
    inline int end() const;

    inline LaVectorInt& ref(const LaGenMatInt &);
    inline LaVectorInt& inject(const LaGenMatInt &);
    inline LaVectorInt& copy(const LaGenMatInt &);
    LaVectorInt& resize(int m)
	{ LaGenMatInt::resize(m, 1); return *this; };

    inline int& operator()(int i);
    inline int& operator()(int i) const ;
    inline LaVectorInt operator()(const LaIndex&);

    inline int& operator()(int i, int j)
	{ return LaGenMatInt::operator()(i,j); }
    inline int& operator()(int i, int j) const
	{ return LaGenMatInt::operator()(i,j); }

    inline LaVectorInt& operator=(int);
    inline LaVectorInt& operator=(const LaGenMatInt&);

    SEXP asSEXP() const;
};

// NOTE: we default to column vectors, since matrices are column
//  oriented.

inline LaVectorInt::LaVectorInt() : LaGenMatInt(0,1) {}
inline LaVectorInt::LaVectorInt(int i) : LaGenMatInt(i,1) {}

// NOTE: one shouldn't be using this method to initalize, but
// it is here so that the constructor can be overloaded with 
// a runtime test.
//
inline LaVectorInt::LaVectorInt(int m, int n) : LaGenMatInt(m,n)
{
    if (!(n==1 || m==1)) throw(LaException("assert failed : n==1 || m==1"));
}

inline LaVectorInt::LaVectorInt(int *d, int n) : 
    LaGenMatInt(d,n,1) {}

inline LaVectorInt::LaVectorInt(int *d, int n, int m) : 
    LaGenMatInt(d,n,m) {}

inline LaVectorInt::LaVectorInt(const LaGenMatInt &G) 
{
        if (!(G.size(0)==1 || G.size(1)==1)) throw(LaException("assert failed : G.size(0)==1 || G.size(1)==1"));

        (*this).ref(G);
}


//note that vectors can be either stored columnwise, or row-wise

// this will handle the 0x0 case as well.

inline int LaVectorInt::size() const 
{ return LaGenMatInt::size(0)*LaGenMatInt::size(1); }

inline int& LaVectorInt::operator()(int i)
{ if (LaGenMatInt::size(0)==1 )
    return LaGenMatInt::operator()(0,i);
  else
    return LaGenMatInt::operator()(i,0);
}

inline int& LaVectorInt::operator()(int i) const
{ if (LaGenMatInt::size(0)==1 )
    return LaGenMatInt::operator()(0,i);
  else
    return LaGenMatInt::operator()(i,0);
}

inline LaVectorInt LaVectorInt::operator()(const LaIndex& I)
{ if (LaGenMatInt::size(0)==1)
    return LaGenMatInt::operator()(LaIndex(0,0),I); 
  else
    return LaGenMatInt::operator()(I,LaIndex(0,0)); 
}


inline LaVectorInt& LaVectorInt::copy(const LaGenMatInt &A)
{
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));   //make sure rhs is a
                                                // a vector.
    LaGenMatInt::copy(A);
    return *this;
}

inline LaVectorInt& LaVectorInt::operator=(const  LaGenMatInt &A)
{
    return inject(A);
}

inline LaVectorInt& LaVectorInt::ref(const LaGenMatInt &A)
{
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));
    LaGenMatInt::ref(A);
    return *this;
}

inline LaVectorInt& LaVectorInt::operator=(int d)
{
    LaGenMatInt::operator=(d);
    return *this;
}

inline LaVectorInt& LaVectorInt::inject(const LaGenMatInt &A)
{
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));
    LaGenMatInt::inject(A);
    return *this;
}
    
inline int LaVectorInt::inc() const
{
  if (LaGenMatInt::size(1)==1 )
    return LaGenMatInt::inc(0);
  else
    return LaGenMatInt::inc(1);
}

inline LaIndex LaVectorInt::index() const
{
  if (LaGenMatInt::size(1)==1 )
    return LaGenMatInt::index(0);
  else
    return LaGenMatInt::index(1);
}

inline int LaVectorInt::start() const
{
  if (LaGenMatInt::size(1)==1 )
    return LaGenMatInt::start(0);
  else
    return LaGenMatInt::start(1);
}

inline int LaVectorInt::end() const
{
  if (LaGenMatInt::size(1)==1 )
    return LaGenMatInt::end(0);
  else
    return LaGenMatInt::end(1);
}

#endif 
// _LA_VECTOR_INT_H_
