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


#ifndef _LA_COL_VECTOR_DOUBLE_H_
#define _LA_COL_VECTOR_DOUBLE_H_

#include "lafnames.h"

#ifndef _LA_GEN_MAT_DOUBLE_H_
#include LA_GEN_MAT_DOUBLE_H
#endif


// a column vector is simply an nx1 matrix, only that it can
// be constructed and accessed by a single dimension

class LaColVectorDouble: public LaGenMatDouble
{
public:

    inline LaColVectorDouble();
    inline LaColVectorDouble(int);
    inline LaColVectorDouble(double*, int);
    inline LaColVectorDouble(const LaGenMatDouble&);

    inline int size() const;
    inline int inc() const;
    inline LaIndex index() const;
    inline int start() const;
    inline int end() const;

    inline LaColVectorDouble& ref(const LaGenMatDouble &);
    LaColVectorDouble& ref(SEXP);
    inline LaColVectorDouble& inject(const LaMatDouble &);
    inline LaColVectorDouble& copy(const LaMatDouble &);
    
    inline double& operator()(int i);
    inline const double& operator()(int i) const ;
    inline LaColVectorDouble operator()(const LaIndex&);

    inline LaColVectorDouble& operator=(const LaMatDouble &A);
    inline LaColVectorDouble& operator=(double d);
    
};


inline LaColVectorDouble::LaColVectorDouble() : LaGenMatDouble(0,1) {}
inline LaColVectorDouble::LaColVectorDouble(int i) : LaGenMatDouble(i,1) {}


inline LaColVectorDouble::LaColVectorDouble(double *d, int m) : 
    LaGenMatDouble(d,m,1) {}


inline LaColVectorDouble::LaColVectorDouble(const LaGenMatDouble& G) : 
        LaGenMatDouble(G)
{
        assert(G.size(1)==1);

}
        
//note that vectors can be either stored columnwise, or row-wise

// this will handle the 0x0 case as well.

inline int LaColVectorDouble::size() const 
{ return LaGenMatDouble::size(0)*LaGenMatDouble::size(1); }

inline double& LaColVectorDouble::operator()(int i)
{ 
    return LaGenMatDouble::operator()(i,0);
}

inline const double& LaColVectorDouble::operator()(int i) const
{ 
    return LaGenMatDouble::operator()(i,0);
}

inline LaColVectorDouble LaColVectorDouble::operator()(const LaIndex& I)
{ 
    return LaGenMatDouble::operator()(I,LaIndex(0,0)); 
}

inline LaColVectorDouble& LaColVectorDouble::operator=(const LaMatDouble &A)
{
    LaGenMatDouble::copy(A);
    return *this;
}

inline LaColVectorDouble& LaColVectorDouble::operator=(double d)
{
    LaGenMatDouble::operator=(d);
    return *this;
}

inline LaColVectorDouble& LaColVectorDouble::copy(const LaMatDouble &A)
{
    assert(A.size(1) == 1);   //make sure A is a column vector.
    LaGenMatDouble::copy(A);
    return *this;
}


inline LaColVectorDouble& LaColVectorDouble::ref(const LaGenMatDouble &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatDouble::ref(A);
    return *this;
}

inline LaColVectorDouble& LaColVectorDouble::inject(const LaMatDouble &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatDouble::inject(A);
    return *this;
}

inline int LaColVectorDouble::inc() const
{
    return LaGenMatDouble::inc(0);
}

inline LaIndex LaColVectorDouble::index() const
{
    return LaGenMatDouble::index(0);
}

inline int LaColVectorDouble::start() const
{
    return LaGenMatDouble::start(0);
}

inline int LaColVectorDouble::end() const
{
    return LaGenMatDouble::end(0);
}

#endif 
    // _LA_COL_VECTOR_DOUBLE_H_
