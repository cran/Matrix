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

#ifndef _LA_ROW_VECTOR_DOUBLE_H_
#define _LA_ROW_VECTOR_DOUBLE_H_

#include "lafnames.h"

#ifndef _LA_GEN_MAT_DOUBLE_H_
#include LA_GEN_MAT_DOUBLE_H
#endif


// a row vector is simply an nx1 matrix, only that it can
// be constructed and accessed by a single dimension

class LaRowVectorDouble: public LaGenMatDouble
{
public:

    inline LaRowVectorDouble();
    inline LaRowVectorDouble(int);
    inline LaRowVectorDouble(double*, int);
    inline LaRowVectorDouble(const LaGenMatDouble&);

    inline int size() const;
    inline int inc() const;
    inline LaIndex index() const;
    inline int start() const;
    inline int end() const;

    inline LaRowVectorDouble& ref(const LaGenMatDouble &);
    LaRowVectorDouble& ref(SEXP);
    inline LaRowVectorDouble& inject(const LaMatDouble &);
    inline LaRowVectorDouble& copy(const LaMatDouble &);
    inline LaRowVectorDouble* clone() const;
    
    inline double& operator()(int i);
    inline const double& operator()(int i) const ;
    inline LaRowVectorDouble operator()(const LaIndex&);

    inline LaRowVectorDouble& operator=(const LaMatDouble &A);
    inline LaRowVectorDouble& operator=(double d);
    
};


inline LaRowVectorDouble::LaRowVectorDouble() : LaGenMatDouble(1,0) {}
inline LaRowVectorDouble::LaRowVectorDouble(int i) : LaGenMatDouble(1, i) {}


inline LaRowVectorDouble::LaRowVectorDouble(double *d, int m) : 
    LaGenMatDouble(d,1,m) {}


inline LaRowVectorDouble::LaRowVectorDouble(const LaGenMatDouble& G) : 
        LaGenMatDouble(G)
{
        if (!(G.size(0)==1)) throw(LaException("assert failed : G.size(0)==1"));

}
        

inline int LaRowVectorDouble::size() const 
{ return LaGenMatDouble::size(0)*LaGenMatDouble::size(1); }

inline double& LaRowVectorDouble::operator()(int i)
{ 
    return LaGenMatDouble::operator()(0,i);
}

inline const double& LaRowVectorDouble::operator()(int i) const
{ 
    return LaGenMatDouble::operator()(0,i);
}

inline LaRowVectorDouble LaRowVectorDouble::operator()(const LaIndex& I)
{ 
    return LaGenMatDouble::operator()(LaIndex(0,0),I); 
}

inline LaRowVectorDouble& LaRowVectorDouble::operator=(const LaMatDouble &A)
{
    LaGenMatDouble::copy(A);
    return *this;
}

inline LaRowVectorDouble& LaRowVectorDouble::operator=(double d)
{
    LaGenMatDouble::operator=(d);
    return *this;
}

inline LaRowVectorDouble& LaRowVectorDouble::copy(const LaMatDouble &A)
{
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));   //make sure rhs is a
                                                // a vector.
    LaGenMatDouble::copy(A);
    return *this;
}

inline LaRowVectorDouble* LaRowVectorDouble::clone() const
{
    LaGenMatDouble* tmp = LaGenMatDouble::clone();
    LaRowVectorDouble* ans = new LaRowVectorDouble();
    ans->ref(*tmp);
    delete tmp;
    return ans;
}

inline LaRowVectorDouble& LaRowVectorDouble::ref(const LaGenMatDouble &A)
{
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));
    LaGenMatDouble::ref(A);
    return *this;
}


inline LaRowVectorDouble& LaRowVectorDouble::inject(const LaMatDouble &A)
{
    if (!(A.size(0) == 1 || A.size(1) == 1)) throw(LaException("assert failed : A.size(0) == 1 || A.size(1) == 1"));
    LaGenMatDouble::inject(A);
    return *this;
}

inline int LaRowVectorDouble::inc() const
{
    return LaGenMatDouble::inc(1);
}

inline LaIndex LaRowVectorDouble::index() const
{
    return LaGenMatDouble::index(1);
}

inline int LaRowVectorDouble::start() const
{
    return LaGenMatDouble::start(1);
}

inline int LaRowVectorDouble::end() const
{
    return LaGenMatDouble::end(1);
}

#endif 
    // _LA_ROW_VECTOR_DOUBLE_H_
