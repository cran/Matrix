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

#ifndef _LA_INDEX_H_
#define _LA_INDEX_H_

// matrix index class.  Note that we name it "LaIndex" to avoid confusion
// with the "index()" string function in C, or other generic Index()
// functions.

#undef Free
#include <iostream>


class LaIndex
{
protected:
    int start_;
    int inc_;
    int end_;
public:
    inline LaIndex()
        : start_(0), inc_(0), end_(0) {}
    inline LaIndex(int i1)
        : start_(i1), inc_(1), end_(i1) {}
    inline LaIndex(int i1, int i2)
        : start_(i1), inc_(1), end_(i2) {}
    inline LaIndex(int i1, int i2, int i3)
        : start_(i1), inc_(i3), end_(i2) {}
    inline LaIndex(const LaIndex &s)
        : start_(s.start_), inc_(s.inc_), end_(s.end_) {}
    ~LaIndex() { }

// ** must have multply defined start(), inc() and end() member functions
// ** for both const and non-const objects because compiler complains in 
// ** LaVector*.h, for example, about assignment to const member. 
// ** (LaVector*.h line 112, 113, 114)

    inline int& start() { return start_;}
    inline const int& start() const { return start_;}
    inline int& inc() { return inc_;}
    inline const int& inc() const { return inc_;}
    inline int& end() { return end_;}
    inline const int& end() const { return end_;}
    inline int length() const { return ((end()-start())/inc() + 1);}
    inline const int null() const{ return (start() == 0 && 
					   inc() == 0  && end() == 0);}
    inline LaIndex& operator()(int i1, int i2){
	start_=i1; inc_=1; end_=i2; return *this;}
    inline LaIndex operator+(int i)
    {
        start_+=i; end_+=i;
        return *this;
    }
    inline LaIndex& operator=(const LaIndex& i){
	start_=i.start_; inc_=1; end_=i.end_; 
	return *this;}
};

inline std::ostream& operator<<(std::ostream& s, const LaIndex i)
{
    s << "(" << i.start() << ":" << i.inc() << ":" << i.end() << ")";
    
    return s;
}

#endif //  LA_INDEX_H_

