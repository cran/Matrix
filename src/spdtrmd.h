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
// Modifications Copyright (C) 2000-2000, 2002 the R Development Core Team

#ifndef _LA_SPD_TRIDIAG_MAT_DOUBLE_H_
#define _LA_SPD_TRIDIAG_MAT_DOUBLE_H_

#include "lafnames.h"
#include LA_VECTOR_DOUBLE_H

class LaSpdTridiagMatDouble
{   
        int size_;
        LaVectorDouble d_;  /* main diag */
        LaVectorDouble dl_; /* lower diag */
        static double outofbounds_; /* return this address, when addresing 
                                        out of bounds */
        static int *info_;        // print matrix info only, not values
                                  //   originally 0, set to 1, and then
                                  //   reset to 0 after use.

    public:

        // constructors / destructor
    
        inline LaSpdTridiagMatDouble();
        inline LaSpdTridiagMatDouble(int N);
        inline LaSpdTridiagMatDouble(const LaSpdTridiagMatDouble &);
        inline ~LaSpdTridiagMatDouble();

        // operators and member functions

            double& operator()(int i, int j);
            double operator()(int i, int j) const;
            LaVectorDouble diag(int); /* 0 main, -1 lower, 1 upper */
        inline LaSpdTridiagMatDouble& ref(LaSpdTridiagMatDouble&); 
        inline LaSpdTridiagMatDouble& copy(const LaSpdTridiagMatDouble&); 
        inline const LaSpdTridiagMatDouble& info() const {
            int *t = info_;
            *t = 1;
            return *this;};
        inline int size() const; 

    friend std::ostream& operator<<(std::ostream&,const LaSpdTridiagMatDouble&);


};

    // constructors

inline LaSpdTridiagMatDouble::LaSpdTridiagMatDouble()
    : size_(0), d_(), dl_()
{
}

inline LaSpdTridiagMatDouble::LaSpdTridiagMatDouble(int N)
    : size_(N), d_(N), dl_(N-1)
{
}
    
inline LaSpdTridiagMatDouble::LaSpdTridiagMatDouble(const LaSpdTridiagMatDouble& td)
    : size_(td.size_), d_(td.d_), dl_(td.dl_)
{
}

    // destructor

inline LaSpdTridiagMatDouble::~LaSpdTridiagMatDouble()
{
}


    // operators and member functions





inline LaSpdTridiagMatDouble& LaSpdTridiagMatDouble::ref(LaSpdTridiagMatDouble&T) 
{
    d_.ref(T.d_);
    dl_.ref(T.dl_); 
    size_ = T.size_;

    return *this;
}


inline LaSpdTridiagMatDouble& LaSpdTridiagMatDouble::copy(const LaSpdTridiagMatDouble&T) 
{
    d_.copy(T.d_);
    dl_.copy(T.dl_);    
    size_ = T.size_;

    return *this;
}

inline int LaSpdTridiagMatDouble::size() const
{
    return size_;
}



#endif 
// _LA_SPD_TRIDIAG_MAT_DOUBLE_H_
