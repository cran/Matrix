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

#ifndef _LA_SPD_BAND_MAT_DOUBLE_H_
#define _LA_SPD_BAND_MAT_DOUBLE_H_

#include LA_GEN_MAT_DOUBLE_H
#include LA_VECTOR_DOUBLE_H

class LaSpdBandMatDouble
{
private:
    LaGenMatDouble data_;	// internal storage.

    int N_;			// N_ is (NxN)
    int kl_;			// kl_ = # subdiags
    static double outofbounds_; // out of range value returned.
    static int *info_;		// print matrix info only, not values
				//   originally 0, set to 1, and then
				//   reset to 0 after use.
public:
				// constructors
    inline LaSpdBandMatDouble();
    inline LaSpdBandMatDouble(int,int);
    inline LaSpdBandMatDouble(const LaSpdBandMatDouble &);
    
				// operators
    inline LaSpdBandMatDouble& operator=(double);
    inline LaSpdBandMatDouble& operator=(const LaSpdBandMatDouble&);
    double& operator()(int,int);
    const double& operator()(int,int) const;
    friend ostream& operator<<(ostream &, const LaSpdBandMatDouble &);

				// member functions
    inline int size(int) const;	// submatrix size
    inline int inc(int d) const; // explicit increment
    inline int gdim(int d) const; // global dimensions
    
    inline LaSpdBandMatDouble& ref(LaSpdBandMatDouble &);
    inline LaSpdBandMatDouble& copy(const LaSpdBandMatDouble &);
    inline double* addr() const { // return address of matrix.
	return data_.addr();}
    inline int ref_count() const { // return ref_count of matrix.
        return data_.ref_count();}
    inline LaIndex index(int d) const { // return indices of matrix.
        return data_.index(d);}
    inline int subdiags() {	// return # of subdiags of matrix.
        return (kl_);}
    inline int shallow() const { // return shallow flag.
        return data_.shallow();}

    inline LaSpdBandMatDouble& resize(const LaSpdBandMatDouble&);

    inline const LaSpdBandMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;};

    inline LaSpdBandMatDouble print_data() const 
	{ cout << data_; return *this;}
				// destructor
    inline ~LaSpdBandMatDouble();
};

				// constructors 
inline LaSpdBandMatDouble::LaSpdBandMatDouble()
    :data_()
{
    N_ = kl_ = 0;
}

inline LaSpdBandMatDouble::LaSpdBandMatDouble(int n,int l)
    :data_(n,2*l+1)
{

    N_ = n;
    kl_ = l;
}

inline LaSpdBandMatDouble::LaSpdBandMatDouble(const LaSpdBandMatDouble &A)
{
    data_.copy(A.data_);
    N_ = A.N_;
    kl_ = A.kl_;
}

				// destructor 
inline LaSpdBandMatDouble::~LaSpdBandMatDouble()
{
}

  
				// member functions and operators
inline LaSpdBandMatDouble& LaSpdBandMatDouble::ref(LaSpdBandMatDouble &ob)
{
    data_.ref(ob.data_);
    N_ = ob.N_;
    kl_ = ob.kl_;

    return *this;
}

inline LaSpdBandMatDouble& LaSpdBandMatDouble::resize(const LaSpdBandMatDouble &ob)
{
    data_.resize(ob.data_);

    return *this;
}


inline LaSpdBandMatDouble& LaSpdBandMatDouble::operator=(const LaSpdBandMatDouble &B)
{
    data_ = B.data_;
    N_ = B.N_;
    kl_ = B.kl_;

    return *this;
}

inline int LaSpdBandMatDouble::size(int d) const
{
    return(data_.size(d));
}

inline int LaSpdBandMatDouble::inc(int d) const
{
    return(data_.inc(d));
}

inline int LaSpdBandMatDouble::gdim(int d) const
{
    return(data_.gdim(d));
}

#endif 
// _LA_SPD_BAND_MAT_DOUBLE_H_

