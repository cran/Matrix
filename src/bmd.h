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


#ifndef _LA_BAND_MAT_DOUBLE_H_
#define _LA_BAND_MAT_DOUBLE_H_

#ifndef _LA_GEN_MAT_DOUBLE_H_
#include LA_GEN_MAT_DOUBLE_H
#endif


#define BOUNDS_CHK
//#define SPARSE_CHK
#ifdef LA_NO_BOUNDS_CHECK
#undef BOUNDS_CHK
#endif
#ifdef LA_NO_SPARSE_CHECK
#undef SPARSE_CHK
#endif

class LaBandMatDouble
{
  LaGenMatDouble data_;		// internal storage.

  int N_;			// N_ is (NxN)
  int kl_;			// kl_ = # subdiags
  int ku_;			// ku_ = # superdiags
  static double outofbounds_;	// value returned if index is out of range.
  static int *info_;		// print matrix info only, not values
				//   originally 0, set to 1, and then
				//   reset to 0 after use.
public:

  // constructors

  inline LaBandMatDouble();
  inline LaBandMatDouble(int,int,int);
  inline LaBandMatDouble(const LaBandMatDouble &);

  friend std::ostream& operator<<(std::ostream &, const LaBandMatDouble &);


  // member functions

  inline double* addr() {	// return address of matrix.
        return data_.addr();}
  inline const double* addr() const {	// return address of matrix.
        return data_.addr();}

  double& operator()(int,int);
  double operator()(int,int) const;
  inline LaBandMatDouble& operator=(const LaBandMatDouble&);
  LaBandMatDouble& operator=(double);

  inline LaBandMatDouble& resize(const LaBandMatDouble&);

  inline int size(int) const;	// submatrix size
  inline int inc(int d) const;	// explicit increment
  inline int gdim(int d) const;	// global dimensions

  inline LaBandMatDouble& ref(LaBandMatDouble &);
        LaBandMatDouble copy(const LaBandMatDouble &);
  inline int ref_count() const { // return ref_count of matrix.
        return data_.ref_count();}
  inline LaIndex index(int d) const { // return indices of matrix.
        return data_.index(d);}
  inline int superdiags() {     // return # of superdiags of matrix.
        return (ku_);}
  inline int superdiags() const { // return # of superdiags of const matrix.
        return (ku_);}
  inline int subdiags() {     // return # of subdiags of matrix.
        return (kl_);}
  inline int subdiags() const {  // return # of subdiags of const matrix.
        return (kl_);}
  inline int shallow() const {	// return shallow flag.
        return data_.shallow();}

  inline const LaBandMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;};

  inline LaBandMatDouble print_data() const 
    { std::cout << data_; return *this;}

  // destructor

  inline ~LaBandMatDouble();
};

  // constructors 

inline LaBandMatDouble::LaBandMatDouble() 
    : data_(), N_(0), kl_(0), ku_(0)
{
}

inline LaBandMatDouble::LaBandMatDouble(int n,int l,int u)
    : data_(2*l+u+1,n), N_(n), kl_(l), ku_(u)
{
}

inline LaBandMatDouble::LaBandMatDouble(const LaBandMatDouble &A)
    : data_(), N_(A.N_), kl_(A.kl_), ku_(A.ku_)
{
}

  // destructor 

inline LaBandMatDouble::~LaBandMatDouble()
{
}

  
  // member functions and operators

inline LaBandMatDouble& LaBandMatDouble::ref(LaBandMatDouble &ob)
{

  data_.ref(ob.data_);
  N_ = ob.N_;
  kl_ = ob.kl_;
  ku_ = ob.ku_;

  return *this;
}

inline LaBandMatDouble& LaBandMatDouble::resize(const LaBandMatDouble &ob)
{

  data_.resize(ob.data_);

  return *this;
}


inline LaBandMatDouble& LaBandMatDouble::operator=(const LaBandMatDouble &B)
{
    data_ = B.data_;
    N_ = B.N_;
    kl_ = B.kl_;
    ku_ = B.ku_;

    return *this;
}


inline int LaBandMatDouble::size(int d) const
{
   return(data_.size(d));
}

inline int LaBandMatDouble::inc(int d) const
{
   return(data_.inc(d));
}

inline int LaBandMatDouble::gdim(int d) const
{
   return(data_.gdim(d));
}

#endif 
    // _LA_BAND_MAT_DOUBLE_H_
