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

//      Lapack++ "Shared" Vector Int Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_INT_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to int*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _VECTOR_INT_H_
#define _VECTOR_INT_H_    

#include <iostream>       // for formatted printing of matrices

typedef  struct {
    int sz;                                        
    int* data;                                       
    int ref_count;
} vrefInt;

class VectorInt
{                                                                      
private:                                                           
    vrefInt *p;
    int *data;			// performance hack, avoid indirection to data.
public:
				// Constructors/Destructors
    VectorInt(int);                             
    VectorInt(int, int);
    VectorInt(int*, int);
    VectorInt(const VectorInt&); 
    ~VectorInt() ;                              
				// Indices and access operations
    inline int& operator[](int);
    inline int operator[](int) const; // read only
    inline int& operator()(int); 
    inline int operator()(int) const; // read only
    inline operator int*(); 
    inline int size() const;
    inline bool null() const;
    int resize(int d);
    inline int ref_count() const; // return the number of ref counts
    inline int* addr();
    inline const int* addr() const;
				//  Assignment
    inline  VectorInt& operator=(const VectorInt&);
            VectorInt& operator=(int);
    inline  VectorInt& ref(const VectorInt &);
            VectorInt& inject(VectorInt&);
            VectorInt& copy(const VectorInt&);
				// I/O
    friend std::ostream&   operator<<(std::ostream&, const VectorInt&);       
};                                                                     
				// operators and member functions
inline bool VectorInt::null() const
{
    return (size() == 0) ;
}

inline int VectorInt::size() const
{
    return   p->sz;
}

inline int VectorInt::ref_count() const
{
    return p->ref_count;
}

inline const int* VectorInt::addr() const
{
    return data;
}

inline int* VectorInt::addr()
{
    return data;
}

inline VectorInt::operator int*() 
{
    return data;
}

inline int& VectorInt::operator()(int i)
{
#ifdef VECTOR_INT_BOUNDS_CHECK
    if (!(0<=i && i<size())) throw(LaException("assert failed : 0<=i && i<size()"));
#endif 
    return data[i];
}

inline int VectorInt::operator()(int i) const
{
#ifdef VECTOR_INT_BOUNDS_CHECK
    if (!(0<=i && i<size())) throw(LaException("assert failed : 0<=i && i<size()"));
#endif
    return data[i];
}

//  *CHANGE*  [] is the same as ()
inline int& VectorInt::operator[](int i)
{
#ifdef VECTOR_INT_BOUNDS_CHECK
    if (!(0<=i && i<size())) throw(LaException("assert failed : 0<=i && i<size()"));
#endif  
    return data[i];
}

//  *CHANGE*  [] is the same as ()
inline int VectorInt::operator[](int i) const
{
#ifdef VECTOR_INT_BOUNDS_CHECK
    if (!(0<=i && i<size())) throw(LaException("assert failed : 0<=i && i<size()"));
#endif  
    return data[i];
}

inline VectorInt& VectorInt::ref(const VectorInt& m)
{
    // always check that the p field has been initialized.
    // Either the lhs or rhs could be a NULL VectorInt...
    
    if (&m != this)		// not really necessary...
    {
	m.p->ref_count++;
	if (--(p->ref_count) == 0) // perform garbage col.
	{
	    delete [] ( p->data);
	    delete p;
	}
	p = m.p;
	data = m.data;
    }
    return *this;
}

inline VectorInt& VectorInt::operator=(const VectorInt& m)
{
    ref(m);
    return *this;
}

#endif 
// _VECTOR_INT_H_

