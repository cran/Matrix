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

//      Lapack++ "Shared" Vector Double Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_DOUBLE_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to double*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _VECTOR_DOUBLE_H_
#define _VECTOR_DOUBLE_H_    

#include <iostream>       // for formatted printing of matrices

typedef  struct {
    int        sz;                                        
    double*    data;                                       
    int        ref_count;
} vrefDouble;

class VectorDouble
{                                                                      
 private:                                                           
    vrefDouble *p;
    double *data;		// performance hack, avoid double
				// indirection to data.
 public:                                                            
				// Constructors/Destructors
    //inline VectorDouble();    // this should behave as VectorDouble(0)
    VectorDouble(int);                             
    VectorDouble(int, double);	// can't be inlined because of 'for'
				// statement.
    VectorDouble(double*, int);
    VectorDouble(const VectorDouble&); 
    ~VectorDouble();
				//  Indices and access operations
    inline double&      operator[](int); 
    inline double       operator[](int) const; // read only
    inline double&      operator()(int); 
    inline double       operator()(int) const; // read only
    inline              operator    double*(); 
    inline int          size() const;
    inline bool         null() const;
    int                 resize(int d);
    inline int          ref_count() const;  // return the number of ref counts
    inline double*      addr();
    inline const double* addr() const;

				//  Assignment
    inline  VectorDouble& operator=(const VectorDouble&);
    VectorDouble& operator=(double);
    inline  VectorDouble& ref(const VectorDouble &);
    VectorDouble& inject(const VectorDouble&);
    VectorDouble& copy(const VectorDouble&);
    
				// I/O
    friend std::ostream&   operator<<(std::ostream&, const VectorDouble&);
};                                                                     

				// operators and member functions
inline bool VectorDouble::null() const
{
    return (size() == 0) ;
}

inline int VectorDouble::size() const
{
    return   p-> sz;
}

inline int VectorDouble::ref_count() const
{
    return p->ref_count;
}

inline const double* VectorDouble::addr() const
{
    return data;
}

inline double* VectorDouble::addr()
{
    return data;
}

inline VectorDouble::operator double*() 
{
    return data;
}

inline double& VectorDouble::operator()(int i)
{
#ifdef VECTOR_DOUBLE_BOUNDS_CHECK
    if (!(0<=i && i<size())) throw(LaException("assert failed : 0<=i && i<size()"));
#endif 
    return data[i];
}

inline double VectorDouble::operator()(int i) const
{
#ifdef VECTOR_DOUBLE_BOUNDS_CHECK
    if (!(0<=i && i<size())) throw(LaException("assert failed : 0<=i && i<size()"));
#endif
    return data[i];
}

//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()
inline double& VectorDouble::operator[](int i)
{
#ifdef VECTOR_DOUBLE_BOUNDS_CHECK
    if (!(0<=i && i<size())) throw(LaException("assert failed : 0<=i && i<size()"));
#endif  
    return data[i];
}

inline double VectorDouble::operator[](int i) const
{
#ifdef VECTOR_DOUBLE_BOUNDS_CHECK
    if (!(0<=i && i<size())) throw(LaException("assert failed : 0<=i && i<size()"));
#endif  
    return data[i];
}

inline VectorDouble& VectorDouble::ref(const VectorDouble& m)
{
				// always check that the p field has
				// been initialized. 
				// Either the lhs or rhs could be a
				// NULL VectorDouble...
    if (&m != this) {		// not really necessary...
	m.p->ref_count++;
	if (--(p->ref_count) == 0) { // perform garbage collection
	    delete [] ( p->data);
	    delete p;
	}
	p = m.p;
	data = m.data;
    }
    return *this;
}

inline VectorDouble& VectorDouble::operator=(const VectorDouble& m)
{
    ref(m);
    return *this;
}

#endif 
// _VECTOR_DOUBLE_H_

