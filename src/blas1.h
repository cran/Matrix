// -*- c++ -*-
// $Id: blas1.h,v 1.4 2000/07/13 02:31:56 bates Exp $

// C++ prototypes for double precision Blas1 routines.

// Based on code labelled
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

// All arguments to Fortran routines must be passed as pointers or
// references.  Arguments that will be used as scalars are passed as
// references.  Arguments that will be used as arrays are passed as
// pointers.  The const keyword modifies those arguments that are
// guaranteed to be unchanged on exit.

// Modifications Copyright (C) 2000-2000 the R Development Core Team

#ifndef _BLAS1_H_
#define _BLAS1_H_

extern "C" {
#include <R_ext/RS.h>		// to define F77_NAME
    
    double F77_NAME(dasum)(const int& n,
			   const double *dx, const int& incx);
    
    void F77_NAME(daxpy)(const int& n, const double& alpha,
			 const double *dx, const int& incx,
			 double *dy, const int& incy);
    
    void F77_NAME(dcopy)(const int& n,
			 const double *dx, const int& incx,
			 double *dy, const int& incy);
    
    double F77_NAME(ddot)(const int& n,
			  const double *dx, const int& incx, 
			  const double *dy, const int& incy);
    
    double F77_NAME(dnrm2)(const int& n,
			   const double *dx, const int& incx); 
    
    void F77_NAME(drot)(const int& n,
			double *dx, const int& incx,
			double *dy, const int& incy,
			const double& c, const double& s);
    
    void F77_NAME(drotg)(const double& a, const double& b,
			 double& c, double& s);
    
    void F77_NAME(dscal)(const int& n, const double& alpha,
			 double *dx, const int& incx);
    
    void F77_NAME(dswap)(const int& n,
			 double *dx, const int& incx,
			 double *dy, const int& incy);
    
    int F77_NAME(idamax)(const int& n, const double *dx, const int& incx);
    
#if defined(LA_COMPLEX_SUPPORT)
    double F77_NAME(zdotc)(doublecomplex *c, const int& n, 
			   const doublecomplex *cx, 
			   const int& incx, const doublecomplex *cy,
			   const int& incy);
    
    double F77_NAME(zdotu)(doublecomplex *c, const int& n, 
			   const doublecomplex *cx, const int& incx, 
			   const doublecomplex *cy, const int& incy);
    
    void F77_NAME(zaxpy)(const int& n, const doublecomplex *da, 
			 const doublecomplex *dx, 
			 const int& incx, doublecomplex *dy, 
			 const int& incy);
    
    void F77_NAME(zcopy)(const int& n, doublecomplex *dx, const int& incx, 
			 doublecomplex *dy, const int& incy);
    
    double  F77_NAME(dzasum)(const int& n, const doublecomplex *dx,
			     const int& incx);
    
    double  F77_NAME(dznrm2)(const int& n, const doublecomplex *dx,
			     const int& incx); 
    
    void F77_NAME(zdscal)(const int& n, const double *da, doublecomplex *dx, 
			  const int& incx);
    
    void F77_NAME(zscal)(const int& n, const doublecomplex *da, doublecomplex *dx, 
			 const int& incx);
    
    int F77_NAME(izamax)(const int& n, const doublecomplex *dx, const int& incx);
    
    void F77_NAME(zswap)(const int& n, doublecomplex *dx, const int& incx, 
			 doublecomplex *dy, int& incy);
#endif
}
#endif

