//
//              LAPACK++ 1.1a Linear Algebra Package 1.1a
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

#include "lafnames.h"
#include LA_VECTOR_DOUBLE_H
#include "blas1++.h"

double Blas_Norm1(const LaVectorDouble &dx)
{
    return F77_CALL(dasum)(dx.size(), &dx(0), dx.inc());
}

void Blas_Add_Mult(LaVectorDouble &dy, double da, const LaVectorDouble &dx) 
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    F77_CALL(daxpy)(dx.size(), da, &dx(0), dx.inc(), &dy(0), dy.inc());
}

void Blas_Mult(LaVectorDouble &dy, double da, LaVectorDouble &dx)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    dy = 0.0;
    F77_CALL(daxpy)(dx.size(), da, &dx(0), dx.inc(), &dy(0), dy.inc());
}

void Blas_Copy(LaVectorDouble &dy, LaVectorDouble &dx)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    F77_CALL(dcopy)(dx.size(), &dx(0), dx.inc(), &dy(0), dy.inc());
}

double Blas_Dot_Prod(const LaVectorDouble &dx, const LaVectorDouble &dy)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    return F77_CALL(ddot)(dx.size(), &dx(0), dx.inc(), &dy(0), dy.inc());
}

double Blas_Norm2(const LaVectorDouble &dx)
{
    return F77_CALL(dnrm2)(dx.size(), &dx(0), dx.inc());
}

void Blas_Apply_Plane_Rot(LaVectorDouble &dx, LaVectorDouble &dy, 
			  double c, double s)
{
    if (!(dx.size() == dy.size())) throw(LaException("assert failed : dx.size() == dy.size()"));
    F77_CALL(drot)(dx.size(), &dx(0), dx.inc(), &dy(0), dy.inc(), c, s);
}

void Blas_Gen_Plane_Rot(double da, double db, double &c, double &s)
{
    F77_CALL(drotg)(da, db, c, s);
}

void Blas_Scale(double da, LaVectorDouble &dx)
{
    F77_CALL(dscal)(dx.size(), da, &dx(0), dx.inc());
}

void Blas_Swap(LaVectorDouble &dx, LaVectorDouble &dy)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    F77_CALL(dswap)(dx.size(), &dx(0), dx.inc(), &dy(0), dy.inc());
}

int Blas_Index_Max(const LaVectorDouble &dx)
{
    // subtract one from index since f77 starts at 1, not 0.
    return F77_CALL(idamax)(dx.size(), &dx(0), dx.inc()) - 1;
}

// Complex Routines
#ifdef LA_COMPLEX_SUPPORT
#include LA_VECTOR_COMPLEX_H

COMPLEX Blas_H_Dot_Prod(const LaVectorComplex &cx, const LaVectorComplex &cy)
{
    if (!(cx.size()==cy.size())) throw(LaException("assert failed : cx.size()==cy.size()"));
    int n = cx.size();
    int incx = cx.inc(), incy = cy.inc();
    COMPLEX tmp;

    F77_CALL(zdotc)(&tmp, &n, &cx(0), &incx, &cy(0), &incy);
    return tmp;
}

COMPLEX Blas_U_Dot_Prod(const LaVectorComplex &cx, const LaVectorComplex &cy)
{
    if (!(cx.size()==cy.size())) throw(LaException("assert failed : cx.size()==cy.size()"));
    int n = cx.size();
    int incx = cx.inc(), incy = cy.inc();
    COMPLEX tmp;

    F77_CALL(zdotu)(&tmp, &n, &cx(0), &incx, &cy(0), &incy);
    return tmp;
}


double Blas_Norm1(const LaVectorComplex &dx)
{
    int n = dx.size();
    int incx = dx.inc();

    return F77_CALL(dzasum)(&n, &dx(0), &incx);
}


void Blas_Add_Mult(LaVectorComplex &dy, COMPLEX da, const LaVectorComplex &dx) 
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    int n = dx.size();
    int incx = dx.inc(), incy = dy.inc();

    F77_CALL(zaxpy)(&n, &da, &dx(0), &incx, &dy(0), &incy);
}

void Blas_Mult(LaVectorComplex &dy, COMPLEX da, LaVectorComplex &dx)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    int n = dx.size();
    int incx = dx.inc(), incy = dy.inc();
    dy = COMPLEX(0.0,0.0);

    F77_CALL(zaxpy)(&n, &da, &dx(0), &incx, &dy(0), &incy);
}


void Blas_Copy(LaVectorComplex &dy, LaVectorComplex &dx)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    int n = dx.size();
    int incx = dx.inc(), incy = dy.inc();

    F77_CALL(zcopy)(&n, &dx(0), &incx, &dy(0), &incy);
}




double Blas_Norm2(const LaVectorComplex &dx)
{
    int n = dx.size();
    int incx = dx.inc();

    return F77_CALL(dznrm2)(&n, &dx(0), &incx);
}



void Blas_Scale(COMPLEX za, LaVectorComplex &dx)
{
    int n = dx.size();
    int incx = dx.inc();

    F77_CALL(zscal)(&n, &za, &dx(0), &incx);
}


void Blas_Swap(LaVectorComplex &dx, LaVectorComplex &dy)
{
    if (!(dx.size()==dy.size())) throw(LaException("assert failed : dx.size()==dy.size()"));
    int n = dx.size();
    int incx = dx.inc(), incy = dy.inc();

    F77_CALL(zswap)(&n, &dx(0), &incx, &dy(0), &incy);
}



int Blas_Index_Max(const LaVectorComplex &dx)
{
    int n = dx.size();
    int incx = dx.inc();

    // subtract one from index since f77 starts at 1, not 0.
    return F77_CALL(izamax)(&n, &dx(0), &incx) - 1;
}

#endif
// LA_COMPLEX_SUPPORT
