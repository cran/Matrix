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

#ifndef _BLAS1_PP_H_
#define _BLAS1_PP_H_

#include "blas1.h"

#ifdef _LA_VECTOR_COMPLEX_H_
COMPLEX Blas_U_Dot_Prod(const LaVectorComplex &cx, const LaVectorComplex &cy);
COMPLEX Blas_H_Dot_Prod(const LaVectorComplex &cx, const LaVectorComplex &cy);
double Blas_Norm1(const LaVectorComplex &dx);
void Blas_Add_Mult(LaVectorComplex &dy, COMPLEX da, const LaVectorComplex &dx);
void Blas_Mult(LaVectorComplex &dy, COMPLEX da, const LaVectorComplex &dx);
void Blas_Copy(LaVectorComplex &dx, LaVectorComplex &dy);
double Blas_Norm1(const LaVectorComplex &dx);
double Blas_Norm2(const LaVectorComplex &dx);
void Blas_Scale(COMPLEX da, LaVectorComplex &dx);
void Blas_Swap(LaVectorComplex &dx, LaVectorComplex &dy);
int Blas_Index_Max(const LaVectorComplex &dx);
#endif


#ifdef _LA_VECTOR_DOUBLE_H_
double Blas_Norm1(const LaVectorDouble &dx);
void Blas_Add_Mult(LaVectorDouble &dy, double da, const LaVectorDouble &dx);
void Blas_Mult(LaVectorDouble &dy, double da, LaVectorDouble &dx);
void Blas_Copy(LaVectorDouble &dx, LaVectorDouble &dy);
double Blas_Dot_Prod(const LaVectorDouble &dx, const LaVectorDouble &dy);
double Blas_Norm2(const LaVectorDouble &dx);
void Blas_Apply_Plane_Rot(LaVectorDouble &dx, LaVectorDouble &dy, 
			  double c, double s);
void Blas_Gen_Plane_Rot(double da, double db, double &c, double &s);
void Blas_Scale(double da, LaVectorDouble &dx);
void Blas_Swap(LaVectorDouble &dx, LaVectorDouble &dy);
int Blas_Index_Max(const LaVectorDouble &dx);
#endif // _LA_VECTOR_DOUBLE_H_
#endif // _BLAS1_PP_H_
