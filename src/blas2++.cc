//
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


#include "lafnames.h"
#include LA_VECTOR_DOUBLE_H
#include LA_SYMM_MAT_DOUBLE_H
#include LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H
#include LA_UPPER_TRIANG_MAT_DOUBLE_H
#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_SPD_MAT_DOUBLE_H
//#include LA_SYMM_BAND_MAT_DOUBLE_H
#include LA_TRIDIAG_MAT_DOUBLE_H

#include "blas2++.h"
void Blas_Mat_Trans_Vec_Mult(LaGenMatDouble &A, LaVectorDouble &dx,
			     LaVectorDouble &dy, double alpha, double beta)
{  
    F77_CALL(dgemv)('T', A.size(0), A.size(1), alpha, &A(0,0),
		   A.gdim(0), &dx(0), dx.inc(), beta, &dy(0),
		   dy.inc());
}

void Blas_Mat_Vec_Mult(LaGenMatDouble &A, LaVectorDouble &dx,
		       LaVectorDouble &dy, double alpha , double beta )
{
    F77_CALL(dgemv)('N', A.size(0), A.size(1), alpha, &A(0,0),
		   A.gdim(0), &dx(0), dx.inc(), beta, &dy(0), dy.inc());
}

//  void Blas_Mat_Vec_Mult(LaSymmBandMatDouble &A, LaVectorDouble &dx, 
//  		       LaVectorDouble &dy, double alpha , double beta )
//  {
//      F77_CALL(dsbmv)('L', A.size(1), A.subdiags(), alpha, &A(0,0),
//  		   A.gdim(0), &dx(0), dx.inc(), beta, &dy(0), dy.inc());
//  }


void Blas_Mat_Vec_Mult(LaSpdMatDouble &AP, LaVectorDouble &dx, 
		       LaVectorDouble &dy, double alpha , double beta )
{
    F77_CALL(dspmv)('L', AP.size(1), alpha, &AP(0,0), &dx(0),
		    dx.inc(), beta, &dy(0), dy.inc());
}

void Blas_Mat_Vec_Mult(LaLowerTriangMatDouble &A, LaVectorDouble &dx)
{
    F77_CALL(dtrmv)('L', 'N', 'N', A.size(1), &A(0,0), A.gdim(0),
		   &dx(0), dx.inc());
}

void Blas_Mat_Vec_Mult(LaUpperTriangMatDouble &A, LaVectorDouble &dx)
{
    F77_CALL(dtrmv)('U', 'N', 'N', A.size(1), &A(0,0), A.gdim(0),
		   &dx(0), dx.inc()); 
}

void Blas_Mat_Vec_Mult(LaUnitLowerTriangMatDouble &A, LaVectorDouble &dx)
{
    F77_CALL(dtrmv)('L', 'N', 'U', A.size(1), &A(0,0), A.gdim(0),
		   &dx(0), dx.inc()); 
}

void Blas_Mat_Vec_Mult(LaUnitUpperTriangMatDouble &A, LaVectorDouble &dx)
{
    F77_CALL(dtrmv)('U', 'N', 'U', A.size(1), &A(0,0), A.gdim(0),
		   &dx(0), dx.inc());
}

void Blas_Mat_Vec_Solve(LaLowerTriangMatDouble &A, LaVectorDouble &dx)
{
    F77_CALL(dtrsv)('L', 'N', 'N', A.size(1), &A(0,0), A.gdim(0),
		    &dx(0), dx.inc());
}

void Blas_Mat_Vec_Solve(LaUpperTriangMatDouble &A, LaVectorDouble &dx)
{
    F77_CALL(dtrsv)('U', 'N', 'N', A.size(1), &A(0,0), A.gdim(0),
		    &dx(0), dx.inc()); 
}


void Blas_Mat_Vec_Solve(LaUnitLowerTriangMatDouble &A, LaVectorDouble &dx)
{
    F77_CALL(dtrsv)('L', 'N', 'U', A.size(1), &A(0,0), A.gdim(0),
		    &dx(0), dx.inc());
}

void Blas_Mat_Vec_Solve(LaUnitUpperTriangMatDouble &A, LaVectorDouble &dx)
{
    F77_CALL(dtrsv)('U', 'N', 'U', A.size(1), &A(0,0), A.gdim(0),
		    &dx(0), dx.inc());
}

void Blas_R1_Update(LaGenMatDouble &A, LaVectorDouble &dx, 
		    LaVectorDouble &dy, double alpha)
{
    F77_CALL(dger)(A.size(0), A.size(1), alpha, &dx(0), dx.inc(),
		   &dy(0), dy.inc(), &A(0,0), A.gdim(0));
}

void Blas_R1_Update(LaSymmMatDouble &A, LaVectorDouble &dx, double alpha)
{
    F77_CALL(dsyr)('L', A.size(1), alpha, &dx(0), dx.inc(), &A(0,0),
		   A.gdim(0));
}

void Blas_R1_Update(LaSpdMatDouble &AP, LaVectorDouble &dx, double alpha )
{
    F77_CALL(dspr)('L', AP.size(1), alpha, &dx(0), dx.inc(), &AP(0,0));
}

void Blas_R2_Update(LaSymmMatDouble &A, LaVectorDouble &dx, 
                LaVectorDouble &dy, double alpha )
{
    F77_CALL(dsyr2)('L', A.size(1), alpha, &dx(0), dx.inc(), &dy(0),
		    dy.inc(), &A(0,0), A.gdim(0));
}

void Blas_R2_Update(LaSpdMatDouble &AP, LaVectorDouble &dx, 
		    LaVectorDouble &dy, double alpha )
{
    F77_CALL(dspr2)('L', AP.size(1), alpha, &dx(0), dx.inc(), &dy(0),
		    dy.inc(), &AP(0,0));
}

