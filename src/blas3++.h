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


#ifndef _BLAS3_PP_H_
#define _BLAS3_PP_H_

#include "blas3.h"


#ifdef _LA_GEN_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha = 1.0, double beta = 0.0);

void Blas_Mat_Trans_Mat_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha = 1.0, double beta = 0.0);

void Blas_Mat_Mat_Trans_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha = 1.0, double beta = 0.0);



#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_

void Blas_Mat_Mat_Solve(LaUnitLowerTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha = 1.0);

#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUnitUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

void Blas_Mat_Mat_Solve(LaUnitUpperTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha = 1.0);

#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

void Blas_Mat_Mat_Solve(LaLowerTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha = 1.0);
#endif


#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

void Blas_Mat_Mat_Solve(LaUpperTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha = 1.0);
#endif


#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaUnitLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaUnitUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);
#endif

#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUnitLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);
#endif


#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_ 
void Blas_Mat_Trans_Mat_Solve(LaUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaSymmMatDouble &A, LaGenMatDouble &B, 
            LaGenMatDouble &C, double alpha = 1.0, double beta = 0.0);

void Blas_R1_Update(LaSymmMatDouble &C, LaGenMatDouble &A,
            double alpha = 1.0, double beta = 0.0);
void Blas_R2_Update(LaSymmMatDouble &C, LaGenMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0, double beta = 0.0);
#endif



#endif
    /* _LA_GEN_MAT_DOUBLE_H_ */


#endif 
    // _BLAS3_PP_H_
            
