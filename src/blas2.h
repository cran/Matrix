// -*- c++ -*-
// $Id: blas2.h,v 1.5 2000/07/12 22:54:57 bates Exp $

// C++ prototypes for double precision Blas2 routines.

// Based on code labelled
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

// All arguments to Fortran routines must be passed as pointers or
// references.  Arguments that will be used as scalars are passed as
// references.  Arguments that will be used as arrays are passed as
// pointers.  The const keyword modifies those arguments that are
// guaranteed to be unchanged on exit.

// Modifications Copyright (C) 2000-2000 the R Development Core Team

#ifndef _BLAS2_H_
#define _BLAS2_H_

extern "C" {
#include <R_ext/RS.h>		// to define F77_NAME
    
    // DGEMV - perform one of the matrix-vector operations
    // y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y, 
    void F77_NAME(dgemv)(const char& trans, const int& m, const int& n,
			 const double& alpha,
			 const double *a, const int& lda,
			 const double *x, const int& incx,
			 const double& beta,
			 double *y, const int& incy); 

    // DGBMV - perform one of the matrix-vector operations
    // y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,
    void F77_NAME(dgbmv)(const char& trans, const int& m, const int& n,
			 const int& kl,const int& ku,
			 const double& alpha,
			 const double *a, const int& lda,
			 const double *x, const int& incx,
			 const double& beta,
			 double *y, const int& incy);
    
    // DSBMV - perform the matrix-vector operation
    // y := alpha*A*x + beta*y,
    void F77_NAME(dsbmv)(const char& uplo, const int& n, const int& k,
			 const double& alpha,
			 const double *a, const int& lda,
			 const double *x, const int& incx,
			 const double& beta,
			 double *y, const int& incy);

    // DSPMV - perform the matrix-vector operation
    // y := alpha*A*x + beta*y,
    void F77_NAME(dspmv)(const char& uplo, const int& n,
			 const double& alpha,
			 const double *ap,
			 const double *x, const int& incx,
			 const double& beta,
			 double *y, const int& incy);

    // DSYMV - perform the matrix-vector operation
    //  y := alpha*A*x + beta*y,
    void F77_NAME(dsymv)(const char& uplo, const int& n,
			 const double& alpha,
			 const double *a, const int& lda,
			 const double *x, const int& incx,
			 const double& beta, double *y, const int& incy);

    // DTBMV - perform one of the matrix-vector operations
    // x := A*x, or x := A'*x,
    void F77_NAME(dtbmv)(const char& uplo, const char& trans,
			 const char& diag, const int& n, const int& k,
			 const double *a, const int& lda,
			 double *x, const int& incx); 

    // DTPMV - perform one of the matrix-vector operations
    // x := A*x, or x := A'*x,
    void F77_NAME(dtpmv)(const char& uplo, const char& trans,
			 const char& diag, const int& n,
			 const double *ap,
			 double *x, const int& incx);

    // DTRMV - perform one of the matrix-vector operations 
    // x := A*x, or x := A'*x,
    void F77_NAME(dtrmv)(const char& uplo, const char& trans,
			 const char& diag, const int& n,
			 const double *a, const int& lda,
			 double *x, const int& incx);

    // DTBSV - solve one of the systems of equations   A*x = b, or
    // A'*x = b,
    void F77_NAME(dtbsv)(const char& uplo, const char& trans,
			 const char& diag, const int& n, const int& k,
			 const double *a, const int& lda,
			 double *x, const int& incx); 

    // DTPSV - solve one of the systems of equations
    // A*x = b, or A'*x = b,
    void F77_NAME(dtpsv)(const char& uplo, const char& trans,
			 const char& diag, const int& n,
			 const double *ap,
			 double *x, const int& incx);

    // DTRSV - solve one of the systems of equations   A*x = b, or
    // A'*x = b,
    void F77_NAME(dtrsv)(const char& uplo, const char& trans,
			 const char& diag, const int& n,
			 const double *a, const int& lda,
			 double *x, const int& incx); 

    // DGER - perform the rank 1 operation   A := alpha*x*y' + A,
    void F77_NAME(dger)(const int& m, const int& n, const double& alpha,
			double *x, const int& incx,
			double *y, const int& incy,
			double *a, const int& lda);

    // DSYR - perform the symmetric rank 1 operation   A := alpha*x*x' + A,
    void F77_NAME(dsyr)(const char& uplo, const int& n,
			const double& alpha,
			const double *x, const int& incx,
			double *a, const int& lda);

    // DSPR - perform the symmetric rank 1 operation   A := alpha*x*x' + A,
    void F77_NAME(dspr)(const char& uplo, const int& n,
			const double& alpha,
			const double *x, const int& incx,
			double *ap);

    // DSYR2 - perform the symmetric rank 2 operation
    // A := alpha*x*y' + alpha*y*x' + A,
    void F77_NAME(dsyr2)(const char& uplo, const int& n,
			 const double& alpha,
			 const double *x, const int& incx,
			 const double *y, const int& incy,
			 double *a, const int& lda); 

    // DSPR2 - perform the symmetric rank 2 operation
    // A := alpha*x*y' + alpha*y*x' + A, 
    void F77_NAME(dspr2)(const char& uplo, const int& n,
			 const double& alpha,
			 const double *x, const int& incx,
			 const double *y, const int& incy,
			 double *ap);

}

#endif
// _BLAS2_H_
