// -*- c++ -*-
// $Id: blas3.h,v 1.2 2000/07/12 22:54:57 bates Exp $

// C++ prototypes for double precision Blas3 routines.

// Based on code labelled
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

// All arguments to Fortran routines must be passed as pointers or
// references.  Arguments that will be used as scalars are passed as
// references.  Arguments that will be used as arrays are passed as
// pointers.  The const keyword modifies those arguments that are
// guaranteed to be unchanged on exit.

// Modifications Copyright (C) 2000-2000 the R Development Core Team

extern "C" {
#include <R_ext/RS.h>		// to define F77_NAME

    // DGEMM - perform one of the matrix-matrix operations
    // C := alpha*op( A )*op( B ) + beta*C, 
    void F77_NAME(dgemm)(const char& transa, const char& transb,
			 const int& m, const int& n,
			 const int& k, const double& alpha,
			 const double *a, const int& lda,
			 const double *b, const int& ldb,
			 const double& beta, double *c, const int& ldc);

    // DTRSM - solve one of the matrix equations   op( A );*X =
    // alpha*B, or X*op( A ); = alpha*B, 
    void F77_NAME(dtrsm)(const char& side, const char& uplo,
			 const char& transa, const char& diag,
			 const int& m, const int& n,
			 const double& alpha,
			 const double *a, const int& lda,
			 double *b, const int& ldb);

    // DTRMM - perform one of the matrix-matrix operations   B :=
    // alpha*op( A )*B, or B := alpha*B*op( A ),
    void F77_NAME(dtrmm)(const char& side, const char& uplo,
			 const char& transa, const char& diag,
			 const int& m, const int& n,
			 const double& alpha,
			 const double *a, const int& lda,
			 double *b, const int& ldb); 
    
    // DSYMM - perform one of the matrix-matrix operations  
    //  C := alpha*A*B + beta*C,
    void F77_NAME(dsymm)(const char& side, const char& uplo,
			 const int& m, const int& n, const double& alpha,
			 const double *a, const int& lda,
			 const double *b, const int& ldb,
			 const double& beta, double *c, const int& ldc);

    // DSYRK - perform one of the symmetric rank k operations
    // C := alpha*A*A' + beta*C,
    void F77_NAME(dsyrk)(const char& uplo, const char& trans,
			 const int& n, const int& k,
			 const double& alpha,
			 const double *a, const int& lda,
			 const double& beta, double *c, const int& ldc);

    // DSYR2K - perform one of the symmetric rank 2k operations   C
    // := alpha*A*B' + alpha*B*A' + beta*C,
    void F77_NAME(dsyr2k)(const char& uplo, const char& trans,
			  const int& n, const int& k,
			  const double& alpha,
			  const double *a, const int& lda,
			  const double *b, const int& ldb,
			  const double& beta, double *c, const int& ldc);

}

