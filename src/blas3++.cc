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

#include "blas3++.h"

void Blas_Mat_Mat_Mult(const LaGenMatDouble &A, 
                       const LaGenMatDouble &B, LaGenMatDouble &C, 
                       double alpha, double beta)
{
    F77_CALL(dgemm)('N', 'N', A.size(0), B.size(1), A.size(1), alpha,
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0), beta,
                    C.addr(), C.gdim(0)); 
}

void Blas_Mat_Trans_Mat_Mult(const LaGenMatDouble &A, 
                             const LaGenMatDouble &B, LaGenMatDouble &C, 
                             double alpha, double beta)
{
    F77_CALL(dgemm)('T', 'N', A.size(0), B.size(1), A.size(1), alpha,
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0), beta,
                    C.addr(), C.gdim(0));
}

void Blas_Mat_Mat_Trans_Mult(const LaGenMatDouble &A, 
                             const LaGenMatDouble &B, LaGenMatDouble &C, 
                             double alpha, double beta)
{
    F77_CALL(dgemm)('N', 'T', A.size(0), B.size(1), A.size(1), alpha,
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0), beta,
                    C.addr(), C.gdim(0)); 
}

#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Solve(LaUnitLowerTriangMatDouble &A, 
                        LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrsm)('L', 'L', 'N', 'U', B.size(0), B.size(1), alpha, 
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUnitUpperTriangMatDouble &A,
                       LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrmm)('L', 'U', 'N', 'U', B.size(0), B.size(1), alpha,
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}

void Blas_Mat_Mat_Solve(LaUnitUpperTriangMatDouble &A, 
                        LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrsm)('L', 'U', 'N', 'U', B.size(0), B.size(1), alpha, 
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}
#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaLowerTriangMatDouble &A,
                       LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrmm)('L', 'L', 'N', 'N', B.size(0), B.size(1), alpha,
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}

void Blas_Mat_Mat_Solve(LaLowerTriangMatDouble &A, 
                        LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrsm)('L', 'L', 'N', 'N', B.size(0), B.size(1), alpha, 
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}
#endif


#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUpperTriangMatDouble &A,
                       LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrmm)('L', 'U', 'N', 'N', B.size(0), B.size(1), alpha,
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}

void Blas_Mat_Mat_Solve(LaUpperTriangMatDouble &A, 
                        LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrsm)('L', 'U', 'N', 'N', B.size(0), B.size(1), alpha, 
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}
#endif

#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaUnitLowerTriangMatDouble &A,
                              LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrsm)('L', 'L', 'T', 'U', B.size(0), B.size(1), alpha,
                    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaUnitUpperTriangMatDouble &A,
			      LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrsm)('L', 'U', 'T', 'U', B.size(0), B.size(1), alpha,
		    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}
#endif

#ifdef LA_LOWER_TRIANG_MAT_DOUBLE_H
void Blas_Mat_Mat_Mult(LaUnitLowerTriangMatDouble &A,
		       LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrmm)('L', 'L', 'N', 'U', B.size(0), B.size(1), alpha,
		    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}

void Blas_Mat_Trans_Mat_Solve(LaLowerTriangMatDouble &A,
			      LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrsm)('L', 'L', 'T', 'N', B.size(0), B.size(1), alpha,
		    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}
#endif


#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H 
void Blas_Mat_Trans_Mat_Solve(LaUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha)
{
    F77_CALL(dtrsm)('L', 'U', 'T', 'N', B.size(0), B.size(1), alpha,
		    A.addr(), A.gdim(0), B.addr(), B.gdim(0));
}
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaSymmMatDouble &A, LaGenMatDouble &B, 
		       LaGenMatDouble &C, double alpha,
		       double beta)
{
    F77_CALL(dsymm)('L', 'L', C.size(0), C.size(1), alpha, A.addr(),
		    A.gdim(0), B.addr(), B.gdim(0), beta, C.addr(),
		    C.gdim(0));
}

void Blas_R1_Update(LaSymmMatDouble &C, LaGenMatDouble &A,
		    double alpha, double beta)
{
    F77_CALL(dsyrk)('L', 'N', C.size(0), A.size(1), alpha, A.addr(),
		    A.gdim(0), beta, C.addr(), C.gdim(0));
}


void Blas_R2_Update(LaSymmMatDouble &C, LaGenMatDouble &A,
		    LaGenMatDouble &B, double alpha,
		    double beta)
{
    F77_CALL(dsyr2k)('L', 'N', C.size(0), A.size(1), alpha, A.addr(),
		     A.gdim(0), B.addr(), B.gdim(0), beta, C.addr(),
		     C.gdim(0));
}
#endif
