#include "cscBlocked.h"

SEXP cscBlocked_validate(SEXP x)
{
    SEXP pp = GET_SLOT(x, Matrix_pSym),
	ip = GET_SLOT(x, Matrix_iSym),
	xp = GET_SLOT(x, Matrix_xSym),
	dp = getAttrib(xp, R_DimSymbol);
    int *pv = INTEGER(pp),
	*dim = INTEGER(dp),
	ncol = length(pp) - 1;
    int nnz = pv[ncol];

    if (!(isReal(xp) && isArray(xp)))
	return ScalarString(mkChar("slot x should be a real array"));
    if (length(dp) != 3)
	return ScalarString(mkChar("slot x should be a 3-dimensional array"));
    if (length(ip) != nnz)
	return ScalarString(mkChar("length of slot i does not matck last element of slot p"));
    if (dim[2] != nnz)
	return ScalarString(mkChar("third dimension of slot x does not match number of nonzeros"));
    return ScalarLogical(1);
}

/** 
 * Perform one of the matrix operations 
 *  C := alpha*op(A)*B + beta*C
 * or
 *  C := alpha*B*op(A) + beta*C
 * where A is a compressed, sparse, blocked matrix and
 * B and C are dense matrices.
 * 
 * @param side 'L' or 'l' for left, otherwise right
 * @param transa 'T' or 't' for transpose.
 * @param m number of rows in C
 * @param n number of columns in C
 * @param k number of columns in B if side = 'L', otherwise
 *        number of rows in B.
 * @param alpha
 * @param nr number of rows per block of A
 * @param nc number of columns per block of A
 * @param ap vector of column pointers in A
 * @param ai vector of row indices in A
 * @param ax contents of non-zero blocks of A
 * @param b matrix to be multiplied
 * @param ldb leading dimension of b as declared in the calling
 *        routine
 * @param beta scalar multiplier of c
 * @param c product matrix to be modified
 * @param ldc leading dimension of c as declared in the calling
 *        routine
 */
void
cscBlocked_mm(char side, char transa, int m, int n, int k,
	      double alpha, int nr, int nc,
	      const int ap[], const int ai[],
	      const double ax[],
	      const double b[], int ldb,
	      double beta, double c[], int ldc)
{
    int j, kk, lside = (side == 'L' || side == 'l');
    int ncb, nrb, sz = nr * nc, tra = (transa == 'T' || transa == 't');

    if (nr < 1 || nc < 1 || m < 0 || n < 0 || k < 0)
	error("improper dims m=%d, n=%d, k=%d, nr=%d, nc=%d",
		  m, n, k, nr, nc);
    if (ldc < n) error("incompatible dims n=%d, ldc=%d", n, ldc);
    if (lside) {
	if (ldb < k)
	    error("incompatible L dims k=%d, ldb=%d, n=%d, nr=%d, nc=%d",
		  k, ldb, n, nr, nc);
	if (tra) {
	    if (m % nc || k % nr)
		error("incompatible LT dims m=%d, k = %d, nr=%d, nc=%d",
		      m, k, nr, nc);
	    nrb = k/nr; ncb = m/nc;
	} else {
	    if (m % nr || k % nc)
		error("incompatible LN dims m=%d, k = %d, nr=%d, nc=%d",
		      m, k, nr, nc);
	    nrb = m/nr; ncb = k/nc;
	}
	for (j = 0; j < ncb; j++) {
	    int j2 = ap[j + 1];
	    for (kk = ap[j]; kk < j2; kk++) {
		int ii = ai[kk];
		if (ii < 0 || ii >= nrb)
		    error("improper row index ii=%d, nrb=%d", ii, nrb);
		if (tra) {
		    F77_CALL(dgemm)("T", "N", &nc, &n, &nr,
				    &alpha, ax + j*sz, &nr,
				    b + ii * nr, &ldb,
				    &beta, c + j * nc, &ldc);
		} else {
		    F77_CALL(dgemm)("N", "N", &nr, &n, &nc,
				    &alpha, ax + j * sz, &nr,
				    b + j*nc, &ldb,
				    &beta, c + ii * nr, &ldc);
		}
	    }
	}
    } else {
	error("Call to cscBlocked_mm must have side == 'L'");
    }
}

/** 
 * Invert a triangular sparse blocked matrix.  This is not done in
 * place because the number of non-zero blocks in A-inverse can be
 * different than the number of non-zero blocks in A.
 * 
 * @param upper 'U' indicates upper triangular, 'L' lower
 * @param unit 'U' indicates unit diagonal, 'N' non-unit
 * @param n number of columns of blocks in A and A-inverse
 * @param nr number of rows per block in A and A-inverse
 * @param nc number of columns per block in A and A-inverse
 * @param ap vector of column pointers for A
 * @param ai vector of row indices for A
 * @param ax contents of the non-zero blocks of A
 * @param aip vector of column pointers for A-inverse
 * @param aii vector of row indices for A-inverse
 * @param aix contents of the non-zero blocks of A-inverse
 */
void
cscBlocked_tri(char upper, char unit, int n, int nr, int nc,
	       const int ap[], const int ai[], const double ax[],
	       int aip[], int aii[], double aix[])
{
    int iup = (upper == 'U' || upper == 'u') ? 1 : 0,
	iunit = (unit == 'U' || unit == 'u') ? 1 : 0;
    
    if (!iunit)
	error("Code for non-unit triangular matrices not yet written");
    if (ap[n] > 0)
	error("Code for non-trivial unit inverse not yet written");
    else {
	if (aip[n] == 0) return;
	error ("Structure of A and A-inverse does not agree");
    }
}
