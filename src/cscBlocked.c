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
    int /* iup = (upper == 'U' || upper == 'u') ? 1 : 0, */
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

/** 
 * Search for the element in a compressed sparse matrix at a given row and column
 * 
 * @param row row index
 * @param col column index
 * @param cp column pointers
 * @param ci row indices
 * 
 * @return index of element in ci, if it exists, else -1
 */
static R_INLINE
int Ind(int row, int col, const int cp[], const int ci[])
{
    int i, i2 = cp[col + 1];
    for (i = cp[col]; i < i2; i++)
	if (ci[i] == row) return i;
    return -1;
}

/** 
 * Perform one of the matrix operations 
 *  C := alpha*A*A' + beta*C,
 * or
 *  C := alpha*A'*A + beta*C,
 * where A is a compressed, sparse, blocked matrix and
 * C is a compressed, sparse, symmetric blocked matrix.
 * 
 * @param uplo 'U' or 'u' for upper triangular storage, else lower.
 * @param trans 'T' or 't' for transpose.
 * @param alpha scalar multiplier of outer product
 * @param A compressed sparse blocked matrix
 * @param beta scalar multiplier of c
 * @param C compressed sparse blocked symmetric matrix to be updated
 */
SEXP cscBlocked_syrk(SEXP uplo, SEXP trans, SEXP alpha, SEXP A, SEXP beta, SEXP C)
{
    SEXP axp = GET_SLOT(A, Matrix_xSym),
	app = GET_SLOT(A, Matrix_pSym),
	cxp = GET_SLOT(C, Matrix_xSym),
	cpp = GET_SLOT(C, Matrix_pSym);
    int *adims = INTEGER(getAttrib(axp, R_DimSymbol)),
	*ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*ap = INTEGER(app),
	*cdims = INTEGER(getAttrib(cxp, R_DimSymbol)),
	*ci = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*cp = INTEGER(cpp),
	j, k,
	nca = length(app) - 1,
	ncc = length(cpp) - 1,
	npairs,
	inplace;
    char *ul, *tr;
    double *ax = REAL(axp),
	*cx = REAL(cxp),
	bta = asReal(beta);

    if (length(uplo) != 1 || length(trans) != 1 || length(alpha) != 1 ||
	length(beta) != 1 || !isReal(alpha) || !isReal(beta) ||
	!isString(uplo) || !isString(trans))
	error("uplo and trans must be character scalars, alpha and beta real scalars");
    if (cdims[0] != cdims[1]) error("blocks in C must be square");
    ul = CHAR(STRING_ELT(uplo, 0));
    tr = CHAR(STRING_ELT(trans, 0));
    if (toupper(tr[0]) == 'T')
	error("Code for trans == 'T' not yet written");
/* FIXME: Write the other version */
    if (adims[0] != cdims[0])
	error("Dimension inconsistency in blocks: dim(A)=[%d,,], dim(C)=[%d,,]",
	      adims[0], cdims[0]);
				/* check the row indices */
    for (k = 0; k < adims[2]; k++) {
	int aik = ai[k];
	if (aik < 0 || aik >= ncc)
	    error("Row index %d = %d is out of range [0, %d]",
		  k, ai[k], ncc - 1);
    }
				/* check if C can be updated in place */
    inplace = 1;
    npairs = 0;
    for (j = 0; j < nca; j++) {
	int k, kk, k2 = ap[j+1];
	int nnz = k2 - ap[j];
	npairs += (nnz * (nnz - 1)) / 2;
	for (k = ap[j]; k < k2; k++) {
/* FIXME: This check assumes uplo == 'L' */
	    for (kk = k; kk < k2; kk++) {
		if (Ind(ai[k], ai[kk], cp, ci) < 0) {
		    inplace = 0;
		    break;
		}
	    }
	    if (!inplace) break;
	}
    }
				/* multiply C by beta */
    for (j = 0; j < cdims[0]*cdims[1]*cdims[2]; j++) cx[j] *= bta;
    if (inplace) {
	int scalar = (adims[0] == 1 && adims[1] == 1);
	
	for (j = 0; j < nca; j++) {
	    int k, kk, k2 = ap[j+1];
	    for (k = ap[j]; k < k2; k++) {
		int ii = ai[k];
		for (kk = k; kk < k2; kk++) {
		    int jj = ai[kk], K = Ind(ii, jj, cp, ci);
		    if (scalar) cx[K] += ax[k] * ax[kk];
		    else {
		    }
		}
	    }
	}
    }
    return C;
}

