#include "lgCMatrix.h"

SEXP lgCMatrix_validate(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym);
    int j,
	ncol = length(pslot) - 1,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow,
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    nrow = dims[0];
    if (length(pslot) <= 0)
	return mkString(_("slot p must have length > 0"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return mkString(_("last element of slot p must match length of slot i"));
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
    }
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    if (csc_unsorted_columns(ncol, xp, xi))
	csc_sort_columns(ncol, xp, xi, (double *) NULL);

    return ScalarLogical(1);
}

/* very parallel to csc_to_matrix() in ./dgCMatrix.c */
SEXP lcsc_to_matrix(SEXP x)
{
    SEXP ans, pslot = GET_SLOT(x, Matrix_pSym);
    int j, ncol = length(pslot) - 1,
	nrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	*xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    int *ax;

    ax = LOGICAL(ans = PROTECT(allocMatrix(LGLSXP, nrow, ncol)));
    for (j = 0; j < (nrow * ncol); j++) ax[j] = 0;
    for (j = 0; j < ncol; j++) {
	int ind;
	for (ind = xp[j]; ind < xp[j+1]; ind++)
	    ax[j * nrow + xi[ind]] = 1;
    }
    UNPROTECT(1);
    return ans;
}

#ifdef _NEED_logical_to_csc_FIRST_
/* very parallel to matrix_to_csc() in ./dgCMatrix.c */
SEXP matrix_to_lcsc(SEXP A)
{
    if (!(isMatrix(A) && isLogical(A)))
	error(_("A must be a logical matrix"));
    return logical_to_csc(LOGICAL(A),
			  INTEGER(getAttrib(A, R_DimSymbol)));
}
#endif

/**
 * C := op(A) %*% op(B) + beta ^ C for logical sparse column-oriented matrices
 *
 * @param tra nonzero if A is to be transposed
 * @param trb nonzero if B is to be transposed
 * @param m number of rows in C
 * @param n number of columns in C
 * @param k number of columns in A if tra == 0, otherwise number of
 *          rows in A
 * @param ai vector of row indices of TRUE elements in A
 * @param ap column pointers for A
 * @param bi vector of row indices of TRUE elements in B
 * @param bp column pointers for B
 * @param beta if non-zero existing TRUE elements in C are retained
 * @param ciP SEXP whose INTEGER part is the column indices of TRUE
 * elements in C (not used if beta == 0).
 * @param cp column pointers for C
 *
 * @return SEXP whose INTEGER part is the column indices of TRUE
 * elements in the product.  Note that the contents of cp may be modified.
 */
SEXP Matrix_lgClgCmm(int tra, int trb, int m, int n, int k,
		     const int ai[], const int ap[],
		     const int bi[], const int bp[],
		     int beta, SEXP CIP, int cp[])
{
    int cnnz = cp[n], extra = 0;
    int *ci, i, j, prot = 0;	/* prot is the number of PROTECTs to UNPROTECT */

    if (beta) {
	ci = INTEGER(CIP);
    } else {			/* blank the C matrix */
	for (j = 0; j <= n; j++) cp[j] = 0;
	cnnz = 0;
	ci = (int *) NULL;
    }

    if (tra) {			/* replace ai and ap by els for transpose */
	int nz = ap[m];
	int *Ai = Calloc(nz, int),
	    *aj = expand_cmprPt(m, ap, Calloc(nz, int)),
	    *Ap = Calloc(k + 1, int);

	triplet_to_col(m, k, nz, aj, ai, (double *) NULL,
		       Ap, Ai, (double *) NULL);
	Free(aj);
	ai = Ai; ap = Ap;
    }

    if (trb) {			/* replace bi and bp by els for transpose */
	int nz = bp[k];
	int *Bi = Calloc(nz, int),
	    *bj = expand_cmprPt(k, bp, Calloc(nz, int)),
	    *Bp = Calloc(n + 1, int);

	triplet_to_col(k, n, nz, bj, bi, (double *) NULL,
		       Bp, Bi, (double *) NULL);
	Free(bj);
	bi = Bi; bp = Bp;
    }

    for (j = 0; j < n; j++) { /* col index for B and C */
	int ii, ii2 = bp[j + 1];
	for (ii = bp[j]; ii < ii2; ii++) { /* index into bi */
	    int jj = bi[ii]; /* row index of B; col index of A */
	    int i, i2 = ap[jj + 1]; /* index into ai */
	    for (i = ap[jj]; i < i2; i++)
		if (check_csc_index(cp, ci, ai[i], j, -1) < 0) extra++;
	}
    }

    if (extra) {
	int ntot = cnnz + extra;
	int *Cp = Calloc(n + 1, int),
	    *Ti = Calloc(ntot, int),
	    *rwInd = Calloc(m, int), /* indicator of TRUE in column j */
	    pos = 0;

	Cp[0] = 0;
	for (j = 0; j < n; j++) {
	    int ii, ii2 = bp[j + 1];

	    AZERO(rwInd, m);	/* initialize column j of C */
	    for (i = cp[j]; i < cp[j+1]; i++) rwInd[ci[i]] = 1;

	    Cp[j + 1] = Cp[j];
	    for (ii = bp[j]; ii < ii2; ii++) { /* index into bi */
		int jj = bi[ii]; /* row index of B; col index of A */
		int i, i2 = ap[jj + 1]; /* index into ai */
		for (i = ap[jj]; i < i2; i++) rwInd[ai[i]] = 1;
	    }
	    for (i = 0; i < m; i++)
		if (rwInd[i]) {Cp[j + 1]++; Ti[pos++] = i;}
	}
	PROTECT(CIP = allocVector(INTSXP, Cp[n])); prot++;
	Memcpy(INTEGER(CIP), Ti, Cp[n]);
	Memcpy(cp, Cp, n + 1);
	Free(Cp); Free(Ti); Free(rwInd);
    }

    if (tra) {Free(ai); Free(ap);}
    if (trb) {Free(bi); Free(bp);}
    UNPROTECT(prot);
    return CIP;
}

SEXP lgCMatrix_lgCMatrix_mm(SEXP a, SEXP b)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lgCMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	*cdims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    int k = adims[1], m = adims[0], n = bdims[1];
    int *cp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1));

    if (bdims[0] != k)
	error(_("Matrices are not conformable for multiplication"));
    cdims[0] = m; cdims[1] = n;
    SET_SLOT(ans, Matrix_iSym,
	     Matrix_lgClgCmm(0, 0, m, n, k,
			     INTEGER(GET_SLOT(a, Matrix_iSym)),
			     INTEGER(GET_SLOT(a, Matrix_pSym)),
			     INTEGER(GET_SLOT(b, Matrix_iSym)),
			     INTEGER(GET_SLOT(b, Matrix_pSym)),
			     0, (SEXP) NULL, cp));
    UNPROTECT(1);
    return ans;
}

SEXP lgCMatrix_trans(SEXP x)
{
    SEXP xi = GET_SLOT(x, Matrix_iSym);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lgCMatrix")));
    int *adims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	*xdims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nz = length(xi);
    int *xj = Calloc(nz, int);
    SEXP adn = ALLOC_SLOT(ans, Matrix_DimNamesSym, VECSXP, 2),
	xdn = GET_SLOT(x, Matrix_DimNamesSym);

    adims[1] = xdims[0]; adims[0] = xdims[1];
    SET_VECTOR_ELT(adn, 0, VECTOR_ELT(xdn, 1));
    SET_VECTOR_ELT(adn, 1, VECTOR_ELT(xdn, 0));
    triplet_to_col(adims[0], adims[1], nz,
		   expand_cmprPt(xdims[1], INTEGER(GET_SLOT(x, Matrix_pSym)), xj),
		   INTEGER(xi), (double *) NULL,
		   INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP,  adims[1] + 1)),
		   INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP,  nz)),
		   (double *) NULL);
    Free(xj);
    UNPROTECT(1);
    return ans;
}

/**
 * Replace C by AA' + beta*C or A'A + beta*C
 *
 * @param up Indicator of upper/lower triangle in the symmetric sparse matrix
 * @param tra Transpose, in the sense of dsyrk.  That is, tra TRUE indicates A'A
 * @param n size of the product matrix
 * @param k number of columns in A if tra is FALSE, otherwise the number of rows
 * @param ai row indices for A
 * @param ap column pointers for A
 * @param beta TRUE if existing elements in C are to be preserved
 * @param CIP SEXP whose INTEGER part is the row indices of C (not used if beta is FALSE)
 * @param cp column pointers for C
 *
 * @return SEXP whose INTEGER part is the updated row indices of C
 */
SEXP Matrix_lgCsyrk(int up, int tra, int n, int k, const int ai[], const int ap[],
		    int beta, SEXP CIP, int cp[])
{
    int extra = 0, i, ii, j, prot = 0;
    int *ci, cnnz = cp[n];

    if (beta) {
	ci = INTEGER(CIP);
    } else {			/* blank the C matrix */
	for (j = 0; j <= n; j++) cp[j] = 0;
	cnnz = 0;
	ci = (int *) NULL;
    }

    if (tra) {			/* replace ai and ap by els for transpose */
	int nz = ap[n];
	int *Ai = Calloc(nz, int),
	    *aj = expand_cmprPt(n, ap, Calloc(nz, int)),
	    *Ap = Calloc(k + 1, int);

	triplet_to_col(n, k, nz, aj, ai, (double *) NULL,
		       Ap, Ai, (double *) NULL);
	Free(aj);
	ai = Ai; ap = Ap;
    }

    for (j = 0; j < k; j++) {
	int i2 = ap[j + 1];
	for (i = ap[j]; i < i2; i++) {
	    int r1 = ai[i];
	    if (r1 < 0 || r1 >= n)
		error(_("row %d not in row range [0,%d]"), r1, n - 1);
	    for (ii = i; ii < i2; ii++) {
		int r2 = ai[ii];
		if (r2 < 0 || r2 >= n)
		    error(_("row %d not in row range [0,%d]"), r2, n - 1);
		if (check_csc_index(cp, ci, up?r1:r2, up?r2:r1, -1) < 0)
		    extra++;
	    }
	}
    }

    if (extra) {
	int ntot = cnnz + extra;
	int *Ti = Memcpy(Calloc(ntot, int), ci, cnnz),
	    *Tj = expand_cmprPt(n, cp, Calloc(ntot, int)),
	    *Ci = Calloc(ntot, int),
	    pos = cnnz;

	for (j = 0; j < k; j++) {
	    int i2 = ap[j + 1];
	    for (i = ap[j]; i < i2; i++) {
		int r1 = ai[i];
		for (ii = i; ii < i2; ii++) {
		    int r2 = ai[ii];
		    int row = up ? r1 : r2, col = up ? r2 : r1;
		    if (r2 < r1) error("[j,i,ii,r1,r2] = [%d,%d,%d,%d,%d]",
				       j,i,ii,r1,r2);
		    if (check_csc_index(cp, ci, row, col, -1) < 0) {
			Ti[pos] = row;
			Tj[pos] = col;
			pos++;
		    }
		}
	    }
	}

	triplet_to_col(n, n, pos, Ti, Tj, (double *) NULL,
		       cp, Ci, (double *) NULL);
	PROTECT(CIP = allocVector(INTSXP, cp[n])); prot++;
	Memcpy(INTEGER(CIP), Ci, cp[n]);
	Free(Ti); Free(Tj); Free(Ci);
    }

    if (tra) {Free(ai); Free(ap);}
    UNPROTECT(prot);
    return CIP;
}

/**
 * Create the cross-product or transpose cross-product of a logical
 * sparse matrix in column-oriented compressed storage mode.
 *
 * @param x Pointer to a lgCMatrix
 * @param trans logical indicator of transpose, in the sense of dsyrk.
 * That is, trans == TRUE is used for crossprod.
 * @param C
 *
 * @return An lsCMatrix of the form if(trans) X'X else XX'
 */
SEXP lgCMatrix_crossprod(SEXP x, SEXP trans, SEXP C)
{
    int tra = asLogical(trans);
    int *adims, *xdims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int k = xdims[tra ? 0 : 1], n = xdims[tra ? 1 : 0];

    if (C == R_NilValue) {
	SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lsCMatrix")));

	adims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
	adims[0] = adims[1] = n;
	SET_SLOT(ans, Matrix_uploSym, mkString("U"));
	SET_SLOT(ans, Matrix_iSym,
		 Matrix_lgCsyrk(1, tra, n, k,
				INTEGER(GET_SLOT(x, Matrix_iSym)),
				INTEGER(GET_SLOT(x, Matrix_pSym)),
				0, R_NilValue,
				INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1))));
	UNPROTECT(1);
	return ans;
    }
    adims = INTEGER(GET_SLOT(C, Matrix_DimSym));
    if (adims[0] != n || adims[1] != n)
	error(_("Dimensions of x and y are not compatible for crossprod"));
    SET_SLOT(C, Matrix_iSym,
	     Matrix_lgCsyrk(uplo_P(C)[0] == 'U',
			    tra, n, k,
			    INTEGER(GET_SLOT(x, Matrix_iSym)),
			    INTEGER(GET_SLOT(x, Matrix_pSym)),
			    1, GET_SLOT(C, Matrix_iSym),
			    INTEGER(GET_SLOT(C, Matrix_pSym))));
    return C;
}

/**
 * Special-purpose function that returns a permutation of the columns
 * of a lgTMatrix for which nrow(x) > ncol(x).  The ordering puts
 * columns with fewer entries on the left.  Once a column has been
 * moved to the left the rows in where that column is TRUE are removed
 * from the counts.
 *
 * @param x Pointer to an lgTMatrix object
 *
 * @return 0-based permutation vector for the columns of x
 */
SEXP lgCMatrix_picky_column(SEXP x)
{
    int *xdims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int *xi = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*xp = INTEGER(GET_SLOT(x, Matrix_pSym)),
	m = xdims[0], n = xdims[1];
    SEXP ans = PROTECT(allocVector(INTSXP, n));
    int *actr = Calloc(m, int),
	*actc = Calloc(n, int),
	cj, i, j, mincount, minloc = -1, pos;

    for (i = 0; i < m; i++) actr[i] = 1;
    mincount = m + 1;
    for (j = 0; j < n; j++) {
	cj = xp[j + 1] - xp[j];
	actc[j] = 1;
	if (cj < mincount) {
	    mincount = cj;
	    minloc = j;
	}
    }

    pos = 0;
    while (pos < n) {
	INTEGER(ans)[pos++] = minloc;
	actc[minloc] = 0;
	for (i = xp[minloc]; i < xp[minloc + 1]; i++) actr[xi[i]] = 0;
	mincount = m + 1;
	for (j = 0; j < n; j++) {
	    if (actc[j]) {
		cj = 0;
		for (i = xp[j]; i < xp[j + 1]; i++) {
		    if (actr[xi[i]]) cj++;
		    if (cj < mincount) {
			mincount = cj;
			minloc = j;
		    }
		}
	    }
	}
    }

    Free(actr); Free(actc);
    UNPROTECT(1);
    return ans;
}
