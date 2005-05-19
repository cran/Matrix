				/* LDL' factorization of a logical,
				 * sparse column-oriented matrix */
#include "lCholCMatrix.h"

SEXP lCholCMatrix_validate(SEXP x)
{
    int i, n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    SEXP perm = GET_SLOT(x, Matrix_permSym),
	Parent = GET_SLOT(x, Matrix_ParentSym);

    if (length(perm) != n)
	return mkString(_("slot perm must have length n"));
    if (length(Parent) != n)
	return mkString(_("slot Parent must have length n"));
    if (!R_ldl_valid_perm(n, INTEGER(perm))) 
	return mkString(_("slot perm is not a valid 0-based permutation"));
    for (i = 0; i < n; i++) {
	int pari = INTEGER(Parent)[i];
	if (pari < -1 || pari > (n - 1))
	    return mkString(_("an element of the Parent array is not in range [-1,n-1]"));
    }

    return ScalarLogical(1);
}

/** 
 * Create the structure of the inverse of L from the LDL' factorization.
 * 
 * @param x Pointer to a lCholCMatrix object
 * 
 * @return An ltCMatrix object representing L^{-1}
 */
SEXP lCholCMatrix_solve(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("ltCMatrix")));
    SEXP Parent = GET_SLOT(x, Matrix_ParentSym);
    int i, n = length(Parent), pari, pos;
    int *ai, *ap, *nzc = Calloc(n, int), ntot;

    ntot = 0;
    for (i = n - 1; i >= 0; i--) { /* count non-zeros in each column */
	pari = INTEGER(Parent)[i];
	ntot += (nzc[i] = (pari >= 0) ? 1 + nzc[pari] : 0);
    }

    slot_dup(ans, x, Matrix_DimSym);
    slot_dup(ans, x, Matrix_DimNamesSym);
    slot_dup(ans, x, Matrix_uploSym);
    slot_dup(ans, x, Matrix_diagSym);
    ai = INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, ntot));
    ap = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1));
    ap[0] = pos = 0;
    for (i = 0; i < n; i++) {
	ap[i + 1] = ap[i] + nzc[i];
	for (pari = INTEGER(Parent)[i]; pari >= 0; pari = INTEGER(Parent)[pari])
	    ai[pos++] = pari;
    }

    Free(nzc);
    UNPROTECT(1);
    return ans;
}

/** 
 * Solve  one  of the matrix equations  op(A)*C = B, or * C*op(A) = B
 * where A is an lCholCMatrix object and B and C are lgCMatrix 
 * objects.
 *
 * An early check is done to see if A is the identity, in which case
 * BIP is returned and bp is unmodified.
 * 
 * @param side LFT or RGT
 * @param transa TRN or NTR
 * @param m number of rows in B and C
 * @param n number of columns in B and C
 * @param Parent Parent array of A
 * @param BIP pointer to the row indices for B
 * @param bp array of column pointers for B (may be overwritten)
 *
 * @return Pointer to the updated row indices for C.  The column
 * pointers bp are overwritten with cp.
 */
SEXP
lCholClgCsm(enum CBLAS_SIDE side, enum CBLAS_TRANSPOSE transa, int m,
	    int n, const int Parent[], SEXP BIP, int bp[])
{
    int *bi = INTEGER(BIP), bnz, extra, ident, j, pari, pos;
    int nca = (transa == TRN) ? n : m;

    ident = 1;
    for (j = 0; j < nca; j++)
	if (Parent[j] >= 0) {
	    ident = 0;
	    break;
	}
    if (ident) return BIP;

    extra = 0;
    if (side == LFT) {
	if (transa == TRN) {
	    error(_("code not yet written"));
	    return R_NilValue;
	} else {
	    int *Tci, *Tj, *Ti, ntot;
	    for (j = 0; j < n; j++) {
		int ii, ii2 = bp[j + 1];
		for (ii = bp[j]; ii < ii2; ii++)
		    for (pari = Parent[bi[ii]]; pari >= 0;
			 pari = Parent[pari]) extra++;
	    }

	    bnz = bp[n];
	    ntot = bnz + extra;
	    Ti = Memcpy(Calloc(ntot, int), bi, bnz);
	    Tj = expand_cmprPt(n, bp, Calloc(ntot, int));
	    Tci = Calloc(ntot, int);

	    pos = bnz;
	    for (j = 0; j < n; j++) {
		int ii, ii2 = bp[j + 1];
		for (ii = bp[j]; ii < ii2; ii++)
		    for (pari = Parent[bi[ii]]; pari >= 0;
			 pari = Parent[pari]) {
			    Ti[pos] = pari;
			    Tj[pos] = j;
			    pos++;
			}
	    }
	    triplet_to_col(m, n, ntot, Ti, Tj, (double *) NULL,
			   bp, Tci, (double *) NULL);

	    bnz = bp[n];
	    BIP = PROTECT(allocVector(INTSXP, bnz));
	    Memcpy(INTEGER(BIP), Tci, bnz);
	    
	    Free(Tci); Free(Ti); Free(Tj);
	    UNPROTECT(1);
	    return BIP;
	}
    } else {
	if (transa == TRN) {
	    int *Tci, *Tj, *Ti, ntot;
	    for (j = 0; j < n; j++) {
		int ii, ii2 = bp[j + 1];
		for (pari = Parent[j]; pari >= 0; pari = Parent[pari])
		    for (ii = bp[j]; ii < ii2; ii++) extra++;
	    }

	    bnz = bp[n];
	    ntot = bnz + extra;
	    Ti = Memcpy(Calloc(ntot, int), bi, bnz);
	    Tj = expand_cmprPt(n, bp, Calloc(ntot, int));
	    Tci = Calloc(ntot, int);

	    pos = bnz;
	    for (j = 0; j < n; j++) {
		int ii, ii2 = bp[j + 1];
		for (pari = Parent[j]; pari >= 0; pari = Parent[pari]) {
		    for (ii = bp[j]; ii < ii2; ii++) {
			Ti[pos] = bi[ii];
			Tj[pos] = pari;
			pos++;
		    }
		}
	    }
	    triplet_to_col(m, n, ntot, Ti, Tj, (double *) NULL,
			   bp, Tci, (double *) NULL);
	    bnz = bp[n];
	    BIP = PROTECT(allocVector(INTSXP, bnz));
	    Memcpy(INTEGER(BIP), Tci, bnz);
	    
	    Free(Tci); Free(Ti); Free(Tj);
	    UNPROTECT(1);
	    return BIP;
	} else {
	    error(_("code not yet written"));
	    return R_NilValue;
	}
    }
}

SEXP lCholCMatrix_lgCMatrix_solve(SEXP a, SEXP b)
{
    SEXP ans = PROTECT(duplicate(b));
    int n = INTEGER(GET_SLOT(a, Matrix_DimSym))[0],
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym));

    if (n != bdims[0])
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      n, n, bdims[0], bdims[1]);
    SET_SLOT(ans, Matrix_iSym,
	     lCholClgCsm(LFT, NTR, bdims[0], bdims[1],
			 INTEGER(GET_SLOT(a, Matrix_ParentSym)),
			 GET_SLOT(ans, Matrix_iSym),
			 INTEGER(GET_SLOT(ans, Matrix_pSym))));
    UNPROTECT(1);
    return ans;
}
			 
	     
