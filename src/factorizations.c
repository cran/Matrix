#include <Rmath.h> /* logspace_add, logspace_sub */
#include "factorizations.h"

SEXP denseLU_expand(SEXP obj)
{
    /* A = P L U   <=>   P' A = L U   where ... 
       
       A -> [m,n]
       P -> [m,m], permutation
       L -> if m <= n then [m,n] else [m,m], lower trapezoidal, unit diagonal
       U -> if m >= n then [m,n] else [n,n], upper trapezoidal
       
       square L,U given as dtrMatrix with appropriate 'uplo', 'diag' slots,
       non-square L,U given as dgeMatrix
    */
    
    const char *nms[] = {"L", "U", "P", ""};
    PROTECT_INDEX pidA, pidB;
    SEXP res = PROTECT(Rf_mkNamed(VECSXP, nms)),
	P = PROTECT(NEW_OBJECT_OF_CLASS("pMatrix")),
	dim, x;
    PROTECT_WITH_INDEX(dim = GET_SLOT(obj, Matrix_DimSym), &pidA);
    PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pidB);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n, j;
    
    if (m == n) {
	SEXP L = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	    U = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	    uplo = PROTECT(mkString("L")),
	    diag = PROTECT(mkString("U"));
	SET_SLOT(L, Matrix_DimSym, dim);
	SET_SLOT(U, Matrix_DimSym, dim);
	SET_SLOT(P, Matrix_DimSym, dim);
	SET_SLOT(L, Matrix_uploSym, uplo);
	SET_SLOT(L, Matrix_diagSym, diag);
	SET_SLOT(L, Matrix_xSym, x);
	SET_SLOT(U, Matrix_xSym, x);
	SET_VECTOR_ELT(res, 0, L);
	SET_VECTOR_ELT(res, 1, U);
	UNPROTECT(4); /* diag, uplo, U, L */
    } else {
	SEXP G = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),
	    T = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	    y = PROTECT(allocVector(REALSXP, (R_xlen_t) r * r));
	REPROTECT(x = duplicate(x), pidB);
	double *px = REAL(x), *py = REAL(y);
	int whichT = (m < n) ? 0 : 1;
	
	SET_SLOT(G, Matrix_DimSym, dim);
	REPROTECT(dim = allocVector(INTSXP, 2), pidA);
	pdim = INTEGER(dim);
	pdim[0] = pdim[1] = r;
	SET_SLOT(T, Matrix_DimSym, dim);
	REPROTECT(dim = allocVector(INTSXP, 2), pidA);
	pdim = INTEGER(dim);
	pdim[0] = pdim[1] = m;
	SET_SLOT(P, Matrix_DimSym, dim);
	
	if (whichT == 0) {
            /* G is upper trapezoidal, T is unit lower triangular */
	    SEXP uplo = PROTECT(mkString("L")),
		diag = PROTECT(mkString("U"));
	    SET_SLOT(T, Matrix_uploSym, uplo);
	    SET_SLOT(T, Matrix_diagSym, diag);
	    UNPROTECT(2); /* diag, uplo */

	    Memcpy(py, px, (size_t) m * m);
	    ddense_unpacked_make_triangular(px, m, n, 'U', 'N');
	} else {
            /* G is unit lower trapezoidal, T is upper triangular */
	    double *tmp = px;
	    for (j = 0; j < n; ++j, px += m, py += r)
		Memcpy(py, px, j+1);
	    ddense_unpacked_make_triangular(tmp, m, n, 'L', 'U');
	}
	SET_SLOT(G, Matrix_xSym, x);
	SET_SLOT(T, Matrix_xSym, y);
	
	SET_VECTOR_ELT(res, !whichT, G);
	SET_VECTOR_ELT(res,  whichT, T);
	UNPROTECT(3); /* y, T, G */
    }

    SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
	perm = PROTECT(allocVector(INTSXP, m));
    int *ppivot = INTEGER(pivot), *pperm = INTEGER(perm), *pinvperm, pos, tmp;
    Calloc_or_Alloca_TO(pinvperm, m, int);

    for (j = 0; j < m; ++j) /* initialize inverse permutation */
	pinvperm[j] = j;
    for (j = 0; j < r; ++j) { /* generate inverse permutation */
	pos = ppivot[j] - 1;
	if (pos != j) {
	    tmp = pinvperm[j];
	    pinvperm[j] = pinvperm[pos];
	    pinvperm[pos] = tmp;
	}
    }
    for (j = 0; j < m; ++j) /* invert inverse permutation (0->1-based) */
	pperm[pinvperm[j]] = j + 1;
    Free_FROM(pinvperm, m);

    SET_SLOT(P, Matrix_permSym, perm);
    SET_VECTOR_ELT(res, 2, P);
    UNPROTECT(6); /* perm, pivot, x, dim, P, res */
    return res;
}

SEXP BunchKaufman_expand(SEXP obj)
{
    SEXP P_ = PROTECT(NEW_OBJECT_OF_CLASS("pMatrix")),
	T_ = PROTECT(NEW_OBJECT_OF_CLASS("dtCMatrix")),
	D_ = PROTECT(NEW_OBJECT_OF_CLASS("dsCMatrix")),
	dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int i, j, s, n = INTEGER(dim)[0];
    R_xlen_t n1a = (R_xlen_t) n + 1;
    if (n > 0) {
	SET_SLOT(P_, Matrix_DimSym, dim);
	SET_SLOT(T_, Matrix_DimSym, dim);
	SET_SLOT(D_, Matrix_DimSym, dim);
    }
    UNPROTECT(1); /* dim */

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
    if (!upper) {
	SET_SLOT(T_, Matrix_uploSym, uplo);
	SET_SLOT(D_, Matrix_uploSym, uplo);
    }
    UNPROTECT(1); /* uplo */
    
    SEXP diag = PROTECT(mkString("U"));
    SET_SLOT(T_, Matrix_diagSym, diag);
    UNPROTECT(1); /* diag */
    
    SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
	D_p = PROTECT(allocVector(INTSXP, n1a));
    int *ppivot = INTEGER(pivot), *D_pp = INTEGER(D_p),
	b = n, dp = (upper) ? 1 : 2;
    D_pp[0] = 0;
    j = 0;
    while (j < n) {
	if (ppivot[j] > 0) {
	    D_pp[j+1] = D_pp[j] + 1;
	    j += 1;
	} else {
	    D_pp[j+1] = D_pp[j] + dp;
	    D_pp[j+2] = D_pp[j] + 3;
	    j += 2;
	    --b;
	}
    }
    SET_SLOT(D_, Matrix_pSym, D_p);
    UNPROTECT(1); /* D_p */

    SEXP P, P_perm, T, T_p, T_i, T_x,
	D_i = PROTECT(allocVector(INTSXP, D_pp[n])),
	D_x = PROTECT(allocVector(REALSXP, D_pp[n])),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *P_pperm, *T_pp, *T_pi, *D_pi = INTEGER(D_i);
    double *T_px, *D_px = REAL(D_x), *px = REAL(x);

    int unpacked = (double) n * n <= R_XLEN_T_MAX &&
	(R_xlen_t) n * n == XLENGTH(x);

    R_xlen_t len = (R_xlen_t) 2 * b + 1, k = (upper) ? len - 1 : 0;
    SEXP res = PROTECT(allocVector(VECSXP, len));

    j = 0;
    while (b--) {
	s = (ppivot[j] > 0) ? 1 : 2;
	dp = (upper) ? j : n - j - s;
	
	PROTECT(P = duplicate(P_));
	PROTECT(P_perm = allocVector(INTSXP, n));
	PROTECT(T = duplicate(T_));
	PROTECT(T_p = allocVector(INTSXP, n1a));
	PROTECT(T_i = allocVector(INTSXP, (R_xlen_t) s * dp));
	PROTECT(T_x = allocVector(REALSXP, (R_xlen_t) s * dp));
	
	P_pperm = INTEGER(P_perm);
	T_pp = INTEGER(T_p);
	T_pi = INTEGER(T_i);
	T_px = REAL(T_x);
	T_pp[0] = 0;
	
	for (i = 0; i < j; ++i) {
	    T_pp[i+1] = 0;
	    P_pperm[i] = i + 1;
	}
	for (i = j; i < j+s; ++i) {
	    T_pp[i+1] = T_pp[i] + dp;
	    P_pperm[i] = i + 1;
	}
	for (i = j+s; i < n; ++i) {
	    T_pp[i+1] = T_pp[i];
	    P_pperm[i] = i + 1;
	}
	
	if (s == 1) {
	    P_pperm[j] = ppivot[j];
	    P_pperm[ppivot[j]-1] = j + 1;
	} else if (upper) {
	    P_pperm[j] = -ppivot[j];
	    P_pperm[-ppivot[j]-1] = j + 1;
	} else {
	    P_pperm[j+1] = -ppivot[j];
	    P_pperm[-ppivot[j]-1] = j + 2;
	}

	if (upper) {
	    for (i = 0; i < j; ++i) {
		*(T_pi++) = i;
		*(T_px++) = *(px++);
	    }
	    *(D_pi++) = j;
	    *(D_px++) = *(px++);
	    ++j;
	    if (unpacked)
		px += n - j;
	    if (s == 2) {
		for (i = 0; i < j-1; ++i) {
		    *(T_pi++) = i;
		    *(T_px++) = *(px++);
		}
		*(D_pi++) = j - 1;
		*(D_pi++) = j;
		*(D_px++) = *(px++);
		*(D_px++) = *(px++);
		++j;
		if (unpacked)
		    px += n - j;
	    }
	} else {
	    if (s == 2) {
		*(D_pi++) = j;
		*(D_pi++) = j + 1;
		*(D_px++) = *(px++);
		*(D_px++) = *(px++);
		for (i = j+2; i < n; ++i) {
		    *(T_pi++) = i;
		    *(T_px++) = *(px++);
		}
		++j;
		if (unpacked)
		    px += j;
	    }
	    *(D_pi++) = j;
	    *(D_px++) = *(px++);
	    for (i = j+1; i < n; ++i) {
		*(T_pi++) = i;
		*(T_px++) = *(px++);
	    }
	    ++j;
	    if (unpacked)
		px += j;
	}

	SET_SLOT(P, Matrix_permSym, P_perm);
	SET_SLOT(T, Matrix_pSym, T_p);
	SET_SLOT(T, Matrix_iSym, T_i);
	SET_SLOT(T, Matrix_xSym, T_x);

	if (upper) {
	    SET_VECTOR_ELT(res, k-1, P);
	    SET_VECTOR_ELT(res, k  , T);
	    k -= 2;
	} else {
	    SET_VECTOR_ELT(res, k  , P);
	    SET_VECTOR_ELT(res, k+1, T);
	    k += 2;
	}
	UNPROTECT(6); /* T_x, T_i, T_p, T, P_perm, P */
    }
    
    SET_SLOT(D_, Matrix_iSym, D_i);
    SET_SLOT(D_, Matrix_xSym, D_x);
    SET_VECTOR_ELT(res, k, D_);

    UNPROTECT(8); /* res, x, D_x, D_i, pivot, D_, T_, P_ */ 
    return res;
}

SEXP denseLU_determinant(SEXP obj, SEXP logarithm)
{
    /* MJ: unfortunately, we do not retain the 'info' given by 'dgetrf'
           ... if we knew info>0, then we could return 0/-Inf "fast" 
	   as base R does, or check for NaN and return NaN in that case 
	   (the "right" thing to do)
    */
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("determinant of non-square matrix is undefined"));
    UNPROTECT(1); /* dim */
    int givelog = asLogical(logarithm) != 0, sign = 1;
    double modulus = (givelog) ? 0.0 : 1.0; /* result for n == 0 */
    if (n > 0) {
	SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
	    x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int j, *ppivot = INTEGER(pivot);
	R_xlen_t n1a = (R_xlen_t) n + 1;
	double *px = REAL(x);
	
	if (givelog) {
	    for (j = 0; j < n; ++j, px += n1a, ++ppivot) {
		if (*px < 0) {
		    modulus += log(-(*px));
		    if (*ppivot == j + 1)
			sign = -sign;
		} else {
		    /* incl. 0, NaN cases */
		    modulus += log(*px);
		    if (*ppivot != j + 1)
			sign = -sign;
		}
	    }
	} else {
	    for (j = 0; j < n; ++j, px += n1a, ++ppivot) {
		modulus *= *px;
		if (*ppivot != j + 1)
		    sign = -sign;
	    }
	    if (modulus < 0.0) {
		modulus = -modulus;
		sign = -sign;
	    }
	}
	UNPROTECT(2); /* x, pivot */
    }
    return as_det_obj(modulus, givelog, sign);
}

SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm)
{
    /* MJ: unfortunately, we do not retain the 'info' given by 'ds[yp]trf'
           ... if we knew info>0, then we could return 0/-Inf "fast" 
	   as base R does, or check for NaN and return NaN in that case 
	   (the "right" thing to do)
    */
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    int givelog = asLogical(logarithm) != 0, sign = 1;
    double modulus = (givelog) ? 0.0 : 1.0; /* result for n == 0 */
    if (n > 0) {
	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	UNPROTECT(1); /* uplo */

	SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
	    x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int j = 0, *ppivot = INTEGER(pivot);
	R_xlen_t n1a = (R_xlen_t) n + 1;
	double *px = REAL(x), a, b, c;

	int unpacked = (double) n * n <= R_XLEN_T_MAX &&
	    (R_xlen_t) n * n == XLENGTH(x);
	
	if (givelog) {
	    double logab, logcc;
	    while (j < n) {
		if (ppivot[j] > 0) {
		    if (*px < 0) {
			modulus += log(-(*px));
			sign = -sign;
		    } else {
			/* incl. 0, NaN cases */
			modulus += log(*px);
		    }
		    px += (unpacked) ? n1a : ((upper) ? j + 2 : n - j);
		    j += 1;
		} else {
		    a = *px;
		    if (upper) {
			px += (unpacked) ? n1a : j + 2;
			b = *px;
			c = *(px - 1);
			px += (unpacked) ? n1a : j + 3;
		    } else {
			c = *(px + 1);
			px += (unpacked) ? n1a : n - j;
			b = *px;
			px += (unpacked) ? n1a : n - j - 1;
		    }
		    logab = log((a < 0.0) ? -a : a) + log((b < 0.0) ? -b : b);
		    logcc = 2.0 * log((c < 0.0) ? -c : c);
		    if ((a < 0.0) != (b < 0.0)) {
			/* det = ab - cc = -(abs(ab) + cc) < 0 */
			modulus += logspace_add(logab, logcc);
			sign = -sign;
		    } else if (logab < logcc) {
			/* det = ab - cc = -(cc - ab) < 0 */
			modulus += logspace_sub(logcc, logab);
			sign = -sign;
		    } else {
			/* det = ab - cc > 0 */
			modulus += logspace_sub(logab, logcc);
		    }
		    j += 2;
		}
	    }
	} else {
	    while (j < n) {
		if (ppivot[j] > 0) {
		    modulus *= *px;
		    px += (unpacked) ? n1a : ((upper) ? j + 2 : n - j);
		    j += 1;
		} else {
		    a = *px;
		    if (upper) {
			px += (unpacked) ? n1a : j + 2;
			b = *px;
			c = *(px - 1);
			px += (unpacked) ? n1a : j + 3;
		    } else {
			c = *(px + 1);
			px += (unpacked) ? n1a : n - j;
			b = *px;
			px += (unpacked) ? n1a : n - j - 1;
		    }
		    modulus *= a * b - c * c;
		    j += 2;
		}
	    }
	    if (modulus < 0.0) {
		modulus = -modulus;
		sign = -sign;
	    }
	}
	UNPROTECT(2); /* x, pivot */
    }
    return as_det_obj(modulus, givelog, sign);
}

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0

SEXP LU_validate(SEXP obj)
{
    /* NB: 'Dim' already checked by MatrixFactorization_validate() */
    
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (!isReal(x))
	return mkString(_("'x' slot is not of type \"double\""));
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if (XLENGTH(x) != pdim[0] * (double) pdim[1])
	return mkString(_("length of 'x' slot is not prod(Dim)"));
    return DimNames_validate(obj, pdim);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer denseLU_expand() */
#if 0

SEXP LU_expand(SEXP x)
{
    const char *nms[] = {"L", "U", "P", ""};
    // x[,] is  m x n    (using LAPACK dgetrf notation)
    SEXP L, U, P, val = PROTECT(Rf_mkNamed(VECSXP, nms)),
	lux = GET_SLOT(x, Matrix_xSym),
	dd = GET_SLOT(x, Matrix_DimSym);
    int *iperm, *perm, *pivot = INTEGER(GET_SLOT(x, Matrix_permSym)),
	*dim = INTEGER(dd), m = dim[0], n = dim[1], nn = m, i;
    size_t m_ = (size_t) m; // to prevent integer (multiplication) overflow
    Rboolean is_sq = (n == m), L_is_tri = TRUE, U_is_tri = TRUE;

    // nn :=  min(n,m) ==  length(pivot[])
    if(!is_sq) {
	if(n < m) { // "long"
	    nn = n;
	    L_is_tri = FALSE;
	} else { // m < n : "wide"
	    U_is_tri = FALSE;
	}
    }

    SET_VECTOR_ELT(val, 0, NEW_OBJECT_OF_CLASS(L_is_tri ? "dtrMatrix":"dgeMatrix"));
    SET_VECTOR_ELT(val, 1, NEW_OBJECT_OF_CLASS(U_is_tri ? "dtrMatrix":"dgeMatrix"));
    SET_VECTOR_ELT(val, 2, NEW_OBJECT_OF_CLASS("pMatrix"));
    L = VECTOR_ELT(val, 0);
    U = VECTOR_ELT(val, 1);
    P = VECTOR_ELT(val, 2);
    if(is_sq || !L_is_tri) {
	SET_SLOT(L, Matrix_xSym, duplicate(lux));
	SET_SLOT(L, Matrix_DimSym, duplicate(dd));
    } else { // !is_sq && L_is_tri -- m < n -- "wide" -- L is  m x m
	size_t m2 = m_ * m;
	double *Lx = REAL(ALLOC_SLOT(L, Matrix_xSym, REALSXP, m2));
	int *dL = INTEGER(ALLOC_SLOT(L, Matrix_DimSym, INTSXP, 2));
	dL[0] = dL[1] = m;
	// fill lower-diagonal (non-{0,1}) part -- remainder by ddense_unpacked_make_*() below:
	Memcpy(Lx, REAL(lux), m2);
    }
    if(is_sq || !U_is_tri) {
	SET_SLOT(U, Matrix_xSym, duplicate(lux));
	SET_SLOT(U, Matrix_DimSym, duplicate(dd));
    } else { // !is_sq && U_is_tri -- m > n -- "long" -- U is  n x n
	double *Ux = REAL(ALLOC_SLOT(U, Matrix_xSym, REALSXP, ((size_t) n) * n)),
	       *xx = REAL(lux);
	int *dU = INTEGER(ALLOC_SLOT(U, Matrix_DimSym, INTSXP, 2));
	dU[0] = dU[1] = n;
	/* fill upper-diagonal (non-0) part -- remainder by ddense_unpacked_make_*() below:
	 * this is more complicated than in the L case, as the x / lux part we need
	 * is  *not*  continguous:  Memcpy(Ux, REAL(lux), n * n); -- is  WRONG */
	for (size_t j = 0; j < n; j++) {
	    Memcpy(Ux+j*n, xx+j*m, j+1);
	    // for (i = 0; i <= j; i++)
	    //   Ux[i + j*n] = xx[i + j*m];
	}
    }
    if(L_is_tri) {
	SET_SLOT(L, Matrix_uploSym, mkString("L"));
	SET_SLOT(L, Matrix_diagSym, mkString("U"));
    }
    // fill the upper right part with 0  *and* the diagonal with 1
    ddense_unpacked_make_triangular(REAL(GET_SLOT(L, Matrix_xSym)),
				    m, (is_sq || !L_is_tri) ? n : m, 'L', 'U');

    if(U_is_tri) {
	SET_SLOT(U, Matrix_uploSym, mkString("U"));
	SET_SLOT(U, Matrix_diagSym, mkString("N"));
	
    }
    // fill the lower left part with 0
    ddense_unpacked_make_triangular(REAL(GET_SLOT(U, Matrix_xSym)),
				    (is_sq || !U_is_tri) ? m : n, n, 'U', 'N');
    
    SET_SLOT(P, Matrix_DimSym, duplicate(dd));
    if(!is_sq) // m != n -- P is  m x m
	INTEGER(GET_SLOT(P, Matrix_DimSym))[1] = m;
    perm = INTEGER(ALLOC_SLOT(P, Matrix_permSym, INTSXP, m));
    Calloc_or_Alloca_TO(iperm, m, int);

    for (i = 0; i < m; i++) iperm[i] = i + 1; /* initialize permutation*/
    for (i = 0; i < nn; i++) {	/* generate inverse permutation */
	int newp = pivot[i] - 1;
	if (newp != i) { // swap
	    int tmp = iperm[i]; iperm[i] = iperm[newp]; iperm[newp] = tmp;
	}
    }
    // invert the inverse
    for (i = 0; i < m; i++) perm[iperm[i] - 1] = i + 1;

    Free_FROM(iperm, m);
    UNPROTECT(1);
    return val;
}

#endif /* MJ */
