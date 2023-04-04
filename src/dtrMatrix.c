/* double (precision) TRiangular Matrices */

#include "dtrMatrix.h"

double get_norm_dtr(SEXP obj, const char *typstr)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(obj, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *pdim = INTEGER(dim);
    double *px = REAL(x), norm, *work = NULL;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
    if (typstr[0] == 'I')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlantr)(typstr, ul, di, pdim, pdim + 1, px, pdim,
			    work FCONE FCONE FCONE);

    UNPROTECT(4);
    return norm;
}

SEXP dtrMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dtr(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dtrMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(obj, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));

    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_rcond_type(CHAR(type));
    
    int *pdim = INTEGER(dim), info;
    double *px = REAL(x), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
    F77_CALL(dtrcon)(typstr, ul, di, pdim, px, pdim, &rcond,
		     (double *) R_alloc((size_t) 3 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE FCONE FCONE);

    UNPROTECT(5);
    return ScalarReal(rcond);
}

SEXP dtrMatrix_solve(SEXP a)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	dim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(a, Matrix_diagSym)),
	x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(a, Matrix_xSym), &pid);
    REPROTECT(x = duplicate(x), pid);
    
    SET_SLOT(val, Matrix_DimSym, dim);
    set_reversed_DimNames(val, dimnames);
    SET_SLOT(val, Matrix_uploSym, uplo);
    SET_SLOT(val, Matrix_diagSym, diag);
    SET_SLOT(val, Matrix_xSym, x);
    
    int *pdim = INTEGER(dim), info;
    double *px = REAL(x);
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
    F77_CALL(dtrtri)(ul, di, pdim, px, pdim, &info FCONE FCONE);

    UNPROTECT(6);
    return val;
}

SEXP dtrMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0)),
	adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	bdim = PROTECT(GET_SLOT(val, Matrix_DimSym));
    int *padim = INTEGER(adim), *pbdim = INTEGER(bdim);
    
    if (padim[0] != pbdim[0] || padim[0] < 1 || pbdim[1] < 1)
	error(_("dimensions of system to be solved are inconsistent"));
    
    SEXP uplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(a, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(a, Matrix_xSym)),
	y = PROTECT(GET_SLOT(val, Matrix_xSym));
    
    double *px = REAL(x), *py = REAL(y), one = 1.0;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));

    F77_CALL(dtrsm)("L", ul, "N", di, pbdim, pbdim + 1, &one,
		    px, pbdim, py, pbdim FCONE FCONE FCONE FCONE);

    UNPROTECT(7);
    return val;
}

SEXP dtrMatrix_chol2inv(SEXP a)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dpoMatrix")),
	dim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
	x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(a, Matrix_xSym), &pid);
    REPROTECT(x = duplicate(x), pid);

    SET_SLOT(val, Matrix_DimSym, dim);
    SET_SLOT(val, Matrix_DimNamesSym, dimnames);
    SET_SLOT(val, Matrix_uploSym, uplo);
    SET_SLOT(val, Matrix_xSym, x);

    int *pdim = INTEGER(dim), info;
    double *px = REAL(x);
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    F77_CALL(dpotri)(ul, pdim, px, pdim, &info FCONE);

    UNPROTECT(5);
    return val;
}

/** Matrix products of dense triangular Matrices
 *
 * @param a triangular matrix of class "dtrMatrix"
 * @param b  ( ditto )
 * @param right logical, if true, compute b %*% a,  else  a %*% b
 * @param trans logical, if true, "transpose a", i.e., use t(a), otherwise a
 *
 * @return the matrix product, one of   a %*% b, t(a) %*% b,  b %*% a, or  b %*% t(a)
 *      depending on (right, trans) =    (F, F)    (F, T)      (T, F)        (T, T)
 */
SEXP dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans)
{
    /* called from "%*%" : (x,y, FALSE,FALSE),
             crossprod() : (x,y, FALSE, TRUE) , and
	     tcrossprod(): (y,x, TRUE , TRUE)
     * 	     -
     * TWO cases : (1) result is triangular  <=> uplo's "match" (i.e., non-equal iff trans)
     * ===         (2) result is "general"
     */
    SEXP val,/* = in case (2):  dense_as_general(b, 'd', 2, 0); */
	d_a = GET_SLOT(a, Matrix_DimSym),
	uplo_a = GET_SLOT(a, Matrix_uploSym),  diag_a = GET_SLOT(a, Matrix_diagSym),
	uplo_b = GET_SLOT(b, Matrix_uploSym),  diag_b = GET_SLOT(b, Matrix_diagSym);
    int rt = asLogical(right);
    int tr = asLogical(trans);
    int *adims = INTEGER(d_a), n = adims[0];
    double *valx = (double *) NULL /*Wall*/;
    const char
	*uplo_a_ch = CHAR(STRING_ELT(uplo_a, 0)), /* = uplo_P(a) */
	*diag_a_ch = CHAR(STRING_ELT(diag_a, 0)), /* = diag_P(a) */
	*uplo_b_ch = CHAR(STRING_ELT(uplo_b, 0)), /* = uplo_P(b) */
	*diag_b_ch = CHAR(STRING_ELT(diag_b, 0)); /* = diag_P(b) */
    Rboolean same_uplo = (*uplo_a_ch == *uplo_b_ch),
	matching_uplo = tr ? (!same_uplo) : same_uplo,
	uDiag_b = /* -Wall: */ FALSE;

    if (INTEGER(GET_SLOT(b, Matrix_DimSym))[0] != n)
	/* validity checking already "assures" square matrices ... */
	error(_("dimension mismatch in matrix multiplication of \"dtrMatrix\": %d != %d"),
	      n, INTEGER(GET_SLOT(b, Matrix_DimSym))[0]);
    if(matching_uplo) {
	/* ==> result is triangular -- "dtrMatrix" ! */
	R_xlen_t sz = n * (R_xlen_t) n, np1 = n+1;
	val = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix"));
	SET_SLOT(val, Matrix_uploSym, duplicate(uplo_b));
	SET_SLOT(val, Matrix_DimSym,  duplicate(d_a));
	set_DimNames(val, GET_SLOT(b, Matrix_DimNamesSym));
	valx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, sz));
	Memcpy(valx, REAL(GET_SLOT(b, Matrix_xSym)), sz);
	if((uDiag_b = (*diag_b_ch == 'U'))) {
	    /* unit-diagonal b - may contain garbage in diagonal */
	    for (int i = 0; i < n; i++)
		valx[i * np1] = 1.;
	}
    } else { /* different "uplo" ==> result is "dgeMatrix" ! */
	val = PROTECT(dense_as_general(b, 'd', 2, 0));
	SEXP
	    dn_a = GET_SLOT( a , Matrix_DimNamesSym),
	    dn   = GET_SLOT(val, Matrix_DimNamesSym);
	/* matrix product   a %*% b, t(a) %*% b,  b %*% a, or  b %*% t(a)
	 * (right, trans) =  (F, F)    (F, T)      (T, F)        (T, T)
	 *   set:from_a   =   0:0       0:1         1:1           1:0
	 */
	SET_VECTOR_ELT(dn, rt ? 1 : 0, VECTOR_ELT(dn_a, (rt+tr) % 2));
    }
    if (n >= 1) {
	double alpha = 1.;
	/* Level 3 BLAS - DTRMM(): Compute one of the matrix multiplication operations
	 * B := alpha*op( A )*B ["L"], or B := alpha*B*op( A ) ["R"],
	 *	where trans_A determines  op(A):=  A   "N"one  or
	 *				  op(A):= t(A) "T"ransposed */
	F77_CALL(dtrmm)(rt ? "R" : "L", uplo_a_ch,
			/*trans_A = */ tr ? "T" : "N", diag_a_ch, &n, &n, &alpha,
			REAL(GET_SLOT(a,   Matrix_xSym)), adims,
			REAL(GET_SLOT(val, Matrix_xSym)),
			&n FCONE FCONE FCONE FCONE);
    }
    if(matching_uplo) {
	/* set "other triangle" to 0 */
	if (tr)
	    ddense_unpacked_make_triangular(valx, n, n, *uplo_b_ch, *diag_b_ch);
	else
	    ddense_unpacked_make_triangular(valx, n, n, *uplo_a_ch, *diag_a_ch);
	if(*diag_a_ch == 'U' && uDiag_b) /* result remains uni-diagonal */
	    SET_SLOT(val, Matrix_diagSym, duplicate(diag_a));
    }
    UNPROTECT(1);
    return val;
}

// to be used for all three: '%*%', crossprod() and tcrossprod()
/** Matrix products  dense triangular Matrices o  <matrix>
 *
 * @param a triangular matrix of class "dtrMatrix"
 * @param b a <matrix> or <any-denseMatrix>
 * @param right logical, if true, compute b %*% a,  else  a %*% b
 * @param trans logical, if true, "transpose a", i.e., use t(a), otherwise a
 *
 * @return the matrix product, one of   a %*% b, t(a) %*% b,  b %*% a, or  b %*% t(a)
 *      depending on (right, trans) =    (F, F)    (F, T)      (T, F)        (T, T)
 */
SEXP dtrMatrix_matrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans)
{
    /* called from "%*%", crossprod() and tcrossprod() in  ../R/products.R
     *
     * Because 'a' must be square, the size of the answer 'val',
     * is the same as the size of 'b' */
    SEXP val = PROTECT(dense_as_general(b, 'd', 2, 0));
    int rt = asLogical(right); /* if(rt), compute b %*% op(a),  else  op(a) %*% b */
    int tr = asLogical(trans);/* if true, use t(a) */
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int m = bdims[0], n = bdims[1];
    double one = 1.;

    if (adims[0] != adims[1])
	error(_("dtrMatrix must be square"));
    if ((rt && adims[0] != n) || (!rt && adims[1] != m))
	error(_("Matrices are not conformable for multiplication"));
    if (m >= 1 && n >= 1) {
	// Level 3 BLAS - DTRMM() --> see call further below
	F77_CALL(dtrmm)(rt ? "R" : "L", uplo_P(a),
			/*trans_A = */ tr ? "T" : "N",
			diag_P(a), &m, &n, &one,
			REAL(GET_SLOT(a, Matrix_xSym)), adims,
			REAL(GET_SLOT(val, Matrix_xSym)),
			&m FCONE FCONE FCONE FCONE);
    }

    SEXP
	dn_a = GET_SLOT( a,  Matrix_DimNamesSym),
	dn   = GET_SLOT(val, Matrix_DimNamesSym);
    /* matrix product   a %*% b, t(a) %*% b,  b %*% a, or  b %*% t(a)
     * (right, trans) =  (F, F)    (F, T)      (T, F)        (T, T)
     *   set:from_a   =   0:0       0:1         1:1           1:0
     */
    SET_VECTOR_ELT(dn, rt ? 1 : 0, VECTOR_ELT(dn_a, (rt+tr) % 2));

    UNPROTECT(1);
    return val;
}

SEXP dtrMatrix_addDiag(SEXP x, SEXP d) {
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    SEXP ret = PROTECT(duplicate(x)),
	r_x = GET_SLOT(ret, Matrix_xSym);
    double *dv = REAL(d), *rv = REAL(r_x);

    if ('U' == diag_P(x)[0])
	error(_("cannot add diag() as long as 'diag = \"U\"'"));
    for (int i = 0; i < n; i++) rv[i * (n + 1)] += dv[i];

    UNPROTECT(1);
    return ret;
}

/* MJ: no longer needed ... prefer more general unpackedMatrix_diag_[gs]et() */
#if 0

#define GET_trMatrix_Diag(_C_TYPE_, _SEXPTYPE_, _SEXP_, _ONE_)		\
    int i, n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];			\
    SEXP x_x = GET_SLOT(x, Matrix_xSym);				\
									\
    SEXP ret = PROTECT(allocVector(_SEXPTYPE_, n));			\
    _C_TYPE_ *rv = _SEXP_(ret),						\
	     *xv = _SEXP_(x_x);						\
									\
    if ('U' == diag_P(x)[0]) {						\
	for (i = 0; i < n; i++) rv[i] = _ONE_;				\
    } else {								\
	for (i = 0; i < n; i++) rv[i] = xv[i * (n + 1)];		\
    }									\
    UNPROTECT(1);							\
    return ret


SEXP dtrMatrix_getDiag(SEXP x) {
    GET_trMatrix_Diag(double, REALSXP, REAL, 1.);
}

SEXP ltrMatrix_getDiag(SEXP x) {
    GET_trMatrix_Diag(  int, LGLSXP, LOGICAL, 1);
}

#define SET_trMatrix_Diag(_C_TYPE_, _SEXP_)				\
    if ('U' == diag_P(x)[0])						\
	error(_("cannot set diag() as long as 'diag = \"U\"'"));	\
			    /* careful to recycle RHS value: */		\
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];			\
    int l_d = LENGTH(d); Rboolean d_full = (l_d == n);			\
    if (!d_full && l_d != 1)						\
	error(_("replacement diagonal has wrong length"));		\
    SEXP ret = PROTECT(duplicate(x)),					\
	r_x = GET_SLOT(ret, Matrix_xSym);				\
    _C_TYPE_ *dv = _SEXP_(d),						\
	     *rv = _SEXP_(r_x);						\
									\
    if(d_full) for (int i = 0; i < n; i++)				\
	rv[i * (n + 1)] = dv[i];					\
    else for (int i = 0; i < n; i++)					\
	rv[i * (n + 1)] = *dv;						\
									\
    UNPROTECT(1);							\
    return ret

SEXP dtrMatrix_setDiag(SEXP x, SEXP d) {
    SET_trMatrix_Diag(double, REAL);
}

SEXP ltrMatrix_setDiag(SEXP x, SEXP d) {
    SET_trMatrix_Diag(  int, LOGICAL);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0

SEXP dtrMatrix_as_dtpMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dtpMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_pack(
	REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, (n*(n+1))/2)),
	REAL(GET_SLOT(from, Matrix_xSym)),
	n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	*CHAR(STRING_ELT(diag, 0)) == 'N' ? NUN : UNT);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0

SEXP dtrMatrix_as_matrix(SEXP from, SEXP keep_dimnames)
{
    int *Dim = INTEGER(GET_SLOT(from, Matrix_DimSym));
    int m = Dim[0], n = Dim[1];
    SEXP val = PROTECT(allocMatrix(REALSXP, m, n));
    ddense_unpacked_make_triangular(Memcpy(REAL(val),
					   REAL(GET_SLOT(from, Matrix_xSym)),
					   m * n),
				    from);
    if(asLogical(keep_dimnames))
	setAttrib(val, R_DimNamesSymbol, GET_SLOT(from, Matrix_DimNamesSym));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */
