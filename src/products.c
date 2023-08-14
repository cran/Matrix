#include "products.h"
#include "chm_common.h"
#include "coerce.h"

/* Given a denseMatrix, diagonalMatrix, matrix, or vector,
   return a newly allocated dgeMatrix with newly allocated 'x'
   and 'Dimnames' slots and an empty 'factors' slot; the 'Dim'
   slot and elements of the 'Dimnames' slots need not be newly
   allocated.

   FIXME: Refactor the products and avoid such complexity altogether.
*/
static SEXP asdge(SEXP from, int transpose_if_vector)
{
	static const char *valid[] = {
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, VALID_DIAGONAL, "" };
	int ivalid = R_check_class_etc(from, valid);

	SEXP to;
	if (ivalid < 0)
		PROTECT(to = matrix_as_dense(from, "dge", '\0', '\0',
		                             transpose_if_vector, 1));
	else {
		const char *cl = valid[ivalid];
		if (cl[0] == 'd' && cl[1] == 'g' && cl[2] == 'e') {
			PROTECT(to = NEW_OBJECT_OF_CLASS("dgeMatrix"));
			SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
				dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym)),
				x = PROTECT(GET_SLOT(from, Matrix_xSym));
			PROTECT(x = duplicate(x));
			SET_SLOT(to, Matrix_DimSym, dim);
			SET_SLOT(to, Matrix_DimNamesSym, dimnames);
			SET_SLOT(to, Matrix_xSym, x);
			UNPROTECT(4);
		} else if (cl[0] == 'd') {
			if (cl[1] == 'd' && cl[2] == 'i')
				PROTECT(to = diagonal_as_dense(from, cl, 'g', 0, '\0'));
			else {
				PROTECT(to = dense_as_general(from, cl, 1));
				SEXP factors = PROTECT(allocVector(VECSXP, 0));
				SET_SLOT(to, Matrix_factorSym, factors);
				UNPROTECT(1);
			}
		} else {
			char cl_[] = "d..Matrix";
			cl_[1] = cl[1]; cl_[2] = cl[2];
			if (cl[1] == 'd' && cl[2] == 'i') {
				to = diagonal_as_kind(from, cl, 'd');
				PROTECT(to);
				to = diagonal_as_dense(to, cl_, 'g', 0, '\0');
			} else {
				to = dense_as_kind(from, cl, 'd');
				PROTECT(to);
				to = dense_as_general(to, cl_, 0);
			}
			UNPROTECT(1);
			PROTECT(to);
		}
	}

	SEXP dn0 = PROTECT(GET_SLOT(to, Matrix_DimNamesSym)),
		dn1 = PROTECT(allocVector(VECSXP, 2));
	for (int i = 0; i < 2; ++i)
		SET_VECTOR_ELT(dn1, i, VECTOR_ELT(dn0, i));
	SET_SLOT(to, Matrix_DimNamesSym, dn1);

	UNPROTECT(3);
	return to;
}

SEXP dgeMatrix_crossprod(SEXP x, SEXP trans)
{
#define DGE_CROSS_1							\
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x) */		\
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dpoMatrix")),		\
	vDnms = PROTECT(ALLOC_SLOT(val, Matrix_DimNamesSym, VECSXP, 2)),\
	nms  = VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym), tr ? 0 : 1);	\
    int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),			\
	*vDims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));	\
    int k = tr ? Dims[1] : Dims[0],					\
	n = tr ? Dims[0] : Dims[1];					\
    R_xlen_t n_ = n, n2 = n_ * n_;					\
    double *vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n2)),	\
	one = 1.0, zero = 0.0;						\
    									\
    Memzero(vx, n2);							\
    SET_SLOT(val, Matrix_uploSym, mkString("U"));			\
    ALLOC_SLOT(val, Matrix_factorSym, VECSXP, 0);			\
    vDims[0] = vDims[1] = n;						\
    SET_VECTOR_ELT(vDnms, 0, duplicate(nms));				\
    SET_VECTOR_ELT(vDnms, 1, duplicate(nms))

#define DGE_CROSS_DO(_X_X_)					\
    if(n)							\
	F77_CALL(dsyrk)("U", tr ? "N" : "T", &n, &k, &one,	\
			_X_X_, Dims, &zero, vx, &n FCONE FCONE);	\
    UNPROTECT(2);						\
    return val

    DGE_CROSS_1;
    DGE_CROSS_DO(REAL(GET_SLOT(x, Matrix_xSym)));
}

static double *gematrix_real_x(SEXP x, int nn)
{
    if(class_P(x)[0] == 'd') // <<- FIXME: use R_check_class_etc(x, valid)  !!!
	return REAL(GET_SLOT(x, Matrix_xSym));
#ifdef _potentically_more_efficient_but_not_working
    // else : 'l' or 'n' (for now !!)
    int *xi = INTEGER(GET_SLOT(x, Matrix_xSym));
    double *x_x;
    Matrix_Calloc(x_x, nn, double);
    for(int i=0; i < nn; i++)
	x_x[i] = (double) xi[i];

    // FIXME: this is not possible either; the *caller* would have to R_Free(.)
    Matrix_Free(x_x, nn);
#else
    // ideally should be PROTECT()ed ==> make sure R does not run gc() now!
    double *x_x = REAL(coerceVector(GET_SLOT(x, Matrix_xSym), REALSXP));
#endif
    return x_x;
}

//!  As  dgeMatrix_crossprod(), but x can be [dln]geMatrix
static SEXP _geMatrix_crossprod(SEXP x, SEXP trans)
{
    DGE_CROSS_1;
    double *x_x = gematrix_real_x(x, k * n_);
    DGE_CROSS_DO(x_x);
}

SEXP geMatrix_crossprod(SEXP x, SEXP trans)
{
    SEXP y = PROTECT(asdge(x, 0)),
	val = _geMatrix_crossprod(y, trans);
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
#define DGE_DGE_CROSS_1							\
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x,y) */		\
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),		\
	dn = PROTECT(allocVector(VECSXP, 2));				\
    int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),			\
	*yDims = INTEGER(GET_SLOT(y, Matrix_DimSym)),			\
	*vDims;								\
    int m  = xDims[!tr],  n = yDims[!tr];/* -> result dim */		\
    int xd = xDims[ tr], yd = yDims[ tr];/* the conformable dims */	\
    double one = 1.0, zero = 0.0;					\
									\
    if (xd != yd)							\
	error(_("Dimensions of x and y are not compatible for %s"),	\
	      tr ? "tcrossprod" : "crossprod");				\
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));		\
    /* establish dimnames */						\
    SET_VECTOR_ELT(dn, 0,						\
		   duplicate(VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym), \
					tr ? 0 : 1)));			\
    SET_VECTOR_ELT(dn, 1,						\
		   duplicate(VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym), \
					tr ? 0 : 1)));			\
    SET_SLOT(val, Matrix_DimNamesSym, dn);				\
    vDims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));		\
    vDims[0] = m; vDims[1] = n;						\
    double *v = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * (R_xlen_t) n))

#define DGE_DGE_CROSS_DO(_X_X_, _Y_Y_)					\
    if (xd > 0 && n > 0 && m > 0)					\
	F77_CALL(dgemm)(tr ? "N" : "T", tr ? "T" : "N", &m, &n, &xd, &one, \
			_X_X_, xDims,					\
			_Y_Y_, yDims, &zero, v, &m FCONE FCONE);	\
    else								\
	Memzero(v,  m * (R_xlen_t) n);					\
    UNPROTECT(2);							\
    return val

    DGE_DGE_CROSS_1;
    DGE_DGE_CROSS_DO(REAL(GET_SLOT(x, Matrix_xSym)),
		     REAL(GET_SLOT(y, Matrix_xSym)));
}

//!  As  dgeMatrix_dgeMatrix_crossprod(), but x and y can be [dln]geMatrix
static SEXP _geMatrix__geMatrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
    DGE_DGE_CROSS_1;

    double *x_x = gematrix_real_x(x, m * (R_xlen_t) xd);
    double *y_x = gematrix_real_x(y, n * (R_xlen_t) yd);

    DGE_DGE_CROSS_DO(x_x, y_x);
}
#undef DGE_DGE_CROSS_1
#undef DGE_DGE_CROSS_DO

SEXP geMatrix_geMatrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
    SEXP gx = PROTECT(asdge(x, 0)),
	gy = PROTECT(asdge(y, 0)),
	val = _geMatrix__geMatrix_crossprod(gx, gy, trans);
    UNPROTECT(2);
    return val;
}

SEXP dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
#define DGE_MAT_CROSS_1							\
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x,y) */		\
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),		\
	dn = PROTECT(allocVector(VECSXP, 2)),				\
	yDnms = R_NilValue, yD;						\
    int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),			\
	*yDims, *vDims, nprot = 2;					\
    int m  = xDims[!tr],						\
	xd = xDims[ tr];						\
    double one = 1.0, zero = 0.0;					\
    Rboolean y_has_dimNames;						\
									\
    if (!isReal(y)) {							\
	if(isInteger(y) || isLogical(y)) {				\
	    y = PROTECT(coerceVector(y, REALSXP));			\
	    nprot++;							\
	}								\
	else								\
	    error(_("Argument y must be numeric, integer or logical"));	\
    }									\
    if(isMatrix(y)) {							\
	yDims = INTEGER(getAttrib(y, R_DimSymbol));			\
	yDnms = getAttrib(y, R_DimNamesSymbol);				\
	y_has_dimNames = yDnms != R_NilValue;				\
    } else { /* ! matrix */ 						\
	yDims = INTEGER(yD = PROTECT(allocVector(INTSXP, 2))); nprot++;	\
	if(xDims[0] == 1) {						\
             /* "new" (2014-10-10): "be tolerant" as for R 3.2.0*/ 	\
	    yDims[0] = 1;						\
	    yDims[1] = LENGTH(y);					\
	} else {							\
	    yDims[0] = LENGTH(y);					\
	    yDims[1] = 1;						\
	}								\
	y_has_dimNames = FALSE;						\
    }									\
    int  n = yDims[!tr],/* (m,n) -> result dim */			\
	yd = yDims[ tr];/* (xd,yd): the conformable dims */		\
    if (xd != yd)							\
	error(_("Dimensions of x and y are not compatible for %s"),	\
	      tr ? "tcrossprod" : "crossprod");				\
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));		\
    vDims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));		\
    vDims[0] = m; vDims[1] = n;						\
    /* establish dimnames */						\
    SET_VECTOR_ELT(dn, 0,						\
		   duplicate(VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym), \
					tr ? 0 : 1)));			\
    if(y_has_dimNames)							\
	SET_VECTOR_ELT(dn, 1,						\
		       duplicate(VECTOR_ELT(yDnms, tr ? 0 : 1)));	\
    SET_SLOT(val, Matrix_DimNamesSym, dn);				\
									\
    double *v = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP,  m * (R_xlen_t) n))

#define DGE_MAT_CROSS_DO(_X_X_)						\
    if (xd > 0 && n > 0 && m > 0)					\
	F77_CALL(dgemm)(tr ? "N" : "T", tr ? "T" : "N", &m, &n, &xd, &one, \
			_X_X_, xDims, REAL(y), yDims, 			\
			&zero, v, &m FCONE FCONE);			\
    else								\
	Memzero(v,  m * (R_xlen_t) n);					\
    UNPROTECT(nprot);							\
    return val

    DGE_MAT_CROSS_1;
    DGE_MAT_CROSS_DO(REAL(GET_SLOT(x, Matrix_xSym)));
}

//! as dgeMatrix_matrix_crossprod() but x can be  [dln]geMatrix
static SEXP _geMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
    DGE_MAT_CROSS_1;

    double *x_x = gematrix_real_x(x, m * (R_xlen_t) xd);

    DGE_MAT_CROSS_DO(x_x);
}

SEXP geMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans) {
    SEXP dx = PROTECT(asdge(x, 0)),
	val = _geMatrix_matrix_crossprod(dx, y, trans);
    UNPROTECT(1);
    return val;
}

//  right = TRUE:  %*%  is called as  *(y, x, right=TRUE)
SEXP dgeMatrix_matrix_mm(SEXP a, SEXP bP, SEXP right)
{
#define DGE_MAT_MM_1(N_PROT)						\
    SEXP val= PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),		\
	 dn = PROTECT(allocVector(VECSXP, 2));				\
    int nprot = N_PROT + 2,						\
	*adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),			\
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),			\
	*cdims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2)),	\
	Rt = asLogical(right), m, k, n;					\
    double one = 1., zero = 0.;						\
									\
    if (Rt) { /* b %*% a : (m x k) (k x n) -> (m x n) */		\
	m = bdims[0]; k = bdims[1]; n = adims[1];			\
	if (adims[0] != k)						\
	    error(_("Matrices are not conformable for multiplication")); \
    } else {  /* a %*% b : (m x k) (k x n) -> (m x n) */		\
	m = adims[0]; k = adims[1]; n = bdims[1];			\
	if (bdims[0] != k)						\
	    error(_("Matrices are not conformable for multiplication")); \
    }									\
									\
    cdims[0] = m; cdims[1] = n;						\
    /* establish dimnames */						\
    SET_VECTOR_ELT(dn, 0, duplicate(					\
		       VECTOR_ELT(GET_SLOT(Rt ? b : a,			\
					   Matrix_DimNamesSym), 0)));	\
    SET_VECTOR_ELT(dn, 1,						\
		   duplicate(						\
		       VECTOR_ELT(GET_SLOT(Rt ? a : b,			\
					   Matrix_DimNamesSym), 1)));	\
    SET_SLOT(val, Matrix_DimNamesSym, dn);				\
    double *v = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * (R_xlen_t) n))

#define DGE_MAT_MM_DO(_A_X_, _B_X_)					\
    if (m < 1 || n < 1 || k < 1) {/* zero extent matrices should work */ \
	Memzero(v, m * (R_xlen_t) n);					\
    } else {								\
	if (Rt) { /* b %*% a  */					\
	    F77_CALL(dgemm) ("N", "N", &m, &n, &k, &one,		\
			     _B_X_, &m, _A_X_, &k, &zero, v, &m FCONE FCONE); \
	} else {  /* a %*% b  */					\
	    F77_CALL(dgemm) ("N", "N", &m, &n, &k, &one,		\
			     _A_X_, &m,	_B_X_, &k, &zero, v, &m FCONE FCONE); \
	}								\
    }									\
    UNPROTECT(nprot);							\
    return val

    SEXP b = PROTECT(asdge(bP, 0));
    DGE_MAT_MM_1(1);
    DGE_MAT_MM_DO(REAL(GET_SLOT(a, Matrix_xSym)),
                  REAL(GET_SLOT(b, Matrix_xSym)));
}

//! as dgeMatrix_matrix_mm() but a can be  [dln]geMatrix
static SEXP _geMatrix_matrix_mm(SEXP a, SEXP b, SEXP right)
{
    DGE_MAT_MM_1(0);
    double *a_x = gematrix_real_x(a, k * (R_xlen_t)(Rt ? n : m));
    double *b_x = gematrix_real_x(b, k * (R_xlen_t)(Rt ? m : n));
    DGE_MAT_MM_DO(a_x, b_x);
}

//! %*% -- generalized from dge to *ge():
SEXP geMatrix_matrix_mm(SEXP a, SEXP b, SEXP right) {
    SEXP
	da = PROTECT(asdge(a, 0)),
	db = PROTECT(asdge(b, 0)),
	val = _geMatrix_matrix_mm(da, db, right);
    UNPROTECT(2);
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
    SEXP val,/* = in case (2):  asdge(b, 0); */
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
	val = PROTECT(asdge(b, 0));
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
    SEXP val = PROTECT(asdge(b, 0));
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

SEXP dtpMatrix_matrix_mm(SEXP x, SEXP y, SEXP right, SEXP trans)
{
    SEXP val = PROTECT(asdge(y, 0));
    int rt = asLogical(right); // if(rt), compute b %*% op(a), else op(a) %*% b
    int tr = asLogical(trans); // if(tr), op(a) = t(a), else op(a) = a
    /* Since 'x' is square (n x n ),   dim(x %*% y) = dim(y) */
    int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDim = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int m = yDim[0], n = yDim[1];
    int ione = 1;
    const char *uplo = uplo_P(x), *diag = diag_P(x);
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)),
	*vx = REAL(GET_SLOT(val, Matrix_xSym));

    if (yDim[0] != xDim[1])
    if ((rt && xDim[0] != n) || (!rt && xDim[1] != m))
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      xDim[0], xDim[1], yDim[0], yDim[1]);
    if (m < 1 || n < 1) {
/* 	error(_("Matrices with zero extents cannot be multiplied")); */
    } else /* BLAS */
	// go via BLAS 2  dtpmv(.); there is no dtpmm in Lapack!
	if(rt) {
	    error(_("right=TRUE is not yet implemented __ FIXME"));
	} else {
	    for (int j = 0; j < n; j++) // X %*% y[,j]
		F77_CALL(dtpmv)(uplo, /*trans = */ tr ? "T" : "N",
				diag, yDim, xx,
				vx + j * (size_t) m, &ione FCONE FCONE FCONE);
	}
    UNPROTECT(1);
    return val;
}

/* FIXME: This function should be removed and a 'right' argument added to
 * dtpMatrix_matrix_mm -- also to be more parallel to ./dtrMatrix.c code */
SEXP dgeMatrix_dtpMatrix_mm(SEXP x, SEXP y)
{
    SEXP val = PROTECT(duplicate(x));
    /* Since 'y' is square (n x n ),   dim(x %*% y) = dim(x) */
    int *xDim = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDim = INTEGER(GET_SLOT(y, Matrix_DimSym));
    const char *uplo = uplo_P(y), *diag = diag_P(y);
    double *yx = REAL(GET_SLOT(y, Matrix_xSym)),
 	*vx = REAL(GET_SLOT(val, Matrix_xSym));

    if (yDim[0] != xDim[1])
	error(_("Dimensions of a (%d,%d) and b (%d,%d) do not conform"),
	      xDim[0], xDim[1], yDim[0], yDim[1]);
    for (int i = 0; i < xDim[0]; i++)/* val[i,] := Y' %*% x[i,]  */
	F77_CALL(dtpmv)(uplo, "T", diag, yDim, yx,
			vx + i, /* incr = */ xDim FCONE FCONE FCONE);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_matrix_mm(SEXP a, SEXP b)
{
    SEXP val = PROTECT(asdge(b, 0));
    int *bdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int i, ione = 1, n = bdims[0], nrhs = bdims[1];
    R_xlen_t nn = n * (R_xlen_t) nrhs;
    const char *uplo = uplo_P(a);
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)), one = 1., zero = 0.,
	*vx = REAL(GET_SLOT(val, Matrix_xSym)), *bx;

    if (bdims[0] != n)
	error(_("Matrices are not conformable for multiplication"));
    if (nrhs >= 1 && n >= 1) {
	Matrix_Calloc(bx, nn, double);
	Memcpy(bx, vx, nn);

	R_xlen_t in;
	for (i = 0, in = 0; i < nrhs; i++, in += n) { // in := i * n (w/o overflow!)
	    F77_CALL(dspmv)(uplo, &n, &one, ax, bx + in, &ione,
			    &zero, vx + in, &ione FCONE);
	}

	Matrix_Free(bx, nn);
    }
    UNPROTECT(1);
    return val;
}

SEXP dsyMatrix_matrix_mm(SEXP a, SEXP b, SEXP right)
{
    SEXP val = PROTECT(asdge(b, 0));// incl. dimnames
    int rt = asLogical(right); /* if(rt), compute b %*% a,  else  a %*% b */
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym)),
	m = bdims[0], n = bdims[1];

    if ((rt && n != adims[0]) || (!rt && m != adims[0]))
	error(_("Matrices are not conformable for multiplication"));

    double one = 1., zero = 0.;
    R_xlen_t mn = m * (R_xlen_t)n;
    double *bcp, *vx = REAL(GET_SLOT(val, Matrix_xSym));
    Matrix_Calloc(bcp, mn, double);
    Memcpy(bcp, vx, mn);

    if (m >=1 && n >= 1)
	F77_CALL(dsymm)(rt ? "R" :"L", uplo_P(a), &m, &n, &one,
			REAL(GET_SLOT(a, Matrix_xSym)), adims, bcp,
			&m, &zero, vx, &m FCONE FCONE);
    // add dimnames:
    int nd = rt ?
	1 : // v <- b %*% a : rownames(v) == rownames(b)  are already there
	0;  // v <- a %*% b : colnames(v) == colnames(b)  are already there
    SEXP nms = PROTECT(VECTOR_ELT(get_symmetrized_DimNames(a, -1), nd));
    SET_VECTOR_ELT(GET_SLOT(val, Matrix_DimNamesSym), nd, nms);
    Matrix_Free(bcp, mn);
    UNPROTECT(2);
    return val;
}

/**
 * All (dense * sparse)  Matrix products and cross products
 *
 *   f( f(<Csparse>)  %*%  f(<dense>) )   where  f ()  is either t () [tranpose] or the identity.
 *
 * @param a CsparseMatrix  (n x m)
 * @param b numeric vector, matrix, or denseMatrix (m x k) or (k x m)  if `trans` is '2' or 'B'
 * @param trans character.
 *        = " " : nothing transposed {apart from a}
 *        = "2" : "transpose 2nd arg": use  t(b) instead of b (= 2nd argument)
 *        = "c" : "transpose c":       Return  t(c) instead of c
 *        = "B" : "transpose both":    use t(b) and return t(c) instead of c
 * NB: For "2", "c", "B", need to transpose a *dense* matrix, B or C --> chm_transpose_dense()
 *
 * @return a dense matrix, the matrix product c = g(a,b) :
 *
 *                                                Condition (R)   Condition (C)
 *   R notation            Math notation          cross   trans   t.a t.b t.ans
 *   ~~~~~~~~~~~~~~~~~     ~~~~~~~~~~~~~~~~~~     ~~~~~~~~~~~~~   ~~~~~~~~~~~~~
 *   c <-   a %*%   b      C :=      A B            .       " "    .   .   .
 *   c <-   a %*% t(b)     C :=      A B'           .       "2"    .   |   .
 *   c <- t(a %*%   b)     C := (A B)'  = B'A'      .	    "c"    .   .   |
 *   c <- t(a %*% t(b))    C := (A B')' = B A'      .	    "B"    .   |   |
 *
 *   c <-   t(a) %*%   b   C :=      A'B           TRUE	    " "    |   .   .
 *   c <-   t(a) %*% t(b)  C :=      A'B'          TRUE	    "2"    |   |   .
 *   c <- t(t(a) %*%   b)  C := (A'B)'  = B'A      TRUE	    "c"    |   .   |
 *   c <- t(t(a) %*% t(b)) C := (A'B')' = B A      TRUE	    "B"    |   |   |
 */
SEXP Csp_dense_products(SEXP a, SEXP b,
                        Rboolean trans_a, Rboolean trans_b, Rboolean trans_ans)
{
	static const char *valid[] = { VALID_CSPARSE, "" };
	int ivalid = R_check_class_etc(a, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(a, __func__);
	const char *cl = valid[ivalid];

	if (cl[0] == 'n') {
#if 0
		warning(_("cholmod_sdmult() not yet implemented for pattern matrices -> coercing to double"));
#endif
		a = sparse_as_kind(a, cl, 'd');
	}
	PROTECT(a);

	CHM_SP cha = AS_CHM_SP(a);
	R_CheckStack();
	size_t
		a_nc = trans_a ? cha->nrow : cha->ncol,
		a_nr = trans_a ? cha->ncol : cha->nrow;

	Rboolean
		b_is_vector = !(IS_S4_OBJECT(b) || isMatrix(b)),
		b_transpose_if_vector = b_is_vector && XLENGTH(b) != a_nc;
	if (b_is_vector)
		trans_b = FALSE; /* don't transpose twice! */
	PROTECT(b = asdge(b, b_transpose_if_vector));

	CHM_DN chb = AS_CHM_DN(b);
	R_CheckStack();
	if (trans_b) {
		CHM_DN chbt = cholmod_allocate_dense(chb->ncol, chb->nrow, chb->ncol,
		                                     chb->xtype, &c);
		chm_transpose_dense(chbt, chb);
		chb = chbt;
	}
	size_t b_nc = chb->ncol;

	CHM_DN chc = cholmod_allocate_dense(a_nr, b_nc, a_nr, chb->xtype, &c);
	double one[] = {1.0, 0.0}, zero[] = {0.0, 0.0};
	cholmod_sdmult(cha, trans_a, one, zero, chb, chc, &c);

	SEXP dna = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
		dnb = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		dnc = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dnc, (trans_ans) ? 1 : 0, VECTOR_ELT(dna, (trans_a) ? 1 : 0));
	SET_VECTOR_ELT(dnc, (trans_ans) ? 0 : 1, VECTOR_ELT(dnb, (trans_b) ? 0 : 1));

	if (trans_b)
		cholmod_free_dense(&chb, &c);
	SEXP ans = chm_dense_to_SEXP(chc, 1, 0, dnc, trans_ans);
	UNPROTECT(5);
	return ans;
}

/** @brief  A %*% B  - for matrices of class CsparseMatrix (R package "Matrix")
 *
 * @param a
 * @param b
 * @param bool_arith
 *
 * @return
 *
 * NOTA BENE:  cholmod_ssmult(A,B, ...) ->  ./CHOLMOD/MatrixOps/cholmod_ssmult.c
 * ---------  computes a patter*n* matrix __always_ when
 * *one* of A or B is pattern*n*, because of this (line 73-74):
 * ---------------------------------------------------------------------------
 * values = values &&
 * (A->xtype != CHOLMOD_PATTERN) && (B->xtype != CHOLMOD_PATTERN) ;
 * ---------------------------------------------------------------------------
 * ==> Often need to copy the patter*n* to a *l*ogical matrix first !!!
 */
SEXP Csparse_Csparse_prod(SEXP a, SEXP b, SEXP boolArith)
{
	static const char *valid[] = { VALID_CSPARSE, "" };
	int ivalid = R_check_class_etc(a, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(a, __func__);
	const char *acl = valid[ivalid];
	ivalid = R_check_class_etc(b, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(b, __func__);
	const char *bcl = valid[ivalid];

	int doBool = asLogical(boolArith);
	if (doBool == NA_LOGICAL)
		doBool = (acl[0] == 'n' && bcl[0] == 'n');
	if (doBool) {
		if (acl[0] != 'n')
			a = sparse_as_kind(a, acl, 'n');
		PROTECT(a);
		if (bcl[0] != 'n')
			b = sparse_as_kind(b, bcl, 'n');
		PROTECT(b);
	} else {
		if (acl[0] != 'd')
			a = sparse_as_kind(a, acl, 'd');
		PROTECT(a);
		if (bcl[0] != 'd')
			b = sparse_as_kind(b, bcl, 'd');
		PROTECT(b);
	}

	CHM_SP cha = AS_CHM_SP(a), chb = AS_CHM_SP(b), chc;
	R_CheckStack();
	chc = cholmod_ssmult(cha, chb, 0, !doBool, 1, &c);

	char ul = '\0', di = '\0';
	if (acl[1] == 't' && bcl[1] == 't') {
		SEXP auplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
			buplo = PROTECT(GET_SLOT(b, Matrix_uploSym));
		char aul = *CHAR(STRING_ELT(auplo, 0)),
			bul = *CHAR(STRING_ELT(buplo, 0));
		if (aul == bul) {
			ul = aul;
			di = 'N';
			SEXP adiag = PROTECT(GET_SLOT(a, Matrix_diagSym)),
				bdiag = PROTECT(GET_SLOT(b, Matrix_diagSym));
			char adi = *CHAR(STRING_ELT(adiag, 0)),
				bdi = *CHAR(STRING_ELT(bdiag, 0));
			if (adi != 'N' && bdi != 'N') {
				di = 'U';
				chm_diagN2U(chc, (ul == 'U') ? 1 : -1, 0);
			}
			UNPROTECT(2);
		}
		UNPROTECT(2);
	}

	SEXP
		dna = PROTECT((acl[1] != 's')
		              ? GET_SLOT(a, Matrix_DimNamesSym)
		              : get_symmetrized_DimNames(a, -1)),
		dnb = PROTECT((bcl[1] != 's')
		              ? GET_SLOT(b, Matrix_DimNamesSym)
		              : get_symmetrized_DimNames(b, -1)),
		dnc = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dnc, 0, VECTOR_ELT(dna, 0));
	SET_VECTOR_ELT(dnc, 1, VECTOR_ELT(dnb, 1));

	SEXP ans = chm_sparse_to_SEXP(chc,
	                              1,
	                              (ul == '\0') ? 0 : ((ul == 'U') ? 1 : -1),
	                              0,
	                              (di == '\0') ? "" : ((di == 'N') ? "N" : "U"),
	                              dnc);
	UNPROTECT(5);
	return ans;
}

/** @brief [t]crossprod (<Csparse>, <Csparse>)
 *
 * @param a a "CsparseMatrix" object
 * @param b a "CsparseMatrix" object
 * @param trans trans = FALSE:  crossprod(a,b)
 *              trans = TRUE : tcrossprod(a,b)
 * @param bool_arith logical (TRUE / NA / FALSE): Should boolean arithmetic be used.
 *
 * @return a CsparseMatrix, the (t)cross product of a and b.
 */
SEXP Csparse_Csparse_crossprod(SEXP a, SEXP b, SEXP trans, SEXP boolArith)
{
	static const char *valid[] = { VALID_CSPARSE, "" };
	int ivalid = R_check_class_etc(a, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(a, __func__);
	const char *acl = valid[ivalid];
	ivalid = R_check_class_etc(b, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(b, __func__);
	const char *bcl = valid[ivalid];

	int doTrans = asLogical(trans), doBool = asLogical(boolArith);
	if (doBool == NA_LOGICAL)
		doBool = (acl[0] == 'n' && bcl[0] == 'n');
	if (doBool) {
		if (acl[0] != 'n')
			a = sparse_as_kind(a, acl, 'n');
		PROTECT(a);
		if (bcl[0] != 'n')
			b = sparse_as_kind(b, bcl, 'n');
		PROTECT(b);
	} else {
		if (acl[0] != 'd')
			a = sparse_as_kind(a, acl, 'd');
		PROTECT(a);
		if (bcl[0] != 'd')
			b = sparse_as_kind(b, bcl, 'd');
		PROTECT(b);
	}

	CHM_SP cha = AS_CHM_SP(a), chb = AS_CHM_SP(b), chc;
	R_CheckStack();
	if (doTrans)
		chb = cholmod_transpose(chb, chb->xtype, &c);
	else
		cha = cholmod_transpose(cha, cha->xtype, &c);
	chc = cholmod_ssmult(cha, chb, 0, !doBool, 1, &c);
	if (doTrans)
		cholmod_free_sparse(&chb, &c);
	else
		cholmod_free_sparse(&cha, &c);

	char ul = '\0', di = '\0';
	if (acl[1] == 't' && bcl[1] == 't') {
		SEXP auplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
			buplo = PROTECT(GET_SLOT(b, Matrix_uploSym));
		char aul = *CHAR(STRING_ELT(auplo, 0)),
			bul = *CHAR(STRING_ELT(buplo, 0));
		if (aul != bul) {
			ul = (doTrans) ? aul : bul;
			di = 'N';
			SEXP adiag = PROTECT(GET_SLOT(a, Matrix_diagSym)),
				bdiag = PROTECT(GET_SLOT(b, Matrix_diagSym));
			char adi = *CHAR(STRING_ELT(adiag, 0)),
				bdi = *CHAR(STRING_ELT(bdiag, 0));
			if (adi != 'N' && bdi != 'N') {
				di = 'U';
				chm_diagN2U(chc, (ul == 'U') ? 1 : -1, 0);
			}
			UNPROTECT(2);
		}
		UNPROTECT(2);
	}

	SEXP
		dna = PROTECT((acl[1] != 's')
		              ? GET_SLOT(a, Matrix_DimNamesSym)
		              : get_symmetrized_DimNames(a, -1)),
		dnb = PROTECT((bcl[1] != 's')
		              ? GET_SLOT(b, Matrix_DimNamesSym)
		              : get_symmetrized_DimNames(b, -1)),
		dnc = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dnc, 0, VECTOR_ELT(dna, (doTrans) ? 0 : 1));
	SET_VECTOR_ELT(dnc, 1, VECTOR_ELT(dnb, (doTrans) ? 0 : 1));

	SEXP ans = chm_sparse_to_SEXP(chc,
	                              1,
	                              (ul == '\0') ? 0 : ((ul == 'U') ? 1 : -1),
	                              0,
	                              (di == '\0') ? "" : ((di == 'N') ? "N" : "U"),
	                              dnc);
	UNPROTECT(5);
	return ans;
}

/** @brief Computes   x'x  or  x x' -- *also* for Tsparse
 *  see Csparse_Csparse_crossprod above for  x'y and x y'
 */
SEXP Csparse_crossprod(SEXP x, SEXP trans, SEXP boolArith)
{
	static const char *valid[] = { VALID_CSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	const char *cl = valid[ivalid];

	int doTrans = asLogical(trans), doBool = asLogical(boolArith);
	if (doBool == NA_LOGICAL)
		doBool = cl[0] == 'n';
	if (doBool) {
		if (cl[0] != 'n')
			x = sparse_as_kind(x, cl, 'n');
	} else {
		if (cl[0] != 'd')
			x = sparse_as_kind(x, cl, 'd');
	}
	PROTECT(x);

	int doFree = 0;
	CHM_SP chx, chc;
	if (cl[2] != 'T') {
		chx = AS_CHM_SP(x);
		R_CheckStack();
	} else {
		/* defined in ./sparse.c : */
		SEXP sparse_diag_U2N(SEXP, const char *);
		x = sparse_diag_U2N(x, cl); /* work around as_cholmod_triplet (?) */
		UNPROTECT(1);
		PROTECT(x);
		CHM_TR tmp = AS_CHM_TR__(x);
		R_CheckStack();
		chx = cholmod_triplet_to_sparse(tmp, tmp->nnz, &c);
		doFree = 1;
	}
	if (!doTrans) {
		CHM_SP tmp = cholmod_transpose(chx, chx->xtype, &c);
		if (doFree)
			cholmod_free_sparse(&chx, &c);
		else doFree = 1;
		chx = tmp;
	}
	if (chx->stype != 0) {
		CHM_SP tmp = cholmod_copy(chx, 0, chx->xtype, &c);
		if (doFree)
			cholmod_free_sparse(&chx, &c);
		else doFree = 1;
		chx = tmp;
	}
	chc = cholmod_aat(chx, (int *) NULL, 0, chx->xtype, &c);
	if (doFree)
		cholmod_free_sparse(&chx, &c);
	chc->stype = 1;

	SEXP
		dnx = PROTECT((cl[1] != 's')
		              ? GET_SLOT(x, Matrix_DimNamesSym)
		              : get_symmetrized_DimNames(x, -1)),
		dnc = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dnc, 0, VECTOR_ELT(dnx, (doTrans) ? 0 : 1));
	SET_VECTOR_ELT(dnc, 1, VECTOR_ELT(dnx, (doTrans) ? 0 : 1));

	SEXP ans = chm_sparse_to_SEXP(chc, 1, 0, 0, "", dnc);
	UNPROTECT(3);
	return ans;
}

SEXP Csparse_dense_prod(SEXP a, SEXP b, SEXP trans)
{
    return
	Csp_dense_products(a, b,
		/* trans_a = */ FALSE,
		/* trans_b   = */ (*CHAR(asChar(trans)) == '2' || *CHAR(asChar(trans)) == 'B'),
		/* trans_ans = */ (*CHAR(asChar(trans)) == 'c' || *CHAR(asChar(trans)) == 'B'));
}

SEXP Csparse_dense_crossprod(SEXP a, SEXP b, SEXP trans)
{
    return
	Csp_dense_products(a, b,
		/* trans_a = */ TRUE,
		/* trans_b   = */ (*CHAR(asChar(trans)) == '2' || *CHAR(asChar(trans)) == 'B'),
		/* trans_ans = */ (*CHAR(asChar(trans)) == 'c' || *CHAR(asChar(trans)) == 'B'));
}
