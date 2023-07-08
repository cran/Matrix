/** @file Csparse.c
 * The "CsparseMatrix" class from R package Matrix:
 *
 * Sparse matrices in compressed column-oriented form
 */
#include "Csparse.h"
#include "chm_common.h"
#include "cs_utils.h" /* -> ./cs.h  for cs_dmperm() */

#define _t_Csparse_validate
#include "t_Csparse_validate.c"

#define _t_Csparse_sort
#include "t_Csparse_validate.c"

// R: .validateCsparse(x, sort.if.needed = FALSE) :
SEXP Csparse_validate2(SEXP x, SEXP maybe_modify) {
    return Csparse_validate_(x, asLogical(maybe_modify));
}

// R: Matrix:::.sortCsparse(x) :
SEXP Csparse_sort (SEXP x) {
   int ok = Csparse_sort_2(x, TRUE); // modifying x directly
   if(!ok) warning(_("Csparse_sort(x): x is not a valid (apart from sorting) CsparseMatrix"));
   return x;
}

/** "Cheap" C version of  Csparse_validate() - *not* sorting : */
Rboolean isValid_Csparse(SEXP x)
{
    /* NB: we do *NOT* check a potential 'x' slot here, at all */
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym);
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)), j,
	nrow = dims[0],
	ncol = dims[1],
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    if (length(pslot) != dims[1] + 1)
	return FALSE;
    if (xp[0] != 0)
	return FALSE;
    if (length(islot) < xp[ncol]) /* allow larger slots from over-allocation!*/
	return FALSE;
    for (j = 0; j < xp[ncol]; j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return FALSE;
    }
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j + 1])
	    return FALSE;
    }
    return TRUE;
}

/** @brief Csparse_drop(x, tol):  drop entries with absolute value < tol, i.e,
 *  at least all "explicit" zeros. */
SEXP Csparse_drop(SEXP x, SEXP tol)
{
    const char *cl = class_P(x);
    /* dtCMatrix, etc; [1] = the second character =?= 't' for triangular */
    int tr = (cl[1] == 't'); // FIXME - rather  R_check_class_etc(..)
    CHM_SP chx = AS_CHM_SP__(x);
    CHM_SP ans = cholmod_copy(chx, chx->stype, chx->xtype, &c);
    double dtol = asReal(tol);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    if(!cholmod_drop(dtol, ans, &c))
	error(_("cholmod_drop() failed"));
   return chm_sparse_to_SEXP(ans, 1,
			      tr ? ((*uplo_P(x) == 'U') ? 1 : -1) : 0,
			      Rkind, tr ? diag_P(x) : "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

/** @brief Horizontal Concatenation -  cbind( <Csparse>,  <Csparse>)
 */
SEXP Csparse_horzcat(SEXP x, SEXP y)
{
#define CSPARSE_CAT(_KIND_)						\
    CHM_SP chx = AS_CHM_SP__(x), chy = AS_CHM_SP__(y);			\
    R_CheckStack();							\
    void* chx_x = chx->x;						\
    void* chx_z = chx->z;						\
    void* chy_x = chy->x;						\
    void* chy_z = chy->z;						\
    int Rk_x = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : x_pattern, \
	Rk_y = (chy->xtype != CHOLMOD_PATTERN) ? Real_kind(y) : x_pattern, Rkind; \
    if(Rk_x == x_pattern || Rk_y == x_pattern) { /* at least one of them is patter"n" */ \
	if(Rk_x == x_pattern && Rk_y == x_pattern) { /* fine */		\
	} else { /* only one is a patter"n"				\
		  * "Bug" in cholmod_horzcat()/vertcat():               \
		  * returns patter"n" matrix if one of them is */	\
	    Rboolean ok;						\
	    if(Rk_x == x_pattern) {					\
		ok = chm_MOD_xtype(CHOLMOD_REAL, chx, &c); Rk_x = 0;	\
	    } else if(Rk_y == x_pattern) {				\
		ok = chm_MOD_xtype(CHOLMOD_REAL, chy, &c); Rk_y = 0;	\
	    } else							\
		error(_("Impossible Rk_x/Rk_y in Csparse_%s(), please report"), _KIND_); \
	    if(!ok)							\
		error(_("chm_MOD_xtype() was not successful in Csparse_%s(), please report"), \
		      _KIND_);						\
	}								\
    }									\
    Rkind = /* logical if both x and y are */ (Rk_x == 1 && Rk_y == 1) ? 1 : 0

    CSPARSE_CAT("horzcat");
    // TODO: currently drops dimnames - and we fix at R level;

    SEXP retval = PROTECT(
	chm_sparse_to_SEXP(cholmod_horzcat(chx, chy, 1, &c),
			   1, 0, Rkind, "", R_NilValue));
/* AS_CHM_SP(x) fills result with points to R-allocated memory but
   chm_MOD_xtype can change ->x and ->z to cholmod_alloc'ed memory.
   The former needs no freeing but the latter does.
   The first 2 arguments to cholmod_free should contain the number
   and size of things being freed, but lying about that is sort of ok. */
#define CSPARSE_CAT_CLEANUP					\
    if (chx_x != chx->x) cholmod_free(0, 0, chx->x, &c);	\
    if (chx_z != chx->z) cholmod_free(0, 0, chx->z, &c);	\
    if (chy_x != chy->x) cholmod_free(0, 0, chy->x, &c);	\
    if (chy_z != chy->z) cholmod_free(0, 0, chy->z, &c);	\
    UNPROTECT(1);

    CSPARSE_CAT_CLEANUP;
    return retval;
}

/** @brief Vertical Concatenation -  rbind( <Csparse>,  <Csparse>)
 */
SEXP Csparse_vertcat(SEXP x, SEXP y)
{
    CSPARSE_CAT("vertcat");
    // TODO: currently drops dimnames - and we fix at R level;

    SEXP retval  = PROTECT(
	chm_sparse_to_SEXP(cholmod_vertcat(chx, chy, 1, &c),
			   1, 0, Rkind, "", R_NilValue));
    CSPARSE_CAT_CLEANUP;
    return retval;
}

#define _d_Csp_
#include "t_Csparse_subassign.c"

#define _l_Csp_
#include "t_Csparse_subassign.c"

#define _i_Csp_
#include "t_Csparse_subassign.c"

#define _n_Csp_
#include "t_Csparse_subassign.c"

#define _z_Csp_
#include "t_Csparse_subassign.c"

SEXP Csparse_MatrixMarket(SEXP x, SEXP fname)
{
    FILE *f = fopen(CHAR(asChar(fname)), "w");

    if (!f)
	error(_("failure to open file \"%s\" for writing"),
	      CHAR(asChar(fname)));
    if (!cholmod_write_sparse(f, AS_CHM_SP(x),
			      (CHM_SP)NULL, (char*) NULL, &c))
	error(_("cholmod_write_sparse returned error code"));
    fclose(f);
    return R_NilValue;
}

// seed will *not* be used unless it's -1 (inverse perm.) or  0 ("no" / identity) perm.
static csd* Csparse_dmperm_raw(SEXP mat, SEXP seed)
{
    mat = PROTECT(duplicate(mat));
    CSP matx = AS_CSP__(mat); /* m x n ; compressed column, *double* 'x' or none */
    int iseed = asInteger(seed);
    R_CheckStack();
    UNPROTECT(1);
    return cs_dmperm(matx, iseed); // -> ./cs.c
}

/* NB:  cs.h  defines the 'csd' struct as  (NB: csi :== int in  Matrix, for now)

   typedef struct cs_dmperm_results    // cs_dmperm or cs_scc output
   {
   csi *p ;        // size m, row permutation
   csi *q ;        // size n, column permutation
   csi *r ;        // size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q)
   csi *s ;        // size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q)
   csi nb ;        // # of blocks in fine dmperm decomposition
   csi rr [5] ;    // coarse row decomposition
   csi cc [5] ;    // coarse column decomposition
   } csd ;
*/

/* MM: should allow to return the full info above
   (Timothy Davis, p.126, explains why it's interesting ..) */

/* Here, return the full *named* list to R */
SEXP Csparse_dmperm(SEXP mat, SEXP seed, SEXP nAns) {
    csd *DMp = Csparse_dmperm_raw(mat, seed);
    if(DMp == NULL) // "failure" in cs_dmperm()
	return(R_NilValue);
    int *dims = INTEGER(GET_SLOT(mat, Matrix_DimSym)),
	m = dims[0],
	n = dims[1],
	n_ans = asInteger(nAns),
	nb = DMp->nb;

    SEXP nms = PROTECT(allocVector(STRSXP, n_ans));
    SEXP ans = PROTECT(allocVector(VECSXP, n_ans));
    R_CheckStack();
    int *ip;
    /* p : */SET_STRING_ELT(nms, 0, mkChar("p"));
             SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, m));
    ip = INTEGER(VECTOR_ELT(ans, 0));
    /* 0-based permutation:
     * Memcpy(ip , (int*)(DMp->p), m); */
    // 1-based permutation:
    for(int i=0; i < m; i++) ip[i] = DMp->p[i] + 1;

    /* q : */SET_STRING_ELT(nms, 1, mkChar("q"));
             SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, n));
    ip = INTEGER(VECTOR_ELT(ans, 1));
    /* 0-based permutation:
     * Memcpy(ip , (int*)(DMp->q), m); */
    // 1-based permutation:
    for(int i=0; i < n; i++) ip[i] = DMp->q[i] + 1;

    if(n_ans > 2) {
      /* r : */  SET_STRING_ELT(nms, 2, mkChar("r"));
		 SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, nb+1));
      Memcpy(INTEGER(VECTOR_ELT(ans, 2)), (int*)(DMp->r),   nb+1);

      /* s : */  SET_STRING_ELT(nms, 3, mkChar("s"));
		 SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, nb+1));
      Memcpy(INTEGER(VECTOR_ELT(ans, 3)), (int*)(DMp->s),   nb+1);
      if(n_ans > 4) {
	/* rr5 :*/ SET_STRING_ELT(nms, 4, mkChar("rr5"));
		   SET_VECTOR_ELT(ans, 4, allocVector(INTSXP, 5));
	Memcpy(INTEGER(VECTOR_ELT(ans, 4)), (int*)(DMp->rr),  5);

	/* cc5 :*/ SET_STRING_ELT(nms, 5, mkChar("cc5"));
		   SET_VECTOR_ELT(ans, 5, allocVector(INTSXP, 5));
	Memcpy(INTEGER(VECTOR_ELT(ans, 5)), (int*)(DMp->cc),  5);
      }
    }
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/**
 * Extract the diagonal entries from *triangular* Csparse matrix  __or__ a
 * cholmod_sparse factor (LDL = TRUE).
 *
 * @param n  dimension of the matrix.
 * @param x_p  'p' (column pointer) slot contents
 * @param x_x  'x' (non-zero entries) slot contents
 * @param perm 'perm' (= permutation vector) slot contents; only used for "diagBack"
 * @param resultKind a (SEXP) string indicating which kind of result is desired.
 *
 * @return  a SEXP, either a (double) number or a length n-vector of diagonal entries
 */
static
SEXP diag_tC_ptr(int n, int *x_p, double *x_x, Rboolean is_U, int *perm,
/*                                ^^^^^^ FIXME[Generalize] to int / ... -- via x_slot_kind ? */
		 SEXP resultKind)
{
    const char* res_ch = CHAR(STRING_ELT(resultKind,0));
    enum diag_kind { diag, diag_backpermuted, trace, prod, sum_log, min, max, range
    } res_kind = ((!strcmp(res_ch, "trace")) ? trace :
		  ((!strcmp(res_ch, "sumLog")) ? sum_log :
		   ((!strcmp(res_ch, "prod")) ? prod :
		    ((!strcmp(res_ch, "min")) ? min :
		     ((!strcmp(res_ch, "max")) ? max :
		      ((!strcmp(res_ch, "range")) ? range :
		       ((!strcmp(res_ch, "diag")) ? diag :
			((!strcmp(res_ch, "diagBack")) ? diag_backpermuted :
			 -1))))))));
    int i, n_x, i_from;
    SEXP ans = PROTECT(allocVector(REALSXP,
/*                                 ^^^^  FIXME[Generalize] */
				   (res_kind == diag ||
				    res_kind == diag_backpermuted) ? n :
				   (res_kind == range ? 2 : 1)));
    double *v = REAL(ans);
/*  ^^^^^^      ^^^^  FIXME[Generalize] */

    i_from = (is_U ? -1 : 0);

#define for_DIAG(v_ASSIGN)					\
    for(i = 0; i < n; i++) {					\
	/* looking at i-th column */				\
	n_x = x_p[i+1] - x_p[i];/* #{entries} in this column */	\
	if( is_U) i_from += n_x;                                \
	v_ASSIGN;						\
	if(!is_U) i_from += n_x;                                \
    }

    /* NOTA BENE: we assume  -- uplo = "L" i.e. lower triangular matrix
     *            for uplo = "U" (makes sense with a "dtCMatrix" !),
     *            should use  x_x[i_from + (n_x - 1)] instead of x_x[i_from],
     *            where n_x = (x_p[i+1] - x_p[i])
     */

    switch(res_kind) {
    case trace: // = sum
	v[0] = 0.;
	for_DIAG(v[0] += x_x[i_from]);
	break;

    case sum_log:
	v[0] = 0.;
	for_DIAG(v[0] += log(x_x[i_from]));
	break;

    case prod:
	v[0] = 1.;
	for_DIAG(v[0] *= x_x[i_from]);
	break;

    case min:
	v[0] = R_PosInf;
	for_DIAG(if(v[0] > x_x[i_from]) v[0] = x_x[i_from]);
	break;

    case max:
	v[0] = R_NegInf;
	for_DIAG(if(v[0] < x_x[i_from]) v[0] = x_x[i_from]);
	break;

    case range:
	v[0] = R_PosInf;
	v[1] = R_NegInf;
	for_DIAG(if(v[0] > x_x[i_from]) v[0] = x_x[i_from];
		 if(v[1] < x_x[i_from]) v[1] = x_x[i_from]);
	break;

    case diag:
	for_DIAG(v[i] = x_x[i_from]);
	break;

    case diag_backpermuted:
	for_DIAG(v[i] = x_x[i_from]);

	warning(_("%s = '%s' (back-permuted) is experimental"),
		"resultKind", "diagBack");
	/* now back_permute : */
	for(i = 0; i < n; i++) {
	    double tmp = v[i]; v[i] = v[perm[i]]; v[perm[i]] = tmp;
	    /*^^^^ FIXME[Generalize] */
	}
	break;

    default: /* -1 from above */
	error(_("diag_tC(): invalid 'resultKind'"));
	/* Wall: */ ans = R_NilValue; v = REAL(ans);
    }

    UNPROTECT(1);
    return ans;
}

/**
 * Extract the diagonal entries from *triangular* Csparse matrix  __or__ a
 * cholmod_sparse factor (LDL = TRUE).
 *
 * @param obj -- now a cholmod_sparse factor or a dtCMatrix
 * @param pslot  'p' (column pointer)   slot of Csparse matrix/factor
 * @param xslot  'x' (non-zero entries) slot of Csparse matrix/factor
 * @param perm_slot  'perm' (= permutation vector) slot of corresponding CHMfactor;
 *		     only used for "diagBack"
 * @param resultKind a (SEXP) string indicating which kind of result is desired.
 *
 * @return  a SEXP, either a (double) number or a length n-vector of diagonal entries
 */
SEXP diag_tC(SEXP obj, SEXP resultKind)
{
    SEXP
	pslot = GET_SLOT(obj, Matrix_pSym),
	xslot = GET_SLOT(obj, Matrix_xSym);
    Rboolean is_U = (R_has_slot(obj, Matrix_uploSym) &&
		     *CHAR(asChar(GET_SLOT(obj, Matrix_uploSym))) == 'U');
    int n = length(pslot) - 1, /* n = ncol(.) = nrow(.) */
	*x_p  = INTEGER(pslot), pp = -1, *perm;
    double *x_x = REAL(xslot);
/*  ^^^^^^        ^^^^ FIXME[Generalize] to INTEGER(.) / LOGICAL(.) / ... xslot !*/

    if(R_has_slot(obj, Matrix_permSym))
	perm = INTEGER(GET_SLOT(obj, Matrix_permSym));
    else perm = &pp;

    return diag_tC_ptr(n, x_p, x_x, is_U, perm, resultKind);
}
