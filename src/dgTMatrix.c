#include <Rinternals.h>
/* for R_LEN... */

#include "dgTMatrix.h"

#include "chm_common.h"
#include "Tsparse.h"

SEXP xTMatrix_validate(SEXP x)
{
    /* Almost everything now in Tsparse_validate ( ./Tsparse.c )
     * *but* the checking of the 'x' slot : */
    if (LENGTH(GET_SLOT(x, Matrix_iSym)) !=
	LENGTH(GET_SLOT(x, Matrix_xSym)))
	return mkString(_("lengths of slots i and x must match"));
    return ScalarLogical(1);
}

static void
d_insert_triplets_in_array(int m, int n, int nnz,
			   const int xi[], const int xj[], const double xx[],
			   /* --> */ double vx[])
{
    // For ( m*n ) > INT_MAX,  we here assume that size_t is using 64-bit !
    size_t m_ = (size_t) m, len = sizeof(double) * m_ * n;
    if(len == sizeof(double) * (double)m_ *n)
	memset(vx, 0, len);
    else { // len did overflow -- this should call memset() several times:
	size_t max_l = (1 << (sizeof(size_t)-1)); // = 2^(N-1)
	max_l += ((long)max_l - 1); // = 2^(N-1) + 2^(N-1) - 1 =  2^N - 1
	double dlen = ((double)m_) * n;
	if(dlen > max_l)
	    error(_("too large matrix: %.0f"), dlen);
	// else :   m * n does fit -- call memset() several times:
	dlen *= sizeof(double);
	memset(vx, 0, max_l);
	double len_done = 0.; // length in bytes // but also need length in double
	while((len_done += max_l) < dlen) {
	    size_t dd = (dlen - len_done < max_l) ? (size_t)(dlen - len_done) : max_l;
	    memset(vx + (int)(len_done/sizeof(double)), 0, dd);
	}
    }
    for (int i = 0; i < nnz; i++) {
	vx[xi[i] + xj[i] * m_] += xx[i]; /* allow redundant entries in x */
    }
}

static void
l_insert_triplets_in_array(int m, int n, int nnz,
			   const int xi[], const int xj[], const int xx[],
			   /* --> */ int vx[])
{
    // For ( m*n ) > INT_MAX,  we here assume that size_t is using 64-bit !
    size_t m_ = (size_t) m, len = sizeof(int) * m_ * n;
    if(len == sizeof(int) * (double)m_ *n)
	memset(vx, 0, len);
    else { // len did overflow -- this should call memset() several times:
	size_t max_l = (1 << (sizeof(size_t)-1)); // = 2^(N-1)
	max_l += ((long)max_l - 1); // = 2^(N-1) + 2^(N-1) - 1 =  2^N - 1
	double dlen = ((double)m_) * n;
	if(dlen > max_l)
	    error(_("too large matrix: %.0f"), dlen);
	// else :   m * n does fit -- call memset() several times:
	dlen *= sizeof(int);
	memset(vx, 0, max_l);
	double len_done = 0.; // length in bytes // but also need length in int
	while((len_done += max_l) < dlen) {
	    size_t dd = (dlen - len_done < max_l) ? (size_t)(dlen - len_done) : max_l;
	    memset(vx + (int)(len_done/sizeof(int)), 0, dd);
	}
    }

    for (int i = 0; i < nnz; i++) {
	size_t ind = xi[i] + xj[i] * m_;
	if(vx[ind] == NA_LOGICAL) {
	    // do nothing: remains NA
	} else if(xx[i] == NA_LOGICAL)
	    vx[ind] = NA_LOGICAL;
	else // "or" :
	    vx[ind] |= xx[i];
    }
}

#define MAKE_gTMatrix_to_geMatrix(_t1_, _SEXPTYPE_, _SEXP_)		\
SEXP _t1_ ## gTMatrix_to_ ## _t1_ ## geMatrix(SEXP x)			\
{									\
    SEXP dd = GET_SLOT(x, Matrix_DimSym),				\
	islot = GET_SLOT(x, Matrix_iSym),				\
	ans = PROTECT(NEW_OBJECT_OF_CLASS(#_t1_ "geMatrix"));	\
									\
    int *dims = INTEGER(dd),						\
	m = dims[0],							\
	n = dims[1];							\
    double len = m * (double)n;						\
									\
    if (len > R_XLEN_T_MAX)						\
	error(_("Cannot coerce to too large *geMatrix with %.0f entries"), \
              len);							\
									\
    SET_SLOT(ans, Matrix_factorSym, allocVector(VECSXP, 0));		\
    SET_SLOT(ans, Matrix_DimSym, duplicate(dd));			\
    SET_DimNames(ans, x);						\
    SET_SLOT(ans, Matrix_xSym, allocVector(_SEXPTYPE_, (R_xlen_t)len));	\
    _t1_ ## _insert_triplets_in_array(m, n, length(islot),		\
				      INTEGER(islot),			\
				      INTEGER(GET_SLOT(x, Matrix_jSym)),\
				      _SEXP_(GET_SLOT(x, Matrix_xSym)),	\
				      _SEXP_(GET_SLOT(ans, Matrix_xSym))); \
    UNPROTECT(1);							\
    return ans;								\
}

MAKE_gTMatrix_to_geMatrix(d, REALSXP, REAL)

MAKE_gTMatrix_to_geMatrix(l, LGLSXP, LOGICAL)

#undef MAKE_gTMatrix_to_geMatrix

#define MAKE_gTMatrix_to_matrix(_t1_, _SEXPTYPE_, _SEXP_)		\
SEXP _t1_ ## gTMatrix_to_matrix(SEXP x)					\
{									\
    SEXP dd = GET_SLOT(x, Matrix_DimSym),				\
	dn = GET_SLOT(x, Matrix_DimNamesSym),				\
	islot = GET_SLOT(x, Matrix_iSym);				\
    int m = INTEGER(dd)[0],						\
	n = INTEGER(dd)[1];						\
    SEXP ans = PROTECT(allocMatrix(_SEXPTYPE_, m, n));			\
    if(VECTOR_ELT(dn, 0) != R_NilValue || VECTOR_ELT(dn, 1) != R_NilValue) \
	/* matrix() with non-trivial dimnames */			\
	setAttrib(ans, R_DimNamesSymbol, duplicate(dn));		\
    _t1_ ## _insert_triplets_in_array(m, n, length(islot),		\
				      INTEGER(islot),			\
				      INTEGER(GET_SLOT(x, Matrix_jSym)),\
				      _SEXP_(GET_SLOT(x, Matrix_xSym)),	\
				      _SEXP_(ans));			\
    UNPROTECT(1);							\
    return ans;								\
}

MAKE_gTMatrix_to_matrix(d, REALSXP, REAL)

MAKE_gTMatrix_to_matrix(l, LGLSXP, LOGICAL)

#undef MAKE_gTMatrix_to_matrix
