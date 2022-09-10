#ifndef MATRIX_MUTILS_H
#define MATRIX_MUTILS_H

#undef Matrix_with_SPQR

#undef HAVE_PROPER_IMATRIX
#undef HAVE_PROPER_ZMATRIX

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h> /* C99 for int64_t */
#include <ctype.h>
#include <R.h> /* includes <Rconfig.h> */
#include <Rversion.h>
#include <Rinternals.h>
#include <R_ext/RS.h> /* for Memzero() */

/* NB: For 'USE_FC_LEN_T' and 'FCONE' (for LTO),
   the "includer" will #include "Lapack-etc.h" */

/* From <Defn.h> : */
/* 'alloca' is neither C99 nor POSIX */
#ifdef __GNUC__
/* This covers GNU, Clang and Intel compilers */
/* #undef needed in case some other header, e.g. malloc.h, already did this */
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
# ifdef HAVE_ALLOCA_H
/* This covers native compilers on Solaris and AIX */
#  include <alloca.h>
# endif
/* It might have been defined via some other standard header, e.g. stdlib.h */
# if !HAVE_DECL_ALLOCA
extern void *alloca(size_t);
# endif
#endif

#define Matrix_Calloc_Threshold 10000 /* R uses same cutoff in several places */

#define Alloca(_N_, _CTYPE_)					\
    (_CTYPE_ *) alloca((size_t) (_N_) * sizeof(_CTYPE_))
    
#define Calloc_or_Alloca_TO(_VAR_, _N_, _CTYPE_)		\
    do {							\
	if (_N_ >= Matrix_Calloc_Threshold) {			\
	    _VAR_ = R_Calloc(_N_, _CTYPE_);			\
	} else {						\
	    _VAR_ = Alloca(_N_, _CTYPE_);			\
	    R_CheckStack();					\
	}							\
    } while (0)
#define Free_FROM(_VAR_, _N_)					\
    do {							\
	if (_N_ >= Matrix_Calloc_Threshold) {			\
	    R_Free(_VAR_);					\
	}							\
    } while (0)

#ifdef ENABLE_NLS
# include <libintl.h>
# define _(String) dgettext("Matrix", String)
#else
# define _(String) (String)
/* <libintl.h> tests N == 1, _not_ N > 1 */
# define dngettext(Domain, String, StringP, N) ((N == 1) ? String : StringP)
#endif

#ifndef LONG_VECTOR_SUPPORT
/* Notably for R <= 2.15.x : */
# define XLENGTH(x) LENGTH(x)
# if R_VERSION < R_Version(2,16,0)
typedef int R_xlen_t;
# endif
#endif

#define Matrix_ErrorBufferSize 4096
#define Matrix_SupportingCachedMethods

/* Previously from <Rdefines.h> : */
#ifndef GET_SLOT
# define GET_SLOT(x, what)        R_do_slot(x, what)
# define SET_SLOT(x, what, value) R_do_slot_assign(x, what, value)
# define MAKE_CLASS(what)	  R_do_MAKE_CLASS(what)
# define NEW_OBJECT(class_def)	  R_do_new_object(class_def)
#endif
/* A safe NEW_OBJECT(MAKE_CLASS(what)) : */
SEXP NEW_OBJECT_OF_CLASS(const char* what);
    
/* Often used symbols (typically for slots), initialized in ./init.c */
extern
#include "Syms.h"

/* Often used numbers, initialized in ./init.c */
extern Rcomplex Matrix_zzero, Matrix_zone;

/* Duplicate the slot with name given by 'sym' from 'src' to 'dest' */
#define slot_dup(dest, src, sym)			\
    SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))
#define slot_dup_if_has(dest, src, sym)				\
    if (R_has_slot(src, sym))					\
	SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))
#define slot_dup_if_not_null(dest, src, sym)			\
    if (!isNull(GET_SLOT(src, sym)))				\
	SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))

#define class_P(x) CHAR(asChar(getAttrib(x, R_ClassSymbol)))
#define uplo_P(x)  CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0))
#define Uplo_P(x)  (R_has_slot(x, Matrix_uploSym) ? uplo_P(x) : " ")
#define diag_P(x)  CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0))
#define Diag_P(x)  (R_has_slot(x, Matrix_diagSym) ? diag_P(x) : " ")

#define MAXOF(x, y) ((x < y) ? y : x)
#define MINOF(x, y) ((x < y) ? x : y)

#define ISNA_LOGICAL(_X_) ((_X_) == NA_LOGICAL)
#define ISNA_INTEGER(_X_) ((_X_) == NA_INTEGER)
#define ISNA_REAL(_X_)    ISNAN(_X_)
#define ISNA_COMPLEX(_X_) (ISNAN((_X_).r) || ISNAN((_X_).i))

#define ISNZ_LOGICAL(_X_) ((_X_) != 0)
#define ISNZ_INTEGER(_X_) ((_X_) != 0)
#define ISNZ_REAL(_X_)    ((_X_) != 0.0)
#define ISNZ_COMPLEX(_X_) ((_X_).r != 0.0 || (_X_).i != 0.0)

#define STRICTLY_ISNZ_LOGICAL(_X_) (!ISNA_LOGICAL(_X_) && ISNZ_LOGICAL(_X_))
#define STRICTLY_ISNZ_INTEGER(_X_) (!ISNA_INTEGER(_X_) && ISNZ_INTEGER(_X_))
#define STRICTLY_ISNZ_REAL(_X_)    (!ISNA_REAL(_X_)    && ISNZ_REAL(_X_))
#define STRICTLY_ISNZ_COMPLEX(_X_) (!ISNA_COMPLEX(_X_) && ISNZ_COMPLEX(_X_))

/* int i, j, m, n; R_xlen_t n2; */
#define PM_AR21_UP(i, j) ((i) + (j) + ((R_xlen_t) (j) * ((j) - 1)) / 2)
#define PM_AR21_LO(i, j, n2) ((i) + ((j) * ((n2) - (j) - 1)) / 2)
#define PM_LENGTH(n) (n + ((R_xlen_t) (n) * ((n) - 1)) / 2)

#define TRIVIAL_DIMNAMES(_DIMNAMES_)			\
    (isNull(VECTOR_ELT(_DIMNAMES_, 0)) &&		\
     isNull(VECTOR_ELT(_DIMNAMES_, 1)) &&		\
     isNull(getAttrib(_DIMNAMES_, R_NamesSymbol)))

#define ERROR_INVALID_CLASS(_CLASS_, _METHOD_)				\
    error(_("invalid class \"%s\" to '%s()'"),				\
	  _CLASS_, _METHOD_)

#define ERROR_INVALID_TYPE(_WHAT_, _SEXPTYPE_, _METHOD_)		\
    error(_("%s of invalid type \"%s\" in '%s()'"),			\
	  _WHAT_, type2char(_SEXPTYPE_), _METHOD_)

/* Zero an array, but note Memzero() which might be FASTER 
   and uses R_SIZE_T (== size_t for C) */
#define AZERO4(_X_, _N_, _ZERO_, _ITYPE_)				\
    do {								\
	for (_ITYPE_ _I_ = 0, _LEN_ = (_N_); _I_ < _LEN_; ++_I_) {	\
	    (_X_)[_I_] = _ZERO_;					\
	}								\
    } while (0)
#define  AZERO(_X_, _N_, _ZERO_) AZERO4(_X_, _N_, 0.0, R_xlen_t)
#define AZEROs(_X_, _N_, _ZERO_) AZERO4(_X_, _N_, 0.0,   size_t)

/* enum constants from cblas.h and some short forms */
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};
#define RMJ CblasRowMajor
#define CMJ CblasColMajor
#define NTR CblasNoTrans
#define TRN CblasTrans
#define CTR CblasConjTrans
#define UPP CblasUpper
#define LOW CblasLower
#define NUN CblasNonUnit
#define UNT CblasUnit
#define LFT CblasLeft
#define RGT CblasRight

enum dense_enum { ddense, ldense, ndense };

/* Define this to be "Cholmod compatible" to some degree */
enum x_slot_kind {
    x_unknown = -2, /* NA */
    x_pattern = -1, /* n */
    x_double  = 0,  /* d */
    x_logical = 1,  /* l */
    x_integer = 2,  /* i */
    x_complex = 3}; /* z */
/* FIXME: use 'x_slot_kind' instead of 'int' 
   everywhere 'Real_(KIND2?|kind_?)' is used */

#define Real_kind_(_x_)							\
    (isReal(_x_) ? x_double : (isLogical(_x_) ? x_logical : x_pattern))

/* Requires 'x' slot, i.e., not for "..nMatrix"
   FIXME? via R_has_slot(obj, name) */
#define Real_kind(_x_)				\
    (Real_kind_(GET_SLOT(_x_, Matrix_xSym)))
    
/* Should also work for "matrix" matrices : */
#define Real_KIND(_x_)						\
    (IS_S4_OBJECT(_x_) ? Real_kind(_x_) : Real_kind_(_x_))
    
/* This one gives 'x_double' also for integer "matrix" :*/
#define Real_KIND2(_x_)						\
    (IS_S4_OBJECT(_x_) ? Real_kind(_x_) :			\
     (isLogical(_x_) ? x_logical : x_double))
    
#define DECLARE_AND_GET_X_SLOT(__C_TYPE, __SEXP)	\
    __C_TYPE *xx = __SEXP(GET_SLOT(x, Matrix_xSym))

SEXP Dim_validate(SEXP dim, const char* domain);
SEXP R_Dim_validate(SEXP dim);

#ifdef Matrix_SupportingCachedMethods
SEXP R_Dim_validate_old(SEXP obj, SEXP domain);
#endif

SEXP DimNames_validate(SEXP dimnames, int pdim[]);
SEXP R_DimNames_validate(SEXP dimnames, SEXP dim);
    
#ifdef Matrix_SupportingCachedMethods
SEXP R_DimNames_validate_old(SEXP obj);
#endif

SEXP string_scalar_validate(SEXP s, char *valid, char *nm);

SEXP Matrix_validate(SEXP obj);
SEXP compMatrix_validate(SEXP obj);
SEXP symmetricMatrix_validate(SEXP obj);
SEXP triangularMatrix_validate(SEXP obj);
SEXP diagonalMatrix_validate(SEXP obj);
SEXP unpackedMatrix_validate(SEXP obj);
SEXP packedMatrix_validate(SEXP obj);
SEXP dMatrix_validate(SEXP obj);
SEXP lMatrix_validate(SEXP obj);
SEXP ndenseMatrix_validate(SEXP obj);
SEXP iMatrix_validate(SEXP obj);
SEXP zMatrix_validate(SEXP obj);
SEXP indMatrix_validate(SEXP obj);
SEXP pMatrix_validate(SEXP obj);

SEXP R_DimNames_fixup(SEXP dn);

Rboolean DimNames_is_symmetric(SEXP dn);
SEXP R_DimNames_is_symmetric(SEXP dn);
    
void symmDN(SEXP dest, SEXP src, int J);
SEXP R_symmDN(SEXP dn);
SEXP get_symmetrized_DimNames(SEXP obj, int J);
void set_symmetrized_DimNames(SEXP obj, SEXP dn, int J);

void revDN(SEXP dest, SEXP src);
SEXP R_revDN(SEXP dn);
SEXP get_reversed_DimNames(SEXP obj);
void set_reversed_DimNames(SEXP obj, SEXP dn);

void set_DimNames(SEXP obj, SEXP dn);

SEXP get_factor(SEXP obj, char *nm);
void set_factor(SEXP obj, char *nm, SEXP val);
SEXP R_set_factor(SEXP obj, SEXP val, SEXP nm, SEXP warn);
SEXP R_empty_factors(SEXP obj, SEXP warn);

#define PACK(_PREFIX_, _CTYPE_)						\
void _PREFIX_ ## dense_pack(_CTYPE_ *dest, const _CTYPE_ *src, int n,	\
			    char uplo, char diag)
PACK(d, double);
PACK(i, int);
PACK(z, Rcomplex);
#undef PACK

#define UNPACK(_PREFIX_, _CTYPE_)					\
void _PREFIX_ ## dense_unpack(_CTYPE_ *dest, const _CTYPE_ *src, int n, \
			      char uplo, char diag)
UNPACK(d, double);
UNPACK(i, int);
UNPACK(z, Rcomplex);
#undef UNPACK

#define UNPACKED_MAKE_SYMMETRIC(_PREFIX_, _CTYPE_)			\
void _PREFIX_ ## dense_unpacked_make_symmetric(_CTYPE_ *x, int n, char uplo)
UNPACKED_MAKE_SYMMETRIC(d, double);
UNPACKED_MAKE_SYMMETRIC(i, int);
UNPACKED_MAKE_SYMMETRIC(z, Rcomplex);
#undef UNPACKED_MAKE_SYMMETRIC

#define UNPACKED_MAKE_TRIANGULAR(_PREFIX_, _CTYPE_)			\
void _PREFIX_ ## dense_unpacked_make_triangular(_CTYPE_ *x, int m, int n, \
						char uplo, char diag)
UNPACKED_MAKE_TRIANGULAR(d, double);
UNPACKED_MAKE_TRIANGULAR(i, int);
UNPACKED_MAKE_TRIANGULAR(z, Rcomplex);
#undef UNPACKED_MAKE_TRIANGULAR

#define UNPACKED_MAKE_BANDED(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_unpacked_make_banded(_CTYPE_ *x,			\
					    int m, int n, int a, int b,	\
					    char diag)
UNPACKED_MAKE_BANDED(d, double);
UNPACKED_MAKE_BANDED(i, int);
UNPACKED_MAKE_BANDED(z, Rcomplex);
#undef UNPACKED_MAKE_BANDED

#define PACKED_MAKE_BANDED(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_packed_make_banded(_CTYPE_ *x,			\
					  int n, int a, int b,		\
					  char uplo, char diag)
PACKED_MAKE_BANDED(d, double);
PACKED_MAKE_BANDED(i, int);
PACKED_MAKE_BANDED(z, Rcomplex);
#undef PACKED_MAKE_BANDED

#define UNPACKED_COPY_DIAGONAL(_PREFIX_, _CTYPE_)			\
void _PREFIX_ ## dense_unpacked_copy_diagonal(_CTYPE_ *dest,	        \
					      const _CTYPE_ *src,	\
					      int n, R_xlen_t len,	\
					      char uplo, char diag)
UNPACKED_COPY_DIAGONAL(d, double);
UNPACKED_COPY_DIAGONAL(i, int);
UNPACKED_COPY_DIAGONAL(z, Rcomplex);
#undef UNPACKED_COPY_DIAGONAL

#define PACKED_COPY_DIAGONAL(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_packed_copy_diagonal(_CTYPE_ *dest,		\
					    const _CTYPE_ *src,		\
					    int n, R_xlen_t len,	\
					    char uplo_dest,		\
					    char uplo_src,		\
					    char diag)
PACKED_COPY_DIAGONAL(d, double);
PACKED_COPY_DIAGONAL(i, int);
PACKED_COPY_DIAGONAL(z, Rcomplex);
#undef PACKED_COPY_DIAGONAL

#define UNPACKED_IS_SYMMETRIC(_PREFIX_, _CTYPE_)			\
Rboolean _PREFIX_ ## dense_unpacked_is_symmetric(const _CTYPE_ *x, int n)
UNPACKED_IS_SYMMETRIC(d, double);
UNPACKED_IS_SYMMETRIC(l, int);
UNPACKED_IS_SYMMETRIC(n, int);
UNPACKED_IS_SYMMETRIC(i, int);
UNPACKED_IS_SYMMETRIC(z, Rcomplex);
#undef UNPACKED_IS_SYMMETRIC
    
#define UNPACKED_IS_TRIANGULAR(_PREFIX_, _CTYPE_)			\
Rboolean _PREFIX_ ## dense_unpacked_is_triangular(const _CTYPE_ *x, int n, \
						  char uplo)
UNPACKED_IS_TRIANGULAR(d, double);
UNPACKED_IS_TRIANGULAR(i, int);
UNPACKED_IS_TRIANGULAR(z, Rcomplex);
#undef UNPACKED_IS_TRIANGULAR

#define UNPACKED_IS_DIAGONAL(_PREFIX_, _CTYPE_)				\
Rboolean _PREFIX_ ## dense_unpacked_is_diagonal(const _CTYPE_ *x, int n)
UNPACKED_IS_DIAGONAL(d, double);
UNPACKED_IS_DIAGONAL(i, int);
UNPACKED_IS_DIAGONAL(z, Rcomplex);
#undef UNPACKED_IS_DIAGONAL
    
#define PACKED_IS_DIAGONAL(_PREFIX_, _CTYPE_)			       \
Rboolean _PREFIX_ ## dense_packed_is_diagonal(const _CTYPE_ *x, int n, \
					      char uplo)
PACKED_IS_DIAGONAL(d, double);
PACKED_IS_DIAGONAL(i, int);
PACKED_IS_DIAGONAL(z, Rcomplex);
#undef PACKED_IS_DIAGONAL

#define PACKED_TRANSPOSE(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_packed_transpose(_CTYPE_ *dest, const _CTYPE_ *src, \
					int n, char uplo)
PACKED_TRANSPOSE(d, double);
PACKED_TRANSPOSE(i, int);
PACKED_TRANSPOSE(z, Rcomplex);
#undef PACKED_TRANSPOSE

SEXP packed_transpose(SEXP x, int n, char uplo);
SEXP unpacked_force(SEXP x, int n, char uplo, char diag);
    
char type2kind(SEXPTYPE type);
SEXPTYPE kind2type(char kind);
size_t kind2size(char kind);

char Matrix_kind(SEXP obj, int i2d);
SEXP R_Matrix_kind(SEXP obj, SEXP i2d);
char Matrix_shape(SEXP obj);
SEXP R_Matrix_shape(SEXP obj);
char Matrix_repr(SEXP obj);
SEXP R_Matrix_repr(SEXP obj);

SEXP matrix_as_dense(SEXP from, const char *code, char uplo, char diag,
		     int new, int transpose_if_vector);
SEXP R_matrix_as_dense(SEXP from, SEXP code, SEXP uplo, SEXP diag);

SEXP dense_as_general(SEXP from, char kind, int new, int transpose_if_vector);
SEXP R_dense_as_general(SEXP from, SEXP kind);

SEXP R_index_triangle(SEXP n_, SEXP upper_, SEXP diag_, SEXP packed_);
SEXP R_index_diagonal(SEXP n_, SEXP upper_,             SEXP packed_);

SEXP R_nnz(SEXP x, SEXP countNA, SEXP nnzmax);

void conjugate(SEXP x);
void zeroRe(SEXP x);
void zeroIm(SEXP x);
void na2one(SEXP x);
    
Rboolean equal_string_vectors(SEXP s1, SEXP s2, int n);
R_xlen_t strmatch(char *nm, SEXP s);
SEXP append_to_named_list(SEXP x, char *nm, SEXP val);

char La_norm_type(const char *typstr);
char La_rcond_type(const char *typstr);
SEXP as_det_obj(double mod, int log, int sign);

/* MJ: no longer needed ... prefer more general (un)?packedMatrix_diag_[gs]et() */
#if 0
void d_packed_getDiag(double *dest, SEXP x, int n);
void l_packed_getDiag(   int *dest, SEXP x, int n);
SEXP d_packed_setDiag(double *diag, int l_d, SEXP x, int n);
SEXP l_packed_setDiag(   int *diag, int l_d, SEXP x, int n);
void tr_d_packed_getDiag(double *dest, SEXP x, int n);
void tr_l_packed_getDiag(   int *dest, SEXP x, int n);
SEXP tr_d_packed_setDiag(double *diag, int l_d, SEXP x, int n);
SEXP tr_l_packed_setDiag(   int *diag, int l_d, SEXP x, int n);
/* were unused, not replaced: */
SEXP d_packed_addDiag(double *diag, int l_d, SEXP x, int n); 
SEXP tr_d_packed_addDiag(double *diag, int l_d, SEXP x, int n);
#endif /* MJ */

#if 0 /* unused */
double get_double_by_name(SEXP obj, char *nm);
SEXP set_double_by_name(SEXP obj, double val, char *nm);
SEXP dgCMatrix_set_Dim(SEXP x, int nrow);
SEXP new_dgeMatrix(int nrow, int ncol);
#endif /* unused */

SEXP Matrix_expand_pointers(SEXP pP);
SEXP m_encodeInd (SEXP ij,        SEXP di, SEXP orig_1, SEXP chk_bnds);
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds);
SEXP Mmatrix(SEXP args);

SEXP R_rbind2_vector(SEXP a, SEXP b);
SEXP R_all0(SEXP x);
SEXP R_any0(SEXP x);

/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.  The validity of this function
 * depends on SET_SLOT not duplicating val when NAMED(val) == 0.  If
 * this behavior changes then ALLOC_SLOT must use SET_SLOT followed by
 * GET_SLOT to ensure that the value returned is indeed the SEXP in
 * the slot.
 * NOTE:  GET_SLOT(x, what)        :== R_do_slot       (x, what)
 * ----   SET_SLOT(x, what, value) :== R_do_slot_assign(x, what, value)
 * and the R_do_slot* are in src/main/attrib.c
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, R_xlen_t length)
{
    SEXP val = allocVector(type, length);

    SET_SLOT(obj, nm, val);
    return val;
}

/**
 * Expand compressed pointers in the array mp into a full set of indices
 * in the array mj.
 *
 * @param ncol number of columns (or rows)
 * @param mp column pointer vector of length ncol + 1
 * @param mj vector of length mp[ncol] to hold the result
 *
 * @return mj
 */
static R_INLINE
int* expand_cmprPt(int ncol, const int mp[], int mj[])
{
    int j;
    for (j = 0; j < ncol; j++) {
	int j2 = mp[j+1], jj;
	for (jj = mp[j]; jj < j2; jj++) mj[jj] = j;
    }
    return mj;
}

/**
 * Check if slot(obj, "x") contains any NA (or NaN).
 *
 * @param obj   a 'Matrix' object with a (double precision) 'x' slot.
 *
 * @return Rboolean :== any(is.na(slot(obj, "x") )
 */
static R_INLINE
Rboolean any_NA_in_x(SEXP obj)
{
    double *x = REAL(GET_SLOT(obj, Matrix_xSym));
    int i, n = LENGTH(GET_SLOT(obj, Matrix_xSym));
    for(i=0; i < n; i++)
	if(ISNAN(x[i])) return TRUE;
    /* else */
    return FALSE;
}

/** Inverse Permutation
 * C version of   .inv.perm.R <- function(p) { p[p] <- seq_along(p) ; p }
 */
static R_INLINE
SEXP inv_permutation(SEXP p_, SEXP zero_p, SEXP zero_res)
{
    int np = 1;
    if(!isInteger(p_)) {p_ = PROTECT(coerceVector(p_, INTSXP)); np++; }
    int *p = INTEGER(p_), n = LENGTH(p_);
    SEXP val = PROTECT(allocVector(INTSXP, n));
    int *v = INTEGER(val), p_0 = asLogical(zero_p), r_0 = asLogical(zero_res);
    if(!p_0) v--; // ==> use 1-based indices
    // shorter (but not 100% sure if ok: is LHS always eval'ed *before* RHS ?) :
    // for(int i=0; i < n; ) v[p[i]] = ++i;
    for(int i=0; i < n; ) {
	int j = p[i]; v[j] = (r_0) ? i++ : ++i;
    }
    UNPROTECT(np);
    return val;
}

/**
 * Return the 0-based index of a string match in a vector of strings
 * terminated by an empty string.  Returns -1 for no match.
 * Is  __cheap__ :  __not__ looking at superclasses --> better use  R_check_class_etc(obj, *)
 *
 * @param x string to match
 * @param valid vector of possible matches terminated by an empty string
 *
 * @return index of match or -1 for no match
 */
static R_INLINE
int Matrix_check_class_(char *x, const char **valid)
{
    int ans = 0;
    while (strlen(valid[ans]) > 0)
	if (strcmp(x, valid[ans]) == 0)
	    return ans;
	else
	    ++ans;
    return -1;
}

static R_INLINE
int Matrix_check_class(SEXP x, const char **valid)
{
    return Matrix_check_class_((char *) class_P(x), valid);
}

/**
 * These are the ones "everyone" should use -- is() versions, also looking
 * at super classes:

 * They now use R(semi_API) from  Rinternals.h :
 * int R_check_class_and_super(SEXP x, const char **valid, SEXP rho);
 * int R_check_class_etc      (SEXP x, const char **valid);

 * R_check_class_etc      (x, v)      basically does  rho <- .classEnv(x)  and then calls
 * R_check_class_and_super(x, v, rho)
 */
// No longer:
#ifdef DEPRECATED_Matrix_check_class_
# define Matrix_check_class_etc R_check_class_etc
# define Matrix_check_class_and_super R_check_class_and_super
#endif

// Keep centralized --- *and* in sync with ../inst/include/Matrix.h :
#define MATRIX_VALID_ge_dense			\
        "dmatrix", "dgeMatrix",			\
	"lmatrix", "lgeMatrix",			\
	"nmatrix", "ngeMatrix",			\
	"zmatrix", "zgeMatrix"

/* NB:  ddiMatrix & ldiMatrix are part of VALID_ddense / VALID_ldense
 * --   even though they are no longer "denseMatrix" formally.
*/
#define MATRIX_VALID_ddense					\
		    "dgeMatrix", "dtrMatrix",			\
		    "dsyMatrix", "dpoMatrix", "ddiMatrix",	\
		    "dtpMatrix", "dspMatrix", "dppMatrix",	\
		    /* sub classes of those above:*/		\
		    /* dtr */ "Cholesky", "LDL", "BunchKaufman",\
		    /* dtp */ "pCholesky", "pBunchKaufman",	\
		    /* dpo */ "corMatrix"

#define MATRIX_VALID_ldense			\
		    "lgeMatrix", "ltrMatrix",	\
		    "lsyMatrix", "ldiMatrix",	\
		    "ltpMatrix", "lspMatrix"

#define MATRIX_VALID_ndense			\
		    "ngeMatrix", "ntrMatrix",	\
		    "nsyMatrix",		\
		    "ntpMatrix", "nspMatrix"

#define MATRIX_VALID_dCsparse			\
 "dgCMatrix", "dsCMatrix", "dtCMatrix"
#define MATRIX_VALID_nCsparse			\
 "ngCMatrix", "nsCMatrix", "ntCMatrix"

#define MATRIX_VALID_Csparse			\
    MATRIX_VALID_dCsparse,			\
 "lgCMatrix", "lsCMatrix", "ltCMatrix",		\
    MATRIX_VALID_nCsparse,			\
 "zgCMatrix", "zsCMatrix", "ztCMatrix"

#define MATRIX_VALID_Tsparse			\
 "dgTMatrix", "dsTMatrix", "dtTMatrix",		\
 "lgTMatrix", "lsTMatrix", "ltTMatrix",		\
 "ngTMatrix", "nsTMatrix", "ntTMatrix",		\
 "zgTMatrix", "zsTMatrix", "ztTMatrix"

#define MATRIX_VALID_Rsparse			\
 "dgRMatrix", "dsRMatrix", "dtRMatrix",		\
 "lgRMatrix", "lsRMatrix", "ltRMatrix",		\
 "ngRMatrix", "nsRMatrix", "ntRMatrix",		\
 "zgRMatrix", "zsRMatrix", "ztRMatrix"

#define MATRIX_VALID_tri_Csparse		\
   "dtCMatrix", "ltCMatrix", "ntCMatrix", "ztCMatrix"

#define MATRIX_VALID_sym_Csparse		\
   "dsCMatrix", "lsCMatrix", "nsCMatrix", "zsCMatrix"

#ifdef __UN_USED__
#define MATRIX_VALID_tri_sparse			\
 "dtCMatrix",  "dtTMatrix", "dtRMatrix",	\
 "ltCMatrix",  "ltTMatrix", "ltRMatrix",	\
 "ntCMatrix",  "ntTMatrix", "ntRMatrix",	\
 "ztCMatrix",  "ztTMatrix", "ztRMatrix"

#define MATRIX_VALID_tri_dense			\
 "dtrMatrix",  "dtpMatrix"			\
 "ltrMatrix",  "ltpMatrix"			\
 "ntrMatrix",  "ntpMatrix"			\
 "ztrMatrix",  "ztpMatrix"
#endif

#define MATRIX_VALID_CHMfactor "dCHMsuper", "dCHMsimpl", "nCHMsuper", "nCHMsimpl"

/** Accessing  *sparseVectors :  fast (and recycling)  v[i] for v = ?sparseVector:
 * -> ./sparseVector.c  -> ./t_sparseVector.c :
 */
// Type_ans sparseVector_sub(int64_t i, int nnz_v, int* v_i, Type_ans* v_x, int len_v):

/* Define all of
 *  dsparseVector_sub(....)
 *  isparseVector_sub(....)
 *  lsparseVector_sub(....)
 *  nsparseVector_sub(....)
 *  zsparseVector_sub(....)
 */
#define _dspV_
#include "t_sparseVector.c"

#define _ispV_
#include "t_sparseVector.c"

#define _lspV_
#include "t_sparseVector.c"

#define _nspV_
#include "t_sparseVector.c"

#define _zspV_
#include "t_sparseVector.c"


#ifdef __cplusplus
}
#endif

#endif /* MATRIX_MUTILS_H_ */
