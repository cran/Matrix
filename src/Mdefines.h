#ifndef MATRIX_DEFINES_H
#define MATRIX_DEFINES_H

/* R has the same in several places : */
#define Matrix_CallocThreshold 10000
#define Matrix_ErrorBufferSize  4096

#define Matrix_SupportingCachedMethods

#undef Matrix_with_SPQR
#undef HAVE_PROPER_IMATRIX
#undef HAVE_PROPER_ZMATRIX

/* Copy and paste from WRE : */
#ifdef ENABLE_NLS
# include <libintl.h>
# define _(String) dgettext("Matrix", String)
#else
# define _(String) (String)
/* <libintl.h> tests N == 1, _not_ N > 1 */
# define dngettext(Domain, String, StringP, N) ((N == 1) ? String : StringP)
#endif

/* Copy and paste from Defn.h : */
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

#define Alloca(_N_, _CTYPE_)					\
    (_CTYPE_ *) alloca((size_t) (_N_) * sizeof(_CTYPE_))
    
#define Calloc_or_Alloca_TO(_VAR_, _N_, _CTYPE_)		\
    do {							\
	if (_N_ >= Matrix_CallocThreshold)			\
	    _VAR_ = R_Calloc(_N_, _CTYPE_);			\
	else {							\
	    _VAR_ = Alloca(_N_, _CTYPE_);			\
	    R_CheckStack();					\
	    /* Memzero(_VAR_, _N_); */				\
	}							\
    } while (0)

#define Free_FROM(_VAR_, _N_)					\
    do {							\
	if (_N_ >= Matrix_CallocThreshold)			\
	    R_Free(_VAR_);					\
    } while (0)

/* To zero an array ... however, note Memzero(), which calls memset()
   and so can be faster in the range of R_SIZE_T (an alias for size_t in C)
*/
#define AZERO(_X_, _N_, _ZERO_, _CTYPE_)				\
    do {								\
	for (_CTYPE_ _I_ = 0, _LEN_ = (_N_); _I_ < _LEN_; ++_I_)	\
	    (_X_)[_I_] = _ZERO_;					\
    } while (0)

/* Copy and paste from now-deprecated Rdefines.h : */
#ifndef R_DEFINES_H
# define GET_SLOT(x, what)        R_do_slot(x, what)
# define SET_SLOT(x, what, value) R_do_slot_assign(x, what, value)
# define MAKE_CLASS(what)	  R_do_MAKE_CLASS(what)
# define NEW_OBJECT(class_def)	  R_do_new_object(class_def)
#endif
#define HAS_SLOT(obj, name)       R_has_slot(obj, name)

/* Often used symbols, defined in ./init.c */
extern
#include "Syms.h"

/* Often used numbers, defined in ./init.c */
extern
Rcomplex Matrix_zzero, Matrix_zone; /* 0+0i, 1+0i */

/* To become deprecated ... defensive code should PROTECT() more */
#define class_P(x) CHAR(asChar(getAttrib(x, R_ClassSymbol)))
#define  uplo_P(x) CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0))
#define  Uplo_P(x) (R_has_slot(x, Matrix_uploSym) ? uplo_P(x) : " ")
#define  diag_P(x) CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0))
#define  Diag_P(x) (R_has_slot(x, Matrix_diagSym) ? diag_P(x) : " ")

/* Ditto */
#define slot_dup(dest, src, sym)			\
    SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))
#define slot_dup_if_has(dest, src, sym)				\
    if (R_has_slot(src, sym))					\
	SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))
#define slot_dup_if_not_null(dest, src, sym)			\
    if (!isNull(GET_SLOT(src, sym)))				\
	SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))

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

#define PM_AR21_UP(i, j) ((i) + (j) + ((R_xlen_t) (j) * ((j) - 1)) / 2)
#define PM_AR21_LO(i, j, n2) ((i) + ((j) * ((n2) - (j) - 1)) / 2)
#define PM_LENGTH(n) (n + ((R_xlen_t) (n) * ((n) - 1)) / 2)

#define ERROR_INVALID_CLASS(_X_, _METHOD_)			\
    do {							\
	SEXP class = PROTECT(getAttrib(_X_, R_ClassSymbol));	\
	if (TYPEOF(class) == STRSXP && LENGTH(class) > 0)	\
	    error(_("invalid class \"%s\" to '%s()'"),		\
		  CHAR(STRING_ELT(class, 0)), _METHOD_);	\
	else							\
	    error(_("unclassed \"%s\" to '%s()'"),		\
		  type2char(TYPEOF(_X_)), _METHOD_);		\
	UNPROTECT(1);						\
    } while (0)

#define ERROR_INVALID_TYPE(_WHAT_, _SEXPTYPE_, _METHOD_)	\
    error(_("%s of invalid type \"%s\" in '%s()'"),		\
	  _WHAT_, type2char(_SEXPTYPE_), _METHOD_)

/* For C-level isTriangular() : */
#define RETURN_TRUE_OF_KIND(_KIND_)			\
    do {						\
	SEXP ans = PROTECT(allocVector(LGLSXP, 1)),	\
	    val = PROTECT(mkString(_KIND_));		\
	static SEXP sym = NULL;				\
	if (!sym)					\
	    sym = install("kind");			\
	LOGICAL(ans)[0] = 1;				\
	setAttrib(ans, sym, val);			\
	UNPROTECT(2); /* val, ans */			\
	return ans;					\
    } while (0)

/* Define this to be "Cholmod-compatible" to some degree : */
enum x_slot_kind {
    x_unknown = -2, /* NA */
    x_pattern = -1, /* n */
    x_double  = 0,  /* d */
    x_logical = 1,  /* l */
    x_integer = 2,  /* i */
    x_complex = 3}; /* z */
/* FIXME: use 'x_slot_kind' instead of 'int' 
   everywhere that Real_(KIND2?|kind_?) is used ...
*/

#define Real_kind_(_x_)							\
    (isReal(_x_) ? x_double : (isLogical(_x_) ? x_logical : x_pattern))

/* Requires 'x' slot, hence not for nsparseMatrix or indMatrix : */
#define Real_kind(_x_)				\
    (Real_kind_(GET_SLOT(_x_, Matrix_xSym)))


/* ==== NO LONGER USED ============================================== */

#if 0

enum dense_enum { ddense, ldense, ndense };

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

/* This one also works for traditional matrices : */
#define Real_KIND(_x_)						\
    (IS_S4_OBJECT(_x_) ? Real_kind(_x_) : Real_kind_(_x_))
    
/* This one gives 'x_double' also for integer matrices : */
#define Real_KIND2(_x_)						\
    (IS_S4_OBJECT(_x_) ? Real_kind(_x_) :			\
     (isLogical(_x_) ? x_logical : x_double))

#define DECLARE_AND_GET_X_SLOT(__C_TYPE, __SEXP)	\
    __C_TYPE *xx = __SEXP(GET_SLOT(x, Matrix_xSym))

#endif


/* ==== CLASS LISTS ================================================= */
/* Keep synchronized with ../inst/include/Matrix.h !                  */

#define MATRIX_VALID_ge_dense			\
    "dmatrix", "dgeMatrix",			\
    "lmatrix", "lgeMatrix",			\
    "nmatrix", "ngeMatrix",			\
    "zmatrix", "zgeMatrix"

/* NB: includes ddiMatrix which is no longer formally denseMatrix */
#define MATRIX_VALID_ddense				\
    "dgeMatrix", "dtrMatrix",				\
    "dsyMatrix", "dpoMatrix", "ddiMatrix",		\
    "dtpMatrix", "dspMatrix", "dppMatrix",		\
    /* subclasses of the above : */			\
    /* dtr */ "Cholesky", "LDL", "BunchKaufman",	\
    /* dtp */ "pCholesky", "pBunchKaufman",		\
    /* dpo */ "corMatrix"

/* NB: includes ldiMatrix which is no longer formally denseMatrix */
#define MATRIX_VALID_ldense			\
    "lgeMatrix",				\
    "ltrMatrix", "lsyMatrix", "ldiMatrix",	\
    "ltpMatrix", "lspMatrix"

#define MATRIX_VALID_ndense			\
    "ngeMatrix",				\
    "ntrMatrix", "nsyMatrix",			\
    "ntpMatrix", "nspMatrix"

#define MATRIX_VALID_dCsparse			\
    "dgCMatrix", "dsCMatrix", "dtCMatrix"
#define MATRIX_VALID_nCsparse			\
    "ngCMatrix", "nsCMatrix", "ntCMatrix"

#define MATRIX_VALID_Csparse			\
    MATRIX_VALID_dCsparse,			\
    "lgCMatrix", "lsCMatrix", "ltCMatrix",	\
    MATRIX_VALID_nCsparse,			\
    "zgCMatrix", "zsCMatrix", "ztCMatrix"

#define MATRIX_VALID_Tsparse			\
    "dgTMatrix", "dsTMatrix", "dtTMatrix",	\
    "lgTMatrix", "lsTMatrix", "ltTMatrix",	\
    "ngTMatrix", "nsTMatrix", "ntTMatrix",	\
    "zgTMatrix", "zsTMatrix", "ztTMatrix"

#define MATRIX_VALID_Rsparse			\
    "dgRMatrix", "dsRMatrix", "dtRMatrix",	\
    "lgRMatrix", "lsRMatrix", "ltRMatrix",	\
    "ngRMatrix", "nsRMatrix", "ntRMatrix",	\
    "zgRMatrix", "zsRMatrix", "ztRMatrix"

#define MATRIX_VALID_tri_Csparse			\
    "dtCMatrix", "ltCMatrix", "ntCMatrix", "ztCMatrix"

#define MATRIX_VALID_sym_Csparse			\
    "dsCMatrix", "lsCMatrix", "nsCMatrix", "zsCMatrix"

#define MATRIX_VALID_CHMfactor				\
    "dCHMsuper", "dCHMsimpl", "nCHMsuper", "nCHMsimpl"

#endif /* MATRIX_DEFINES_H */
