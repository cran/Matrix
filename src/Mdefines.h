#ifndef MATRIX_MDEFINES_H
#define MATRIX_MDEFINES_H

#define Matrix_Domain          "Matrix"
#define Matrix_CallocThreshold    10000
#define Matrix_ErrorBufferSize     4096

/* Eventually these will no longer be needed : */
#undef Matrix_WithSPQR

/* NB: system headers should come before R headers */

#ifdef __GLIBC__
/* ensure that strdup() and others are declared when string.h is included : */
# define _POSIX_C_SOURCE 200809L
#endif

#include <string.h>
#include <stdint.h>
#include <limits.h>

#ifdef INT_FAST64_MAX
typedef int_fast64_t Matrix_int_fast64_t;
# define MATRIX_INT_FAST64_MAX INT_FAST64_MAX
#else
typedef    long long Matrix_int_fast64_t;
# define MATRIX_INT_FAST64_MAX      LLONG_MAX
#endif

#ifndef STRICT_R_HEADERS
# define STRICT_R_HEADERS
#endif

#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>

/* Copy and paste from WRE : */
#ifdef ENABLE_NLS
# include <libintl.h>
# define _(String) dgettext(Matrix_Domain, String)
#else
# define _(String) (String)
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

#define Matrix_Calloc(_VAR_, _N_, _CTYPE_) \
do { \
	if (_N_ >= Matrix_CallocThreshold) \
		_VAR_ = R_Calloc(_N_, _CTYPE_); \
	else { \
		_VAR_ = (_CTYPE_ *) alloca((size_t) (_N_) * sizeof(_CTYPE_)); \
		R_CheckStack(); \
		memset(_VAR_, 0, (size_t) (_N_) * sizeof(_CTYPE_)); \
	} \
} while (0)

#define Matrix_Free(_VAR_, _N_) \
do { \
	if (_N_ >= Matrix_CallocThreshold) \
		R_Free(_VAR_); \
} while (0)

/* Copy and paste from now-deprecated Rdefines.h : */
#ifndef R_DEFINES_H
# define GET_SLOT(x, what)        R_do_slot(x, what)
# define SET_SLOT(x, what, value) R_do_slot_assign(x, what, value)
# define MAKE_CLASS(what)         R_do_MAKE_CLASS(what)
# define NEW_OBJECT(class_def)    R_do_new_object(class_def)
#endif
#define HAS_SLOT(obj, name)       R_has_slot(obj, name)

/* Often used symbols, defined in ./init.c */
extern
#include "Syms.h"

/* Often used numbers, defined in ./init.c */
extern
Rcomplex Matrix_zzero, Matrix_zone, Matrix_zna; /* 0+0i, 1+0i, NA+NAi */

/* To become deprecated ... defensive code should PROTECT() more */
#define class_P(x) CHAR(asChar(getAttrib(x, R_ClassSymbol)))
#define  uplo_P(x) CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0))
#define  Uplo_P(x) (R_has_slot(x, Matrix_uploSym) ? uplo_P(x) : " ")
#define  diag_P(x) CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0))
#define  Diag_P(x) (R_has_slot(x, Matrix_diagSym) ? diag_P(x) : " ")

/* Ditto */
#define slot_dup(dest, src, sym) \
	SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))
#define slot_dup_if_has(dest, src, sym) \
	if (R_has_slot(src, sym)) \
		SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))
#define slot_dup_if_not_null(dest, src, sym) \
	if (!isNull(GET_SLOT(src, sym))) \
		SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))

#define MAXOF(x, y) ((x < y) ? y : x)
#define MINOF(x, y) ((x < y) ? x : y)
#define  FIRSTOF(x, y) (x)
#define SECONDOF(x, y) (y)

#define ISNA_LOGICAL(_X_) ((_X_) == NA_LOGICAL)
#define ISNA_INTEGER(_X_) ((_X_) == NA_INTEGER)
#define ISNA_REAL(_X_)    (ISNAN(_X_))
#define ISNA_COMPLEX(_X_) (ISNAN((_X_).r) || ISNAN((_X_).i))

#define ISNZ_LOGICAL(_X_) ((_X_) != 0)
#define ISNZ_INTEGER(_X_) ((_X_) != 0)
#define ISNZ_REAL(_X_)    ((_X_) != 0.0)
#define ISNZ_COMPLEX(_X_) ((_X_).r != 0.0 || (_X_).i != 0.0)

#define STRICTLY_ISNZ_LOGICAL(_X_) \
	(!ISNA_LOGICAL(_X_) && ISNZ_LOGICAL(_X_))
#define STRICTLY_ISNZ_INTEGER(_X_) \
	(!ISNA_INTEGER(_X_) && ISNZ_INTEGER(_X_))
#define STRICTLY_ISNZ_REAL(_X_) \
	(!ISNA_REAL(   _X_) && ISNZ_REAL(   _X_))
#define STRICTLY_ISNZ_COMPLEX(_X_) \
	(!ISNA_COMPLEX(_X_) && ISNZ_COMPLEX(_X_))

#define INCREMENT_PATTERN(_X_, _Y_) \
	do { \
		_X_ = 1; \
	} while (0)
#define INCREMENT_LOGICAL(_X_, _Y_) \
	do { \
		if (_Y_ == NA_LOGICAL) { \
			if (_X_ == 0) \
				_X_ = NA_LOGICAL; \
		} else if (_Y_ != 0) \
			_X_ = 1; \
	} while (0)
#define INCREMENT_INTEGER(_X_, _Y_) \
	do { \
		if (_X_ != NA_INTEGER) { \
			if (_Y_ == NA_INTEGER) \
				_X_ = NA_INTEGER; \
			else if ((_Y_ < 0) \
					 ? (_X_ <= INT_MIN - _Y_) \
					 : (_X_ >  INT_MAX - _Y_)) { \
				warning(_("NAs produced by integer overflow")); \
				_X_ = NA_INTEGER; \
			} else \
				_X_ += _Y_; \
		} \
	} while (0)
#define INCREMENT_REAL(_X_, _Y_) \
	do { \
		_X_ += _Y_; \
	} while (0)
#define INCREMENT_COMPLEX(_X_, _Y_) \
	do { \
		_X_.r += _Y_.r; \
		_X_.i += _Y_.i; \
	} while (0)

#define PM_AR21_UP(i, j) \
	((R_xlen_t) ((i) + ((Matrix_int_fast64_t) (j) * (       (j) + 1)) / 2))
#define PM_AR21_LO(i, j, m2) \
	((R_xlen_t) ((i) + ((Matrix_int_fast64_t) (j) * ((m2) - (j) - 1)) / 2))
#define PM_LENGTH(m) \
	((R_xlen_t) ((m) + ((Matrix_int_fast64_t) (m) * (       (m) - 1)) / 2))

#define SHOW(...) __VA_ARGS__
#define HIDE(...)

#define ERROR_INVALID_TYPE(_X_, _FUNC_) \
	error(_("invalid type \"%s\" in %s()"), \
	      type2char(TYPEOF(_X_)), _FUNC_)

#define ERROR_INVALID_CLASS(_X_, _FUNC_) \
do { \
	if (!OBJECT(_X_)) \
		ERROR_INVALID_TYPE(_X_, _FUNC_); \
	else { \
		SEXP class = PROTECT(getAttrib(_X_, R_ClassSymbol)); \
		error(_("invalid class \"%s\" in %s()"), \
		      CHAR(STRING_ELT(class, 0)), _FUNC_); \
		UNPROTECT(1); \
	} \
} while (0)

/* For C-level isTriangular() : */
#define RETURN_TRUE_OF_KIND(_KIND_) \
do { \
	SEXP ans = PROTECT(allocVector(LGLSXP, 1)), \
		val = PROTECT(mkString(_KIND_)); \
	static SEXP sym = NULL; \
	if (!sym) \
		sym = install("kind"); \
	LOGICAL(ans)[0] = 1; \
	setAttrib(ans, sym, val); \
	UNPROTECT(2); /* val, ans */ \
	return ans; \
} while (0)

/* Define this to be CHOLMOD-compatible to some degree : */
enum x_slot_kind {
	x_unknown = -2,  /* NA */
	x_pattern = -1,  /*  n */
	x_double  =  0,  /*  d */
	x_logical =  1,  /*  l */
	x_integer =  2,  /*  i */
	x_complex =  3}; /*  z */

#define Real_kind_(_x_) \
	(isReal(_x_) ? x_double : (isLogical(_x_) ? x_logical : x_pattern))

/* Requires 'x' slot, hence not for nsparseMatrix or indMatrix : */
#define Real_kind(_x_) \
	(Real_kind_(GET_SLOT(_x_, Matrix_xSym)))


/* ==== CLASS LISTS ================================================= */
/* Keep synchronized with ../inst/include/Matrix.h !                  */

/* dpoMatrix->dsyMatrix, etc. */
#define VALID_NONVIRTUAL_SHIFT(i, p2ind) \
	((i >= 5) ? 0 : ((i >= 4) ? p2ind != 0 : ((i >= 2) ? 12 : 14)))

#define VALID_NONVIRTUAL_MATRIX \
/*  0 */ "dpoMatrix", "dppMatrix", \
/*  2 */ "corMatrix", "pcorMatrix", \
/*  4 */ "pMatrix", "indMatrix", \
/*  6 */ "dgCMatrix", "dgRMatrix", "dgTMatrix", "dgeMatrix", "ddiMatrix", \
/* 11 */ "dsCMatrix", "dsRMatrix", "dsTMatrix", "dsyMatrix", "dspMatrix", \
/* 16 */ "dtCMatrix", "dtRMatrix", "dtTMatrix", "dtrMatrix", "dtpMatrix", \
/* 21 */ "lgCMatrix", "lgRMatrix", "lgTMatrix", "lgeMatrix", "ldiMatrix", \
/* 26 */ "lsCMatrix", "lsRMatrix", "lsTMatrix", "lsyMatrix", "lspMatrix", \
/* 31 */ "ltCMatrix", "ltRMatrix", "ltTMatrix", "ltrMatrix", "ltpMatrix", \
/* 36 */ "ngCMatrix", "ngRMatrix", "ngTMatrix", "ngeMatrix", "ndiMatrix", \
/* 41 */ "nsCMatrix", "nsRMatrix", "nsTMatrix", "nsyMatrix", "nspMatrix", \
/* 46 */ "ntCMatrix", "ntRMatrix", "ntTMatrix", "ntrMatrix", "ntpMatrix", \
/* 51 */ "igCMatrix", "igRMatrix", "igTMatrix", "igeMatrix", "idiMatrix", \
/* 56 */ "isCMatrix", "isRMatrix", "isTMatrix", "isyMatrix", "ispMatrix", \
/* 61 */ "itCMatrix", "itRMatrix", "itTMatrix", "itrMatrix", "itpMatrix", \
/* 66 */ "zgCMatrix", "zgRMatrix", "zgTMatrix", "zgeMatrix", "zdiMatrix", \
/* 71 */ "zsCMatrix", "zsRMatrix", "zsTMatrix", "zsyMatrix", "zspMatrix", \
/* 76 */ "ztCMatrix", "ztRMatrix", "ztTMatrix", "ztrMatrix", "ztpMatrix"

#define VALID_NONVIRTUAL_VECTOR \
/* 81 */ "dsparseVector", "lsparseVector", "nsparseVector", \
         "isparseVector", "zsparseVector"

#define VALID_NONVIRTUAL VALID_NONVIRTUAL_MATRIX, VALID_NONVIRTUAL_VECTOR

#define VALID_NDENSE \
"ngeMatrix", "ntrMatrix", "nsyMatrix", "ntpMatrix", "nspMatrix"

#define VALID_LDENSE \
"lgeMatrix", "ltrMatrix", "lsyMatrix", "ltpMatrix", "lspMatrix"

#define VALID_DDENSE \
"dgeMatrix", "dtrMatrix", "dsyMatrix", "dtpMatrix", "dspMatrix"

#define VALID_NSPARSE \
"ngCMatrix", "ngRMatrix", "ngTMatrix", \
"ntCMatrix", "ntRMatrix", "ntTMatrix", \
"nsCMatrix", "nsRMatrix", "nsTMatrix"

#define VALID_LSPARSE \
"lgCMatrix", "lgRMatrix", "lgTMatrix", \
"ltCMatrix", "ltRMatrix", "ltTMatrix", \
"lsCMatrix", "lsRMatrix", "lsTMatrix"

#define VALID_DSPARSE \
"dgCMatrix", "dgRMatrix", "dgTMatrix", \
"dtCMatrix", "dtRMatrix", "dtTMatrix", \
"dsCMatrix", "dsRMatrix", "dsTMatrix"

#define VALID_CSPARSE \
"dgCMatrix", "dtCMatrix", "dsCMatrix", \
"lgCMatrix", "ltCMatrix", "lsCMatrix", \
"ngCMatrix", "ntCMatrix", "nsCMatrix"

#define VALID_RSPARSE \
"dgRMatrix", "dtRMatrix", "dsRMatrix", \
"lgRMatrix", "ltRMatrix", "lsRMatrix", \
"ngRMatrix", "ntRMatrix", "nsRMatrix"

#define VALID_TSPARSE \
"dgTMatrix", "dtTMatrix", "dsTMatrix", \
"lgTMatrix", "ltTMatrix", "lsTMatrix", \
"ngTMatrix", "ntTMatrix", "nsTMatrix"

#define VALID_DIAGONAL \
"ddiMatrix", "ldiMatrix"

/* Older ones : */

#define MATRIX_VALID_ge_dense \
"dmatrix", "dgeMatrix", \
"lmatrix", "lgeMatrix", \
"nmatrix", "ngeMatrix", \
"zmatrix", "zgeMatrix"

#define MATRIX_VALID_ddense \
"dgeMatrix", "dtrMatrix", \
"dsyMatrix", "dpoMatrix", "ddiMatrix", \
"dtpMatrix", "dspMatrix", "dppMatrix", \
/* subclasses of the above : */ \
/* dtr */ "Cholesky", "LDL", "BunchKaufman", \
/* dtp */ "pCholesky", "pBunchKaufman", \
/* dpo */ "corMatrix"

#define MATRIX_VALID_ldense \
"lgeMatrix", \
"ltrMatrix", "lsyMatrix", "ldiMatrix", \
"ltpMatrix", "lspMatrix"

#define MATRIX_VALID_ndense \
"ngeMatrix", \
"ntrMatrix", "nsyMatrix", \
"ntpMatrix", "nspMatrix"

#define MATRIX_VALID_dCsparse \
"dgCMatrix", "dsCMatrix", "dtCMatrix"
#define MATRIX_VALID_nCsparse \
"ngCMatrix", "nsCMatrix", "ntCMatrix"

#define MATRIX_VALID_Csparse \
MATRIX_VALID_dCsparse, \
"lgCMatrix", "lsCMatrix", "ltCMatrix", \
MATRIX_VALID_nCsparse, \
"zgCMatrix", "zsCMatrix", "ztCMatrix"

#define MATRIX_VALID_Tsparse \
"dgTMatrix", "dsTMatrix", "dtTMatrix", \
"lgTMatrix", "lsTMatrix", "ltTMatrix", \
"ngTMatrix", "nsTMatrix", "ntTMatrix", \
"zgTMatrix", "zsTMatrix", "ztTMatrix"

#define MATRIX_VALID_Rsparse \
"dgRMatrix", "dsRMatrix", "dtRMatrix", \
"lgRMatrix", "lsRMatrix", "ltRMatrix", \
"ngRMatrix", "nsRMatrix", "ntRMatrix", \
"zgRMatrix", "zsRMatrix", "ztRMatrix"

#define MATRIX_VALID_tri_Csparse \
"dtCMatrix", "ltCMatrix", "ntCMatrix", "ztCMatrix"

#define MATRIX_VALID_sym_Csparse \
"dsCMatrix", "lsCMatrix", "nsCMatrix", "zsCMatrix"

#define MATRIX_VALID_CHMfactor \
"dCHMsuper", "dCHMsimpl", "nCHMsuper", "nCHMsimpl"

#endif /* MATRIX_MDEFINES_H */
