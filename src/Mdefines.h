#ifndef MATRIX_MDEFINES_H
#define MATRIX_MDEFINES_H

#include "version.h"

#define Matrix_Domain "Matrix"
#define Matrix_CallocThreshold 8192
#define Matrix_ErrorBufferSize 4096

/* NB: system headers should come before R headers */

#ifdef __GLIBC__
/* ensure that strdup() and others are declared when string.h is included : */
# define _POSIX_C_SOURCE 200809L
#endif

#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <float.h>

#ifdef INT_FAST64_MAX
typedef int_fast64_t Matrix_int_fast64_t;
# define MATRIX_INT_FAST64_MIN INT_FAST64_MIN
# define MATRIX_INT_FAST64_MAX INT_FAST64_MAX
#else
typedef    long long Matrix_int_fast64_t;
# define MATRIX_INT_FAST64_MIN      LLONG_MIN
# define MATRIX_INT_FAST64_MAX      LLONG_MAX
#endif

#ifndef STRICT_R_HEADERS
# define STRICT_R_HEADERS
#endif

#include <R.h>
#include <Rinternals.h>

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
#endif

/* Often used symbols, defined in ./init.c */
extern
#include "Msymbols.h"

/* Often used numbers, defined in ./init.c */
extern
Rcomplex Matrix_zzero, Matrix_zone, Matrix_zna; /* 0+0i, 1+0i, NA+NAi */

#define    MINOF(x, y) ((x < y) ? x : y)
#define    MAXOF(x, y) ((x < y) ? y : x)
#define  FIRSTOF(x, y) (x)
#define SECONDOF(x, y) (y)

#define ISNA_PATTERN(_X_) (0)
#define ISNA_LOGICAL(_X_) ((_X_) == NA_LOGICAL)
#define ISNA_INTEGER(_X_) ((_X_) == NA_INTEGER)
#define ISNA_REAL(_X_)    (ISNAN(_X_))
#define ISNA_COMPLEX(_X_) (ISNAN((_X_).r) || ISNAN((_X_).i))

#define ISNZ_PATTERN(_X_) ((_X_) != 0)
#define ISNZ_LOGICAL(_X_) ((_X_) != 0)
#define ISNZ_INTEGER(_X_) ((_X_) != 0)
#define ISNZ_REAL(_X_)    ((_X_) != 0.0)
#define ISNZ_COMPLEX(_X_) ((_X_).r != 0.0 || (_X_).i != 0.0)

#define STRICTLY_ISNZ_PATTERN(_X_) \
	(                      ISNZ_PATTERN(_X_))
#define STRICTLY_ISNZ_LOGICAL(_X_) \
	(!ISNA_LOGICAL(_X_) && ISNZ_LOGICAL(_X_))
#define STRICTLY_ISNZ_INTEGER(_X_) \
	(!ISNA_INTEGER(_X_) && ISNZ_INTEGER(_X_))
#define STRICTLY_ISNZ_REAL(_X_) \
	(!ISNA_REAL(   _X_) && ISNZ_REAL(   _X_))
#define STRICTLY_ISNZ_COMPLEX(_X_) \
	(!ISNA_COMPLEX(_X_) && ISNZ_COMPLEX(_X_))

#define NOTREAL_PATTERN(_X_) 0
#define NOTREAL_LOGICAL(_X_) 0
#define NOTREAL_INTEGER(_X_) 0
#define NOTREAL_REAL(_X_)    0
#define NOTREAL_COMPLEX(_X_) (_X_.i != 0.0)

#define NOTCONJ_PATTERN(_X_, _Y_) \
	((_X_ != 0) != (_Y_ != 0))
#define NOTCONJ_LOGICAL(_X_, _Y_) \
	(_X_ != _Y_)
#define NOTCONJ_INTEGER(_X_, _Y_) \
	(_X_ != _Y_)
#define NOTCONJ_REAL(_X_, _Y_) \
	((ISNAN(_X_)) ? !ISNAN(_Y_) : ISNAN(_Y_) || _X_ != _Y_)
#define NOTCONJ_COMPLEX(_X_, _Y_) \
	(((ISNAN(_X_.r)) ? !ISNAN(_Y_.r) : ISNAN(_Y_.r) || _X_.r !=  _Y_.r) || \
	 ((ISNAN(_X_.i)) ? !ISNAN(_Y_.i) : ISNAN(_Y_.i) || _X_.r != -_Y_.r))

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

#define ASSIGN_REAL(_X_, _Y_) \
	do { _X_   = _Y_  ;                } while (0)
#define ASSIGN_COMPLEX(_X_, _Y_) \
	do { _X_.r = _Y_.r; _X_.i = _Y_.i; } while (0)

#define SCALE1_REAL(_X_, _A_) \
	do { _X_   *= _A_;               } while (0)
#define SCALE1_COMPLEX(_X_, _A_) \
	do { _X_.r *= _A_; _X_.i *= _A_; } while (0)

#define SCALE2_REAL(_X_, _A_) \
	do { _X_   /= _A_;               } while (0)
#define SCALE2_COMPLEX(_X_, _A_) \
	do { _X_.r /= _A_; _X_.i /= _A_; } while (0)

#define PACKED_AR21_UP(i, j) \
	((R_xlen_t) ((i) + ((Matrix_int_fast64_t) (j) * (       (j) + 1)) / 2))
#define PACKED_AR21_LO(i, j, m2) \
	((R_xlen_t) ((i) + ((Matrix_int_fast64_t) (j) * ((m2) - (j) - 1)) / 2))
#define PACKED_LENGTH(m) \
	((R_xlen_t) ((m) + ((Matrix_int_fast64_t) (m) * (       (m) - 1)) / 2))

#define SHOW(...) __VA_ARGS__
#define HIDE(...)

#define ERROR_INVALID_TYPE(_X_, _FUNC_) \
	error(_("invalid type \"%s\" in '%s'"), \
	      type2char(TYPEOF(_X_)), _FUNC_)

#define ERROR_INVALID_CLASS(_X_, _FUNC_) \
do { \
	if (!OBJECT(_X_)) \
		ERROR_INVALID_TYPE(_X_, _FUNC_); \
	else { \
		SEXP class = PROTECT(getAttrib(_X_, R_ClassSymbol)); \
		error(_("invalid class \"%s\" in '%s'"), \
		      CHAR(STRING_ELT(class, 0)), _FUNC_); \
		UNPROTECT(1); \
	} \
} while (0)

#define VALID_NONVIRTUAL_MATRIX \
/*  0 */ "dpoMatrix", "dppMatrix", \
/*  2 */ "corMatrix", "copMatrix", \
/*  4 */   "pMatrix", "indMatrix", \
/*  6 */ "ngCMatrix", "ngRMatrix", "ngTMatrix", "ngeMatrix", "ndiMatrix", \
/* 11 */ "nsCMatrix", "nsRMatrix", "nsTMatrix", "nsyMatrix", "nspMatrix", \
/* 16 */ "ntCMatrix", "ntRMatrix", "ntTMatrix", "ntrMatrix", "ntpMatrix", \
/* 21 */ "lgCMatrix", "lgRMatrix", "lgTMatrix", "lgeMatrix", "ldiMatrix", \
/* 26 */ "lsCMatrix", "lsRMatrix", "lsTMatrix", "lsyMatrix", "lspMatrix", \
/* 31 */ "ltCMatrix", "ltRMatrix", "ltTMatrix", "ltrMatrix", "ltpMatrix", \
/* 36 */ "igCMatrix", "igRMatrix", "igTMatrix", "igeMatrix", "idiMatrix", \
/* 41 */ "isCMatrix", "isRMatrix", "isTMatrix", "isyMatrix", "ispMatrix", \
/* 46 */ "itCMatrix", "itRMatrix", "itTMatrix", "itrMatrix", "itpMatrix", \
/* 51 */ "dgCMatrix", "dgRMatrix", "dgTMatrix", "dgeMatrix", "ddiMatrix", \
/* 56 */ "dsCMatrix", "dsRMatrix", "dsTMatrix", "dsyMatrix", "dspMatrix", \
/* 61 */ "dtCMatrix", "dtRMatrix", "dtTMatrix", "dtrMatrix", "dtpMatrix", \
/* 66 */ "zgCMatrix", "zgRMatrix", "zgTMatrix", "zgeMatrix", "zdiMatrix", \
/* 71 */ "zsCMatrix", "zsRMatrix", "zsTMatrix", "zsyMatrix", "zspMatrix", \
/* 76 */ "ztCMatrix", "ztRMatrix", "ztTMatrix", "ztrMatrix", "ztpMatrix"

#define VALID_NONVIRTUAL_VECTOR \
/* 81 */ "nsparseVector", "lsparseVector", "isparseVector", \
         "dsparseVector", "zsparseVector"

#define VALID_NONVIRTUAL VALID_NONVIRTUAL_MATRIX, VALID_NONVIRTUAL_VECTOR

/* dpoMatrix->dsyMatrix, etc. */
#define VALID_NONVIRTUAL_SHIFT(i, pToInd) \
	((i >= 5) ? 0 : ((i >= 4) ? pToInd != 0 : ((i >= 2) ? 57 : 59)))

#define VALID_DENSE \
"ngeMatrix", "nsyMatrix", "nspMatrix", "ntrMatrix", "ntpMatrix", \
"lgeMatrix", "lsyMatrix", "lspMatrix", "ltrMatrix", "ltpMatrix", \
"igeMatrix", "isyMatrix", "ispMatrix", "itrMatrix", "itpMatrix", \
"dgeMatrix", "dsyMatrix", "dspMatrix", "dtrMatrix", "dtpMatrix", \
"zgeMatrix", "zsyMatrix", "zspMatrix", "ztrMatrix", "ztpMatrix"

#define VALID_CSPARSE \
"ngCMatrix", "nsCMatrix", "ntCMatrix", \
"lgCMatrix", "lsCMatrix", "ltCMatrix", \
"igCMatrix", "isCMatrix", "itCMatrix", \
"dgCMatrix", "dsCMatrix", "dtCMatrix", \
"zgCMatrix", "zsCMatrix", "ztCMatrix"

#define VALID_RSPARSE \
"ngRMatrix", "nsRMatrix", "ntRMatrix", \
"lgRMatrix", "lsRMatrix", "ltRMatrix", \
"igRMatrix", "isRMatrix", "itRMatrix", \
"dgRMatrix", "dsRMatrix", "dtRMatrix", \
"zgRMatrix", "zsRMatrix", "ztRMatrix"

#define VALID_TSPARSE \
"ngTMatrix", "nsTMatrix", "ntTMatrix", \
"lgTMatrix", "lsTMatrix", "ltTMatrix", \
"igTMatrix", "isTMatrix", "itTMatrix", \
"dgTMatrix", "dsTMatrix", "dtTMatrix", \
"zgTMatrix", "zsTMatrix", "ztTMatrix"

#define VALID_DIAGONAL \
"ndiMatrix", "ldiMatrix", "idiMatrix", "ddiMatrix", "zdiMatrix"


/* What we want declared "everywhere" : */

#include "utils.h"

SEXP newObject(const char *);
void validObject(SEXP, const char *);

char typeToKind(SEXPTYPE);
SEXPTYPE kindToType(char);
size_t kindToSize(char);

int DimNames_is_trivial(SEXP);
int DimNames_is_symmetric(SEXP);

void symDN(SEXP, SEXP, int);
void revDN(SEXP, SEXP);

SEXP get_symmetrized_DimNames(SEXP, int);
SEXP get_reversed_DimNames(SEXP);

void set_symmetrized_DimNames(SEXP, SEXP, int);
void set_reversed_DimNames(SEXP, SEXP);

#endif /* MATRIX_MDEFINES_H */
