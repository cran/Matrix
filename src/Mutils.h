#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#ifndef STRICT_R_HEADERS
# define STRICT_R_HEADERS
#endif

#ifdef __GLIBC__
/* to ensure that strdup() and others are declared
   when string.h is included with R.h (WRE) :
*/
# define _POSIX_C_SOURCE 200809L
#endif

/* NB: system headers must come before R headers */

#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>

#include "Mdefines.h"
#include "Minlines.h"

#ifdef __cplusplus
extern "C" {
/* NB: this block must not include system or R headers */
#endif

SEXP NEW_OBJECT_OF_CLASS(const char* what);

Rboolean DimNames_is_trivial(SEXP dn);
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

SEXP get_factor(SEXP obj, const char *nm);
void set_factor(SEXP obj, const char *nm, SEXP val);
SEXP R_set_factor(SEXP obj, SEXP val, SEXP nm, SEXP warn);
SEXP R_empty_factors(SEXP obj, SEXP warn);

char type2kind(SEXPTYPE type);
SEXPTYPE kind2type(char kind);
size_t kind2size(char kind);

char Matrix_kind(SEXP obj, int i2d);
SEXP R_Matrix_kind(SEXP obj, SEXP i2d);
char Matrix_shape(SEXP obj);
SEXP R_Matrix_shape(SEXP obj);
char Matrix_repr(SEXP obj);
SEXP R_Matrix_repr(SEXP obj);

SEXP R_index_triangle(SEXP n_, SEXP upper_, SEXP diag_, SEXP packed_);
SEXP R_index_diagonal(SEXP n_, SEXP upper_,             SEXP packed_);

SEXP R_nnz(SEXP x, SEXP countNA, SEXP nnzmax);

void conjugate(SEXP x);
void zeroRe(SEXP x);
void zeroIm(SEXP x);
void na2one(SEXP x);

SEXP v2spV(SEXP from);

Rboolean equal_string_vectors(SEXP s1, SEXP s2, int n);
SEXP append_to_named_list(SEXP x, const char *nm, SEXP val);

char La_norm_type(const char *typstr);
char La_rcond_type(const char *typstr);
SEXP as_det_obj(double mod, int log, int sign);

SEXP Matrix_expand_pointers(SEXP pP);
SEXP m_encodeInd (SEXP ij,        SEXP di, SEXP orig_1, SEXP chk_bnds);
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds);
SEXP Mmatrix(SEXP args);

SEXP R_rbind2_vector(SEXP a, SEXP b);
SEXP R_all0(SEXP x);
SEXP R_any0(SEXP x);

    
/* ================================================================== */
/* Defined elsewhere but used in a few places, hence "exported" here: */
/* ================================================================== */
    
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

#define UNPACKED_MAKE_TRIANGULAR(_PREFIX_, _CTYPE_)			\
void _PREFIX_ ## dense_unpacked_make_triangular(_CTYPE_ *x, int m, int n, \
						char uplo, char diag)
UNPACKED_MAKE_TRIANGULAR(d, double);
UNPACKED_MAKE_TRIANGULAR(i, int);
UNPACKED_MAKE_TRIANGULAR(z, Rcomplex);
#undef UNPACKED_MAKE_TRIANGULAR

#define UNPACKED_MAKE_SYMMETRIC(_PREFIX_, _CTYPE_)			\
void _PREFIX_ ## dense_unpacked_make_symmetric(_CTYPE_ *x, int n, char uplo)
UNPACKED_MAKE_SYMMETRIC(d, double);
UNPACKED_MAKE_SYMMETRIC(i, int);
UNPACKED_MAKE_SYMMETRIC(z, Rcomplex);
#undef UNPACKED_MAKE_SYMMETRIC

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

SEXP unpacked_force(SEXP x, int n, char uplo, char diag);
SEXP packed_transpose(SEXP x, int n, char uplo);

SEXP dense_as_general(SEXP from, char kind, int new, int transpose_if_vector);

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_UTILS_H */
