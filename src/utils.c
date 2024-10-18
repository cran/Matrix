#include <stdio.h> /* vsnprintf */
#include "Mdefines.h"
#include "utils.h"

#if R_VERSION < R_Version(4, 5, 0)
int ANY_ATTRIB(SEXP x)
{
	return ATTRIB(x) != R_NilValue;
}

void CLEAR_ATTRIB(SEXP x)
{
	SET_ATTRIB(x, R_NilValue);
	SET_OBJECT(x, 0);
	UNSET_S4_OBJECT(x);
	return;
}
#endif

/* memset() but passing length and size rather than their product
   which can overflow size_t ... hence _safer_ than Memzero()
*/
void *Matrix_memset(void *dest, int ch, R_xlen_t length, size_t size)
{
	if (dest && length > 0 && size > 0) {

		char *dest_ = (char *) dest;
		size_t N = SIZE_MAX / size;

#if (SIZE_MAX < R_XLEN_T_MAX)
		R_xlen_t S_M = (R_xlen_t) SIZE_MAX;
		if (length <= S_M) {
#endif

			/* 'length' is representable as size_t : */

			size_t n = (size_t) length;
			if (n <= N)
				memset(dest_, ch, n * size);
			else {
				size_t d = N * size;
				while (n > N) {
					memset(dest_, ch, d);
					dest_ += d;
					n -= d;
				}
				memset(dest_, ch, n * size);
			}

#if (SIZE_MAX < R_XLEN_T_MAX)
		} else {

			/* 'length' would overflow size_t : */

			size_t n, d = N * size;
			while (length > S_M) {
				n = SIZE_MAX;
				while (n > N) {
					memset(dest_, ch, d);
					dest_ += d;
					n -= d;
				}
				memset(dest_, ch, n * size);
				length -= S_M;
			}
			n = (size_t) length;
			while (n > N) {
				memset(dest_, ch, d);
				dest_ += d;
				n -= d;
			}
			memset(dest_, ch, n * size);

		}
#endif

	}

	return dest;
}

/* memcpy() but passing length and size rather than their product
   which can overflow size_t ... hence _safer_ than Memcpy()
*/
void *Matrix_memcpy(void *dest, const void *src, R_xlen_t length, size_t size)
{
	if (dest && src && length > 0 && size > 0) {

		char *dest_ = (char *) dest;
		const char *src_ = (const char *) src;

		size_t N = SIZE_MAX / size;

#if (SIZE_MAX < R_XLEN_T_MAX)
		R_xlen_t S_M = (R_xlen_t) SIZE_MAX;
		if (length <= S_M) {
#endif

			/* 'length' is representable as size_t : */

			size_t n = (size_t) length;
			if (n <= N)
				memcpy(dest_, src_, n * size);
			else {
				size_t d = N * size;
				while (n > N) {
					memcpy(dest_, src_, d);
					dest_ += d;
					src_ += d;
					n -= d;
				}
				memcpy(dest_, src_, n * size);
			}

#if (SIZE_MAX < R_XLEN_T_MAX)
		} else {

			/* 'length' would overflow size_t : */

			size_t n, d = N * size;
			while (length > S_M) {
				n = SIZE_MAX;
				while (n > N) {
					memcpy(dest_, src_, d);
					dest_ += d;
					src_ += d;
					n -= d;
				}
				memcpy(dest_, src_, n * size);
				length -= S_M;
			}
			n = (size_t) length;
			while (n > N) {
				memcpy(dest_, src_, d);
				dest_ += d;
				n -= d;
			}
			memcpy(dest_, src_, n * size);

		}
#endif

	}

	return dest;
}

char *Matrix_sprintf(const char *format, ...)
{
	char *buf = R_alloc(Matrix_ErrorBufferSize, sizeof(char));
	va_list args;
	va_start(args, format);
	vsnprintf(buf, Matrix_ErrorBufferSize, format, args);
	va_end(args);
	return buf;
}

int equal_character_vectors(SEXP s1, SEXP s2, int n)
{
	/* FIXME? not distinguishing between NA_STRING and "NA" */
	for (int i = 0; i < n; ++i)
		if (strcmp(CHAR(STRING_ELT(s1, i)), CHAR(STRING_ELT(s2, i))) != 0)
			return 0;
	return 1;
}

void conjugate(SEXP x)
{
	Rcomplex *px = COMPLEX(x);
	R_xlen_t nx = XLENGTH(x);
	while (nx--) {
		(*px).i = -(*px).i;
		++px;
	}
	return;
}

void zeroRe(SEXP x)
{
	Rcomplex *px = COMPLEX(x);
	R_xlen_t nx = XLENGTH(x);
	while (nx--) {
		(*px).r = 0.0;
		++px;
	}
	return;
}

void zeroIm(SEXP x)
{
	Rcomplex *px = COMPLEX(x);
	R_xlen_t nx = XLENGTH(x);
	while (nx--) {
		(*px).i = 0.0;
		++px;
	}
	return;
}

void naToOne(SEXP x)
{
	R_xlen_t i, n = XLENGTH(x);
	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *px = LOGICAL(x);
		for (i = 0; i < n; ++i, ++px)
			if (*px == NA_LOGICAL)
				*px = 1;
		break;
	}
	case INTSXP:
	{
		int *px = INTEGER(x);
		for (i = 0; i < n; ++i, ++px)
			if (*px == NA_INTEGER)
				*px = 1;
		break;
	}
	case REALSXP:
	{
		double *px = REAL(x);
		for (i = 0; i < n; ++i, ++px)
			if (ISNAN(*px))
				*px = 1.0;
		break;
	}
	case CPLXSXP:
	{
		Rcomplex *px = COMPLEX(x);
		for (i = 0; i < n; ++i, ++px)
			if (ISNAN((*px).r) || ISNAN((*px).i))
				*px = Matrix_zone;
		break;
	}
	default:
		ERROR_INVALID_TYPE(x, __func__);
		break;
	}
	return;
}
