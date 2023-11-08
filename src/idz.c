#include "Mdefines.h"
#include "idz.h"

#define IDZ \
TEMPLATE(i,      int,          0  ,         1  ) \
TEMPLATE(d,   double,          0.0,         1.0) \
TEMPLATE(z, Rcomplex, Matrix_zzero, Matrix_zone)

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
static void _PREFIX_ ## \
swap(int n, _CTYPE_ *x, int incx, _CTYPE_ *y, int incy) \
{ \
	_CTYPE_ tmp; \
	while (n--) { \
		tmp = *x; \
		*x = *y; \
		*y = tmp; \
		x += incx; \
		y += incy; \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
static void _PREFIX_ ## \
syswapr(char uplo, int n, _CTYPE_ *x, int k0, int k1) \
{ \
	_CTYPE_ *x0 = x + (R_xlen_t) k0 * n, *x1 = x + (R_xlen_t) k1 * n, \
		tmp; \
	if (uplo == 'U') { \
		_PREFIX_ ## swap(k0, x0, 1, x1, 1); \
		tmp = x0[k0]; \
		x0[k0] = x1[k1]; \
		x1[k1] = tmp; \
		_PREFIX_ ## swap(k1 - k0 - 1, x0 + k0 + n, n, x1 + k0 + 1, 1); \
		_PREFIX_ ## swap(n  - k1 - 1, x1 + k0 + n, n, x1 + k1 + n, n); \
	} else { \
		_PREFIX_ ## swap(k0, x + k0, n, x + k1, n); \
		tmp = x0[k0]; \
		x0[k0] = x1[k1]; \
		x1[k1] = tmp; \
		_PREFIX_ ## swap(k1 - k0 - 1, x0 + k0 + 1, 1, x0 + k1 + n, n); \
		_PREFIX_ ## swap(n  - k1 - 1, x0 + k1 + 1, 1, x1 + k1 + 1, 1); \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
rowperm2(_CTYPE_ *x, int m, int n, int *p, int off, int invert) \
{ \
	int i, k0, k1; \
	for (i = 0; i < m; ++i) \
		p[i] = -(p[i] - off + 1); \
	if (!invert) { \
		for (i = 0; i < m; ++i) { \
			if (p[i] > 0) \
				continue; \
			k0 = i; \
			p[k0] = -p[k0]; \
			k1 = p[k0] - 1; \
			while (p[k1] < 0) { \
				_PREFIX_ ## swap(n, x + k0, m, x + k1, m); \
				k0 = k1; \
				p[k0] = -p[k0]; \
				k1 = p[k0] - 1; \
			} \
		} \
	} else { \
		for (i = 0; i < m; ++i) { \
			if (p[i] > 0) \
				continue; \
			k0 = i; \
			p[k0] = -p[k0]; \
			k1 = p[k0] - 1; \
			while (k1 != k0) { \
				_PREFIX_ ## swap(n, x + k0, m, x + k1, m); \
				p[k1] = -p[k1]; \
				k1 = p[k1] - 1; \
			} \
		} \
	} \
	for (i = 0; i < m; ++i) \
		p[i] = p[i] + off - 1; \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
symperm2(_CTYPE_ *x, int n, char uplo, int *p, int off, int invert) \
{ \
	int i, k0, k1; \
	for (i = 0; i < n; ++i) \
		p[i] = -(p[i] - off + 1); \
	if (!invert) { \
		for (i = 0; i < n; ++i) { \
			if (p[i] > 0) \
				continue; \
			k0 = i; \
			p[k0] = -p[k0]; \
			k1 = p[k0] - 1; \
			while (p[k1] < 0) { \
				if (k0 < k1) \
					_PREFIX_ ## syswapr(uplo, n, x, k0, k1); \
				else \
					_PREFIX_ ## syswapr(uplo, n, x, k1, k0); \
				k0 = k1; \
				p[k0] = -p[k0]; \
				k1 = p[k0] - 1; \
			} \
		} \
	} else { \
		for (i = 0; i < n; ++i) { \
			if (p[i] > 0) \
				continue; \
			k0 = i; \
			p[k0] = -p[k0]; \
			k1 = p[k0] - 1; \
			while (k1 != k0) { \
				if (k0 < k1) \
					_PREFIX_ ## syswapr(uplo, n, x, k0, k1); \
				else \
					_PREFIX_ ## syswapr(uplo, n, x, k1, k0); \
				p[k1] = -p[k1]; \
				k1 = p[k1] - 1; \
			} \
		} \
	} \
	for (i = 0; i < n; ++i) \
		p[i] = p[i] + off - 1; \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
pack2(_CTYPE_ *dest, const _CTYPE_ *src, int n, char uplo, char diag) \
{ \
	int i, j; \
	R_xlen_t dpos = 0, spos = 0; \
	if (uplo == 'U') { \
		for (j = 0; j < n; spos += n-(++j)) \
			for (i = 0; i <= j; ++i) \
				dest[dpos++] = src[spos++]; \
		if (diag != 'N') { \
			dpos = 0; \
			for (j = 0; j < n; dpos += (++j)+1) \
				dest[dpos] = _ONE_; \
		} \
	} else { \
		for (j = 0; j < n; spos += (++j)) \
			for (i = j; i < n; ++i) \
				dest[dpos++] = src[spos++]; \
		if (diag != 'N') { \
			dpos = 0; \
			for (j = 0; j < n; dpos += n-(j++)) \
				dest[dpos] = _ONE_; \
		} \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
unpack1(_CTYPE_ *dest, const _CTYPE_ *src, int n, char uplo, char diag) \
{ \
	int i, j; \
	R_xlen_t dpos = 0, spos = 0; \
	if (uplo == 'U') { \
		for (j = 0; j < n; dpos += n-(++j)) \
			for (i = 0; i <= j; ++i) \
				dest[dpos++] = src[spos++]; \
	} else { \
		for (j = 0; j < n; dpos += (++j)) \
			for (i = j; i <  n; ++i) \
				dest[dpos++] = src[spos++]; \
	} \
	if (diag != 'N') { \
		dpos = 0; \
		R_xlen_t n1a = (R_xlen_t) n + 1; \
		for (j = 0; j < n; ++j, dpos += n1a) \
			dest[dpos] = _ONE_; \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
transpose2(_CTYPE_ *dest, const _CTYPE_ *src, int m, int n) \
{ \
	R_xlen_t mn1s = (R_xlen_t) m * n - 1; \
	int i, j; \
	for (j = 0; j < m; ++j, src -= mn1s) \
		for (i = 0; i < n; ++i, src += m) \
			*(dest++) = *src; \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
transpose1(_CTYPE_ *dest, const _CTYPE_ *src, int n, char uplo) \
{ \
	int i, j; \
	if (uplo == 'U') { \
		for (j = 0; j < n; ++j) \
			for (i = j; i < n; ++i) \
				*(dest++) = *(src + PACKED_AR21_UP(j, i)); \
	} else { \
		R_xlen_t n2 = (R_xlen_t) n * 2; \
		for (j = 0; j < n; ++j) \
			for (i = 0; i <= j; ++i) \
				*(dest++) = *(src + PACKED_AR21_LO(j, i, n2)); \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define ASSIGN_JJ_i(_X_)
#define ASSIGN_JJ_d(_X_)
#define ASSIGN_JJ_z(_X_) \
	_X_.i = 0.0
#define ASSIGN_JI_i(_X_, _Y_) \
	_X_ = _Y_
#define ASSIGN_JI_d(_X_, _Y_) \
	_X_ = _Y_
#define ASSIGN_JI_z(_X_, _Y_) \
	do { \
		_X_.r =  _Y_.r; \
		_X_.i = -_Y_.i; \
	} while (0)

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
syforce2(_CTYPE_ *x, int n, char uplo) \
{ \
	_CTYPE_ *y = x; \
	int i, j; \
	if (uplo == 'U') { \
		for (j = 0; j < n; ++j) { \
			ASSIGN_JJ_ ## _PREFIX_((*x)); \
			x += 1; \
			y += n; \
			for (i = j + 1; i < n; ++i) { \
				ASSIGN_JI_ ## _PREFIX_((*x), (*y)); \
				x += 1; \
				y += n; \
			} \
			x = y = x + j + 1; \
		} \
	} else { \
		for (j = 0; j < n; ++j) { \
			ASSIGN_JJ_ ## _PREFIX_((*y)); \
			x += 1; \
			y += n; \
			for (i = j + 1; i < n; ++i) { \
				ASSIGN_JI_ ## _PREFIX_((*y), (*x)); \
				x += 1; \
				y += n; \
			} \
			x = y = x + j + 1; \
		} \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#undef ASSIGN_JJ_i
#undef ASSIGN_JJ_d
#undef ASSIGN_JJ_z
#undef ASSIGN_JI_i
#undef ASSIGN_JI_d
#undef ASSIGN_JI_z

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
trforce2(_CTYPE_ *x, int m, int n, char uplo, char diag) \
{ \
	_CTYPE_ *y = x; \
	int i, j, r = (m < n) ? m : n; \
	if (uplo == 'U') { \
		for (j = 0; j < r; ++j) { \
			for (i = j + 1; i < m; ++i) \
				*(++x) = _ZERO_; \
			x += j + 2; \
		} \
	} else { \
		for (j = 0; j < r; ++j) { \
			for (i = 0; i < j; ++i) \
				*(x++) = _ZERO_; \
			x += m - j; \
		} \
		for (j = r; j < n; ++j) \
			for (i = 0; i < m; ++i) \
				*(x++) = _ZERO_; \
	} \
	if (diag != 'N') { \
		R_xlen_t m1a = (R_xlen_t) m + 1; \
		for (j = 0; j < r; ++j) { \
			*y = _ONE_; \
			y += m1a; \
		} \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
band2(_CTYPE_ *x, int m, int n, int a, int b, char diag) \
{ \
	if (m == 0 || n == 0) \
		return; \
	if (a > b || a >= n || b <= -m) { \
		Matrix_memset(x, 0, (R_xlen_t) m * n, sizeof(_CTYPE_)); \
		return; \
	} \
	if (a <= -m) a = 1-m; \
	if (b >=  n) b = n-1; \
 \
	int i, j, i0, i1, \
		j0 = (a < 0) ? 0 : a, \
		j1 = (b < n-m) ? m+b : n; \
 \
	if (j0 > 0) { \
		R_xlen_t dx = (R_xlen_t) m * j0; \
		Matrix_memset(x, 0, dx, sizeof(_CTYPE_)); \
		x += dx; \
	} \
	for (j = j0; j < j1; ++j, x += m) { \
		i0 = j - b; \
		i1 = j - a + 1; \
		for (i = 0; i < i0; ++i) \
			*(x + i) = _ZERO_; \
		for (i = i1; i < m; ++i) \
			*(x + i) = _ZERO_; \
	} \
	if (j1 < n) \
		Matrix_memset(x, 0, (R_xlen_t) m * (n - j1), sizeof(_CTYPE_)); \
	if (diag != 'N' && a <= 0 && b >= 0) { \
		x -= m * (R_xlen_t) j; \
		R_xlen_t m1a = (R_xlen_t) m + 1; \
		for (j = 0; j < n; ++j, x += m1a) \
			*x = _ONE_; \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
band1(_CTYPE_ *x, int n, int a, int b, char uplo, char diag) \
{ \
	if (n == 0) \
		return; \
	if (a > b || a >= n || b <= -n) { \
		Matrix_memset(x, 0, PACKED_LENGTH(n), sizeof(_CTYPE_)); \
		return; \
	} \
	if (uplo == 'U') { \
		if (a <   0) a = 0; \
		if (b >=  n) b = n-1; \
	} else { \
		if (b >   0) b = 0; \
		if (a <= -n) a = 1-n; \
	} \
 \
	int i, j, i0, i1, \
		j0 = (a < 0) ? 0 : a, \
		j1 = (b < 0) ? n+b : n; \
 \
	if (uplo == 'U') { \
		if (j0 > 0) { \
			R_xlen_t dx; \
			Matrix_memset(x, 0, dx = PACKED_LENGTH(j0), \
			              sizeof(_CTYPE_)); \
			x += dx; \
		} \
		for (j = j0; j < j1; x += (++j)) { \
			i0 = j - b; \
			i1 = j - a + 1; \
			for (i = 0; i < i0; ++i) \
				*(x + i) = _ZERO_; \
			for (i = i1; i <= j; ++i) \
				*(x + i) = _ZERO_; \
		} \
		if (j1 < n) \
			Matrix_memset(x, 0, PACKED_LENGTH(n) - PACKED_LENGTH(j1), \
			              sizeof(_CTYPE_)); \
		if (diag != 'N' && a == 0) { \
			x -= PACKED_LENGTH(j); \
			for (j = 0; j < n; x += (++j)+1) \
				*x = _ONE_; \
		} \
	} else { \
		if (j0 > 0) { \
			R_xlen_t dx = PACKED_LENGTH(n) - PACKED_LENGTH(j0); \
			Matrix_memset(x, 0, dx, sizeof(_CTYPE_)); \
			x += dx; \
		} \
		for (j = j0; j < j1; x += n-(j++)) { \
			i0 = j - b; \
			i1 = j - a + 1; \
			for (i = j; i < i0; ++i) \
				*(x + i - j) = _ZERO_; \
			for (i = i1; i < n; ++i) \
				*(x + i - j) = _ZERO_; \
		} \
		if (j1 < n) \
			Matrix_memset(x, 0, PACKED_LENGTH(n - j1), \
			              sizeof(_CTYPE_)); \
		if (diag != 'N' && b == 0) { \
			x -= PACKED_LENGTH(n) - PACKED_LENGTH(j); \
			for (j = 0; j < n; x += n-(j++)) \
				*x = _ONE_; \
		} \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
dcpy2(_CTYPE_ *dest, const _CTYPE_ *src, int n, R_xlen_t length, \
      char uplo, char diag) \
{ \
	int j; \
	R_xlen_t n1a = (R_xlen_t) n + 1; \
	if (diag != 'N') { \
		for (j = 0; j < n; ++j, dest += n1a) \
			*dest = _ONE_; \
	} else if (length == n) { \
		/* copying from diagonalMatrix */ \
		for (j = 0; j < n; ++j, dest += n1a, ++src) \
			*dest = *src; \
	} else if (length == (n * n1a) / 2) { \
		/* copying from packedMatrix */ \
		if (uplo == 'U') { \
			for (j = 0; j < n; dest += n1a, src += (++j)+1) \
				*dest = *src; \
		} else { \
			for (j = 0; j < n; dest += n1a, src += n-(j++)) \
				*dest = *src; \
		} \
	} else if (length == (R_xlen_t) n * n) { \
		/* copying from square unpackedMatrix */ \
		for (j = 0; j < n; ++j, dest += n1a, src += n1a) \
			*dest = *src; \
	} else { \
		error(_("incompatible '%s' and '%s' in '%s'"), \
		      "n", "length", __func__); \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## \
dcpy1(_CTYPE_ *dest, const _CTYPE_ *src, int n, R_xlen_t length, \
      char uplo_dest, char uplo_src, char diag) \
{ \
	int j; \
	if (diag != 'N') { \
		if (uplo_dest == 'U') { \
			for (j = 0; j < n; dest += (++j)+1) \
				*dest = _ONE_; \
		} else { \
			for (j = 0; j < n; dest += n-(j++)) \
				*dest = _ONE_; \
		} \
	} else if (length == n) { \
		/* copying from diagonalMatrix */ \
		if (uplo_dest == 'U') { \
			for (j = 0; j < n; dest += (++j)+1, ++src) \
				*dest = *src; \
		} else { \
			for (j = 0; j < n; dest += n-(j++), ++src) \
				*dest = *src; \
		} \
	} else if (length == PACKED_LENGTH(n)) { \
		/* copying from packedMatrix */ \
		if (uplo_dest == 'U') { \
			if (uplo_src == 'U') { \
				for (j = 0; j < n; src += (++j)+1, dest += j+1) \
					*dest = *src; \
			} else { \
				for (j = 0; j < n; src += n-j, dest += (++j)+1) \
					*dest = *src; \
			} \
		} else { \
			if (uplo_src == 'U') { \
				for (j = 0; j < n; dest += n-(j++), src += j+1) \
					*dest = *src; \
			} else { \
				for (j = 0; j < n; dest += n-j, src += n-(j++)) \
					*dest = *src; \
			} \
		} \
	} else if (length == (R_xlen_t) n * n) { \
		/* copying from square unpackedMatrix */ \
		R_xlen_t n1a = (R_xlen_t) n + 1; \
		if (uplo_dest == 'U') { \
			for (j = 0; j < n; dest += (++j)+1, src += n1a) \
				*dest = *src; \
		} else { \
			for (j = 0; j < n; dest += n-(j++), src += n1a) \
				*dest = *src; \
		} \
	} else { \
		error(_("incompatible '%s' and '%s' in '%s'"), \
		      "n", "length", __func__); \
	} \
	return; \
}
IDZ
#undef TEMPLATE

#undef IDZ
