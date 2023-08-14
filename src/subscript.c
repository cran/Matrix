#include "subscript.h"

#define F_X( _X_)  (_X_)
#define F_ND(_X_) ((_X_) ? 1 : 0)
#define F_NS(_X_)          1

#define AR21_UP(i, j, m) i + j + (        j * (    j - 1)) / 2
#define AR21_LO(i, j, m) i +     (j * m + j * (m - j - 1)) / 2

static SEXP unpackedMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{

#define SUB1_START(_SEXPTYPE_) \
	SEXPTYPE typ = _SEXPTYPE_; \
	R_xlen_t l = 0, len = XLENGTH(w); \
	SEXP res = allocVector(typ, len); \
	if (len == 0) \
		return res; \
	PROTECT(res); \
	 \
	SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	Matrix_int_fast64_t mn64 = (Matrix_int_fast64_t) m * n; \
	UNPROTECT(1);

#define SUB1_START_EXTRA(_SEXPTYPE_) \
	SUB1_START(_SEXPTYPE_); \
	 \
	int ge = cl[1] == 'g', tr = cl[1] == 't', upper = 1, nonunit = 1; \
	if (!ge) { \
		SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym)); \
		upper = *CHAR(STRING_ELT(uplo, 0)) == 'U'; \
		UNPROTECT(1); \
		if (tr) { \
			SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym)); \
			nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N'; \
			UNPROTECT(1); \
		} \
	}

	SUB1_START_EXTRA(kind2type(cl[0]));

#define SUB1_CASES(_SUB1_N_, _SUB1_X_, _F_N_, _F_X_) \
	do { \
		switch (cl[0]) { \
		case 'n': \
			_SUB1_N_(int, LOGICAL, NA_LOGICAL, 0, 1, _F_N_); \
			break; \
		case 'l': \
			_SUB1_X_(int, LOGICAL, NA_LOGICAL, 0, 1, _F_X_); \
			break; \
		case 'i': \
			_SUB1_X_(int, INTEGER, NA_INTEGER, 0, 1, _F_X_); \
			break; \
		case 'd': \
			_SUB1_X_(double, REAL, NA_REAL, 0.0, 1.0, _F_X_); \
			break; \
		case 'z': \
			_SUB1_X_(Rcomplex, COMPLEX, \
			         Matrix_zna, Matrix_zzero, Matrix_zone, _F_X_); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define SUB1_N(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		_CTYPE_ *pres = _PTR_(res); \
		if (TYPEOF(w) == INTSXP) { \
			int *pw = INTEGER(w); \
			if (mn64 >= INT_MAX) { \
				/* index is never out of bounds */ \
				SUB1_LOOP((pw[l] == NA_INTEGER), \
				          _NA_, _ZERO_, _ONE_, _F_, Matrix_int_fast64_t); \
			} else { \
				int mn = m * n; \
				SUB1_LOOP((pw[l] == NA_INTEGER || pw[l] > mn), \
				          _NA_, _ZERO_, _ONE_, _F_, int); \
			} \
		} else { \
			double *pw = REAL(w); \
			if (mn64 >= 0x1.0p+53) \
				/* m*n may not be exactly representable as double */ \
				/* but it does not exceed INT_MAX * INT_MAX       */ \
				SUB1_LOOP((ISNAN(pw[l]) || pw[l] >= 0x1.0p+62 || \
				           (Matrix_int_fast64_t) pw[l] > mn64), \
				          _NA_, _ZERO_, _ONE_, _F_, Matrix_int_fast64_t); \
			else { \
				double mn1a = (double) m * n + 1.0; \
				SUB1_LOOP((ISNAN(pw[l]) || pw[l] >= mn1a), \
				          _NA_, _ZERO_, _ONE_, _F_, Matrix_int_fast64_t); \
			} \
		} \
	} while (0)

#define SUB1_X(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		PROTECT(x = GET_SLOT(x, Matrix_xSym)); \
		_CTYPE_ *px = _PTR_(x); \
		SUB1_N(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_, _F_); \
		UNPROTECT(1); \
	} while (0)

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_, _INT_) \
	do { \
		_INT_ index, i_, j_; \
		for (l = 0; l < len; ++l) { \
			if (_NA_SUBSCRIPT_) \
				pres[l] = _NA_; \
			else { \
				index = (_INT_) pw[l] - 1; \
				if (ge) \
					pres[l] = _F_(px[index]); \
				else { \
					i_ = index % m; \
					j_ = index / m; \
					if (tr) { \
						if ((upper) ? i_ > j_ : i_ < j_) \
							pres[l] = _ZERO_; \
						else if (!nonunit && i_ == j_) \
							pres[l] = _ONE_; \
						else \
							pres[l] = _F_(px[index]); \
					} else { \
						if ((upper) ? i_ > j_ : i_ < j_) \
							pres[l] = _F_(px[i_ * m + j_]); \
						else \
							pres[l] = _F_(px[index]); \
					} \
				} \
			} \
		} \
	} while (0)

	SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);

#undef SUB1_LOOP

	UNPROTECT(1);
	return res;
}

static SEXP packedMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
	SUB1_START_EXTRA(kind2type(cl[0]));

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_, _INT_) \
	do { \
		_INT_ index, i_, j_; \
		for (l = 0; l < len; ++l) { \
			if (_NA_SUBSCRIPT_) \
				pres[l] = _NA_; \
			else { \
				index = (_INT_) pw[l] - 1; \
				i_ = index % m; \
				j_ = index / m; \
				if (tr) { \
					if (upper) { \
						if (i_ > j_) \
							pres[l] = _ZERO_; \
						else if (!nonunit && i_ == j_) \
							pres[l] = _ONE_; \
						else \
							pres[l] = _F_(px[AR21_UP(i_, j_, m)]); \
					} else { \
						if (i_ < j_) \
							pres[l] = _ZERO_; \
						else if (!nonunit && i_ == j_) \
							pres[l] = _ONE_; \
						else \
							pres[l] = _F_(px[AR21_LO(i_, j_, m)]); \
					} \
				} else { \
					if (upper) { \
						if (i_ > j_) \
							pres[l] = _F_(px[AR21_UP(j_, i_, m)]); \
						else \
							pres[l] = _F_(px[AR21_UP(i_, j_, m)]); \
					} else { \
						if (i_ < j_) \
							pres[l] = _F_(px[AR21_LO(j_, i_, m)]); \
						else \
							pres[l] = _F_(px[AR21_LO(i_, j_, m)]); \
					} \
				} \
			} \
		} \
	} while (0)

	SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);

#undef SUB1_LOOP

	UNPROTECT(1);
	return res;
}

static SEXP CsparseMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
	SUB1_START(kind2type(cl[0]));

	SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
		i = PROTECT(GET_SLOT(x, Matrix_iSym));
	int *pp = INTEGER(p), *pi = INTEGER(i), i_, j_, j, k = 0, kend;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_, _INT_) \
	do { \
		_INT_ index; \
		if (_NA_SUBSCRIPT_) \
			j = -1; \
		else { \
			index = (_INT_) pw[l] - 1; \
			i_ =     (int) (index % m); \
			j_ = j = (int) (index / m); \
		} \
		while (j >= 0) { \
			k = pp[j]; \
			kend = pp[j + 1]; \
			while (k < kend && j_ == j) { \
				if (pi[k] < i_) \
					++k; \
				else { \
					if (pi[k] > i_) \
						pres[l] = _ZERO_; \
					else \
						pres[l] = _F_(px[k]); \
					++l; \
					if (l == len || _NA_SUBSCRIPT_) \
						j_ = -1; \
					else { \
						index = (_INT_) pw[l] - 1; \
						i_ = (int) (index % m); \
						j_ = (int) (index / m); \
					} \
				} \
			} \
			while (j_ == j) { \
				pres[l] = _ZERO_; \
				++l; \
				if (l == len || _NA_SUBSCRIPT_) \
					j_ = -1; \
				else { \
					index = (_INT_) pw[l] - 1; \
					i_ = (int) (index % m); \
					j_ = (int) (index / m); \
				} \
			} \
			j = j_; \
		} \
		while (l < len) { \
			pres[l] = _NA_; \
			++l; \
		} \
	} while (0)

	SUB1_CASES(SUB1_N, SUB1_X, F_NS, F_X);

#undef SUB1_LOOP

	UNPROTECT(3);
	return res;
}

static SEXP RsparseMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
	SUB1_START(kind2type(cl[0]));

	SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
		j = PROTECT(GET_SLOT(x, Matrix_jSym));
	int *pp = INTEGER(p), *pj = INTEGER(j), i, i_, j_, k = 0, kend;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_, _INT_) \
	do { \
		_INT_ index; \
		if (_NA_SUBSCRIPT_) \
			i = -1; \
		else { \
			index = (_INT_) pw[l] - 1; \
			i_ = i = (int) (index % m); \
			j_ =     (int) (index / m); \
		} \
		while (i >= 0) { \
			k = pp[i]; \
			kend = pp[i + 1]; \
			while (k < kend && i_ == i) { \
				if (pj[k] < j_) \
					++k; \
				else { \
					if (pj[k] > j_) \
						pres[l] = _ZERO_; \
					else \
						pres[l] = _F_(px[k]); \
					++l; \
					if (l == len || _NA_SUBSCRIPT_) \
						i_ = -1; \
					else { \
						index = (_INT_) pw[l] - 1; \
						i_ = (int) (index % m); \
						j_ = (int) (index / m); \
					} \
				} \
			} \
			while (i_ == i) { \
				pres[l] = _ZERO_; \
				++l; \
				if (l == len || _NA_SUBSCRIPT_) \
					i_ = -1; \
				else { \
					index = (_INT_) pw[l] - 1; \
					i_ = (int) (index % m); \
					j_ = (int) (index / m); \
				} \
			} \
			i = i_; \
		} \
		while (l < len) { \
			pres[l] = _NA_; \
			++l; \
		} \
	} while (0)

	SUB1_CASES(SUB1_N, SUB1_X, F_NS, F_X);

#undef SUB1_LOOP

	UNPROTECT(3);
	return res;
}

static SEXP diagonalMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
	SUB1_START(kind2type(cl[0]));

	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	int nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	UNPROTECT(1);

	Matrix_int_fast64_t index, n1a = (Matrix_int_fast64_t) n + 1;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_, _INT_) \
	do { \
		for (l = 0; l < len; ++l) { \
			if (_NA_SUBSCRIPT_) \
				pres[l] = _NA_; \
			else if ((index = (Matrix_int_fast64_t) pw[l] - 1) % n1a != 0) \
				pres[l] = _ZERO_; \
			else if (!nonunit) \
				pres[l] = _ONE_; \
			else \
				pres[l] = _F_(px[index / n1a]); \
		} \
	} while (0)

	SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);

#undef SUB1_LOOP

	UNPROTECT(1);
	return res;
}

static SEXP indMatrix_subscript_1ary(SEXP x, SEXP w)
{
	SUB1_START(LGLSXP);

	SEXP perm = PROTECT(GET_SLOT(x, Matrix_permSym));
	int *pperm = INTEGER(perm), i_, j_;

	SEXP margin = PROTECT(GET_SLOT(x, Matrix_marginSym));
	int mg = INTEGER(margin)[0] - 1;
	UNPROTECT(1);

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_, _INT_) \
	do { \
		_INT_ index; \
		for (l = 0; l < len; ++l) { \
			if (_NA_SUBSCRIPT_) \
				pres[l] = _NA_; \
			else { \
				index = (_INT_) pw[l] - 1; \
				i_ = (int) (index % m); \
				j_ = (int) (index / m); \
				if (!mg) \
					pres[l] = j_ == pperm[i_] - 1; \
				else \
					pres[l] = i_ == pperm[j_] - 1; \
			} \
		} \
	} while (0)

	SUB1_N(int, LOGICAL, NA_LOGICAL, 0, 1, );

#undef SUB1_LOOP
#undef SUB1_N
#undef SUB1_START
#undef SUB1_START_EXTRA

	UNPROTECT(2);
	return res;
}

/* x[i] with 'i' of type "integer" or "double" {i >= 1 or NA} */
SEXP R_subscript_1ary(SEXP x, SEXP i)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
	const char *cl = valid[ivalid];
	validObject(x, cl);

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
		return unpackedMatrix_subscript_1ary(x, i, cl);
	case 'p':
		return   packedMatrix_subscript_1ary(x, i, cl);

	/* NB: for [CRT], the caller must preprocess 'x' and 'i';
	   symmetric and unit triangular 'x' are not handled specially,
	   and it is assumed for speed that 'i' is sorted by row [R]
	   or column [CT] with NA last
	*/

	case 'C':
		return  CsparseMatrix_subscript_1ary(x, i, cl);
	case 'R':
		return  RsparseMatrix_subscript_1ary(x, i, cl);
	case 'T':
	{
		char cl_[] = "..CMatrix";
		cl_[0] = cl[0];
		cl_[1] = cl[1];

		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);

		x = sparse_as_Csparse(x, cl);
		PROTECT(x);
		x = CsparseMatrix_subscript_1ary(x, i, cl_);
		UNPROTECT(1);
		return x;
	}
	case 'i':
		return diagonalMatrix_subscript_1ary(x, i, cl);
	default:
		return      indMatrix_subscript_1ary(x, i);
	}
}

static SEXP unpackedMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{

#define SUB1_START(_SEXPTYPE_) \
	SEXPTYPE typ = _SEXPTYPE_; \
	int l = 0, len = (int) (XLENGTH(w) / 2); \
	SEXP res = allocVector(typ, len); \
	if (len == 0) \
		return res; \
	PROTECT(res); \
	int *pw0 = INTEGER(w), *pw1 = pw0 + len;

#define SUB1_START_EXTRA(_SEXPTYPE_) \
	SUB1_START(_SEXPTYPE_); \
	 \
	SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym)); \
	int m = INTEGER(dim)[0]; \
	UNPROTECT(1); \
	 \
	int ge = cl[1] == 'g', tr = cl[1] == 't', upper = 1, nonunit = 1; \
	if (!ge) { \
		SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym)); \
		upper = *CHAR(STRING_ELT(uplo, 0)) == 'U'; \
		UNPROTECT(1); \
		if (tr) { \
			SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym)); \
			nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N'; \
			UNPROTECT(1); \
		} \
	}

	SUB1_START_EXTRA(kind2type(cl[0]));

	Matrix_int_fast64_t i_, j_;

#define SUB1_N(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		_CTYPE_ *pres = _PTR_(res); \
		SUB1_LOOP((pw0[l] == NA_INTEGER || pw1[l] == NA_INTEGER), \
		          _NA_, _ZERO_, _ONE_, _F_); \
	} while (0)

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		for (l = 0; l < len; ++l) { \
			if (_NA_SUBSCRIPT_) \
				pres[l] = _NA_; \
			else { \
				i_ = pw0[l] - 1; \
				j_ = pw1[l] - 1; \
				if (ge) \
					pres[l] = _F_(px[j_ * m + i_]); \
				else if (tr) { \
					if ((upper) ? i_ > j_ : i_ < j_) \
						pres[l] = _ZERO_; \
					else if (!nonunit && i_ == j_) \
						pres[l] = _ONE_; \
					else \
						pres[l] = _F_(px[j_ * m + i_]); \
				} else { \
					if ((upper) ? i_ > j_ : i_ < j_) \
						pres[l] = _F_(px[i_ * m + j_]); \
					else \
						pres[l] = _F_(px[j_ * m + i_]); \
				} \
			} \
		} \
	} while (0)

	SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);

#undef SUB1_LOOP

	UNPROTECT(1);
	return res;
}

static SEXP packedMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
	SUB1_START_EXTRA(kind2type(cl[0]));

	Matrix_int_fast64_t i_, j_;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		for (l = 0; l < len; ++l) { \
			if (_NA_SUBSCRIPT_) \
				pres[l] = _NA_; \
			else { \
				i_ = pw0[l] - 1; \
				j_ = pw1[l] - 1; \
				if (tr) { \
					if (upper) { \
						if (i_ > j_) \
							pres[l] = _ZERO_; \
						else if (!nonunit && i_ == j_) \
							pres[l] = _ONE_; \
						else \
							pres[l] = _F_(px[AR21_UP(i_, j_, m)]); \
					} else { \
						if (i_ < j_) \
							pres[l] = _ZERO_; \
						else if (!nonunit && i_ == j_) \
							pres[l] = _ONE_; \
						else \
							pres[l] = _F_(px[AR21_LO(i_, j_, m)]); \
					} \
				} else { \
					if (upper) { \
						if (i_ > j_) \
							pres[l] = _F_(px[AR21_UP(j_, i_, m)]); \
						else \
							pres[l] = _F_(px[AR21_UP(i_, j_, m)]); \
					} else { \
						if (i_ < j_) \
							pres[l] = _F_(px[AR21_LO(j_, i_, m)]); \
						else \
							pres[l] = _F_(px[AR21_LO(i_, j_, m)]); \
					} \
				} \
			} \
		} \
	} while (0)

	SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);

#undef SUB1_LOOP

	UNPROTECT(1);
	return res;
}

static SEXP CsparseMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
	SUB1_START(kind2type(cl[0]));

	SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
		i = PROTECT(GET_SLOT(x, Matrix_iSym));
	int *pp = INTEGER(p), *pi = INTEGER(i), i_, j_, j, k = 0, kend;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		if (_NA_SUBSCRIPT_) \
			j = -1; \
		else { \
			i_ =     pw0[l] - 1; \
			j_ = j = pw1[l] - 1; \
		} \
		while (j >= 0) { \
			k = pp[j]; \
			kend = pp[j + 1]; \
			while (k < kend && j_ == j) { \
				if (pi[k] < i_) \
					++k; \
				else { \
					if (pi[k] > i_) \
						pres[l] = _ZERO_; \
					else \
						pres[l] = _F_(px[k]); \
					++l; \
					if (l == len || _NA_SUBSCRIPT_) \
						j_ = -1; \
					else { \
						i_ = pw0[l] - 1; \
						j_ = pw1[l] - 1; \
					} \
				} \
			} \
			while (j_ == j) { \
				pres[l] = _ZERO_; \
				++l; \
				if (l == len || _NA_SUBSCRIPT_) \
					j_ = -1; \
				else { \
					i_ = pw0[l] - 1; \
					j_ = pw1[l] - 1; \
				} \
			} \
			j = j_; \
		} \
		while (l < len) { \
			pres[l] = _NA_; \
			++l; \
		} \
	} while (0)

	SUB1_CASES(SUB1_N, SUB1_X, F_NS, F_X);

#undef SUB1_LOOP

	UNPROTECT(3);
	return res;
}

static SEXP RsparseMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
	SUB1_START(kind2type(cl[0]));

	SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
		j = PROTECT(GET_SLOT(x, Matrix_jSym));
	int *pp = INTEGER(p), *pj = INTEGER(j), i, i_, j_, k = 0, kend;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		if (_NA_SUBSCRIPT_) \
			i = -1; \
		else { \
			i_ = i = pw0[l] - 1; \
			j_ =     pw1[l] - 1; \
		} \
		while (i >= 0) { \
			k = pp[i]; \
			kend = pp[i + 1]; \
			while (k < kend && i_ == i) { \
				if (pj[k] < j_) \
					++k; \
				else { \
					if (pj[k] > j_) \
						pres[l] = _ZERO_; \
					else \
						pres[l] = _F_(px[k]); \
					++l; \
					if (l == len || _NA_SUBSCRIPT_) \
						i_ = -1; \
					else { \
						i_ = pw0[l] - 1; \
						j_ = pw1[l] - 1; \
					} \
				} \
			} \
			while (i_ == i) { \
				pres[l] = _ZERO_; \
				++l; \
				if (l == len || _NA_SUBSCRIPT_) \
					i_ = -1; \
				else { \
					i_ = pw0[l] - 1; \
					j_ = pw1[l] - 1; \
				} \
			} \
			i = i_; \
		} \
		while (l < len) { \
			pres[l] = _NA_; \
			++l; \
		} \
	} while (0)

	SUB1_CASES(SUB1_N, SUB1_X, F_NS, F_X);

#undef SUB1_LOOP

	UNPROTECT(3);
	return res;
}

static SEXP diagonalMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
	SUB1_START(kind2type(cl[0]));

	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	int nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	UNPROTECT(1);

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		for (l = 0; l < len; ++l) { \
			if (_NA_SUBSCRIPT_) \
				pres[l] = _NA_; \
			else if (pw0[l] != pw1[l]) \
				pres[l] = _ZERO_; \
			else if (!nonunit) \
				pres[l] = _ONE_; \
			else \
				pres[l] = _F_(px[pw0[l] - 1]); \
		} \
	} while (0)

	SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);

#undef SUB1_LOOP

	UNPROTECT(1);
	return res;
}

static SEXP indMatrix_subscript_1ary_mat(SEXP x, SEXP w)
{
	SUB1_START(LGLSXP);

	SEXP perm = PROTECT(GET_SLOT(x, Matrix_permSym));
	int *pperm = INTEGER(perm);

	SEXP margin = PROTECT(GET_SLOT(x, Matrix_marginSym));
	int mg = INTEGER(margin)[0] - 1;
	UNPROTECT(1);

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_) \
	do { \
		for (l = 0; l < len; ++l) { \
			if (_NA_SUBSCRIPT_) \
				pres[l] = _NA_; \
			else if (!mg) \
				pres[l] = pw1[l] == pperm[pw0[l] - 1]; \
			else \
				pres[l] = pw0[l] == pperm[pw1[l] - 1]; \
		} \
	} while (0)

	SUB1_N(int, LOGICAL, NA_LOGICAL, 0, 1, );

#undef SUB1_LOOP
#undef SUB1_N
#undef SUB1_X
#undef SUB1_CASES
#undef SUB1_START
#undef SUB1_START_EXTRA

	UNPROTECT(2);
	return res;
}

/* x[i] with 'i' of type "integer" and dimensions c(.,2)
   {i[,1] in 1:m or NA, i[,2] in 1:n or NA}
*/
SEXP R_subscript_1ary_mat(SEXP x, SEXP i)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 1);
	const char *cl = valid[ivalid];
	validObject(x, cl);

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
		return unpackedMatrix_subscript_1ary_mat(x, i, cl);
	case 'p':
		return   packedMatrix_subscript_1ary_mat(x, i, cl);

	/* NB: for [CRT], the caller must preprocess 'x' and 'i';
	   symmetric and unit triangular 'x' are not handled specially,
	   and it is assumed for speed that 'i' is sorted by row [R]
	   or column [CT] with NA (incl. out-of-bounds indices) last
	*/

	case 'C':
		return  CsparseMatrix_subscript_1ary_mat(x, i, cl);
	case 'R':
		return  RsparseMatrix_subscript_1ary_mat(x, i, cl);
	case 'T':
	{
		char cl_[] = "..CMatrix";
		cl_[0] = cl[0];
		cl_[1] = cl[1];

		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);

		x = sparse_as_Csparse(x, cl);
		PROTECT(x);
		x = CsparseMatrix_subscript_1ary_mat(x, i, cl_);
		UNPROTECT(1);
		return x;
	}
	case 'i':
		return diagonalMatrix_subscript_1ary_mat(x, i, cl);
	default:
		return      indMatrix_subscript_1ary_mat(x, i);
	}
}

static int keep_tr(int *pi, int *pj, int n, int upper, int nonunit, int checkNA)
{
	int k, ident = memcmp(pi, pj, n * sizeof(int)) == 0;
	if (checkNA) {
		if (ident) {
			for (k = 0; k < n; ++k)
				if (pi[k] == NA_INTEGER)
					return 0;
		} else {
			for (k = 0; k < n; ++k)
				if (pi[k] == NA_INTEGER || pj[k] == NA_INTEGER)
					return 0;
		}
	}
	int r = (upper) ? 1 : -1;
	if (ident) {
		/* triangular iff monotone; unit diagonal is preserved */
		if (n >= 2) {
			if (pi[0] == pi[1])
				return 0;
			else if (pi[0] < pi[1]) {
				for (k = 2; k < n; ++k)
					if (pi[k - 1] >= pi[k])
						return 0;
				/* up->up, lo->lo */
			} else {
				for (k = 2; k < n; ++k)
					if (pi[k - 1] <= pi[k])
						return 0;
				/* up->lo, lo->up */
				r = -r;
			}
		}
		if (!nonunit)
			r *= 2;
		return r;
	} else {
		/* brute force ... */
		int ki, kj, j_;
		if (upper) {
			for (kj = 0; kj < n; ++kj)
				for (ki = kj + 1, j_ = pj[kj]; ki < n; ++ki)
					if (pi[ki] <= j_)
						goto LO;
			/* up->up */
			return  r;
		LO:
			for (kj = 0; kj < n; ++kj)
				for (ki = 0, j_ = pj[kj]; ki < kj; ++ki)
					if (pi[ki] <= j_)
						return 0;
			/* up->lo */
			return -r;
		} else {
			for (kj = 0; kj < n; ++kj)
				for (ki = 0, j_ = pj[kj]; ki < kj; ++ki)
					if (pi[ki] >= j_)
						goto UP;
			/* lo->lo */
			return  r;
		UP:
			for (kj = 0; kj < n; ++kj)
				for (ki = kj + 1, j_ = pj[kj]; ki < n; ++ki)
					if (pi[ki] >= j_)
						return 0;
			/* lo->up */
			return -r;
		}
	}
}

static int keep_sy(int *pi, int *pj, int n, int upper, int checkNA)
{
	if (memcmp(pi, pj, n * sizeof(int)) != 0)
		return 0;
	int k, r = (upper) ? 1 : -1;
	if (checkNA) {
		for (k = 0; k < n; ++k)
			if (pi[k] == NA_INTEGER)
				return r;
	}
	if (n >= 2) {
		/* triangular iff monotone */
		if (pi[0] == pi[1])
			return r;
		else if (pi[0] < pi[1]) {
			for (k = 2; k < n; ++k)
				if (pi[k - 1] >= pi[k])
					return r;
			/* up->up, lo->lo */
		} else {
			for (k = 2; k < n; ++k)
				if (pi[k - 1] <= pi[k])
					return r;
			/* up->lo, lo->up */
			r = -r;
		}
	}
	return 2 * r;
}

static int keep_di(int *pi, int *pj, int n, int nonunit, int checkNA, int lwork)
{
	int k, ident = memcmp(pi, pj, n * sizeof(int)) == 0;
	if (checkNA) {
		if (ident) {
			for (k = 0; k < n; ++k)
				if (pi[k] == NA_INTEGER)
					return 0;
		} else {
			for (k = 0; k < n; ++k)
				if (pi[k] == NA_INTEGER || pj[k] == NA_INTEGER)
					return 0;
		}
	}
	if (ident) {
		/* diagonal iff no duplicates; unit diagonal is preserved */
		char *work;
		Matrix_Calloc(work, lwork, char);
		--work;
		for (k = 0; k < n; ++k) {
			if (work[pi[k]])
				return 0;
			work[pi[k]] = 1;
		}
		++work;
		Matrix_Free(work, lwork);
		return (nonunit) ? 1 : 2;
	} else {
		/* brute force ... */
		int ki, kj, j_;
		for (kj = 0; kj < n; ++kj) {
			j_ = pj[kj];
			for (ki = 0; ki < kj; ++ki)
				if (pi[ki] == j_)
					return 0;
			for (ki = kj + 1; ki < n; ++ki)
				if (pi[ki] == j_)
					return 0;
		}
		return 1;
	}
}

static void sort_cr(SEXP obj, const char *cl)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim),
		m = pdim[(cl[2] == 'C') ? 0 : 1],
		n = pdim[(cl[2] == 'C') ? 1 : 0],
		r = (m < n) ? n : m;
	UNPROTECT(1); /* dim */

	SEXP iSym = (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym,
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, iSym));
	int *pp = INTEGER(p), *pi = INTEGER(i);

	int i_, j_, k, kend, nnz = pp[n], *workA, *workB, *workC;
	size_t lwork = (size_t) m + 1 + r + nnz;
	Matrix_Calloc(workA, lwork, int);
	workB = workA + m + 1;
	workC = workB + r;

#define SORT_LOOP(_MASK_) \
	do { \
		++workA; \
		for (k = 0; k < nnz; ++k) \
			++workA[pi[k]]; \
		--workA; \
		 \
		for (i_ = 1; i_ < m; ++i_) \
			workA[i_] = workB[i_] = workB[i_ - 1] + workA[i_]; \
		workA[m] = nnz; \
		 \
		++pp; \
		k = 0; \
		for (j_ = 0; j_ < n; ++j_) { \
			kend = pp[j_]; \
			while (k < kend) { \
				i_ = pi[k]; \
				workC[workB[i_]] = j_; \
				_MASK_(workD[workB[i_]] = px[k]); \
				++workB[i_]; \
				++k; \
			} \
		} \
		--pp; \
		 \
		for (j_ = 0; j_ < n; ++j_) \
			workB[j_] = pp[j_]; \
		 \
		++workA; \
		k = 0; \
		for (i_ = 0; i_ < m; ++i_) { \
			kend = workA[i_]; \
			while (k < kend) { \
				j_ = workC[k]; \
				pi[workB[j_]] = i_; \
				_MASK_(px[workB[j_]] = workD[k]); \
				++workB[j_]; \
				++k; \
			} \
		} \
		--workA; \
	} while (0)

#define SORT(_CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *workD; \
		Matrix_Calloc(workD, nnz,	_CTYPE_); \
		SORT_LOOP(SHOW); \
		Matrix_Free(workD, nnz); \
	} while (0)

	if (cl[0] == 'n')
		SORT_LOOP(HIDE);
	else {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		switch (cl[0]) {
		case 'l':
			SORT(int, LOGICAL);
			break;
		case 'i':
			SORT(int, INTEGER);
			break;
		case 'd':
			SORT(double, REAL);
			break;
		case 'z':
			SORT(Rcomplex, COMPLEX);
			break;
		default:
			break;
		}
		UNPROTECT(1); /* x */
	}

#undef SORT_LOOP
#undef SORT

	Matrix_Free(workA, lwork);
	UNPROTECT(2); /* i, p */
	return;
}

#define XIJ_GE(    _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	*(_X_ + _J_ * _M_ + _I_)

#define XIJ_TR_U_N(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ <= _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : _ZERO_)

#define XIJ_TR_U_U(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ < _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_TR_L_N(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ >= _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : _ZERO_)

#define XIJ_TR_L_U(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ > _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_TP_U_N(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ <= _J_) \
	 ? *(_X_ + AR21_UP(_I_, _J_, _M_)) \
	 : _ZERO_)

#define XIJ_TP_U_U(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ < _J_) \
	 ? *(_X_ + AR21_UP(_I_, _J_, _M_)) \
	 : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_TP_L_N(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ >= _J_) \
	 ? *(_X_ + AR21_LO(_I_, _J_, _M_)) \
	 : _ZERO_)

#define XIJ_TP_L_U(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ > _J_) \
	 ? *(_X_ + AR21_LO(_I_, _J_, _M_)) \
	 : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_SY_U(  _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ <= _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : XIJ_GE(_X_, _J_, _I_, _M_, _ZERO_, _ONE_))

#define XIJ_SY_L(  _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ >= _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : XIJ_GE(_X_, _J_, _I_, _M_, _ZERO_, _ONE_))

#define XIJ_SP_U(  _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ <= _J_) \
	 ? *(_X_ + AR21_UP(_I_, _J_, _M_)) \
	 : *(_X_ + AR21_UP(_J_, _I_, _M_)))

#define XIJ_SP_L(  _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ >= _J_) \
	 ? *(_X_ + AR21_LO(_I_, _J_, _M_)) \
	 : *(_X_ + AR21_LO(_J_, _I_, _M_)))

static SEXP unpackedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
                                          const char *cl)
{

#define SUB2_START \
	SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	UNPROTECT(1); /* dim */ \
	 \
	int ki, kj, \
		mi = i == R_NilValue, \
		mj = j == R_NilValue, \
		ni = (mi) ? m : LENGTH(i), \
		nj = (mj) ? n : LENGTH(j), \
		*pi = (mi) ? NULL : INTEGER(i), \
		*pj = (mj) ? NULL : INTEGER(j);

#define SUB2_START_EXTRA(_E_, _R_, _Y_, _DENSE_) \
	SUB2_START; \
	 \
	int upper = 1, nonunit = 1, keep = 0; \
	SEXP uplo, diag; \
	if (cl[1] != 'g') { \
		PROTECT(uplo = GET_SLOT(x, Matrix_uploSym)); \
		upper = *CHAR(STRING_ELT(uplo, 0)) == 'U'; \
		UNPROTECT(1); /* uplo */ \
		if (cl[1] == 't') { \
			PROTECT(diag = GET_SLOT(x, Matrix_diagSym)); \
			nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N'; \
			UNPROTECT(1); /* diag */ \
		} \
	} \
	 \
	char cl_[] = "...Matrix"; \
	cl_[0] = cl[0]; \
	cl_[1] = 'g'; \
	cl_[2] = _E_; \
	if (cl[1] != 'g' && !(mi || mj) && ni == nj) { \
		if (cl[1] == 't') { \
			keep = keep_tr(pi, pj, ni, upper, nonunit, _DENSE_); \
			if (keep != 0) { \
				cl_[1] = 't'; \
				cl_[2] = _R_; \
			} \
		} else { \
			keep = keep_sy(pi, pj, ni, upper, 0); \
			if ((_DENSE_) ? keep != 0 : keep < -1 || keep > 1) { \
				cl_[1] = 's'; \
				cl_[2] = _Y_; \
			} \
		} \
	} \
	SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl_)); \
	 \
	PROTECT(dim = GET_SLOT(res, Matrix_DimSym)); \
	pdim = INTEGER(dim); \
	pdim[0] = ni; \
	pdim[1] = nj; \
	UNPROTECT(1); /* dim */ \
	 \
	if ((cl[1] != 's') ? keep < 0 : keep < -1) { \
		PROTECT(uplo = GET_SLOT(res, Matrix_uploSym)); \
		SEXP uplo_ = PROTECT(mkChar("L")); \
		SET_STRING_ELT(uplo, 0, uplo_); \
		UNPROTECT(2); /* uplo_, uplo */ \
	} \
	if (cl[1] == 't' && (keep < -1 || keep > 1)) { \
		PROTECT(diag = GET_SLOT(res, Matrix_diagSym)); \
		SEXP diag_ = PROTECT(mkChar("U")); \
		SET_STRING_ELT(diag, 0, diag_); \
		UNPROTECT(2); /* diag_, diag */ \
	}

	SUB2_START_EXTRA('e', 'r', 'y', 1);

	double ninj = (double) ni * nj;
	if (ninj > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) ninj));

	int i_, j_;
	Matrix_int_fast64_t m_ = m;

#define SUB2_CASES(_SUB2_) \
	do { \
		switch (cl[0]) { \
		case 'n': \
		case 'l': \
			_SUB2_(int, LOGICAL, NA_LOGICAL, 0, 1); \
			break; \
		case 'i': \
			_SUB2_(int, INTEGER, NA_INTEGER, 0, 1); \
			break; \
		case 'd': \
			_SUB2_(double, REAL, NA_REAL, 0.0, 1.0); \
			break; \
		case 'z': \
			_SUB2_(Rcomplex, COMPLEX, \
			       Matrix_zna, Matrix_zzero, Matrix_zone); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (cl_[1] == 'g') { \
			if (cl[1] == 'g') \
				SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
				          XIJ_GE, , , _NA_, _ZERO_, _ONE_); \
			else if (cl[1] == 't') { \
				if (upper) { \
					if (nonunit) \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TR_U_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TR_U_U, , , _NA_, _ZERO_, _ONE_); \
				} else { \
					if (nonunit) \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TR_L_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TR_L_U, , , _NA_, _ZERO_, _ONE_); \
				} \
			} else { \
				if (upper) \
					SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
					          XIJ_SY_U, , , _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
					          XIJ_SY_L, , , _NA_, _ZERO_, _ONE_); \
			} \
		} else if (cl_[1] == 't') { \
			Matrix_memset(px1, 0, XLENGTH(x1), sizeof(_CTYPE_)); \
			if (upper) { \
				if (nonunit) { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TR_U_N, , px1 += ni - kj - 1, \
						          _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TR_U_N, px1 += kj, , \
						          _NA_, _ZERO_, _ONE_); \
				} else { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TR_U_U, , px1 += ni - kj - 1, \
						          _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TR_U_U, px1 += kj, , \
						          _NA_, _ZERO_, _ONE_); \
				} \
			} else { \
				if (nonunit) { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TR_L_N, , px1 += ni - kj - 1, \
						          _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TR_L_N, px1 += kj, , \
						          _NA_, _ZERO_, _ONE_); \
				} else { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TR_L_U, , px1 += ni - kj - 1, \
						          _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TR_L_U, px1 += kj, , \
						          _NA_, _ZERO_, _ONE_); \
				} \
			} \
		} else { \
			Matrix_memset(px1, 0, XLENGTH(x1), sizeof(_CTYPE_)); \
			if (upper) { \
				if (keep > 0) \
					SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
					          XIJ_SY_U, , px1 += ni - kj - 1, \
					          _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
					          XIJ_SY_U, px1 += kj, , \
					          _NA_, _ZERO_, _ONE_); \
			} else { \
				if (keep > 0) \
					SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
					          XIJ_SY_L, , px1 += ni - kj - 1, \
					          _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
					          XIJ_SY_L, px1 += kj, , \
					          _NA_, _ZERO_, _ONE_); \
			} \
		} \
	} while (0)

#define SUB2_LOOP(_FOR_, _XIJ_, _JUMP1_, _JUMP2_, \
                  _NA_, _ZERO_, _ONE_) \
	do { \
		for (kj = 0; kj < nj; ++kj) { \
			if (mj) \
				j_ = kj; \
			else { \
				j_ = pj[kj]; \
				if (j_ != NA_INTEGER) \
					j_ -= 1; \
				else { \
					_JUMP1_; \
					_FOR_ { \
						*(px1++) = _NA_; \
					} \
					_JUMP2_; \
					continue; \
				} \
			} \
			_JUMP1_; \
			_FOR_ { \
				if (mi) \
					i_ = ki; \
				else { \
					i_ = pi[ki]; \
					if (i_ != NA_INTEGER) \
						i_ -= 1; \
					else { \
						*(px1++) = _NA_; \
						continue; \
					} \
				} \
				*(px1++) = _XIJ_(px0, i_, j_, m_, _ZERO_, _ONE_); \
			} \
			_JUMP2_; \
		} \
	} while (0)

	SUB2_CASES(SUB2);

#undef SUB2

	SET_SLOT(res, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, res */
	return res;
}

static SEXP packedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
                                        const char *cl)
{
	SUB2_START_EXTRA('e', 'p', 'p', 1);

	double ninj = (double) ni * nj,
	ninj_ = (keep) ? 0.5 * (ninj + ni) : ninj;
	if (ninj_ > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) ninj_));

	int i_, j_;
	Matrix_int_fast64_t m_ = m;

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (cl_[1] == 'g') { \
			if (cl[1] == 't') { \
				if (upper) { \
					if (nonunit) \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TP_U_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TP_U_U, , , _NA_, _ZERO_, _ONE_); \
				} else { \
					if (nonunit) \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TP_L_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TP_L_U, , , _NA_, _ZERO_, _ONE_); \
				} \
			} else { \
				if (upper) \
					SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
					          XIJ_SP_U, , , _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
					          XIJ_SP_L, , , _NA_, _ZERO_, _ONE_); \
			} \
		} else if (cl_[1] == 't') { \
			if (upper) { \
				if (nonunit) { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TP_U_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TP_U_N, , , _NA_, _ZERO_, _ONE_); \
				} else { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TP_U_U, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TP_U_U, , , _NA_, _ZERO_, _ONE_); \
				} \
			} else { \
				if (nonunit) { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TP_L_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TP_L_N, , , _NA_, _ZERO_, _ONE_); \
				} else { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TP_L_U, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TP_L_U, , , _NA_, _ZERO_, _ONE_); \
				} \
			} \
		} else { \
			if (upper) { \
				if (keep > 0) \
					SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
					          XIJ_SP_U, , , _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
					          XIJ_SP_U, , , _NA_, _ZERO_, _ONE_); \
			} else { \
				if (keep > 0) \
					SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
					          XIJ_SP_L, , , _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
					          XIJ_SP_L, , , _NA_, _ZERO_, _ONE_); \
			} \
		} \
	} while (0)

	SUB2_CASES(SUB2);

#undef SUB2_LOOP
#undef SUB2

	SET_SLOT(res, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, res */
	return res;
}

static SEXP CsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
                                         const char *cl)
{
	SUB2_START_EXTRA('C', 'C', 'C', 0);

	if (cl[1] != 'g' && cl_[1] == 'g') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		x = sparse_as_general(x, cl);
	}
	PROTECT(x);

	SEXP p0 = PROTECT(GET_SLOT(x, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(x, Matrix_iSym)),
		p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1)),
		i1 = NULL;
	int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
		*pp1 = INTEGER(p1), *pi1 = NULL,
		d, k, kend, doSort = 0;
	Matrix_int_fast64_t nnz = 0;
	*(pp1++) = 0;

#define SUB2_FINISH \
	if (nnz > INT_MAX) \
		error(_("%s too dense for %s; would have more than %s nonzero entries"), \
		      "x[i,j]", "[CR]sparseMatrix", "2^31-1"); \
	 \
	PROTECT(i1 = allocVector(INTSXP, (R_xlen_t) nnz)); \
	pi1 = INTEGER(i1); \
	 \
	if (cl[0] == 'n') \
		SUB2_LOOP(HIDE); \
	else { \
		SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)), \
			x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) nnz)); \
		SUB2_CASES(SUB2); \
		SET_SLOT(res, Matrix_xSym, x1); \
		UNPROTECT(2); /* x1, x0 */ \
	}

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		SUB2_LOOP(SHOW); \
	} while (0)

	if (mi) {

		for (kj = 0; kj < nj; ++kj) {
			nnz += pp0[pj[kj]] - pp0[pj[kj] - 1];
			pp1[kj] = (int) nnz;
		}

#define SUB2_LOOP(_MASK_) \
		do { \
			for (kj = 0; kj < nj; ++kj) { \
				k = pp0[pj[kj] - 1]; \
				kend = pp0[pj[kj]]; \
				d = kend - k; \
				if (d) { \
					Matrix_memcpy(pi1, pi0 + k, d, sizeof(int)); \
					pi1 += d; \
					_MASK_(Matrix_memcpy(px1, px0 + k, d, sizeof(*px1))); \
					_MASK_(px1 += d); \
				} \
			} \
		} while (0)

		SUB2_FINISH;

#undef SUB2_LOOP

	} else {

		int *workA, *workB, *workC;
		size_t lwork = (size_t) m + m + ni;
		Matrix_Calloc(workA, lwork, int);
		workB = workA + m;
		workC = workB + m;

		/* workA[ i] : size of the set { ki : pi[ki] - 1 == i }
		   workB[ i] : smallest ki such that pi[ki] - 1 == i
		   workC[ki] : smallest ki' > ki such that pi[ki'] == pi[ki]
		*/

		int i_, j_, i_prev = m;

		for (ki = ni - 1; ki >= 0; --ki) {
			i_ = pi[ki] - 1;
			++workA[i_];
			workC[ki] = workB[i_];
			workB[i_] = ki;
			if (i_ > i_prev)
				doSort = 1;
			i_prev = i_;
		}

		for (kj = 0; kj < nj; ++kj) {
			j_ = (mj) ? kj : pj[kj] - 1;
			k = pp0[j_];
			kend = pp0[j_ + 1];
			while (k < kend) {
				nnz += workA[pi0[k]];
				++k;
			}
			pp1[kj] = (int) nnz;
		}

#define SUB2_LOOP(_MASK_) \
		do { \
			for (kj = 0; kj < nj; ++kj) { \
				j_ = (mj) ? kj : pj[kj] - 1; \
				k = pp0[j_]; \
				kend = pp0[j_ + 1]; \
				while (k < kend) { \
					i_ = pi0[k]; \
					d = workA[i_]; \
					ki = workB[i_]; \
					while (d--) { \
						*(pi1++) = ki; \
						_MASK_(*(px1++) = px0[k]); \
						ki = workC[ki]; \
					} \
					++k; \
				} \
			} \
		} while (0)

		SUB2_FINISH;

#undef SUB2_LOOP

		Matrix_Free(workA, lwork);

	}

#undef SUB2_FINISH

	SET_SLOT(res, Matrix_pSym, p1);
	SET_SLOT(res, Matrix_iSym, i1);
	UNPROTECT(4); /* i1, p1, i0, p0 */

	if (doSort)
		sort_cr(res, cl);
	if (cl[1] == 's' && (keep == -1 || keep == 1)) {
		/* defined in ./sparse.c : */
		SEXP sparse_force_symmetric(SEXP, const char *, char);
		res = sparse_force_symmetric(res, cl_, (keep == 1) ? 'U' : 'L');
	}
	UNPROTECT(2); /* x, res */
	return res;
}

static SEXP RsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
                                         const char *cl)
{
	SUB2_START_EXTRA('R', 'R', 'R', 0);

	if (cl[1] != 'g' && cl_[1] == 'g') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		x = sparse_as_general(x, cl);
	}
	PROTECT(x);

	SEXP p0 = PROTECT(GET_SLOT(x, Matrix_pSym)),
		j0 = PROTECT(GET_SLOT(x, Matrix_jSym)),
		p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) ni + 1)),
		j1 = NULL;
	int *pp0 = INTEGER(p0), *pj0 = INTEGER(j0),
		*pp1 = INTEGER(p1), *pj1 = NULL,
		d, k, kend, doSort = 0;
	Matrix_int_fast64_t nnz = 0;
	*(pp1++) = 0;

#define SUB2_FINISH \
	if (nnz > INT_MAX) \
		error(_("%s too dense for %s; would have more than %s nonzero entries"), \
		      "x[i,j]", "[CR]sparseMatrix", "2^31-1"); \
	 \
	PROTECT(j1 = allocVector(INTSXP, (R_xlen_t) nnz)); \
	pj1 = INTEGER(j1); \
	 \
	if (cl[0] == 'n') \
		SUB2_LOOP(HIDE); \
	else { \
		SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)), \
			x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) nnz)); \
		SUB2_CASES(SUB2); \
		SET_SLOT(res, Matrix_xSym, x1); \
		UNPROTECT(2); /* x1, x0 */ \
	}

	if (mj) {

		for (ki = 0; ki < ni; ++ki) {
			nnz += pp0[pi[ki]] - pp0[pi[ki] - 1];
			pp1[ki] = (int) nnz;
		}

#define SUB2_LOOP(_MASK_) \
		do { \
			for (ki = 0; ki < ni; ++ki) { \
				k = pp0[pi[ki] - 1]; \
				kend = pp0[pi[ki]]; \
				d = kend - k; \
				if (d) { \
					Matrix_memcpy(pj1, pj0 + k, d, sizeof(int)); \
					pj1 += d; \
					_MASK_(Matrix_memcpy(px1, px0 + k, d, sizeof(*px1))); \
					_MASK_(px1 += d); \
				} \
			} \
		} while (0)

		SUB2_FINISH;

#undef SUB2_LOOP

	} else {

		int *workA, *workB, *workC;
		size_t lwork = (size_t) n + n + nj;
		Matrix_Calloc(workA, lwork, int);
		workB = workA + n;
		workC = workB + n;

		/* workA[ j] : size of the set { kj : pj[kj] - 1 == j }
		   workB[ j] : smallest ki such that pj[kj] - 1 == j
		   workC[kj] : smallest kj' > kj such that pj[kj'] == pj[kj]
		*/

		int i_, j_, j_prev = n;

		for (kj = nj - 1; kj >= 0; --kj) {
			j_ = pj[kj] - 1;
			++workA[j_];
			workC[kj] = workB[j_];
			workB[j_] = kj;
			if (j_ > j_prev)
				doSort = 1;
			j_prev = j_;
		}

		for (ki = 0; ki < ni; ++ki) {
			i_ = (mi) ? ki : pi[ki] - 1;
			k = pp0[i_];
			kend = pp0[i_ + 1];
			while (k < kend) {
				nnz += workA[pj0[k]];
				++k;
			}
			pp1[ki] = (int) nnz;
		}

#define SUB2_LOOP(_MASK_) \
		do { \
			for (ki = 0; ki < ni; ++ki) { \
				i_ = (mi) ? ki : pi[ki] - 1; \
				k = pp0[i_]; \
				kend = pp0[i_ + 1]; \
				while (k < kend) { \
					j_ = pj0[k]; \
					d = workA[j_]; \
					kj = workB[j_]; \
					while (d--) { \
						*(pj1++) = kj; \
						_MASK_(*(px1++) = px0[k]); \
						kj = workC[kj]; \
					} \
					++k; \
				} \
			} \
		} while (0)

		SUB2_FINISH;

#undef SUB2_LOOP

		Matrix_Free(workA, lwork);

	}

#undef SUB2_FINISH
#undef SUB2

	SET_SLOT(res, Matrix_pSym, p1);
	SET_SLOT(res, Matrix_jSym, j1);
	UNPROTECT(4); /* j1, p1, j0, p0 */

	if (doSort)
		sort_cr(res, cl);
	if (cl[1] == 's' && (keep == -1 || keep == 1)) {
		/* defined in ./sparse.c : */
		SEXP sparse_force_symmetric(SEXP, const char *, char);
		res = sparse_force_symmetric(res, cl_, (keep == 1) ? 'U' : 'L');
	}
	UNPROTECT(2); /* x, res */
	return res;
}

static SEXP diagonalMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
                                          const char *cl)
{
	SUB2_START;

	int nonunit = 1, keep = 0;
	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	UNPROTECT(1); /* diag */

	char cl_[] = ".gCMatrix";
	cl_[0] = cl[0];
	if (!(mi || mj) && ni == nj) {
		keep = keep_di(pi, pj, ni, nonunit, 0, m);
		if (keep) {
			cl_[1] = 'd';
			cl_[2] = 'i';
		}
	}
	SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl_));

	PROTECT(dim = GET_SLOT(res, Matrix_DimSym));
	pdim = INTEGER(dim);
	pdim[0] = (int) ni;
	pdim[1] = (int) nj;
	UNPROTECT(1); /* dim */

	if (keep > 1) {

		PROTECT(diag = GET_SLOT(res, Matrix_diagSym));
		SEXP diag_ = PROTECT(mkChar("U"));
		SET_STRING_ELT(diag, 0, diag_);
		UNPROTECT(2); /* diag_, diag */

	} else if (keep) {

		SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
			x1 = PROTECT(allocVector(TYPEOF(x0), ni));
		int j_;

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			while (ni--) \
				*(px1++) = \
					(*(pi++) != (j_ = *(pj++))) \
					? _ZERO_ \
					: ((nonunit) ? px0[j_ - 1] : _ONE_); \
		} while (0)

		SUB2_CASES(SUB2);

#undef SUB2

		SET_SLOT(res, Matrix_xSym, x1);
		UNPROTECT(2); /* x0, x1 */

	} else {

		SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1));
		int *pp1 = INTEGER(p1), j_;
		*(pp1++) = 0;

		for (kj = 0; kj < nj; ++kj) {
			pp1[kj] = 0;
			j_ = (mj) ? kj : pj[kj] - 1;
			if (mi) {
				for (ki = 0; ki < ni; ++ki)
					if (ki == j_)
						++pp1[kj];
			} else {
				for (ki = 0; ki < ni; ++ki)
					if (pi[ki] - 1 == j_)
						++pp1[kj];
			}
			if (pp1[kj] > INT_MAX - pp1[kj - 1])
				error(_("%s too dense for %s; would have more than %s nonzero entries"), \
				      "x[i,j]", "[CR]sparseMatrix", "2^31-1"); \

			pp1[kj] += pp1[kj-1];
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[nj - 1])),
			x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
			x1 = PROTECT(allocVector(TYPEOF(x0), pp1[nj - 1]));
		int *pi1 = INTEGER(i1);

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			for (kj = 0; kj < nj; ++kj) { \
				j_ = (mj) ? kj : pj[kj] - 1; \
				for (ki = 0; ki < ni; ++ki) { \
					if (((mi) ? ki : pi[ki] - 1) == j_) { \
						*(pi1++) = ki; \
						*(px1++) = (nonunit) ? px0[j_] : _ONE_; \
					} \
				} \
			} \
		} while (0)

		SUB2_CASES(SUB2);

#undef SUB2

		SET_SLOT(res, Matrix_pSym, p1);
		SET_SLOT(res, Matrix_iSym, i1);
		SET_SLOT(res, Matrix_xSym, x1);
		UNPROTECT(4); /* x1, x0, i1, p1 */

	}

	UNPROTECT(1); /* res */
	return res;
}

static SEXP indMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
                                     const char *cl)
{
	PROTECT_INDEX pidA;
	PROTECT_WITH_INDEX(x, &pidA);

	PROTECT_INDEX pidB;
	SEXP perm0 = GET_SLOT(x, Matrix_permSym);
	int *pperm0 = INTEGER(perm0);
	PROTECT_WITH_INDEX(perm0, &pidB);

	SEXP margin = PROTECT(GET_SLOT(x, Matrix_marginSym));
	int mg = INTEGER(margin)[0] - 1;

	SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[mg], n = pdim[!mg];
	UNPROTECT(1); /* dim */

	if (mg) {
		SEXP i_tmp = i;
		i = j;
		j = i_tmp;
	}

	int ki, kj,
		mi = i == R_NilValue,
		mj = j == R_NilValue,
		ni = (mi) ? m : LENGTH(i),
		nj = (mj) ? n : LENGTH(j),
		*pi = (mi) ? NULL : INTEGER(i),
		*pj = (mj) ? NULL : INTEGER(j),
		isP = cl[0] == 'p';

	if (!mi) {
		isP = isP && ni == m;
		if (isP) {
			char *work;
			Matrix_Calloc(work, m, char);
			--work; /* now 1-indexed */
			for (ki = 0; ki < ni; ++ki) {
				if (work[pi[ki]]) {
					isP = 0;
					break;
				}
				work[pi[ki]] = 1;
			}
			++work; /* now 0-indexed */
			Matrix_Free(work, m);
		}

		x = NEW_OBJECT_OF_CLASS((isP) ? "pMatrix" : "indMatrix");
		REPROTECT(x, pidA);

		PROTECT(dim = GET_SLOT(x, Matrix_DimSym));
		pdim = INTEGER(dim);
		pdim[ mg] = ni;
		pdim[!mg] = n;
		UNPROTECT(1); /* dim */

		if (mg)
			SET_SLOT(x, Matrix_marginSym, margin);

		SEXP perm1 = PROTECT(allocVector(INTSXP, ni));
		int *pperm1 = INTEGER(perm1);
		--pperm0; /* now 1-indexed */
		for (ki = 0; ki < ni; ++ki)
			pperm1[ki] = pperm0[pi[ki]];
		SET_SLOT(x, Matrix_permSym, perm1);
		UNPROTECT(1); /* perm1 */

		perm0 = perm1;
		pperm0 = pperm1;
		REPROTECT(perm0, pidB);

		m = ni;
	}

	if (!mj) {
		isP = isP && nj == n;
		if (isP) {
			char *work;
			Matrix_Calloc(work, nj, char);
			--work; /* now 1-indexed */
			for (kj = 0; kj < nj; ++kj) {
				if (work[pj[kj]]) {
					isP = 0;
					break;
				}
				work[pj[kj]] = 1;
			}
			++work; /* now 0-indexed */
			Matrix_Free(work, nj);
		}

		x = NEW_OBJECT_OF_CLASS((isP)
		                        ? "pMatrix"
		                        : ((!mg) ? "ngCMatrix" : "ngRMatrix"));
		REPROTECT(x, pidA);

		PROTECT(dim = GET_SLOT(x, Matrix_DimSym));
		pdim = INTEGER(dim);
		pdim[ mg] = m;
		pdim[!mg] = nj;
		UNPROTECT(1); /* dim */

		if (isP) {
			SEXP perm1 = PROTECT(allocVector(INTSXP, nj));
			int *pperm1 = INTEGER(perm1), *work;
			Matrix_Calloc(work, nj, int);
			--work; /* now 1-indexed */
			for (kj = 0; kj < nj; ++kj)
				work[pj[kj]] = kj + 1;
			for (kj = 0; kj < nj; ++kj)
				pperm1[kj] = work[pperm0[kj]];
			++work; /* now 0-indexed */
			Matrix_Free(work, nj);
			SET_SLOT(x, Matrix_permSym, perm1);
			UNPROTECT(1); /* perm1 */
		} else {
			int *workA, *workB, *workC;
			size_t lwork = (size_t) n + n + m;
			Matrix_Calloc(workA, lwork, int);
			workB = workA + n;
			workC = workB + n;
			--workA; /* now 1-indexed */
			--workB; /* now 1-indexed */

			SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1));
			int *pp1 = INTEGER(p1), k, kend;

			/* 1. Compute old column counts in 'workA' */
			for (k = 0; k < m; ++k)
				++workA[pperm0[k]];

			/* 2. Compute new column pointers in 'pp1' */
			*(pp1++) = 0;
			for (kj = 0; kj < nj; ++kj) {
				pp1[kj] = workA[pj[kj]];
				if (pp1[kj] > INT_MAX - pp1[kj - 1])
					error(_("%s too dense for %s; would have more than %s nonzero entries"), \
					      "x[i,j]", "[CR]sparseMatrix", "2^31-1"); \

				pp1[kj] += pp1[kj - 1];
			}

			/* 3. Compute old column pointers in 'workB' and copy to 'workA' */
			workB[1] = 0;
			for (k = 1; k < n; ++k) {
				workB[k + 1] = workB[k] + workA[k];
				workA[k] = workB[k];
			}
			workA[n] = workB[n];

			/* 4. Sort old row indices into 'workC' */
			for (k = 0; k < m; ++k)
				workC[workA[pperm0[k]]++] = k;

			SEXP i1 = PROTECT(allocVector(INTSXP, pp1[nj - 1]));
			int *pi1 = INTEGER(i1), pos;

			/* 5. Copy row indices from 'workC' to 'pi1' */
			k = 0;
			for (kj = 0; kj < nj; ++kj) {
				kend = pp1[kj];
				pos = workB[pj[kj]];
				while (k < kend)
					pi1[k++] = workC[pos++];
			}

			++workA; /* now 0-indexed */
			Matrix_Free(workA, lwork);
			SET_SLOT(x, Matrix_pSym, p1);
			SET_SLOT(x, (!mg) ? Matrix_iSym : Matrix_jSym, i1);
			UNPROTECT(2); /* i1, p1 */
		}

		n = nj;
	}

#undef SUB2_CASES
#undef SUB2_START
#undef SUB2_START_EXTRA

	UNPROTECT(3); /* margin, perm0, x */
	return x;
}

/* x[i,j,drop=FALSE] with 'i' and 'j' of type "integer" and length
   not exceeding 2^31-1 {'i' in 1:m or NA, 'j' in 1:n or NA} ...
   but _not_ handling 'Dimnames'
*/
SEXP R_subscript_2ary(SEXP x, SEXP i, SEXP j)
{
	if (i == R_NilValue && j == R_NilValue)
		return x;

	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	ivalid += VALID_NONVIRTUAL_SHIFT(ivalid, 0);
	const char *cl = valid[ivalid];
	validObject(x, cl);

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
		return unpackedMatrix_subscript_2ary(x, i, j, cl);
	case 'p':
		return   packedMatrix_subscript_2ary(x, i, j, cl);
	default:
		break;
	}

	static SEXP anyNA = NULL;
	if (!anyNA)
		anyNA = install("anyNA");
	SEXP call = PROTECT(lang2(anyNA, R_NilValue)), value;

#define ERROR_IF_ANYNA(_I_) \
	do { \
		if ((_I_) != R_NilValue) { \
			SETCADR(call, _I_); \
			PROTECT(value = eval(call, R_BaseEnv)); \
			if (asLogical(value)) \
				error(_("NA subscripts in %s not supported for '%s' inheriting from %s"), \
					  "x[i,j]", "x", "sparseMatrix"); \
			UNPROTECT(1); \
		} \
	} while (0)

	ERROR_IF_ANYNA(i);
	ERROR_IF_ANYNA(j);

#undef ERROR_IF_ANYNA

	UNPROTECT(1);

	switch (cl[2]) {
	case 'C':
		return  CsparseMatrix_subscript_2ary(x, i, j, cl);
	case 'R':
		return  RsparseMatrix_subscript_2ary(x, i, j, cl);
	case 'T':
	{
		char cl_[] = "..CMatrix";
		cl_[0] = cl[0];
		cl_[1] = cl[1];

		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		SEXP sparse_as_Tsparse(SEXP, const char *);

		x = sparse_as_Csparse(x, cl);
		PROTECT(x);
		x = CsparseMatrix_subscript_2ary(x, i, j, cl_);
		UNPROTECT(1);
		PROTECT(x);
		x = sparse_as_Tsparse(x, valid[R_check_class_etc(x, valid)]);
		UNPROTECT(1);
		return x;
	}
	case 'i':
		return diagonalMatrix_subscript_2ary(x, i, j, cl);
	default:
		return      indMatrix_subscript_2ary(x, i, j, cl);
	}
}
