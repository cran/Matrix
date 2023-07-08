#include "sparse.h"

SEXP sparse_as_dense(SEXP from, int packed)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_as_dense");
	const char *clf = valid[ivalid];
	packed = packed && clf[1] != 'g';

	char clt[] = "...Matrix";
	clt[0] = clf[0];
	clt[1] = clf[1];
	clt[2] = (clf[1] == 'g')
		? 'e' : ((packed) ? 'p' : ((clf[1] == 't') ? 'r' : 'y'));
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
	++nprotect;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	double mn = (double) m * n, dnx = (packed) ? 0.5 * (mn + n) : mn,
	dsize = dnx * kind2size(clt[0]);
	if (dnx > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
	if (packed && clf[2] != 'C' && mn > R_XLEN_T_MAX)
		error(_("coercing n-by-n [RT]sparseMatrix to packedMatrix "
		        "is not supported for n*n exceeding R_XLEN_T_MAX"));
	if (dsize > 0x1.0p+30 /* 1 GiB */)
		warning(_("sparse->dense coercion: "
		          "allocating vector of size %0.1f GiB"),
		        0x1.0p-30 * dsize);
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = '\0', di = '\0';
	if (clf[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		SET_SLOT(to, Matrix_uploSym, uplo);
		ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */
		if (clf[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			SET_SLOT(to, Matrix_diagSym, diag);
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	}

	SEXPTYPE tx = kind2type(clf[0]);
	R_xlen_t nx = (R_xlen_t) dnx;
	SEXP x0 = NULL, x1 = PROTECT(allocVector(tx, nx));
	++nprotect;
	if (clf[0] != 'n') {
		PROTECT(x0 = GET_SLOT(from, Matrix_xSym));
		++nprotect;
	}

	/* It remains to fill 'x' ... */

	if (clf[2] == 'C') {

		SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i = PROTECT(GET_SLOT(from, Matrix_iSym));
		nprotect += 2;
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
		++pp;

#define SAD_C(_VAL_) \
		do { \
			if (!packed) { \
				for (j = 0, k = 0; j < n; ++j, px1 += m) { \
					kend = pp[j]; \
					while (k < kend) { \
						px1[*pi] = _VAL_; \
						++pi; ++k; \
					} \
				} \
			} else if (ul == 'U') { \
				for (j = 0, k = 0; j < n; px1 += (++j)) { \
					kend = pp[j]; \
					while (k < kend) { \
						px1[*pi] = _VAL_; \
						++pi; ++k; \
					} \
				} \
			} else { \
				for (j = 0, k = 0; j < n; px1 += n-(j++)) { \
					kend = pp[j]; \
					while (k < kend) { \
						px1[*pi - j] = _VAL_; \
						++pi; ++k; \
					} \
				} \
			} \
		} while (0)

#define SAD_C_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			Matrix_memset(px1, 0, nx, sizeof(_CTYPE_)); \
			SAD_C(*(px0++)); \
		} while (0)

		if (clf[0] == 'n') {
			int *px1 = LOGICAL(x1);
			Matrix_memset(px1, 0, nx, sizeof(int));
			SAD_C(1);
		} else {
			SPARSE_CASES(tx, SAD_C_X);
		}

#undef SAD_C_X
#undef SAD_C

	} else if (clf[2] == 'R') {

		SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym)),
			j = PROTECT(GET_SLOT(from, Matrix_jSym));
		nprotect += 2;
		int *pp = INTEGER(p), *pj = INTEGER(j), i, k, kend;
		++pp;

#define SAD_R(_VAL_) \
		do { \
			if (!packed) { \
				R_xlen_t m1 = (R_xlen_t) m; \
				for (i = 0, k = 0; i < m; ++i, ++px1) { \
					kend = pp[i]; \
					while (k < kend) { \
						px1[m1 * *pj] = _VAL_; \
						++pj; ++k; \
					} \
				} \
			} else if (ul == 'U') { \
				for (i = 0, k = 0; i < n; ++i) { \
					kend = pp[i]; \
					while (k < kend) { \
						px1[PM_AR21_UP(i, *pj)] = _VAL_; \
						++pj; ++k; \
					} \
				} \
			} else { \
				R_xlen_t n2 = (R_xlen_t) n * 2; \
				for (i = 0, k = 0; i < n; ++i) { \
					kend = pp[i]; \
					while (k < kend) { \
						px1[PM_AR21_LO(i, *pj, n2)] = _VAL_; \
						++pj; ++k; \
					} \
				} \
			} \
		} while (0)

#define SAD_R_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			Matrix_memset(px1, 0, nx, sizeof(_CTYPE_)); \
			SAD_R(*(px0++)); \
		} while (0)

			if (clf[0] == 'n') {
				int *px1 = LOGICAL(x1);
				Matrix_memset(px1, 0, nx, sizeof(int));
				SAD_R(1);
			} else {
				SPARSE_CASES(tx, SAD_R_X);
			}

#undef SAD_R_X
#undef SAD_R

	} else { /* clf[2] == 'T' */

		SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j = PROTECT(GET_SLOT(from, Matrix_jSym));
		nprotect += 2;
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t index, k, nnz = XLENGTH(i);

#define SAD_T(_DO_FOR_) \
		do { \
			if (!packed) { \
				R_xlen_t m1 = (R_xlen_t) m; \
				_DO_FOR_(m1 * *pj + *pi); \
			} else if (ul == 'U') { \
				_DO_FOR_(PM_AR21_UP(*pi, *pj)); \
			} else { \
				R_xlen_t n2 = (R_xlen_t) n * 2; \
				_DO_FOR_(PM_AR21_LO(*pi, *pj, n2)); \
			} \
		} while (0)

		switch (clf[0]) {
		case 'n':
		{
			int *px1 = LOGICAL(x1);
			Matrix_memset(px1, 0, nx, sizeof(int));

#define SAD_T_N(_INDEX_) \
			do { \
				for (k = 0; k < nnz; ++k, ++pi, ++pj) \
					px1[_INDEX_] = 1; \
			} while (0)

			SAD_T(SAD_T_N);
			break;
		}
		case 'l':
		{
			int *px0 = LOGICAL(x0), *px1 = LOGICAL(x1);
			Matrix_memset(px1, 0, nx, sizeof(int));

#define SAD_T_L(_INDEX_) \
			do { \
				for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0) { \
					if (*px0 != 0) { \
						index = _INDEX_; \
						if (*px0 != NA_LOGICAL) \
							px1[index] = 1; \
						else if (px1[index] == 0) \
							px1[index] = NA_LOGICAL; \
					} \
				} \
			} while (0)

			SAD_T(SAD_T_L);
			break;
		}
		case 'i':
		{
			int *px0 = INTEGER(x0), *px1 = INTEGER(x1);
			Matrix_memset(px1, 0, nx, sizeof(int));

/* FIXME: not detecting integer overflow here */
#define SAD_T_I(_INDEX_) \
			do { \
				for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0) \
					px1[_INDEX_] += *px0; \
			} while (0)

			SAD_T(SAD_T_I);
			break;
		}
		case 'd':
		{
			double *px0 = REAL(x0), *px1 = REAL(x1);
			Matrix_memset(px1, 0, nx, sizeof(double));

#define SAD_T_D(_INDEX_) SAD_T_I(_INDEX_)

			SAD_T(SAD_T_D);
			break;
		}
		case 'z':
		{
			Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
			Matrix_memset(px1, 0, nx, sizeof(Rcomplex));
			R_xlen_t index;

#define SAD_T_Z(_INDEX_) \
			do { \
				for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0) { \
					index = _INDEX_; \
					px1[index].r += (*px0).r; \
					px1[index].i += (*px0).i; \
				} \
			} while (0)

			SAD_T(SAD_T_Z);
			break;
		}
		default:
			ERROR_INVALID_TYPE("'x' slot", tx, "sparse_as_dense");
			break;
		}

#undef SAD_T_N
#undef SAD_T_L
#undef SAD_T_I
#undef SAD_T_D
#undef SAD_T_Z
#undef SAD_T

	}

	if (di != '\0' && di != 'N') {

		int j;

#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px1 = _PTR_(x1); \
			if (!packed) { \
				R_xlen_t n1a = (R_xlen_t) n + 1; \
				for (j = 0; j < n; ++j, px1 += n1a) \
					*px1 = _ONE_; \
			} else if (ul == 'U') { \
				for (j = 0; j < n; px1 += (++j)+1) \
					*px1 = _ONE_; \
			} else { \
				for (j = 0; j < n; px1 += n-(j++)) \
					*px1 = _ONE_; \
			} \
		} while (0)

		SPARSE_CASES(tx, SET1);

#undef SET1

	}

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(nprotect);
	return to;
}

/* as(<[CRT]sparseMatrix>,       "denseMatrix") */
/* as(<[CRT]sparseMatrix>, "(un)?packedMatrix") */
SEXP R_sparse_as_dense(SEXP from, SEXP packed)
{
	return sparse_as_dense(from, asLogical(packed));
}

/* as(<[CRT]sparseMatrix>, "matrix") */
SEXP R_sparse_as_matrix(SEXP from)
{
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from = sparse_as_dense(from, 0), &pid);
	REPROTECT(from = dense_as_general(from, '.', -1, 0), pid); /* in-place */
	SEXP to = PROTECT(GET_SLOT(from, Matrix_xSym)),
	dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	setAttrib(to, R_DimSymbol, dim);
	if (!DimNames_is_trivial(dimnames))
	setAttrib(to, R_DimNamesSymbol, dimnames);
	UNPROTECT(4); /* dimnames, dim, to, from */
	return to;
}

/* as(<[CRT]sparseMatrix>, "vector") */
SEXP R_sparse_as_vector(SEXP from)
{
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from = sparse_as_dense(from, 0), &pid);
	REPROTECT(from = dense_as_general(from, '.', -1, 0), pid); /* in-place */
	from = GET_SLOT(from, Matrix_xSym);
	UNPROTECT(1); /* from */
	return from;
}

SEXP sparse_as_kind(SEXP from, char kind, int drop0)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_as_kind");
	const char *clf = valid[ivalid];

	if (kind == '.')
		kind = clf[0];
	SEXPTYPE tx = kind2type(kind); /* validating before doing more */

	PROTECT_INDEX pidA;
	PROTECT_WITH_INDEX(from, &pidA);
	if (drop0 && clf[0] != 'n')
		REPROTECT(from = R_sparse_drop0(from), pidA);
	if (kind == clf[0]) {
		UNPROTECT(1); /* from */
		return from;
	}
	if (clf[2] == 'T' && (clf[0] == 'n' || clf[0] == 'l') &&
	    kind != 'n' && kind != 'l')
		REPROTECT(from = Tsparse_aggregate(from), pidA);

	char clt[] = "...Matrix";
	clt[0] = kind;
	clt[1] = clf[1];
	clt[2] = clf[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	R_xlen_t nx = -1;
	if (clf[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (clf[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}
	if (clf[2] != 'T') {
		SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
		SET_SLOT(to, Matrix_pSym, p);
		if (nx < 1)
			nx = INTEGER(p)[XLENGTH(p) - 1];
		UNPROTECT(1); /* p */
	}
	if (clf[2] != 'R') {
		SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
		SET_SLOT(to, Matrix_iSym, i);
		if (nx < 1)
			nx = XLENGTH(i);
		UNPROTECT(1); /* i */
	}
	if (clf[2] != 'C') {
		SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
		SET_SLOT(to, Matrix_jSym, j);
		UNPROTECT(1); /* j */
	}
	if (clf[0] != 'n' && kind != 'n') {
		SEXP x;
		PROTECT_INDEX pidB;
		PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pidB);
		REPROTECT(x = coerceVector(x, tx), pidB);
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	} else if (clf[0] == 'n') {
		SEXP x = PROTECT(allocVector(tx, nx));

#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px = _PTR_(x); \
			while (nx--) \
				*(px++) = _ONE_; \
		} while (0)

		SPARSE_CASES(tx, SET1);

#undef SET1

		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	}

	UNPROTECT(2); /* to, from */
	return to;
}

/* as(<[CRT]sparseMatrix>, "[nlidz]Matrix") */
SEXP R_sparse_as_kind(SEXP from, SEXP kind, SEXP drop0)
{
	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = *CHAR(kind)) == '\0')
		error(_("invalid 'kind' to 'R_sparse_as_kind()'"));

	return sparse_as_kind(from, kind_, asLogical(drop0));
}

/* as(<[CRT]sparseMatrix>, "generalMatrix") */
SEXP R_sparse_as_general(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_as_general");
	const char *clf = valid[ivalid];
	if (clf[1] == 'g')
		return from;

	char clt[] = ".g.Matrix";
	clt[0] = clf[0];
	clt[2] = clf[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
	++nprotect;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (clf[1] != 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (clf[1] != 's') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */

		if (di == 'N') {
			if (clf[2] != 'T') {
				SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
				SET_SLOT(to, Matrix_pSym, p);
				UNPROTECT(1); /* p */
			}
			if (clf[2] != 'R') {
				SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
				SET_SLOT(to, Matrix_iSym, i);
				UNPROTECT(1); /* i */
			}
			if (clf[2] != 'C') {
				SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
				SET_SLOT(to, Matrix_jSym, j);
				UNPROTECT(1); /* j */
			}
			if (clf[0] != 'n') {
				SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
				SET_SLOT(to, Matrix_xSym, x);
				UNPROTECT(1); /* x */
			}
			UNPROTECT(nprotect);
			return to;
		}
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorSym, factors);
		UNPROTECT(1); /* factors */
	}

	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */

	if (clf[2] != 'T') {

		SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
			i0 = PROTECT(GET_SLOT(from, iSym));
		nprotect += 3;
		int j, k, kend, nnz1,
			*pp0 = INTEGER(p0), *pp1 = INTEGER(p1), *pi0 = INTEGER(i0);
		pp0++;
		*(pp1++) = 0;

		if (clf[1] == 's') {
			Matrix_memset(pp1, 0, n, sizeof(int));
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] != j)
						++pp1[pi0[k]];
					++k;
				}
			}
			for (j = 1; j < n; ++j)
				pp1[j] += pp1[j-1];
			for (j = 0; j < n; ++j)
				pp1[j] += pp0[j];
		} else {
			/* FIXME: not detecting integer overflow here */
			for (j = 0; j < n; ++j)
				pp1[j] = pp0[j] + j + 1;
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1 = pp1[n-1]));
		++nprotect;
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to, iSym, i1);

		SEXP x0 = NULL, x1 = NULL;
		SEXPTYPE tx = NILSXP;
		if (clf[0] != 'n') {
			PROTECT(x0 = GET_SLOT(from, Matrix_xSym));
			PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
			nprotect += 2;
			SET_SLOT(to, Matrix_xSym, x1);
		}

		if (clf[1] == 's') {

			int *pp1_;
			Matrix_Calloc(pp1_, n, int);
			Matrix_memcpy(pp1_, pp1 - 1, n, sizeof(int));

#define SAG_CR(_XASSIGN_IJ_, _XASSIGN_JI_) \
			do { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						pi1[pp1_[j]] = pi0[k]; \
						_XASSIGN_IJ_; /* px1[pp1_[j]] = px0[k]; */ \
						++pp1_[j]; \
						if (pi0[k] != j) { \
							pi1[pp1_[pi0[k]]] = j; \
							_XASSIGN_JI_; /* px1[pp1_[pi0[k]]] = px0[k]; */ \
							++pp1_[pi0[k]]; \
						} \
						++k; \
					} \
				} \
			} while (0)

#define SAG_CR_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
				SAG_CR(px1[pp1_[j]] = px0[k], \
				       px1[pp1_[pi0[k]]] = px0[k]); \
			} while (0)

				if (clf[0] == 'n')
					SAG_CR(, );
				else
					SPARSE_CASES(tx, SAG_CR_X);

#undef SAG_CR_X
#undef SAG_CR

				Matrix_Free(pp1_, n);

		} else {

#define SAG_CR(_XASSIGN_, _XASSIGN_ONE_) \
			do { \
				if (ul == ((clf[2] == 'C') ? 'U' : 'L')) { \
					for (j = 0, k = 0; j < n; ++j) { \
						kend = pp0[j]; \
						while (k < kend) { \
							*(pi1++) = *(pi0++); \
							_XASSIGN_; /* *(px1++) = *(px0++); */ \
							++k; \
						} \
						*(pi1++) = j; \
						_XASSIGN_ONE_; /* *(px1++) = _ONE_; */ \
					} \
				} else { \
					for (j = 0, k = 0; j < n; ++j) { \
						kend = pp0[j]; \
						*(pi1++) = j; \
						_XASSIGN_ONE_; /* *(px1++) = _ONE_; */ \
						while (k < kend) { \
							*(pi1++) = *(pi0++); \
							_XASSIGN_; /* *(px1++) = *(px0++); */ \
							++k; \
						} \
					} \
				} \
			} while (0)

#define SAG_CR_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
				SAG_CR(*(px1++) = *(px0++), *(px1++) = _ONE_); \
			} while (0)

			if (clf[0] == 'n')
				SAG_CR(, );
			else
				SPARSE_CASES(tx, SAG_CR_X);

#undef SAG_CR_X
#undef SAG_CR

		}

	} else { /* clf[2] == 'T' */

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		nprotect += 2;
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = nnz0;

		if (clf[1] == 's') {
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] == pj0[k])
					--nnz1;
			nnz1 += nnz0;
		} else {
			nnz1 += n;
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		nprotect += 2;
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		Matrix_memcpy(pi1, pi0, nnz0, sizeof(int));
		Matrix_memcpy(pj1, pj0, nnz0, sizeof(int));
		pi1 += nnz0;
		pj1 += nnz0;

		SEXP x0 = NULL, x1 = NULL;
		SEXPTYPE tx = NILSXP;
		if (clf[0] != 'n') {
			PROTECT(x0 = GET_SLOT(from, Matrix_xSym));
			PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
			nprotect += 2;
			SET_SLOT(to, Matrix_xSym, x1);
		}

		if (clf[1] == 's') {

#define SAG_T(_XASSIGN_, _XINCR_) \
			do { \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 != *pj0) { \
						*(pi1++) = *pj0; \
						*(pj1++) = *pi0; \
						_XASSIGN_; /* *(px1++) = *px0; */ \
					} \
					++pi0; ++pj0; _XINCR_; /* ++px0; */ \
				} \
			} while (0)

#define SAG_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
				Matrix_memcpy(px1, px0, nnz0, sizeof(_CTYPE_)); \
				px1 += nnz0; \
				SAG_T(*(px1++) = *px0, ++px0); \
			} while (0)

			if (clf[0] == 'n')
				SAG_T(, );
			else
				SPARSE_CASES(tx, SAG_T_X);

#undef SAG_T_X
#undef SAG_T

		} else {

#define SAG_T(_XASSIGN_ONE_) \
			do { \
				int j; \
				for (j = 0; j < n; ++j) { \
					*(pi1++) = *(pj1++) = j; \
					_XASSIGN_ONE_; /* *(px1++) = _ONE_; */ \
				} \
			} while (0)

#define SAG_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
				Matrix_memcpy(px1, px0, nnz0, sizeof(_CTYPE_)); \
				px1 += nnz0; \
				SAG_T(*(px1++) = _ONE_); \
			} while (0)

			if (clf[0] == 'n')
				SAG_T();
			else
				SPARSE_CASES(tx, SAG_T_X);

#undef SAG_T_X
#undef SAG_T

		}

	}

	UNPROTECT(nprotect);
	return to;
}

/* as(<diagonalMatrix>, "[nlidz][gts][CRT]Matrix") */
SEXP R_diagonal_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP drop0)
{
	static const char *valid[] = { VALID_DIAGONAL, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_diagonal_as_sparse");
	const char *clf = valid[ivalid];

	const char *zzz;
	if (TYPEOF(code) != STRSXP || LENGTH(code) < 1 ||
	    (code = STRING_ELT(code, 0)) == NA_STRING ||
	    (zzz = CHAR(code))[0] == '\0' ||
	    (zzz[1] != 'g' && zzz[1] != 't' && zzz[1] != 's') ||
	    (zzz[2] != 'C' && zzz[2] != 'R' && zzz[2] != 'T'))
		error(_("invalid 'code' to 'R_diagonal_as_sparse()'"));
	SEXPTYPE tx = kind2type((zzz[0] == '.' || zzz[0] == 'n') ? clf[0] : zzz[0]);

	char clt[] = "...Matrix";
	clt[0] = (zzz[0] == '.') ? clf[0] : zzz[0];
	clt[1] = zzz[1];
	clt[2] = zzz[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
	++nprotect;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));;
	if (clt[1] != 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (clt[1] != 'g') {
		char ul;
		if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
			error(_("invalid 'uplo' to 'R_diagonal_as_sparse()'"));
		if (ul != 'U') {
			PROTECT(uplo = mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
	}

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (clt[1] == 't' && di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	SEXP p = NULL, i = NULL, x = NULL;
	int k, *pp = NULL, *pi = NULL;
	R_xlen_t n1a;
	PROTECT_INDEX pid;

	if (clt[2] != 'T') {
		n1a = (R_xlen_t) n + 1;
		PROTECT(p = allocVector(INTSXP, n1a));
		++nprotect;
		SET_SLOT(to, Matrix_pSym, p);
		pp = INTEGER(p);
	}

	if (di != 'N') { /* unit diagonal */

		if (clt[2] != 'T') {
			if (clt[1] == 't')
				Matrix_memset(pp, 0, n1a, sizeof(int));
			else
				for (k = 0; k <= n; ++k)
					*(pp++) = k;
		}
		if (clt[1] == 't') {
			UNPROTECT(nprotect);
			return to;
		}
		PROTECT(i = allocVector(INTSXP, n));
		++nprotect;
		pi = INTEGER(i);
		if (clt[0] == 'n') {
			for (k = 0; k < n; ++k)
				*(pi++) = k;
		} else {
			PROTECT(x = allocVector(tx, n));
			++nprotect;

#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
			do { \
				_CTYPE_ *px = _PTR_(x); \
				for (k = 0; k < n; ++k) { \
					*(pi++) = k; \
					*(px++) = _ONE_; \
				} \
			} while (0)

			SPARSE_CASES(tx, SET1);

#undef SET1

		}

	} else { /* non-unit diagonal */

		PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);
		REPROTECT(x = coerceVector(x, tx), pid);
		++nprotect;

		if (clt[0] == 'n' || asLogical(drop0) != 0) { /* drop zeros (if any) */

			int nnz = 0;

#define DROP0_DIAGONAL(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
			do { \
				_CTYPE_ *px = _PTR_(x); \
				if (clt[2] != 'T') { \
					*(pp++) = 0; \
					for (k = 0; k < n; ++k) \
						*(pp++) = (_NZ_(px[k])) ? ++nnz : nnz; \
				} else { \
					for (k = 0; k < n; ++k) \
						if (_NZ_(px[k])) \
							++nnz; \
				} \
				if (nnz == 0) { \
					UNPROTECT(nprotect); \
					return to; \
				} \
				PROTECT(i = allocVector(INTSXP, nnz)); \
				++nprotect; \
				pi = INTEGER(i); \
				if (nnz == n) { \
					for (k = 0; k < n; ++k) \
						*(pi++) = k; \
				} else if (clt[0] == 'n') { \
					for (k = 0; k < n; ++k) \
						if (_NZ_(px[k])) \
							*(pi++) = k; \
				} else { \
					SEXP y = PROTECT(allocVector(tx, nnz)); \
					_CTYPE_ *py = _PTR_(y); \
					for (k = 0; k < n; ++k) { \
						if (_NZ_(px[k])) { \
							*(pi++) = k; \
							*(py++) = px[k]; \
						} \
					} \
					UNPROTECT(1); /* y */ \
					REPROTECT(x = y, pid); \
				} \
			} while (0)

			SPARSE_CASES(tx, DROP0_DIAGONAL);

#undef DROP0_DIAGONAL

		} else { /* _do not_ drop zeros */

			PROTECT(i = allocVector(INTSXP, n));
			++nprotect;
			pi = INTEGER(i);
			if (clt[2] != 'T') {
				*(pp++) = 0;
				for (k = 0; k < n; ++k)
					*(pi++) = *(pp++) = k;
			} else {
				for (k = 0; k < n; ++k)
					*(pi++) = k;
			}

		}

	}

	if (clt[2] != 'R')
		SET_SLOT(to, Matrix_iSym, i);
	if (clt[2] != 'C')
		SET_SLOT(to, Matrix_jSym, i);
	if (clt[0] != 'n')
		SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(nprotect);
	return to;
}

/* as(<diagonalMatrix>, ".(ge|tr|sy|tp|sp)Matrix") */
SEXP R_diagonal_as_dense(SEXP from, SEXP code, SEXP uplo)
{
	static const char *valid[] = { VALID_DIAGONAL, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_diagonal_as_dense");
	const char *clf = valid[ivalid];

	const char *zzz;
	if (TYPEOF(code) != STRSXP || LENGTH(code) < 1 ||
	    (code = STRING_ELT(code, 0)) == NA_STRING ||
	    (zzz = CHAR(code))[0] == '\0' ||
	    (zzz             )[1] == '\0' ||
	    !((zzz[1] == 'g' && (zzz[2] == 'e'                 )) ||
	      (zzz[1] == 't' && (zzz[2] == 'r' || zzz[2] == 'p')) ||
	      (zzz[1] == 's' && (zzz[2] == 'y' || zzz[2] == 'p'))))
		error(_("invalid 'code' to 'R_diagonal_as_dense()'"));
	SEXPTYPE tx = kind2type((zzz[0] == '.') ? clf[0] : zzz[0]);

	char clt[] = "...Matrix";
	clt[0] = (zzz[0] == '.') ? clf[0] : zzz[0];
	clt[1] = zzz[1];
	clt[2] = zzz[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	double nn = (double) n * n,
	nx = (clt[2] == 'p') ? 0.5 * (nn + n) : nn,
	size = nx * kind2size(clt[0]);
	if (nx > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
	if (size > 0x1.0p+30) /* 1 GiB */
		warning(_("sparse->dense coercion: "
		          "allocating vector of size %0.1f GiB"),
		        0x1.0p-30 * size);
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (clt[1] != 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	char ul = '\0';
	if (clt[1] != 'g') {
		if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
			error(_("invalid 'uplo' to 'R_diagonal_as_sparse()'"));
		if (ul != 'U') {
			PROTECT(uplo = mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
	}

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (clt[1] == 't' && di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	SEXP x1 = PROTECT(allocVector(tx, (R_xlen_t) nx));

#define DAD_COPY_DIAGONAL(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		Matrix_memset(_PTR_(x1), 0, (R_xlen_t) nx, sizeof(_CTYPE_)); \
		if (clt[1] != 't' || di == 'N') { \
			SEXP x0; \
			PROTECT_INDEX pid; \
			PROTECT_WITH_INDEX(x0 = GET_SLOT(from, Matrix_xSym), &pid); \
			REPROTECT(x0 = coerceVector(x0, tx), pid); \
			if (clt[2] != 'p') \
				_PREFIX_ ## dense_unpacked_copy_diagonal( \
					_PTR_(x1), _PTR_(x0), n, n,     ul /* unused */, di); \
			else \
				_PREFIX_ ## dense_packed_copy_diagonal( \
					_PTR_(x1), _PTR_(x0), n, n, ul, ul /* unused */, di); \
			UNPROTECT(1); /* x0 */ \
		} \
	} while (0)

	switch (tx) {
	case LGLSXP:
		DAD_COPY_DIAGONAL(i, int, LOGICAL);
		break;
	case INTSXP:
		DAD_COPY_DIAGONAL(i, int, INTEGER);
		break;
	case REALSXP:
		DAD_COPY_DIAGONAL(d, double, REAL);
		break;
	case CPLXSXP:
		DAD_COPY_DIAGONAL(z, Rcomplex, COMPLEX);
		break;
	default:
		break;
	}

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(2); /* x1, to */
	return to;
}

/* as(<diagonalMatrix>, "[lidz]Matrix") */
SEXP R_diagonal_as_kind(SEXP from, SEXP kind)
{
	static const char *valid[] = { VALID_DIAGONAL, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_diagonal_as_kind");
	const char *clf = valid[ivalid];

	char k;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (k = *CHAR(kind)) == '\0')
		error(_("invalid 'kind' to 'R_diagonal_as_kind()'"));
	if (k == '.' || k == clf[0])
		return from;
	if (k == 'n')
		error(_("class ndiMatrix is unimplemented"));
	SEXPTYPE tt = kind2type(k); /* validating before doing more */

	char clt[] = ".diMatrix";
	clt[0] = k;
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	SEXP x;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);
	REPROTECT(x = coerceVector(x, tt), pid);
	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(2); /* x, to */
	return to;
}

/* drop0(<[CRT]sparseMatrix>)
   TODO: support 'tol' argument, to be interpreted as modulus for zMatrix
*/
SEXP R_sparse_drop0(SEXP from)
{
	static const char *valid[] = { VALID_DSPARSE, VALID_LSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_drop0");
	const char *cl = valid[ivalid];

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), p0 = NULL;
	++nprotect;
	SEXPTYPE tx = TYPEOF(x0);
	int *pp0 = NULL;
	R_xlen_t n1a = 0, k, kend, nnz_ = 0, nnz0, nnz1 = 0;

	if (cl[2] != 'T') {
		PROTECT(p0 = GET_SLOT(from, Matrix_pSym));
		++nprotect;
		pp0 = INTEGER(p0);
		n1a = XLENGTH(p0);
		nnz0 = pp0[n1a - 1];
	} else
		nnz0 = XLENGTH(x0);

#define DROP0_START(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0); \
		while (nnz1 < nnz0 && _NZ_(*px0)) { \
			++nnz1; \
			++px0; \
		} \
		if (nnz1 == nnz0) { \
			UNPROTECT(nprotect); \
			return from; \
		} \
		nnz_ = nnz1; \
		for (k = nnz_; k < nnz0; ++k, ++px0) \
			if (_NZ_(*px0)) \
				++nnz1; \
	} while (0)

	SPARSE_CASES(tx, DROP0_START);

#undef DROP0_START

	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));
	++nprotect;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (cl[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (cl[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorSym, factors);
		UNPROTECT(1); /* factors */
	}

	/* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

	SEXP iSym = (cl[2] == 'R') ? Matrix_jSym : Matrix_iSym,
		i0 = PROTECT(GET_SLOT(from, iSym)),
		i1 = PROTECT(allocVector(INTSXP, nnz1)),
		x1 = PROTECT(allocVector(tx, nnz1));
	nprotect += 3;
	int *pi0 = INTEGER(i0),
		*pi1 = INTEGER(i1);

	if (cl[2] != 'T') {

		SEXP p1 = PROTECT(allocVector(INTSXP, n1a));
		++nprotect;
		int *pp1 = INTEGER(p1), j;
		n = (int) n1a - 1;

#define DROP0_END(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			Matrix_memcpy(pi1, pi0, nnz_, sizeof(int)); \
			Matrix_memcpy(px1, px0, nnz_, sizeof(_CTYPE_)); \
			j = 0; \
			while ((kend = pp0[j]) <= nnz_) \
				pp1[j++] = kend; \
			for (k = nnz_; k < kend; ++k) { \
				if (_NZ_(px0[k])) { \
					pi1[nnz_] = pi0[k]; \
					px1[nnz_] = px0[k]; \
					++nnz_; \
				} \
			} \
			pp1[j] = nnz_; \
			while (++j <= n) { \
				kend = pp0[j]; \
				while (k < kend) { \
					if (_NZ_(px0[k])) { \
						pi1[nnz_] = pi0[k]; \
						px1[nnz_] = px0[k]; \
						++nnz_; \
					} \
					++k; \
				} \
				pp1[j] = nnz_; \
			} \
		} while (0)

		SPARSE_CASES(tx, DROP0_END);

#undef DROP0_END

		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_xSym, x1);

	} else {

		SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		nprotect += 2;
		int *pj0 = INTEGER(j0),
			*pj1 = INTEGER(j1);

#define DROP0_END(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			Matrix_memcpy(pi1, pi0, nnz_, sizeof(int)); \
			Matrix_memcpy(pj1, pj0, nnz_, sizeof(int)); \
			Matrix_memcpy(px1, px0, nnz_, sizeof(_CTYPE_)); \
			for (k = nnz_; k < nnz0; ++k) { \
				if (_NZ_(px0[k])) { \
					pi1[nnz_] = pi0[k]; \
					pj1[nnz_] = pj0[k]; \
					px1[nnz_] = px0[k]; \
					++nnz_; \
				} \
			} \
		} while (0)

		SPARSE_CASES(tx, DROP0_END);

#undef DROP0_END

		SET_SLOT(to,        iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		SET_SLOT(to, Matrix_xSym, x1);

	}

	UNPROTECT(nprotect);
	return to;
}

/* band(<[CRT]sparseMatrix>, k1, k2), tri[ul](<[CRT]sparseMatrix>, k) */
/* NB: argument validation more or less copied from R_dense_band() */
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_band");
	const char *clf = valid[ivalid];

	SEXP dim;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(dim = GET_SLOT(from, Matrix_DimSym), &pid);
	++nprotect;
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], a, b;
	if (isNull(k1))
		a = (m > 0) ? 1-m : 0;
	else if ((a = asInteger(k1)) == NA_INTEGER || a < -m || a > n)
		error(_("'k1' must be an integer from -Dim[1] to Dim[2]"));
	if (isNull(k2))
		b = (n > 0) ? n-1 : 0;
	else if ((b = asInteger(k2)) == NA_INTEGER || b < -m || b > n)
		error(_("'k2' must be an integer from -Dim[1] to Dim[2]"));
	else if (b < a)
		error(_("'k1' must be less than or equal to 'k2'"));
	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
	if (a <= 1-m && b >= n-1 && (clf[1] == 't' || m != n || m > 1 || n > 1)) {
		UNPROTECT(nprotect); /* dim */
		return from;
	}

	char ulf = 'U', ult = 'U', di = 'N';
	if (clf[1] != 'g') {
		SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ulf = *CHAR(STRING_ELT(uplo_from, 0));
		UNPROTECT(1); /* uplo_from */
		if (clf[1] == 't') {
			/* Be fast if band contains entire triangle */
			if ((ulf == 'U') ? (a <= 0 && b >= n-1) : (b >= 0 && a <= 1-m)) {
				UNPROTECT(nprotect);
				return from;
			} else if (a <= 0 && b >= 0) {
				SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
				di = *CHAR(STRING_ELT(diag, 0));
				UNPROTECT(1); /* diag */
			}
		}
	}

	/* band(<R>, a, b) is equivalent to t(band(t(<R>), -b, -a)) ! */

	if (clf[2] == 'R') {
		int r;
		r = m; m =  n; n =  r;
		r = a; a = -b; b = -r;
		ulf = (ulf == 'U') ? 'L' : 'U';
		PROTECT(from = tCRsparse_as_RCsparse(from));
		++nprotect;
		if (m != n)
			REPROTECT(dim = GET_SLOT(from, Matrix_DimSym), pid);
	}

	int ge = 0, tr = 0, sy = 0;
	ge = m != n || !((tr = a >= 0 || b <= 0 || clf[1] == 't') ||
	                 (sy = a == -b && clf[1] == 's'));

	char clt[] = "...Matrix";
	clt[0] = clf[0];
	clt[1] = (ge) ? 'g' : ((tr) ? 't' : 's');
	clt[2] = (clf[2] == 'R') ? 'C' : clf[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
	++nprotect;

	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (tr && clf[1] == 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (!ge) {
		ult = (tr && clf[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ulf;
		if (ult != 'U') {
			SEXP uplo_to = PROTECT(mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo_to);
			UNPROTECT(1); /* uplo_to */
		}
		if (di != 'N') {
			SEXP diag = PROTECT(mkString("U"));
			SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}
	}

	/* It remains to set some subset of 'p', 'i', 'j', 'x' ... */

	SEXP p0 = NULL, p1 = NULL, i0 = NULL, j0 = NULL;
	int *pp0 = NULL, *pp1 = NULL, *pi0 = NULL, *pj0 = NULL, d, j;
	R_xlen_t k, kend, nnz0, nnz1;

	PROTECT(i0 = GET_SLOT(from, Matrix_iSym));
	++nprotect;
	pi0 = INTEGER(i0);

	/* Counting number of nonzero elements in band ... */

	nnz1 = 0;
	if (clf[2] != 'T') {

		PROTECT(p0 = GET_SLOT(from, Matrix_pSym));
		PROTECT(p1 = allocVector(INTSXP, (R_xlen_t) n + 1));
		nprotect += 2;

		pp0 = INTEGER(p0);
		pp1 = INTEGER(p1);

		nnz0 = pp0[n];
		pp1[0] = 0;

		++pp0;
		++pp1;

		if (!sy && clf[1] == 's') {
			Matrix_memset(pp1, 0, n, sizeof(int));
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((d = j - pi0[k]) >= a && d <= b)
						++pp1[j];
					if (d != 0 && -d >= a && -d <= b)
						++pp1[pi0[k]];
					++k;
				}
			}
			for (j = 0; j < n; ++j) {
				nnz1 += pp1[j];
				pp1[j] = nnz1;
			}
		} else {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((d = j - pi0[k]) >= a && d <= b)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
		}

		SET_SLOT(to, Matrix_pSym, p1);

	} else {

		PROTECT(j0 = GET_SLOT(from, Matrix_jSym));
		++nprotect;

		pj0 = INTEGER(j0);
		nnz0 = XLENGTH(j0);

		if (!sy && clf[1] == 's') {
			for (k = 0; k < nnz0; ++k) {
				if ((d = pj0[k] - pi0[k]) >= a && d <= b)
					++nnz1;
				if (d != 0 && -d >= a && -d <= b)
					++nnz1;
			}
		} else {
			for (k = 0; k < nnz0; ++k) {
				if ((d = pj0[k] - pi0[k]) >= a && d <= b)
					++nnz1;
			}
		}

	}

	if (nnz1 == nnz0 && (sy || clf[1] != 's')) {
		/* No need to allocate in this case: band has all nonzero elements */
		SET_SLOT(to, Matrix_iSym, i0);
		if (clf[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			SET_SLOT(to, Matrix_xSym, x0);
			UNPROTECT(1); /* x0 */
		}
		if (clf[2] == 'T')
			SET_SLOT(to, Matrix_jSym, j0);
		else if (clf[2] == 'R')
			to = tCRsparse_as_RCsparse(to);
		UNPROTECT(nprotect);
		return to;
	}

	/* Now allocating and filling out slots ... */

	SEXP i1 = NULL, j1 = NULL;
	int *pi1 = NULL, *pj1 = NULL;
	PROTECT(i1 = allocVector(INTSXP, nnz1));
	++nprotect;
	pi1 = INTEGER(i1);
	if (clf[2] == 'T') {
		PROTECT(j1 = allocVector(INTSXP, nnz1));
		++nprotect;
		pj1 = INTEGER(j1);
	}

#define SPARSE_BAND(_XASSIGN_, _XASSIGN_IJ_, _XASSIGN_JI_) \
	do { \
		if (clf[2] != 'T') { \
			if (!sy && clf[1] == 's') { \
				int *pp1_; \
				Matrix_Calloc(pp1_, n, int); \
				Matrix_memcpy(pp1_, pp1 - 1, n, sizeof(int)); \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((d = j - pi0[k]) >= a && d <= b) { \
							pi1[pp1_[j]] = pi0[k]; \
							_XASSIGN_IJ_; /* px1[pp1_[j]] = px0[k]; */ \
							++pp1_[j]; \
						} \
						if (d != 0 && -d >= a && -d <= b) { \
							pi1[pp1_[pi0[k]]] = j; \
							_XASSIGN_JI_; /* px1[pp1_[pi0[k]]] = px0[k]; */ \
							++pp1_[pi0[k]]; \
						} \
						++k; \
					} \
				} \
				Matrix_Free(pp1_, n); \
			} else { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((d = j - pi0[k]) >= a && d <= b) { \
							*(pi1++) = pi0[k]; \
							_XASSIGN_; /* *(px1++) = px0[k]; */ \
						} \
						++k; \
					} \
				} \
			} \
		} else { \
			if (!sy && clf[1] == 's') { \
				for (k = 0; k < nnz0; ++k) { \
					if ((d = pj0[k] - pi0[k]) >= a && d <= b) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						_XASSIGN_; /* *(px1++) = px0[k]; */ \
					} \
					if (d != 0 && -d >= a && -d <= b) { \
						*(pi1++) = pj0[k]; \
						*(pj1++) = pi0[k]; \
						_XASSIGN_; /* *(px1++) = px0[k]; */ \
					} \
				} \
			} else { \
				for (k = 0; k < nnz0; ++k) { \
					if ((d = pj0[k] - pi0[k]) >= a && d <= b) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						_XASSIGN_; /* *(px1++) = px0[k]; */ \
					} \
				} \
			} \
		} \
	} while (0)

#define SPARSE_BAND_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		SPARSE_BAND(*(px1++) = px0[k], \
		            px1[pp1_[j]] = px0[k], \
		            px1[pp1_[pi0[k]]] = px0[k]); \
	} while (0)

	if (clf[0] == 'n')
		SPARSE_BAND(, , );
	else {
		SEXPTYPE tx;
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
			x1 = PROTECT(allocVector(tx = TYPEOF(x0), nnz1));
		SPARSE_CASES(tx, SPARSE_BAND_X);
		SET_SLOT(to, Matrix_xSym, x1);
		UNPROTECT(2);
	}

#undef SPARSE_BAND_X
#undef SPARSE_BAND

	SET_SLOT(to, Matrix_iSym, i1);
	if (clf[2] == 'T')
		SET_SLOT(to, Matrix_jSym, j1);
	else {
		SET_SLOT(to, Matrix_pSym, p1);
		if (clf[2] == 'R')
			to = tCRsparse_as_RCsparse(to);
	}

	UNPROTECT(nprotect);
	return to;
}

/* diag(<[CRT]sparseMatrix>, names) */
SEXP R_sparse_diag_get(SEXP obj, SEXP nms)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, "R_sparse_diag_get");
	const char *cl = valid[ivalid];

	int do_nms = asLogical(nms);
	if (do_nms == NA_LOGICAL)
		error(_("'names' must be TRUE or FALSE"));

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	UNPROTECT(1); /* dim */

	SEXPTYPE type = kind2type(cl[0]);
	SEXP res = PROTECT(allocVector(type, r));
	++nprotect;

	char ul = '\0', di = '\0';
	if (cl[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		if (cl[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	}

	if (di != '\0' && di != 'N') {

		/* .t[CRT]Matrix with unit diagonal */

		int i;

#define DO_ONES(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *pres = _PTR_(res); \
			for (i = 0; i < r; ++i) \
				*(pres++) = _ONE_; \
		} while (0)

		SPARSE_CASES(type, DO_ONES);

#undef DO_ONES

	} else if (cl[2] != 'T') {

		/* ..[CR]Matrix with non-unit diagonal */

		SEXP iSym = (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj, iSym)),
			x = NULL;
		nprotect += 2;
		int j, k, kend, *pp = INTEGER(p) + 1, *pi = INTEGER(i);
		if (cl[0] != 'n') {
			PROTECT(x = GET_SLOT(obj, Matrix_xSym));
			++nprotect;
		}

#define DO_DIAG(_RES_, _VAL_GENERAL_, _VAL_TRAILING_, _VAL_LEADING_, _ZERO_) \
		do { \
			if (ul == '\0') { \
				/* .g[CR]Matrix */ \
				for (j = 0, k = 0; j < r; ++j) { \
					_RES_[j] = _ZERO_; \
					kend = pp[j]; \
					while (k < kend) { \
						if (pi[k] == j) { \
							_RES_[j] = _VAL_GENERAL_ /* px[k] */; \
							k = kend; \
							break; \
						} \
						++k; \
					} \
				} \
			} else if (ul == ((cl[2] == 'C') ? 'U' : 'L')) { \
				/* .[ts][CR]Matrix with "trailing" diagonal */ \
				for (j = 0, k = 0; j < r; ++j) { \
					kend = pp[j]; \
					_RES_[j] = (kend - k > 0 && pi[kend-1] == j \
					            ? _VAL_TRAILING_ /* px[kend-1] */ \
					            : _ZERO_); \
					k = kend; \
				} \
			} else { \
				/* .[ts][CR]Matrix with "leading" diagonal */ \
				for (j = 0, k = 0; j < r; ++j) { \
					kend = pp[j]; \
					_RES_[j] = (kend - k > 0 && pi[k] == j \
					            ? _VAL_LEADING_	/* px[k] */ \
					            : _ZERO_); \
					k = kend; \
				} \
			} \
		} while (0)

#define DO_DIAG_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px = _PTR_(x), *pres = _PTR_(res); \
			DO_DIAG(pres, px[k], px[kend-1], px[k], _ZERO_); \
		} while (0)

		if (cl[0] == 'n') {
			int *pres = LOGICAL(res);
			DO_DIAG(pres, 1, 1, 1, 0);
		} else {
			SPARSE_CASES(type, DO_DIAG_X);
		}

#undef DO_DIAG_X
#undef DO_DIAG

	} else {

		/* ..TMatrix with non-unit diagonal */

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym)),
			x = NULL;
		nprotect += 2;
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, nnz = XLENGTH(i);
		if (cl[0] != 'n') {
			PROTECT(x = GET_SLOT(obj, Matrix_xSym));
			++nprotect;
		}

		switch (cl[0]) {
		case 'n':
		{
			int *pres = LOGICAL(res);
			Matrix_memset(pres, 0, r, sizeof(int));
			for (k = 0; k < nnz; ++k, ++pi, ++pj)
				if (*pi == *pj)
					pres[*pi] = 1;
			break;
		}
		case 'l':
		{
			int *px = LOGICAL(x), *pres = LOGICAL(res);
			Matrix_memset(pres, 0, r, sizeof(int));
			for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px) {
				if (*pi == *pj && *px != 0) {
					if (*px != NA_LOGICAL)
						pres[*pi] = 1;
					else if (pres[*pi] == 0)
						pres[*pi] = NA_LOGICAL;
				}
			}
			break;
		}
		case 'i':
		{
			/* FIXME: not detecting integer overflow here */
			int *px = INTEGER(x), *pres = INTEGER(res);
			Matrix_memset(pres, 0, r, sizeof(int));
			for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px)
				if (*pi == *pj)
					pres[*pi] += *px;
			break;
		}
		case 'd':
		{
			double *px = REAL(x), *pres = REAL(res);
			Matrix_memset(pres, 0, r, sizeof(double));
			for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px)
				if (*pi == *pj)
					pres[*pi] += *px;
			break;
		}
		case 'z':
		{
			Rcomplex *px = COMPLEX(x), *pres = COMPLEX(res);
			Matrix_memset(pres, 0, r, sizeof(Rcomplex));
			for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px) {
				if (*pi == *pj) {
					pres[*pi].r += (*px).r;
					pres[*pi].i += (*px).i;
				}
			}
			break;
		}
		default:
			break;
		}

	}

	if (do_nms) {
		/* NB: The logic here must be adjusted once the validity method
	       for 'symmetricMatrix' enforces symmetric 'Dimnames' */
		SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
			rn = VECTOR_ELT(dn, 0),
			cn = VECTOR_ELT(dn, 1);
		if (isNull(cn)) {
			if (ul != '\0' && di == '\0' && !isNull(rn))
				setAttrib(res, R_NamesSymbol, rn);
		} else {
			if (ul != '\0' && di == '\0')
				setAttrib(res, R_NamesSymbol, cn);
			else if (!isNull(rn) &&
			         (rn == cn || equal_string_vectors(rn, cn, r)))
				setAttrib(res, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

	UNPROTECT(nprotect);
	return res;
}

/* diag(<[CRT]sparseMatrix>) <- value */
SEXP R_sparse_diag_set(SEXP obj, SEXP val)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, "R_sparse_diag_set");
	const char *clf = valid[ivalid];

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	++nprotect;
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;

	SEXPTYPE tv = TYPEOF(val);
	if (tv < LGLSXP || tv > REALSXP)
		/* Upper bound can become CPLXSXP once we have proper zMatrix */
		error(_("replacement diagonal has incompatible type \"%s\""),
		      type2char(tv));

	R_xlen_t nv = XLENGTH(val);
	if (nv != 1 && nv != r)
		error(_("replacement diagonal has wrong length"));

	SEXP x0 = NULL, x1 = NULL;
	SEXPTYPE tx = LGLSXP;
	PROTECT_INDEX pid;
	if (clf[0] != 'n') {
		PROTECT_WITH_INDEX(x0 = GET_SLOT(obj, Matrix_xSym), &pid);
		++nprotect;
		tx = TYPEOF(x0);
	}

	SEXP res;
	if (tv <= tx) {
		PROTECT(val = coerceVector(val, tv = tx));
		PROTECT(res = NEW_OBJECT_OF_CLASS(clf));
		nprotect += 2;
	} else { /* tv > tx */
		/* dMatrix is only possibility until we have proper [iz]Matrix */
		PROTECT(val = coerceVector(val, tv = tx = REALSXP));
		char clt[] = "d..Matrix";
		clt[1] = clf[1];
		clt[2] = clf[2];
		PROTECT(res = NEW_OBJECT_OF_CLASS(clt));
		nprotect += 2;
		if (clf[0] != 'n')
			REPROTECT(x0 = coerceVector(x0, tv), pid);
	}

	if (m != n || n > 0)
		SET_SLOT(res, Matrix_DimSym, dim);

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = 'U', di = 'N';
	if (clf[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */
		if (ul != 'U')
			SET_SLOT(res, Matrix_uploSym, uplo);

		if (clf[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
	}

	if (clf[2] != 'T') {

		SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
			i0 = PROTECT(GET_SLOT(obj, iSym));
		nprotect += 3;
		int *pp0 = INTEGER(p0), *pp1 = INTEGER(p1), *pi0 = INTEGER(i0),
			j, k = 0, kend, nd0 = 0, nd1 = 0, n_ = (clf[2] == 'C') ? n : m;
		pp0++;
		*(pp1++) = 0;

		if (clf[1] == 'g') {
			for (j = 0; j < r; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] >= j) {
						if (pi0[k] == j)
							++nd0;
						k = kend;
						break;
					}
					++k;
				}
				pp1[j] = kend - nd0;
			}
			for (j = r; j < n_; ++j)
				pp1[j] = pp0[j] - nd0;
		} else if (di != 'N') {
			for (j = 0; j < n_; ++j)
				pp1[j] = pp0[j];
		} else if (ul == ((clf[2] == 'C') ? 'U' : 'L')) {
			for (j = 0; j < n_; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[kend-1] == j)
					++nd0;
				k = kend;
				pp1[j] = kend - nd0;
			}
		} else {
			for (j = 0; j < n_; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[k] == j)
					++nd0;
				k = kend;
				pp1[j] = kend - nd0;
			}
		}

#define DO_COUNT(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *pval = _PTR_(val); \
			if (nv != 1) { \
				for (j = 0; j < r; ++j) { \
					if (_NZ_(pval[j])) \
						++nd1; \
					pp1[j] += nd1; \
				} \
				for (j = r; j < n_; ++j) \
					pp1[j] += nd1; \
			} else if (_NZ_(pval[0])) { \
				nd1 = r; \
				for (j = 0; j < r; ++j) \
					pp1[j] += j + 1; \
				for (j = r; j < n_; ++j) \
					pp1[j] += r; \
			} \
		} while (0)

		SPARSE_CASES(tv, DO_COUNT);

#undef DO_COUNT

		if (nd1 - nd0 > INT_MAX - pp0[n_-1])
			error(_("p[length(p)] cannot exceed 2^31-1"));

		int nnz1 = pp1[n_-1];
		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1));
		++nprotect;
		int *pi1 = INTEGER(i1);

		if (clf[0] != 'n' || tv != LGLSXP) {
			PROTECT(x1 = allocVector(tx, nnz1));
			++nprotect;
		}

#define SPARSE_D_S(_X_COPY_, _X_INSERT_, _NZ_VAL_) \
		do { \
			k = 0; \
			for (j = 0; j < r; ++j) { \
				kend = pp0[j]; \
				while (k < kend && pi0[k] < j) { \
					*(pi1++) = pi0[k]; \
					_X_COPY_; \
					++k; \
				} \
				if (k < kend && pi0[k] == j) \
					++k; \
				if (_NZ_VAL_) { \
					*(pi1++) = j; \
					_X_INSERT_; \
				} \
				while (k < kend) { \
					*(pi1++) = pi0[k]; \
					_X_COPY_; \
					++k; \
				} \
			} \
			for (j = r; j < n_; ++j) { \
				kend = pp0[j]; \
				while (k < kend && pi0[k] < j) { \
					*(pi1++) = pi0[k]; \
					_X_COPY_; \
					++k; \
				} \
			} \
		} while (0)

#define SPARSE_D_S_N(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px1 = _PTR_(x1), *pval = _PTR_(val); \
			if (nv != 1) \
				SPARSE_D_S(*(px1++) = _ONE_, *(px1++) = pval[j], \
				           _NZ_(pval[j])); \
			else { \
				int nz = _NZ_(pval[0]); \
				SPARSE_D_S(*(px1++) = _ONE_, *(px1++) = pval[j], \
				           nz); \
			} \
		} while (0)

#define SPARSE_D_S_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1), *pval = _PTR_(val); \
			if (nv != 1) \
				SPARSE_D_S(*(px1++) = px0[k], *(px1++) = pval[j], \
				           _NZ_(pval[j])); \
			else { \
				int nz = _NZ_(pval[0]); \
				SPARSE_D_S(*(px1++) = px0[k], *(px1++) = pval[0], \
				           nz); \
			} \
		} while (0)

#define SPARSE_D_S_END \
		do { \
			if (clf[0] != 'n') \
				SPARSE_CASES(tv, SPARSE_D_S_X); \
			else if (tv != LGLSXP) \
				SPARSE_CASES(tv, SPARSE_D_S_N); \
			else { \
				int *pval = LOGICAL(val); \
				if (nv != 1) \
					SPARSE_D_S(, , pval[j]); \
				else \
					SPARSE_D_S(, , pval[0]); \
			} \
		} while (0)

		SPARSE_D_S_END;

#undef SPARSE_D_S

		SET_SLOT(res, Matrix_pSym, p1);
		SET_SLOT(res,        iSym, i1);

	} else {

		SEXP i0 = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(obj, Matrix_jSym));
		nprotect += 2;
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), j, nd0 = 0, nd1 = 0;
		R_xlen_t k, nnz0 = XLENGTH(i0);

		for (k = 0; k < nnz0; ++k)
			if (pi0[k] == pj0[k])
				++nd0;

#define DO_COUNT(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
		do { \
			_CTYPE_ *pval = _PTR_(val); \
			if (nv != 1) { \
				for (j = 0; j < r; ++j) \
					if (_NZ_(pval[j])) \
						++nd1; \
			} else if (_NZ_(pval[0])) \
				nd1 = r; \
		} while (0)

		SPARSE_CASES(tv, DO_COUNT);

#undef DO_COUNT

		if (nd1 - nd0 > R_XLEN_T_MAX - nnz0)
			error(_("length(i) cannot exceed R_XLEN_T_MAX"));

		R_xlen_t nnz1 = nnz0 + (nd1 - nd0);
		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		nprotect += 2;
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

		if (clf[0] != 'n' || tv != LGLSXP) {
			PROTECT(x1 = allocVector(tx, nnz1));
			++nprotect;
		}

#define SPARSE_D_S(_X_COPY_, _X_INSERT_, _NZ_VAL_) \
		do { \
			for (k = 0; k < nnz0; ++k) { \
				if (pi0[k] != pj0[k]) { \
					*(pi1++) = pi0[k]; \
					*(pj1++) = pj0[k]; \
					_X_COPY_; \
				} \
			} \
			for (j = 0; j < r; ++j) { \
				if (_NZ_VAL_) { \
					*(pi1++) = *(pj1++) = j; \
					_X_INSERT_; \
				} \
			} \
		} while (0)

		SPARSE_D_S_END;

#undef SPARSE_D_S_END
#undef SPARSE_D_S_X
#undef SPARSE_D_S_N
#undef SPARSE_D_S

		SET_SLOT(res, Matrix_iSym, i1);
		SET_SLOT(res, Matrix_jSym, j1);

	}

	if (clf[0] != 'n' || tv != LGLSXP)
		SET_SLOT(res, Matrix_xSym, x1);

	UNPROTECT(nprotect);
	return res;
}

/* diagU2N(<[CRT]sparseMatrix>), parallel to R-level ..diagU2N(),
   though that is more general, working for _all_ Matrix */
SEXP R_sparse_diag_U2N(SEXP obj) {
	if (!HAS_SLOT(obj, Matrix_diagSym))
		return obj;
	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	if (di == 'N')
		return obj;
	SEXP val = PROTECT(ScalarLogical(1)),
		res = PROTECT(R_sparse_diag_set(obj, val));
	UNPROTECT(2); /* res, val */
	return res;
}

/* diagU2N(<[CRT]sparseMatrix>), parallel to R-level ..diagN2U(),
   though that is more general, working for _all_ Matrix */
SEXP R_sparse_diag_N2U(SEXP obj) {
	if (!HAS_SLOT(obj, Matrix_diagSym))
		return obj;
	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	if (di != 'N')
		return obj;
	PROTECT(diag = mkString("U"));
	SEXP res, dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n > 0) {
		SEXP k, uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */
		if (ul == 'U') {
			PROTECT(k = ScalarInteger( 1));
			PROTECT(res = R_sparse_band(obj, k, R_NilValue));
		} else {
			PROTECT(k = ScalarInteger(-1));
			PROTECT(res = R_sparse_band(obj, R_NilValue, k));
		}
		SET_SLOT(res, Matrix_diagSym, diag);
		UNPROTECT(3); /* res, k, diag */
	} else {
		PROTECT(res = duplicate(obj));
		SET_SLOT(res, Matrix_diagSym, diag);
		UNPROTECT(2); /* res,    diag */
	}
	return res;
}

/* t(<[CRT]sparseMatrix>) */
SEXP R_sparse_transpose(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_transpose");
	const char *cl = valid[ivalid];

	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n) {
		UNPROTECT(1); /* dim */
		PROTECT(dim = GET_SLOT(to, Matrix_DimSym));
		pdim = INTEGER(dim);
		pdim[0] = n;
		pdim[1] = m;
	} else if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (cl[1] != 's')
		set_reversed_DimNames(to, dimnames);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (cl[1] != 'g') {
		SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ulf = *CHAR(STRING_ELT(uplo_from, 0));
		UNPROTECT(1); /* uplo_from */
		if (ulf == 'U') {
			SEXP uplo_to = PROTECT(mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo_to);
			UNPROTECT(1); /* uplo_to */
		}
		if (cl[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			char di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		} else {
			SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
			if (LENGTH(factors) > 0)
				SET_SLOT(to, Matrix_factorSym, factors);
			UNPROTECT(1); /* factors */
		}
	}

	/* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

	if (cl[2] == 'T') {
		/* No need to allocate in this case: need only reverse 'i' and 'j' */
		SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j = PROTECT(GET_SLOT(from, Matrix_jSym));
		SET_SLOT(to, Matrix_iSym, j);
		SET_SLOT(to, Matrix_jSym, i);
		UNPROTECT(2); /* j, i */
		if (cl[0] != 'n') {
			SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
			SET_SLOT(to, Matrix_xSym, x);
			UNPROTECT(1); /* x */
		}
		UNPROTECT(1); /* to */
		return to;
	}

	/* Now dealing only with [CR]sparseMatrix ... */

	int m_, n_;
	SEXP iSym;
	if (cl[2] == 'C') {
		m_ = m;
		n_ = n;
		iSym = Matrix_iSym;
	} else {
		m_ = n;
		n_ = m;
		iSym = Matrix_jSym;
	}

	R_xlen_t m1a = (R_xlen_t) m_ + 1;
	SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
		p1 = PROTECT(allocVector(INTSXP, m1a));
	int i, j, k, kend,
		*pp0 = INTEGER(p0),
		*pp1 = INTEGER(p1),
		nnz = pp0[n_];
	SEXP i0 = PROTECT(GET_SLOT(from, iSym)),
		i1 = PROTECT(allocVector(INTSXP, nnz));
	int *pi0 = INTEGER(i0),
		*pi1 = INTEGER(i1);
	++pp0;

	/* Counting number of nonzero elements, by "row" */
	Matrix_memset(pp1, 0, m1a, sizeof(int));
	++pp1;
	for (k = 0; k < nnz; ++k)
		++pp1[pi0[k]];

	/* Computing cumulative sum, in place */
	for (i = 1; i < m_; ++i)
		pp1[i] += pp1[i-1];

	/* Allocating work space */
	int *pp1_;
	Matrix_Calloc(pp1_, m_, int);
	Matrix_memcpy(pp1_, pp1 - 1, m_, sizeof(int));

#define SPARSE_T(_XASSIGN_) \
	do { \
		for (j = 0, k = 0; j < n_; ++j) { \
			kend = pp0[j]; \
			while (k < kend) { \
				i = pi0[k]; \
				pi1[pp1_[i]] = j; \
				_XASSIGN_; /* px1[pp1_[i]] = px0[k] */ \
				++pp1_[i]; \
				++k; \
			} \
		} \
	} while (0)

#define SPARSE_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		SPARSE_T(px1[pp1_[i]] = px0[k]); \
	} while (0);

	if (cl[0] == 'n')
		SPARSE_T();
	else {
		SEXPTYPE tx;
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
			x1 = PROTECT(allocVector(tx = TYPEOF(x0), nnz));
		SPARSE_CASES(tx, SPARSE_T_X);
		SET_SLOT(to, Matrix_xSym, x1);
		UNPROTECT(2); /* x1, x0 */
	}

#undef SPARSE_T_X
#undef SPARSE_T

	Matrix_Free(pp1_, m_);
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to,        iSym, i1);

	UNPROTECT(5); /* i1, i0, p1, p0, to */
	return to;
}

/* forceSymmetric(<[CRT]sparseMatrix>, uplo) */
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo_to)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_force_symmetric");
	const char *clf = valid[ivalid];

	char ulf = 'U', ult = 'U';
	if (clf[1] != 'g') {
		SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ulf = ult = *CHAR(STRING_ELT(uplo_from, 0));
		UNPROTECT(1); /* uplo_from */
	}

	if (!isNull(uplo_to) &&
	    (TYPEOF(uplo_to) != STRSXP || LENGTH(uplo_to) < 1 ||
	    (uplo_to = STRING_ELT(uplo_to, 0)) == NA_STRING ||
	    ((ult = *CHAR(uplo_to)) != 'U' && ult != 'L')))
		error(_("invalid 'uplo' to 'R_sparse_force_symmetric()'"));

	if (clf[1] == 's') {
		/* .s[CRT]Matrix */
		if (ulf == ult)
			return from;
		SEXP to = PROTECT(R_sparse_transpose(from));
		if (clf[0] == 'z') {
			/* Need _conjugate_ transpose */
			SEXP x_to = PROTECT(GET_SLOT(from, Matrix_xSym));
			conjugate(x_to);
			UNPROTECT(1); /* x_to */
		}
		UNPROTECT(1) /* to */;
		return to;
	}

	/* Now handling just .[gt][CRT]Matrix ... */

	char clt[] = ".s.Matrix";
	clt[0] = clf[0];
	clt[2] = clf[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
	++nprotect;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to symmetrize a non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (ult != 'U') {
		PROTECT(uplo_to = mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo_to);
		UNPROTECT(1); /* uplo_to */
	}

	/* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

	char di = 'N';
	if (clf[1] == 't') {
		/* .t[CRT]Matrix */
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */
	}

	if (clf[1] == 't' && di == 'N' && ulf == ult) {

		/* No need to allocate in this case: we have the triangle we want */
		if (clf[2] != 'T') {
			SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
			SET_SLOT(to, Matrix_pSym, p);
			UNPROTECT(1); /* p */
		}
		if (clf[2] != 'R') {
			SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
			SET_SLOT(to, Matrix_iSym, i);
			UNPROTECT(1); /* i */
		}
		if (clf[2] != 'C') {
			SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
			SET_SLOT(to, Matrix_jSym, j);
			UNPROTECT(1); /* j */
		}
		if (clf[0] != 'n') {
			SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
			SET_SLOT(to, Matrix_xSym, x);
			UNPROTECT(1); /* x */
		}
		UNPROTECT(nprotect);
		return to;

	} else if (clf[2] != 'T') {

		/* Symmetrizing square .[gt][CR]Matrix ... */

		SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
			i0 = PROTECT(GET_SLOT(from, iSym));
		nprotect += 3;
		int j, k, kend,
			*pp0 = INTEGER(p0),
			*pp1 = INTEGER(p1),
			*pi0 = INTEGER(i0),
			nnz0 = pp0[n],
			nnz1 = 0;
		pp0++;
		*(pp1++) = 0;

		/* Counting number of nonzero elements in triangle, by "column" ... */

		if (clf[1] == 't') {
			if (di != 'N') {
				/* Have triangular matrix with unit diagonal */
				if (ulf != ult) {
					/* Returning identity matrix */
					for (j = 0; j < n; ++j)
						pp1[j] = ++nnz1;
				} else {
					/* Returning symmetric matrix with unit diagonal */
					for (j = 0; j < n; ++j)
						pp1[j] = ++nnz1 + pp0[j];
					nnz1 += nnz0;
				}
			} else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
				/* Have triangular matrix with non-unit "trailing" diagonal
				   and returning diagonal part */
				for (j = 0; j < n; ++j) {
					if (pp0[j-1] < pp0[j] && pi0[pp0[j]-1] == j)
						++nnz1;
					pp1[j] = nnz1;
				}
			} else {
				/* Have triangular matrix with non-unit "leading" diagonal
				   and returning diagonal part */
				for (j = 0; j < n; ++j) {
					if (pp0[j-1] < pp0[j] && pi0[pp0[j-1]] == j)
						++nnz1;
					pp1[j] = nnz1;
				}
			}
		} else if (ult == ((clf[2] == 'C') ? 'U' : 'L')) {
			/* Have general matrix and returning upper triangle */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] <= j)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
		} else {
			/* Have general matrix and returning lower triangle */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] >= j)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
		}

		/* Now allocating and filling out slots ... */

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1));
		++nprotect;
		int *pi1 = INTEGER(i1);

		SEXP x0 = NULL, x1 = NULL;
		SEXPTYPE tx = NILSXP;
		if (clf[0] != 'n') {
			PROTECT(x0 = GET_SLOT(from, Matrix_xSym));
			PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
			nprotect += 2;
		}

		if (clf[1] == 't') {
			if (di != 'N') {
				/* Have triangular matrix with unit diagonal */
				if (ulf != ult) {
					/* Returning identity matrix */

#define SPARSE_FS(_XASSIGN_) \
					do { \
						for (j = 0; j < n; ++j) { \
							*(pi1++) = j; \
							_XASSIGN_; /* *(px1++) = _ONE_; */ \
						} \
					} while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
					do { \
						_CTYPE_ *px1 = _PTR_(x1); \
						SPARSE_FS(*(px1++) = _ONE_); \
					} while (0)

					if (clf[0] == 'n')
						SPARSE_FS();
					else
						SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS

				} else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
					/* Returning symmetric matrix
					   with unit "trailing" diagonal */

#define SPARSE_FS(_XASSIGN_, _XASSIGN_ONE_) \
					do { \
						for (j = 0, k = 0; j < n; ++j) { \
							kend = pp0[j]; \
							while (k < kend) { \
								*(pi1++) = pi0[k]; \
								_XASSIGN_; /* *(px1++) = px0[k]; */ \
								++k; \
							} \
							*(pi1++) = j; \
							_XASSIGN_ONE_; /* *(px1++) = _ONE_; */ \
						} \
					} while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
					do { \
						_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
						SPARSE_FS(*(px1++) = px0[k], *(px1++) = _ONE_); \
					} while (0)

					if (clf[0] == 'n')
						SPARSE_FS(, );
					else
						SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS

				} else {
					/* Returning symmetric matrix
					   with unit "leading" diagonal */

#define SPARSE_FS(_XASSIGN_, _XASSIGN_ONE_) \
					do { \
						for (j = 0, k = 0; j < n; ++j) { \
							*(pi1++) = j; \
							_XASSIGN_ONE_; /* *(px1++) = _ONE_; */ \
							kend = pp0[j]; \
							while (k < kend) { \
								*(pi1++) = pi0[k]; \
								_XASSIGN_; /* *(px1++) = px0[k]; */ \
								++k; \
							} \
						} \
					} while (0)

					if (clf[0] == 'n')
						SPARSE_FS(, );
					else
						SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS

				}
			} else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
				/* Have triangular matrix with non-unit "trailing" diagonal
				   and returning diagonal part */

#define SPARSE_FS(_XASSIGN_) \
				do { \
					for (j = 0; j < n; ++j) { \
						if (pp0[j-1] < pp0[j] && pi0[pp0[j]-1] == j) { \
							*(pi1++) = j; \
							_XASSIGN_; /* *(px1++) = px0[pp0[j]-1]; */ \
						} \
					} \
				} while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
				do { \
					_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
					SPARSE_FS(*(px1++) = px0[pp0[j]-1]); \
				} while (0)

				if (clf[0] == 'n')
					SPARSE_FS();
				else
					SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS

			} else {
				/* Have triangular matrix with non-unit "leading" diagonal
				   and returning diagonal part */

#define SPARSE_FS(_XASSIGN_) \
				do { \
					for (j = 0; j < n; ++j) { \
						if (pp0[j-1] < pp0[j] && pi0[pp0[j-1]] == j) { \
							*(pi1++) = j; \
							_XASSIGN_; /* *(px1++) = px0[pp0[j-1]]; */ \
						} \
					} \
				} while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
				do { \
					_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
					SPARSE_FS(*(px1++) = px0[pp0[j-1]]); \
				} while (0)

				if (clf[0] == 'n')
					SPARSE_FS();
				else
					SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS

			}
		} else if (ult == ((clf[2] == 'C') ? 'U' : 'L')) {
			/* Have general matrix and returning upper triangle */

#define SPARSE_FS(_XASSIGN_) \
			do { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] <= j) { \
							*(pi1++) = pi0[k]; \
							_XASSIGN_; /* *(px1++) = px0[k]; */ \
						} \
						++k; \
					} \
				} \
			} while (0)

#define SPARSE_FS_X_BASIC(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
				SPARSE_FS(*(px1++) = px0[k]); \
			} while (0)

			if (clf[0] == 'n')
				SPARSE_FS();
			else
				SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS

		} else {
			/* Have general matrix and returning lower triangle */

#define SPARSE_FS(_XASSIGN_) \
			do { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] >= j) { \
							*(pi1++) = pi0[k]; \
							_XASSIGN_; /* *(px1++) = px0[k]; */ \
						} \
						++k; \
					} \
				} \
			} while (0)

			if (clf[0] == 'n')
				SPARSE_FS();
			else
				SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS

		}

		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to,        iSym, i1);
		if (clf[0] != 'n')
			SET_SLOT(to, Matrix_xSym, x1);

	} else {

		/* Symmetrizing square .[gt]TMatrix ... */

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		nprotect += 2;
		int *pi0 = INTEGER(i0),
			*pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = 0;

	/* Counting number of nonzero elements in triangle ... */

		if (clf[1] == 't' && di != 'N')
			nnz1 = (ulf == ult) ? n + nnz0 : n;
		else {
			if (ult == 'U') {
				for (k = 0; k < nnz0; ++k)
					if (pi0[k] <= pj0[k])
						++nnz1;
			} else {
				for (k = 0; k < nnz0; ++k)
					if (pi0[k] >= pj0[k])
						++nnz1;
			}
		}

		/* Now allocating and filling out slots ... */

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		nprotect += 2;
		int *pi1 = INTEGER(i1),
			*pj1 = INTEGER(j1);

		SEXP x0 = NULL, x1 = NULL;
		SEXPTYPE tx = NILSXP;
		if (clf[0] != 'n') {
			PROTECT(x0 = GET_SLOT(from, Matrix_xSym));
			PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
			nprotect += 2;
		}

		if (clf[1] == 't' && di != 'N') {
			if (ulf == ult) {
				Matrix_memcpy(pi1, pi0, nnz0, sizeof(int));
				Matrix_memcpy(pj1, pj0, nnz0, sizeof(int));
				pi1 += nnz0;
				pj1 += nnz0;
			}

#define SPARSE_FS(_XASSIGN_) \
			do { \
				int j; \
				for (j = 0; j < n; ++j) { \
					*(pi1++) = *(pj1++) = j; \
					_XASSIGN_; /* *(px1++) = _ONE_; */ \
				} \
			} while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_, _NZ_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
				if (ulf == ult) { \
					Matrix_memcpy(px1, px0, nnz0, sizeof(_CTYPE_)); \
					px1 += nnz0; \
				} \
				SPARSE_FS(*(px1++) = _ONE_); \
			} while (0)

			if (clf[0] == 'n')
				SPARSE_FS();
			else
				SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS

		} else {

#define SPARSE_FS(_XASSIGN_) \
			do { \
				if (ult == 'U') { \
					for (k = 0; k < nnz0; ++k) { \
						if (pi0[k] <= pj0[k]) { \
							*(pi1++) = pi0[k]; \
							*(pj1++) = pj0[k]; \
							_XASSIGN_; /* *(px1++) = px0[k]; */ \
						} \
					} \
				} else { \
					for (k = 0; k < nnz0; ++k) { \
						if (pi0[k] <= pj0[k]) { \
							*(pi1++) = pi0[k]; \
							*(pj1++) = pj0[k]; \
							_XASSIGN_; /* *(px1++) = px0[k]; */ \
						} \
					} \
				} \
			} while (0)

			if (clf[0] == 'n')
				SPARSE_FS();
			else
				SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS_X_BASIC
#undef SPARSE_FS

		}

		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		if (clf[0] != 'n')
			SET_SLOT(to, Matrix_xSym, x1);

	}

	UNPROTECT(nprotect);
	return to;
}

/* symmpart(<[CRT]sparseMatrix>) */
SEXP R_sparse_symmpart(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_symmpart");
	const char *clf = valid[ivalid];
	if (clf[0] == 'd' && clf[1] == 's')
		return from;

	PROTECT_INDEX pidA;
	PROTECT_WITH_INDEX(from, &pidA);
	++nprotect;

	char clt[] = ".s.Matrix";
	clt[0] = (clf[0] != 'z') ? 'd' : 'z';
	clt[2] = clf[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
	++nprotect;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to get symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (clf[1] != 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (clf[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	} else if (clf[2] == 'R') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	char di = 'N';
	if (clf[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N')
			REPROTECT(from = R_sparse_as_general(from), pidA); /* U->N */
		UNPROTECT(1); /* diag */
	}

	SEXP x0 = NULL, x1 = NULL;
	PROTECT_INDEX pidB;
	if (clf[0] != 'n') {
		PROTECT_WITH_INDEX(x0 = GET_SLOT(from, Matrix_xSym), &pidB);
		++nprotect;
	}

#define SPARSE_SYMMPART_CASES(_SYMMPART_) \
	do { \
		switch (clf[0]) { \
		case 'n': \
		{ \
			double *px1 = REAL(x1); \
			_SYMMPART_(px1[k] = 1.0, px1[k] *= 0.5); \
			break; \
		} \
		case 'l': \
		{ \
			int *px0 = LOGICAL(x0); \
			double *px1 = REAL(x1); \
			_SYMMPART_( \
				px1[k]  = (px0[k] == NA_LOGICAL \
				           ? NA_REAL : ((px0[k] != 0) ? 1.0 : 0.0)), \
				px1[k] *= 0.5); \
			break; \
		} \
		case 'i': \
		{ \
			int *px0 = INTEGER(x0); \
			double *px1 = REAL(x1); \
			_SYMMPART_( \
				px1[k]  = (px0[k] == NA_INTEGER \
						   ? NA_REAL : (double) px0[k]), \
				px1[k] *= 0.5); \
			break; \
		} \
		case 'd': \
		{ \
			double *px0 = REAL(x0), *px1 = REAL(x1); \
			_SYMMPART_(px1[k] = px0[k], px1[k] *= 0.5); \
			break; \
		} \
		case 'z': \
		{ \
			Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1); \
			_SYMMPART_( \
				px1[k] = px0[k], \
				do { px1[k].r *= 0.5; px1[k].i *= 0.5; } while (0)); \
			break; \
		} \
		default: \
			break; \
		} \
	} while (0)

	if (clf[2] != 'T') {

		SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from, iSym));
		nprotect += 2;
		int j, k, kend, *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), nnz = pp0[n];
		++pp0;

		if (clf[1] == 'g') {

			REPROTECT(from = R_sparse_transpose(from), pidA);
			SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
				p0_ = PROTECT(GET_SLOT(from, Matrix_pSym)),
				i0_ = PROTECT(GET_SLOT(from, iSym));
			nprotect += 3;
			int k_, kend_,
				*pp1 = INTEGER(p1),
				*pp0_ = INTEGER(p0_) + 1,
				*pi0_ = INTEGER(i0_);
			*(pp1++) = 0;

			SEXP x0_ = NULL;
			if (clf[0] != 'n') {
				PROTECT(x0_ = GET_SLOT(from, Matrix_xSym));
				++nprotect;
			}

			/* Counting number of nonzero elements in each "column"
			   of result ... */

			for (j = 0, k = 0, k_ = 0; j < n; ++j) {
				pp1[j] = pp1[j-1];
				kend = pp0[j];
				kend_ = pp0_[j];
				while (k < kend) {
					if (pi0[k] > j) {
						k = kend;
						break;
					}
					while (k_ < kend_ && pi0_[k_] < pi0[k]) {
						++pp1[j];
						++k_;
					}
					++pp1[j];
					if (k_ < kend_ && pi0_[k_] == pi0[k])
						++k_;
					++k;
				}
				while (k_ < kend_) {
					if (pi0_[k_] > j) {
						k_ = kend_;
						break;
					}
					++pp1[j];
					++k_;
				}
			}

			SEXP i1 = PROTECT(allocVector(INTSXP, pp1[n-1]));
			++nprotect;
			int *pi1 = INTEGER(i1);

			PROTECT(x1 = allocVector(kind2type(clt[0]), pp1[n-1]));
			++nprotect;

#define CRSPARSE_SYMMPART_GENERAL(_DO_ASSIGN_, \
                                  _DO_ASSIGN_FROM_TRANSPOSE_, \
                                  _DO_INCR_FROM_TRANSPOSE_) \
			do { \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend = pp0[j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (pi0[k] > j) { \
							k = kend; \
							break; \
						} \
						while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
							*pi1 = pi0_[k_]; \
							_DO_ASSIGN_FROM_TRANSPOSE_; \
							++pi1; ++px1; ++k_; \
						} \
						*pi1 = pi0[k]; \
						_DO_ASSIGN_; \
						if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
							_DO_INCR_FROM_TRANSPOSE_; \
							++k_; \
						} \
						++pi1; ++px1; ++k; \
					} \
					while (k_ < kend_) { \
						if (pi0_[k_] > j) { \
							k_ = kend_; \
							break; \
						} \
						*pi1 = pi0_[k_]; \
						_DO_ASSIGN_FROM_TRANSPOSE_; \
						++pi1; ++px1; ++k_; \
					} \
				} \
			} while (0)

			switch (clf[0]) {
			case 'n':
			{
				double *px1 = REAL(x1);
				CRSPARSE_SYMMPART_GENERAL(
					*px1  = 0.5,
					*px1  = 0.5,
					*px1 += 0.5);
				break;
			}
			case 'l':
			{
				int *px0 = LOGICAL(x0), *px0_ = LOGICAL(x0_);
				double *px1 = REAL(x1);
				CRSPARSE_SYMMPART_GENERAL(
					*px1  = (px0[k] == NA_LOGICAL
					         ? NA_REAL : ((px0[k] != 0) ? 0.5 : 0.0)),
					*px1  = (px0_[k_] == NA_LOGICAL
					         ? NA_REAL : ((px0_[k_] != 0) ? 0.5 : 0.0)),
					*px1 += (px0_[k_] == NA_LOGICAL
					         ? NA_REAL : ((px0_[k_] != 0) ? 0.5 : 0.0)));
				break;
			}
			case 'i':
			{
				int *px0 = INTEGER(x0), *px0_ = INTEGER(x0_);
				double *px1 = REAL(x1);
				CRSPARSE_SYMMPART_GENERAL(
					*px1  = (px0[k] == NA_INTEGER
					         ? NA_REAL : 0.5 * (double) px0[k]),
					*px1  = (px0_[k_] == NA_INTEGER
					         ? NA_REAL : 0.5 * (double) px0_[k_]),
					*px1 += (px0_[k_] == NA_INTEGER
					         ? NA_REAL : 0.5 * (double) px0_[k_]));
				break;
			}
			case 'd':
			{
				double *px0 = REAL(x0), *px0_ = REAL(x0_),
					*px1 = REAL(x1);
				CRSPARSE_SYMMPART_GENERAL(
					*px1  = 0.5 * px0[k],
					*px1  = 0.5 * px0_[k_],
					*px1 += 0.5 * px0_[k_]);
				break;
			}
			case 'z':
			{
				Rcomplex *px0 = COMPLEX(x0), *px0_ = COMPLEX(x0_),
					*px1 = COMPLEX(x1);
				CRSPARSE_SYMMPART_GENERAL(
					do {
						(*px1).r  = 0.5 * px0[k].r;
						(*px1).i  = 0.5 * px0[k].i;
					} while (0),
					do {
						(*px1).r  = 0.5 * px0_[k_].r;
						(*px1).i  = 0.5 * px0_[k_].i;
					} while (0),
					do {
						(*px1).r += 0.5 * px0_[k_].r;
						(*px1).i += 0.5 * px0_[k_].i;
					} while (0));
				break;
			}
			default:
				break;
			}

#undef CRSPARSE_SYMMPART_GENERAL

			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to,        iSym, i1);

		} else if (clf[1] == 't') {

			if (clf[0] == clt[0] && di != 'N')
				x1 = x0;
			else {
				PROTECT(x1 = allocVector(kind2type(clt[0]), nnz));
				++nprotect;
			}

#define CRSPARSE_SYMMPART_TRIANGULAR(_DO_ASSIGN_, _DO_HALF_) \
			do { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						_DO_ASSIGN_; \
						if (pi0[k] != j) \
							_DO_HALF_; \
						++k; \
					} \
				} \
			} while (0)


			SPARSE_SYMMPART_CASES(CRSPARSE_SYMMPART_TRIANGULAR);

#undef CRSPARSE_SYMMPART_TRIANGULAR

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

		} else {

#define SPARSE_SYMMPART_CASES_TRIVIAL \
			do { \
				switch (clf[0]) { \
				case 'n': \
				{ \
					PROTECT(x1 = allocVector(REALSXP, nnz)); \
					++nprotect; \
					double *px1 = REAL(x1); \
					while (nnz--) \
						*(px1++) = 1.0; \
					break; \
				} \
				case 'l': \
				case 'i': \
					REPROTECT(x1 = coerceVector(x0, REALSXP), pidB); \
					break; \
				case 'd': \
					x1 = x0; \
					break; \
				case 'z': \
					REPROTECT(x1 = duplicate(x0), pidB); \
					zeroIm(x1); \
					break; \
				default: \
					break; \
				} \
			} while (0)

			SPARSE_SYMMPART_CASES_TRIVIAL;

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

		}

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		nprotect += 2;
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz = XLENGTH(i0);

		PROTECT(x1 = allocVector(kind2type(clt[0]), nnz));
		++nprotect;

		if (clf[1] == 'g') {

			SEXP i1 = PROTECT(allocVector(INTSXP, nnz)),
				j1 = PROTECT(allocVector(INTSXP, nnz));
			nprotect += 2;
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

#define TSPARSE_SYMMPART_GENERAL(_DO_ASSIGN_, _DO_HALF_) \
			do { \
				for (k = 0; k < nnz; ++k) { \
					if (pi0[k] <= pj0[k]) { \
						pi1[k] = pi0[k]; \
						pj1[k] = pj0[k]; \
					} else { \
						pi1[k] = pj0[k]; \
						pj1[k] = pi0[k]; \
					} \
					_DO_ASSIGN_; \
					if (pi0[k] != pj0[k]) \
						_DO_HALF_; \
				} \
			} while (0)

			SPARSE_SYMMPART_CASES(TSPARSE_SYMMPART_GENERAL);

#undef TSPARSE_SYMMPART_GENERAL

			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);

	} else if (clf[1] == 't') {

#define TSPARSE_SYMMPART_TRIANGULAR(_DO_ASSIGN_, _DO_HALF_) \
			do { \
				for (k = 0; k < nnz; ++k) { \
					_DO_ASSIGN_; \
					if (pi0[k] != pj0[k]) \
						_DO_HALF_; \
				} \
			} while (0)

			SPARSE_SYMMPART_CASES(TSPARSE_SYMMPART_TRIANGULAR);

#undef TSPARSE_SYMMPART_TRIANGULAR

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

		} else {

			SPARSE_SYMMPART_CASES_TRIVIAL;

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

		}

	}

#undef SPARSE_SYMMPART_CASES_TRIVIAL
#undef SPARSE_SYMMPART_CASES

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(nprotect);
	return to;
}

/* skewpart(<[CRT]sparseMatrix>) */
SEXP R_sparse_skewpart(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "R_sparse_skewpart");
	const char *clf = valid[ivalid];

	PROTECT_INDEX pidA;
	PROTECT_WITH_INDEX(from, &pidA);
	++nprotect;

	char clt[] = "...Matrix";
	clt[0] = (clf[0] != 'z') ? 'd' : 'z';
	clt[1] = (clf[1] != 's') ? 'g' : 's';
	clt[2] = clf[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
	++nprotect;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to get skew-symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (clf[1] != 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP x0 = NULL, x1 = NULL;
	PROTECT_INDEX pidB;
	if (clf[0] != 'n') {
		PROTECT_WITH_INDEX(x0 = GET_SLOT(from, Matrix_xSym), &pidB);
		++nprotect;
	}

	if (clf[1] == 's') {

		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		if (clf[0] != 'z') {
			/* Skew-symmetric part of symmetric matrix is zero matrix */
			if (clf[2] != 'T') {
				R_xlen_t n1a = (R_xlen_t) n + 1;
				SEXP p1 = PROTECT(allocVector(INTSXP, n1a));
				int *pp1 = INTEGER(p1);
				Matrix_memset(pp1, 0, n1a, sizeof(int));
				SET_SLOT(to, Matrix_pSym, p1);
				UNPROTECT(1); /* p1 */
			}
		} else {
			/* Skew-symmetric part of Hermitian matrix is imaginary part */
			REPROTECT(x1 = duplicate(x0), pidB);
			zeroRe(x1);
			SET_SLOT(to, Matrix_xSym, x1);
			if (clf[2] != 'T') {
				SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym));
				SET_SLOT(to, Matrix_pSym, p0);
				UNPROTECT(1); /* p0 */
			}
			if (clf[2] != 'R') {
				SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym));
				SET_SLOT(to, Matrix_iSym, i0);
				UNPROTECT(1); /* i0 */
			}
			if (clf[2] != 'C') {
				SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
				SET_SLOT(to, Matrix_jSym, j0);
				UNPROTECT(1); /* j0 */
			}
		}

	} else if (clf[2] != 'T') {

		SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
			i0 = PROTECT(GET_SLOT(from, iSym));
		nprotect += 3;
		int j, k, kend, k_, kend_,
			*pp0 = INTEGER(p0),
			*pp1 = INTEGER(p1),
			*pi0 = INTEGER(i0);
		++pp0;

		REPROTECT(from = R_sparse_transpose(from), pidA);

		SEXP p0_ = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0_ = PROTECT(GET_SLOT(from, iSym));
		nprotect += 2;
		int *pp0_ = INTEGER(p0_), *pi0_ = INTEGER(i0_);
		++pp0_;

		int *pp1_;
		Matrix_Calloc(pp1_, n, int);

		SEXP x0_ = NULL;
		if (clf[0] != 'n') {
			PROTECT(x0_ = GET_SLOT(from, Matrix_xSym));
			++nprotect;
		}

		/* Counting number of nonzero elements in each "column"
		   of result ... */

		for (j = 0, k = 0, k_ = 0; j < n; ++j) {
			kend = pp0[j];
			kend_ = pp0_[j];
			while (k < kend) {
				if (pi0[k] >= j) {
					k = kend;
					break;
				}
				while (k_ < kend_ && pi0_[k_] < pi0[k]) {
					++pp1_[j];
					++pp1_[pi0_[k_]];
					++k_;
				}
				++pp1_[j];
				++pp1_[pi0[k]];
				if (k_ < kend_ && pi0_[k_] == pi0[k])
					++k_;
				++k;
			}
			while (k_ < kend_) {
				if (pi0_[k_] >= j) {
					k_ = kend_;
					break;
				}
				++pp1_[j];
				++pp1_[pi0_[k_]];
				++k_;
			}
		}

		*(pp1++) = 0;
		for (j = 0; j < n; ++j) {
			pp1[j] = pp1[j-1] + pp1_[j];
			pp1_[j] = pp1[j-1];
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[n-1]));
		++nprotect;
		int *pi1 = INTEGER(i1);

		PROTECT(x1 = allocVector(kind2type(clt[0]), pp1[n-1]));
		++nprotect;

#define CRSPARSE_SKEWPART_GENERAL(_DO_ASSIGN_, \
                                  _DO_ASSIGN_FROM_TRANSPOSE_, \
                                  _DO_INCR_FROM_TRANSPOSE_, \
                                  _DO_NEGATE_, \
                                  _DO_NEGATE_FROM_TRANSPOSE_) \
		do { \
			for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
				kend = pp0[j]; \
				kend_ = pp0_[j]; \
				while (k < kend) { \
					if (pi0[k] >= j) { \
						k = kend; \
						break; \
					} \
					while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
						pi1[pp1_[j]] = pi0_[k_]; \
						_DO_ASSIGN_FROM_TRANSPOSE_; \
						pi1[pp1_[pi0_[k_]]] = j; \
						_DO_NEGATE_FROM_TRANSPOSE_; \
						++pp1_[j]; \
						++pp1_[pi0_[k_]]; \
						++k_; \
					} \
					pi1[pp1_[j]] = pi0[k]; \
					_DO_ASSIGN_; \
					if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
						_DO_INCR_FROM_TRANSPOSE_; \
						++k_; \
					} \
					pi1[pp1_[pi0[k]]] = j; \
					_DO_NEGATE_; \
					++pp1_[j]; \
					++pp1_[pi0[k]]; \
					++k; \
				} \
				while (k_ < kend_) { \
					if (pi0_[k_] >= j) { \
						k_ = kend_; \
						break; \
					} \
					pi1[pp1_[j]] = pi0_[k_]; \
					_DO_ASSIGN_FROM_TRANSPOSE_; \
					pi1[pp1_[pi0_[k_]]] = j; \
					_DO_NEGATE_FROM_TRANSPOSE_; \
					++pp1_[j]; \
					++pp1_[pi0_[k_]]; \
					++k_; \
				} \
			} \
		} while (0)

		switch (clf[0]) {
		case 'n':
		{
			double *px1 = REAL(x1);
			CRSPARSE_SKEWPART_GENERAL(
				px1[pp1_[j]]         =  0.5,
				px1[pp1_[j]]         = -0.5,
				px1[pp1_[j]]        -=  0.5,
				px1[pp1_[pi0[k]]]    = -px1[pp1_[j]],
				px1[pp1_[pi0_[k_]]]  = -px1[pp1_[j]]);
			break;
		}
		case 'l':
		{
			int *px0 = LOGICAL(x0), *px0_ = LOGICAL(x0_);
			double *px1 = REAL(x1);
			CRSPARSE_SKEWPART_GENERAL(
				px1[pp1_[j]]         = (px0[k]   == NA_LOGICAL
				                        ? NA_REAL
				                        : ((px0[k]   != 0) ?  0.5 : 0.0)),
				px1[pp1_[j]]         = (px0_[k_] == NA_LOGICAL
				                        ? NA_REAL
				                        : ((px0_[k_] != 0) ? -0.5 : 0.0)),
				px1[pp1_[j]]        -= (px0_[k_] == NA_LOGICAL
				                        ? NA_REAL
				                        : ((px0_[k_] != 0) ?  0.5 : 0.0)),
				px1[pp1_[pi0[k]]]    = -px1[pp1_[j]],
				px1[pp1_[pi0_[k_]]]  = -px1[pp1_[j]]);
			break;
		}
		case 'i':
		{
			int *px0 = INTEGER(x0), *px0_ = INTEGER(x0_);
			double *px1 = REAL(x1);
			CRSPARSE_SKEWPART_GENERAL(
				px1[pp1_[j]]         = (px0[k]   == NA_INTEGER
				                        ? NA_REAL :  0.5 * (double) px0[k]),
				px1[pp1_[j]]         = (px0_[k_] == NA_INTEGER
				                        ? NA_REAL : -0.5 * (double) px0_[k_]),
				px1[pp1_[j]]        -= (px0_[k_] == NA_INTEGER
				                        ? NA_REAL :  0.5 * (double) px0_[k_]),
				px1[pp1_[pi0[k]]]    = -px1[pp1_[j]],
				px1[pp1_[pi0_[k_]]]  = -px1[pp1_[j]]);
			break;
		}
		case 'd':
		{
			double *px0 = REAL(x0), *px0_ = REAL(x0_),
				*px1 = REAL(x1);
			CRSPARSE_SKEWPART_GENERAL(
				px1[pp1_[j]]         =  0.5 * px0[k],
				px1[pp1_[j]]         = -0.5 * px0_[k_],
				px1[pp1_[j]]        -=  0.5 * px0_[k_],
				px1[pp1_[pi0[k]]]    = -px1[pp1_[j]],
				px1[pp1_[pi0_[k_]]]  = -px1[pp1_[j]]);
			break;
		}
		case 'z':
		{
			Rcomplex *px0 = COMPLEX(x0), *px0_ = COMPLEX(x0_),
				*px1 = COMPLEX(x1);
			CRSPARSE_SKEWPART_GENERAL(
				do {
					px1[pp1_[j]].r          =  0.5 * px0[k].r;
					px1[pp1_[j]].i          =  0.5 * px0[k].i;
				} while (0),
				do {
					px1[pp1_[j]].r          = -0.5 * px0_[k_].r;
					px1[pp1_[j]].i          = -0.5 * px0_[k_].i;
				} while (0),
				do {
					px1[pp1_[j]].r         -=  0.5 * px0_[k_].r;
					px1[pp1_[j]].i         -=  0.5 * px0_[k_].i;
				} while (0),
				do {
					px1[pp1_[pi0[k]]].r     = -px1[pp1_[j]].r;
					px1[pp1_[pi0[k]]].i     = -px1[pp1_[j]].i;
				} while (0),
				do {
					px1[pp1_[pi0_[k_]]].r   = -px1[pp1_[j]].r;
					px1[pp1_[pi0_[k_]]].i   = -px1[pp1_[j]].i;
				} while (0));
			break;
		}
		default:
			break;
		}

#undef CRSPARSE_SKEWPART_GENERAL

		Matrix_Free(pp1_, n);
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to,        iSym, i1);
		SET_SLOT(to, Matrix_xSym, x1);

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		nprotect += 2;
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = nnz0;

		char di = 'N';
		if (clf[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
		}
		if (di == 'N')
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] == pj0[k])
					--nnz1;
		nnz1 *= 2;

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		nprotect += 2;
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);

		PROTECT(x1 = allocVector(kind2type(clt[0]), nnz1));
		++nprotect;

#define TSPARSE_SKEWPART_GENERAL(_DO_ASSIGN_) \
		do { \
			for (k = 0; k < nnz0; ++k) { \
				if (pi0[k] != pj0[k]) { \
					*(pi1++) = pi0[k]; \
					*(pj1++) = pj0[k]; \
					*(pi1++) = pj0[k]; \
					*(pj1++) = pi0[k]; \
					_DO_ASSIGN_; \
				} \
			} \
		} while (0)

		switch (clf[0]) {
		case 'n':
		{
			double *px1 = REAL(x1);
			TSPARSE_SKEWPART_GENERAL(
				do {
					*(px1++) =  0.5;
					*(px1++) = -0.5;
				} while (0));
			break;
		}
		case 'l':
		{
			int *px0 = LOGICAL(x0);
			double *px1 = REAL(x1);
			TSPARSE_SKEWPART_GENERAL(
				do {
					if (px0[k] == NA_LOGICAL) {
						*(px1++) = NA_REAL;
						*(px1++) = NA_REAL;
					} else if (px0[k] != 0) {
						*(px1++) =  0.5;
						*(px1++) = -0.5;
					} else {
						*(px1++) = 0.0;
						*(px1++) = 0.0;
					}
				} while (0));
			break;
		}
		case 'i':
		{
			int *px0 = INTEGER(x0);
			double *px1 = REAL(x1);
			TSPARSE_SKEWPART_GENERAL(
				do {
					if (px0[k] == NA_INTEGER) {
						*(px1++) = NA_REAL;
						*(px1++) = NA_REAL;
					} else {
						*(px1++) =  0.5 * (double) px0[k];
						*(px1++) = -0.5 * (double) px0[k];
					}
				} while (0));
			break;
		}
		case 'd':
		{
			double *px0 = REAL(x0), *px1 = REAL(x1);
			TSPARSE_SKEWPART_GENERAL(
				do {
					*(px1++) =  0.5 * px0[k];
					*(px1++) = -0.5 * px0[k];
				} while (0));
			break;
		}
		case 'z':
		{
			Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
			TSPARSE_SKEWPART_GENERAL(
				do {
					(*(px1  )).r =  0.5 * px0[k].r;
					(*(px1++)).i =  0.5 * px0[k].i;
					(*(px1  )).r = -0.5 * px0[k].r;
					(*(px1++)).i = -0.5 * px0[k].i;
				} while (0));
			break;
		}
		default:
			break;
		}

#undef TSPARSE_SKEWPART_GENERAL

		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		SET_SLOT(to, Matrix_xSym, x1);

	}

	UNPROTECT(nprotect);
	return to;
}

/* as(<[CR]sparseMatrix>, "TsparseMatrix") */
SEXP CRsparse_as_Tsparse(SEXP from)
{
	static const char *valid[] = { VALID_CSPARSE, VALID_RSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "CRsparse_as_Tsparse");
	const char *clf = valid[ivalid];

	char clt[] = "..TMatrix";
	clt[0] = clf[0];
	clt[1] = clf[1];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (clf[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (clf[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorSym, factors);
		UNPROTECT(1); /* factors */
	}

	SEXP iSym, jSym, p0 = PROTECT(GET_SLOT(from, Matrix_pSym));
	int n_, *pp0 = INTEGER(p0);
	if (clf[2] == 'C') {
		iSym = Matrix_iSym;
		jSym = Matrix_jSym;
		n_ = n;
	} else {
		iSym = Matrix_jSym;
		jSym = Matrix_iSym;
		n_ = m;
	}

	SEXP i0 = PROTECT(GET_SLOT(from, iSym));
	SET_SLOT(to, iSym, i0);
	UNPROTECT(1); /* i0 */

	SEXP j0 = PROTECT(allocVector(INTSXP, pp0[n_]));
	int j, k, kend, *pj0 = INTEGER(j0);
	for (j = 0, k = 0; j < n_; ++j) {
		kend = *(++pp0);
		while (k < kend) {
			*(pj0++) = j;
			++k;
		}
	}
	SET_SLOT(to, jSym, j0);
	UNPROTECT(2); /* j0, p0 */

	if (clf[0] != 'n') {
		SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	}

	UNPROTECT(1); /* to */
	return to;
}

SEXP Tsparse_as_CRsparse(SEXP from, SEXP Csparse)
{
	static const char *valid[] = { VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "Tsparse_as_CRsparse");
	const char *clf = valid[ivalid];

	int doC = asLogical(Csparse) != 0;

	char clt[] = "...Matrix";
	clt[0] = clf[0];
	clt[1] = clf[1];
	clt[2] = (doC) ? 'C' : 'R';
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
	++nprotect;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (clf[1] != 'g') {
	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	if (ul != 'U')
		SET_SLOT(to, Matrix_uploSym, uplo);
	UNPROTECT(1); /* uplo */
	}
	if (clf[1] == 't') {
	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (di != 'N')
	   SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */
	} else {
	SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
	if (LENGTH(factors) > 0)
		SET_SLOT(to, Matrix_factorSym, factors);
	UNPROTECT(1); /* factors */
	}

	SEXP iSym, jSym;
	int m_, n_;
	if (doC) {
		iSym = Matrix_iSym;
		jSym = Matrix_jSym;
		m_ = m;
		n_ = n;
	} else {
		iSym = Matrix_jSym;
		jSym = Matrix_iSym;
		m_ = n;
		n_ = m;
	}

	SEXP i0 = PROTECT(GET_SLOT(from, iSym)), j0 = PROTECT(GET_SLOT(from, jSym));
	nprotect += 2;
	int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), r_ = (m_ < n_) ? n_ : m_;
	R_xlen_t nnz0 = XLENGTH(i0), nnz1 = 0;

	SEXP x0 = NULL, x1 = NULL;
	if (clf[0] != 'n') {
		PROTECT(x0 = GET_SLOT(from, Matrix_xSym));
		++nprotect;
	}

	/* FIXME? we would ideally only throw an error if the number
	   of _unique_ (i,j) pairs exceeds INT_MAX ...
	*/
	if (nnz0 > INT_MAX)
		error(_("unable to coerce from TsparseMatrix to [CR]sparseMatrix"
		        "when length of 'i' slot exceeds 2^31-1"));

	SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n_ + 1)), i1 = NULL;
	++nprotect;
	int *pp1 = INTEGER(p1), *pi1, *pj_, *workA, *workB, *workC, i, j;
	R_xlen_t k, kstart, kend, kend_, w = (R_xlen_t) m_ + r_ + m_;
	*(pp1++) = 0;
	Matrix_Calloc(pj_, nnz0, int);
	Matrix_Calloc(workA, w, int);
	workB = workA + m_;
	workC = workB + r_;

	/* 1. Tabulate column indices in workA[i]

	      workA[.]: unused
	      workB[.]: unused
	      workC[.]: unused
	*/

#define T_AS_CR_1 \
	do { \
		for (k = 0; k < nnz0; ++k) \
			++workA[pi0[k]]; \
	} while (0)

	/* 2. Compute cumulative sum in workA[i], copying to workB[i]

	      workA[i]: number of column indices listed for row i,
	                incl. duplicates
	*/

#define T_AS_CR_2 \
	do { \
		if (r_ > 0) \
			for (i = 1; i < m_; ++i) \
				workA[i] += (workB[i] = workA[i-1]); \
	} while (0)

	/* 3. Group column indices and data by row in pj_[k], px_[k]

	      workA[i]: number of column indices listed for row <= i,
	                incl. duplicates
	      workB[i]: number of column indices listed for row <  i,
	                incl. duplicates
	*/

#define T_AS_CR_3(_XASSIGN_) \
	do { \
		for (k = 0; k < nnz0; ++k) { \
			pj_[workB[pi0[k]]] = pj0[k]; \
			_XASSIGN_; /* px_[workB[pi0[k]]] = px0[k]; */ \
			++workB[pi0[k]]; \
		} \
	} while (0)

	/* 4. Gather _unique_ column indices at the front of each group,
	      aggregating data accordingly; record in workC[i] where the
	      unique column indices stop and the duplicates begin

	  workB[.]: unused
	    pj_[k]: column indices grouped by row, incl. duplicates, unsorted
	    px_[k]: corresponding data
	*/

#define T_AS_CR_4(_XASSIGN_, _XINCR_) \
	do { \
		k = 0; \
		for (j = 0; j < n_; ++j) \
			workB[j] = -1; \
		for (i = 0; i < m_; ++i) { \
			kstart = k; \
			kend_ = k; \
			kend = workA[i]; \
			while (k < kend) { \
				if (workB[pj_[k]] < kstart) { \
					/* Have not yet seen this column index */ \
					workB[pj_[k]] = kend_; \
					pj_[kend_] = pj_[k]; \
					_XASSIGN_; /* px_[kend_] = px_[k]; */ \
					++kend_; \
				} else { \
					/* Have already seen this column index */ \
					_XINCR_; /* px_[workB[pj_[k]]] += px_[k]; */ \
				} \
				++k; \
			} \
			workC[i] = kend_; \
			nnz1 += kend_ - kstart; \
		} \
	} while (0)

	/* 5. Tabulate _unique_ column indices in workB[j]

	      workC[i]: pointer to first non-unique column index in row i
	        pi_[k]: column indices grouped by row, with unique indices in front
	                i.e., in positions workA[i-1] <= k < workC[i]
	        px_[k]: corresponding data, "cumulated" appropriately
	*/

#define T_AS_CR_5 \
	do { \
		k = 0; \
		Matrix_memset(workB, 0, n_, sizeof(int)); \
		for (i = 0; i < m_; ++i) { \
			kend_ = workC[i]; \
			while (k < kend_) { \
				++workB[pj_[k]]; \
				++k; \
			} \
			k = workA[i]; \
		} \
	} while (0)

	/* 6. Compute cumulative sum in pp1[j], copying to workB[j]

	      workB[j]: number of nonzero elements in column j
	*/

#define T_AS_CR_6 \
	do { \
		for (j = 0; j < n_; ++j) { \
			pp1[j] = pp1[j-1] + workB[j]; \
			workB[j] = pp1[j-1]; \
		} \
	} while (0)

	/* 7. Pop unique (i,j) pairs from the unsorted stacks 0 <= i < m
	      onto new stacks 0 <= j < n, which will be sorted

	  workB[j]: number of nonzero elements in columns <  j
	    pp1[j]: number of nonzero elements in columns <= j
	*/

#define T_AS_CR_7(_XASSIGN_) \
	do { \
		k = 0; \
		for (i = 0; i < m_; ++i) { \
			kend_ = workC[i]; \
			while (k < kend_) { \
				pi1[workB[pj_[k]]] = i; \
				_XASSIGN_; /* px1[workB[pj_[k]]] = px_[k]; */ \
				++workB[pj_[k]]; \
				++k; \
			} \
			k = workA[i]; \
		} \
	} while (0)

#define T_AS_CR_N \
	do { \
		T_AS_CR_1; \
		T_AS_CR_2; \
		T_AS_CR_3(); \
		T_AS_CR_4(, ); \
		T_AS_CR_5; \
		T_AS_CR_6; \
		PROTECT(i1 = allocVector(INTSXP, nnz1)); \
		++nprotect; \
		pi1 = INTEGER(i1); \
		T_AS_CR_7(); \
	} while (0)

#define T_AS_CR_X(_CTYPE_, _PTR_, _SEXPTYPE_, _XINCR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1, *px_; \
		Matrix_Calloc(px_, nnz0, _CTYPE_); \
		T_AS_CR_1; \
		T_AS_CR_2; \
		T_AS_CR_3(px_[workB[pi0[k]]] = px0[k]); \
		T_AS_CR_4(px_[kend_] = px_[k], _XINCR_); \
		T_AS_CR_5; \
		T_AS_CR_6; \
		PROTECT(i1 = allocVector(INTSXP, nnz1)); \
		PROTECT(x1 = allocVector(_SEXPTYPE_, nnz1)); \
		nprotect += 2; \
		pi1 = INTEGER(i1); \
		px1 = _PTR_(x1); \
		T_AS_CR_7(px1[workB[pj_[k]]] = px_[k]); \
		Matrix_Free(px_, nnz0); \
	} while (0)

#define T_AS_CR_CASES(_KIND_, _DO_N_, _DO_X_) \
	do { \
		switch (_KIND_) { \
		case 'n': \
			_DO_N_; \
			break; \
		case 'l': \
			_DO_X_(int, LOGICAL, LGLSXP, \
				do { \
					if (px_[k] != 0) { \
						if (px_[k] != NA_LOGICAL) \
							px_[workB[pj_[k]]] = 1; \
						else if (px_[workB[pj_[k]]] == 0) \
							px_[workB[pj_[k]]] = NA_LOGICAL; \
					} \
				} while (0)); \
			break; \
		case 'i': \
			_DO_X_(int, INTEGER, INTSXP, \
			       /* FIXME: not detecting integer overflow here */ \
			       px_[workB[pj_[k]]] += px_[k]); \
			break; \
		case 'd': \
			_DO_X_(double, REAL, REALSXP, \
			       px_[workB[pj_[k]]] += px_[k]); \
			break; \
		case 'z': \
			_DO_X_(Rcomplex, COMPLEX, CPLXSXP, \
				do { \
					px_[workB[pj_[k]]].r += px_[k].r; \
					px_[workB[pj_[k]]].i += px_[k].i; \
				} while (0)); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	T_AS_CR_CASES(clf[0], T_AS_CR_N, T_AS_CR_X);
	Matrix_Free(workA, w);
	Matrix_Free(pj_, nnz0);

	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to,        iSym, i1);
	if (clf[0] != 'n')
		SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(nprotect);
	return to;
}

SEXP Tsparse_aggregate(SEXP from)
{
	static const char *valid[] = { VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "Tsparse_aggregate");
	const char *cl = valid[ivalid];

	/* Need to behave as Tsparse_as_CRsparse(from, FALSE)
	   in order to get aggregated triplets sorted by column */

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
		i0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
		j0 = PROTECT(GET_SLOT(from, Matrix_iSym));
	nprotect += 3;
	int *pdim = INTEGER(dim), *pi0 = INTEGER(i0), *pj0 = INTEGER(j0),
		m_ = pdim[1], n_ = pdim[0], r_ = (m_ < n_) ? n_ : m_;
	R_xlen_t nnz0 = XLENGTH(i0), nnz1 = 0;

	SEXP x0 = NULL, x1 = NULL;
	if (cl[0] != 'n') {
		PROTECT(x0 = GET_SLOT(from, Matrix_xSym));
		++nprotect;
	}

	if (nnz0 > INT_MAX)
		error(_("unable to aggregate TsparseMatrix with 'i' slot "
		        "of length exceeding 2^31-1"));

	SEXP i1 = NULL, j1 = NULL;
	int *pi1, *pj1, *pj_, *workA, *workB, *workC, i, j;
	R_xlen_t k, kstart, kend, kend_, w = (R_xlen_t) m_ + r_ + m_;
	Matrix_Calloc(pj_, nnz0, int);
	Matrix_Calloc(workA, w, int);
	workB = workA + m_;
	workC = workB + r_;

#define SET_TRIPLET(_XASSIGN_) \
	do { \
		k = 0; \
		for (i = 0; i < m_; ++i) { \
			kend_ = workC[i]; \
			while (k < kend_) { \
				*(pi1++) = i; \
				*(pj1++) = pj_[k]; \
				_XASSIGN_; /* *(px1++) = px_[k]; */ \
				++k; \
			} \
			k = workA[i]; \
		} \
	} while (0)

#define T_AGGR_N \
	do { \
		T_AS_CR_1; \
		T_AS_CR_2; \
		T_AS_CR_3(); \
		T_AS_CR_4(, ); \
		if (nnz1 != nnz0) { \
			PROTECT(i1 = allocVector(INTSXP, nnz1)); \
			PROTECT(j1 = allocVector(INTSXP, nnz1)); \
			nprotect += 2; \
			pi1 = INTEGER(i1); \
			pj1 = INTEGER(j1); \
			SET_TRIPLET(); \
		} \
		Matrix_Free(workA, w); \
		Matrix_Free(pj_, nnz0); \
		if (nnz1 == nnz0) { \
			UNPROTECT(nprotect); \
			return from; \
		} \
	} while (0)

#define T_AGGR_X(_CTYPE_, _PTR_, _SEXPTYPE_, _XINCR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1, *px_; \
		Matrix_Calloc(px_, nnz0, _CTYPE_); \
		T_AS_CR_1; \
		T_AS_CR_2; \
		T_AS_CR_3(px_[workB[pi0[k]]] = px0[k]); \
		T_AS_CR_4(px_[kend_] = px_[k], _XINCR_); \
		if (nnz1 != nnz0) { \
			PROTECT(i1 = allocVector(INTSXP, nnz1)); \
			PROTECT(j1 = allocVector(INTSXP, nnz1)); \
			PROTECT(x1 = allocVector(_SEXPTYPE_, nnz1)); \
			nprotect += 3; \
			pi1 = INTEGER(i1); \
			pj1 = INTEGER(j1); \
			px1 = _PTR_(x1); \
			SET_TRIPLET(*(px1++) = px_[k]); \
		} \
		Matrix_Free(workA, w); \
		Matrix_Free(pj_, nnz0); \
		Matrix_Free(px_, nnz0); \
		if (nnz1 == nnz0) { \
			UNPROTECT(nprotect); \
			return from; \
		} \
	} while (0)

	T_AS_CR_CASES(cl[0], T_AGGR_N, T_AGGR_X);

	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));
	++nprotect;

	if (m_ != n_ || n_ > 0)
		SET_SLOT(to, Matrix_DimSym, dim);

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (cl[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (cl[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorSym, factors);
		UNPROTECT(1); /* factors */
	}

	SET_SLOT(to, Matrix_iSym, j1);
	SET_SLOT(to, Matrix_jSym, i1);
	if (cl[0] != 'n')
		SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(nprotect);
	return to;
}

/* as(t(<[CR]sparseMatrix>), "[RC]sparseMatrix") */
SEXP tCRsparse_as_RCsparse(SEXP from)
{
	static const char *valid[] = { VALID_CSPARSE, VALID_RSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, "tCRsparse_as_RCsparse");
	const char *clf = valid[ivalid];

	char clt[] = "...Matrix";
	clt[0] = clf[0];
	clt[1] = clf[1];
	clt[2] = (clf[2] == 'C') ? 'R' : 'C';
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n)	{
		UNPROTECT(1); /* dim */
		PROTECT(dim = GET_SLOT(to, Matrix_DimSym));
		pdim = INTEGER(dim);
		pdim[0] = n;
		pdim[1] = m;
	} else if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (clf[1] != 's')
		set_reversed_DimNames(to, dimnames);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (clf[1] != 'g') {
		SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ulf = *CHAR(STRING_ELT(uplo_from, 0));
		if (ulf == 'U') {
			SEXP uplo_to = PROTECT(mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo_to);
			UNPROTECT(1); /* uplo_to */
		}
		UNPROTECT(1); /* uplo_from */

		if (clf[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			char di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		} else {
			SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
			if (LENGTH(factors) > 0)
				SET_SLOT(to, Matrix_factorSym, factors);
			UNPROTECT(1); /* factors */
		}
	}

	SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
	SET_SLOT(to, Matrix_pSym, p);
	UNPROTECT(1); /* p */
	if (clf[2] == 'C') {
		SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
		SET_SLOT(to, Matrix_jSym, i);
		UNPROTECT(1); /* i */
	} else {
		SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
		SET_SLOT(to, Matrix_iSym, j);
		UNPROTECT(1); /* j */
	}
	if (clf[0] != 'n') {
		SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	}

	UNPROTECT(1); /* to */
	return to;
}

/* isDiagonal(<[CR]sparseMatrix>) */
#define CR_IS_DIAGONAL(_C_, _I_) \
SEXP _C_ ## sparse_is_diagonal(SEXP obj) \
{ \
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	UNPROTECT(1); /* dim */ \
	if (m != n) \
		return ScalarLogical(0); \
	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)); \
	int *pp = INTEGER(p); \
	if (pp[n] > n) { \
		UNPROTECT(1); /* p */ \
		return ScalarLogical(0); \
	} \
	SEXP i = PROTECT(GET_SLOT(obj, Matrix_ ## _I_ ## Sym)); \
	int d, j, *pi = INTEGER(i); \
	Rboolean res = TRUE; \
	for (j = 0; j < n; ++j) { \
		if ((d = pp[j+1] - pp[j]) > 1 || (d == 1 && *(pi++) != j)) { \
			res = FALSE; \
			break; \
		} \
	} \
	UNPROTECT(2); /* i, p */ \
	return ScalarLogical(res); \
}

/* Csparse_is_diagonal() */
CR_IS_DIAGONAL(C, i)
/* Rsparse_is_diagonal() */
CR_IS_DIAGONAL(R, j)

/* isDiagonal(<TsparseMatrix>) */
SEXP Tsparse_is_diagonal(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */
	if (m != n)
		return ScalarLogical(0);
	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	int *pi = INTEGER(i), *pj = INTEGER(j);
	R_xlen_t k, nnz = XLENGTH(i);
	Rboolean res = TRUE;
	for (k = 0; k < nnz; ++k) {
		if (*(pi++) != *(pj++)) {
			res = FALSE;
			break;
		}
	}
	UNPROTECT(2); /* j, i */
	return ScalarLogical(res);
}

/* isTriangular(<.g[CR]Matrix>, upper) */
#define CR_IS_TRIANGULAR(_C_, _I_, _UPPER_, _LOWER_) \
SEXP _C_ ## sparse_is_triangular(SEXP obj, SEXP upper) \
{ \
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	UNPROTECT(1); /* dim */ \
	if (m != n) \
		return ScalarLogical(0); \
	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)), \
		i = PROTECT(GET_SLOT(obj, Matrix_ ## _I_ ## Sym)); \
	int j, k, kend, *pp = INTEGER(p), *pi = INTEGER(i), \
		need_upper = asLogical(upper); \
	Rboolean res = TRUE; \
	++pp; \
	if (need_upper == NA_LOGICAL) { \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if (_LOWER_) \
					goto opposite; \
				++k; \
			} \
		} \
		UNPROTECT(2); /* i, p */ \
		RETURN_TRUE_OF_KIND("U"); \
	opposite: \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if (_UPPER_) { \
					res = FALSE; \
					goto nokind; \
				} \
				++k; \
			} \
		} \
		UNPROTECT(2); /* i, p */ \
		RETURN_TRUE_OF_KIND("L"); \
	} else if (need_upper != 0) { \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if (_LOWER_) { \
					res = FALSE; \
					goto nokind; \
				} \
				++k; \
			} \
		} \
	} else { \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if (_UPPER_) { \
					res = FALSE; \
					goto nokind; \
				} \
				++k; \
			} \
		} \
	} \
nokind: \
	UNPROTECT(2); /* i, p */ \
	return ScalarLogical(res); \
}

/* Csparse_is_triangular() */
CR_IS_TRIANGULAR(C, i, pi[k] < j, pi[k] > j)
/* Rsparse_is_triangular() */
CR_IS_TRIANGULAR(R, j, pi[k] > j, pi[k] < j)

/* isTriangular(<.gTMatrix>, upper) */
SEXP Tsparse_is_triangular(SEXP obj, SEXP upper)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */
	if (m != n)
		return ScalarLogical(0);
	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	int *pi = INTEGER(i), *pj = INTEGER(j), need_upper = asLogical(upper);
	R_xlen_t k, nnz = XLENGTH(i);
	Rboolean res = TRUE;
	if (need_upper == NA_LOGICAL) {
		for (k = 0; k < nnz; ++k)
			if (pi[k] > pj[k])
				goto opposite;
		UNPROTECT(2); /* j, i */
		RETURN_TRUE_OF_KIND("U");
	opposite:
		for (k = 0; k < nnz; ++k)
			if (pi[k] < pj[k]) {
				res = FALSE;
				goto nokind;
			}
		UNPROTECT(2); /* j, i */
		RETURN_TRUE_OF_KIND("L");
	} else if (need_upper != 0) {
		for (k = 0; k < nnz; ++k)
			if (pi[k] > pj[k]) {
				res = FALSE;
				goto nokind;
			}
	} else {
		for (k = 0; k < nnz; ++k)
			if (pi[k] < pj[k]) {
				res = FALSE;
				goto nokind;
			}
	}
nokind:
	UNPROTECT(2); /* j, i */
	return ScalarLogical(res);
}

#define CR_IS_SYMMETRIC_LOOP(_XCOND_) \
	do { \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				if ((i = pi[k]) >= j) { \
					if (i == j) \
						++pp_[j]; \
					k = kend; \
					break; \
				} \
				if (pp_[i] == pp[i] || pi[pp_[i]] != j || (_XCOND_)) { \
					res = FALSE; \
					goto finish; \
				} \
				++pp_[i]; \
				++pp_[j]; \
				++k; \
			} \
		} \
	} while (0)

/* isSymmetric(<.g[CR]Matrix>, tol = 0, checkDN) */
#define CR_IS_SYMMETRIC(_C_, _I_) \
SEXP _C_ ## sparse_is_symmetric(SEXP obj, SEXP checkDN) \
{ \
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), n = pdim[0], s = pdim[1] == n; \
	UNPROTECT(1); /* dim */ \
	if (!s) \
		return ScalarLogical(0); \
	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)); \
	s = asLogical(checkDN) == 0 || DimNames_is_symmetric(dimnames); \
	UNPROTECT(1); /* dimnames */ \
	if (!s) \
		return ScalarLogical(0); \
	SEXP p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)), \
		i0 = PROTECT(GET_SLOT(obj, Matrix_ ## _I_ ## Sym)); \
	int i, j, k, kend, *pp_, *pp = INTEGER(p0), *pi = INTEGER(i0), \
		nprotect = 2; \
	Rboolean res = TRUE; \
	Matrix_Calloc(pp_, n, int); \
	Matrix_memcpy(pp_, pp, n, sizeof(int)); \
	++pp; \
	/* For all X[i,j] in "leading" triangle, */ \
	/* need that X[j,i] exists and X[j,i] == X[i,j] */ \
	if (!HAS_SLOT(obj, Matrix_xSym)) \
		CR_IS_SYMMETRIC_LOOP(0); \
	else { \
		SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)); \
		++nprotect; \
		switch (TYPEOF(x0)) { \
		case LGLSXP: \
		{ \
			int *px = LOGICAL(x0); \
			CR_IS_SYMMETRIC_LOOP( \
				px[pp_[i]] == NA_LOGICAL \
				? (px[k] != NA_LOGICAL) \
				: (px[k] == NA_LOGICAL || px[pp_[i]] != px[k])); \
			break; \
		} \
		case INTSXP: \
		{ \
			int *px = INTEGER(x0); \
			CR_IS_SYMMETRIC_LOOP( \
				px[pp_[i]] == NA_INTEGER \
				? (px[k] != NA_INTEGER) \
				: (px[k] == NA_INTEGER || px[pp_[i]] != px[k])); \
			break; \
		} \
		case REALSXP: \
		{ \
			double *px = REAL(x0); \
			CR_IS_SYMMETRIC_LOOP( \
				ISNAN(px[pp_[i]]) \
				? !ISNAN(px[k]) \
				: (ISNAN(px[k]) || px[pp_[i]] != px[k])); \
			break; \
		} \
		case CPLXSXP: \
		{ \
			Rcomplex *px = COMPLEX(x0); \
			CR_IS_SYMMETRIC_LOOP( \
				ISNAN(px[pp_[i]].r) || ISNAN(px[pp_[i]].i) \
				? !(ISNAN(px[k].r) || ISNAN(px[k].i)) \
				: (ISNAN(px[k].r) || ISNAN(px[k].i) || \
				   px[pp_[i]].r != px[k].r || px[pp_[i]].i != px[k].i)); \
			break; \
		} \
		default: \
			ERROR_INVALID_TYPE( \
				"'x' slot", TYPEOF(x0), "[CR]sparse_is_symmetric"); \
			break; \
		} \
	} \
	/* Need upper, lower triangles to have same number of nonzero elements */ \
	for (j = 0; j < n; ++j) { \
		if (pp_[j] != pp[j]) { \
			res = FALSE; \
			goto finish; \
		} \
	} \
finish: \
	Matrix_Free(pp_, n); \
	UNPROTECT(nprotect); /* x0, i0, p0 */ \
	return ScalarLogical(res); \
}

/* Csparse_is_symmetric() */
/* FIXME: not checking for real diagonal in complex case */
CR_IS_SYMMETRIC(C, i)
/* Rsparse_is_symmetric() */
/* FIXME: not checking for real diagonal in complex case */
CR_IS_SYMMETRIC(R, j)

/* colSums(<CsparseMatrix>), rowSums(<RsparseMatrix>) */
SEXP CRsparse_colSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse)
{
	static const char *valid[] = { VALID_CSPARSE, VALID_RSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, "CRsparse_colSums");
	const char *cl = valid[ivalid];
	if (cl[1] == 's')
		return CRsparse_rowSums(obj, narm, mean, sparse);

	int doSparse = asLogical(sparse) != 0,
		doNaRm = asLogical(narm) != 0,
		doMean = asLogical(mean) != 0,
		doCount = doNaRm && doMean;

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int margin = (cl[2] == 'C') ? 1 : 0,
		*pdim = INTEGER(dim), m = pdim[!margin], n = pdim[margin];
	UNPROTECT(1); /* dim */

	char di = 'N';
	if (cl[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */
		if (doSparse && di != 'N')
			warning(_("sparseResult=TRUE inefficient for unit triangular 'x'"));
	}

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	++nprotect;
	int *pp = INTEGER(p) + 1, j, k, kend, count = m;

	PROTECT_INDEX pid;
	SEXP res, vl = NULL, vi = NULL, vx = NULL;
	int *pvi = NULL;

	if (!doSparse) {
		PROTECT_WITH_INDEX(
			res = allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, n), &pid);
		++nprotect;
	} else {
		int nnz = n;
		if (di == 'N') {
			nnz = 0;
			for (j = 0; j < n; ++j)
				if (pp[j-1] < pp[j])
					++nnz;
		}

		char cl_[] = ".sparseVector";
		cl_[0] = (((cl[0] == 'n' || cl[0] == 'l') && !doMean)
		          ? 'i' : ((cl[0] != 'z') ? 'd' : 'z'));
		PROTECT(res = NEW_OBJECT_OF_CLASS(cl_));
		PROTECT(vl = ScalarInteger(n));
		PROTECT(vi = allocVector(INTSXP, nnz));
		PROTECT_WITH_INDEX(
			vx = allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, nnz), &pid);
		nprotect += 4;
		pvi = INTEGER(vi);
	}

	if (cl[0] == 'n') {
		double *pres = (doSparse) ? REAL(vx) : REAL(res);
		if (!doSparse) {
			int u = (di == 'N') ? 0 : 1;
			for (j = 0; j < n; ++j) {
				*pres = pp[j] - pp[j-1] + u;
				if (doMean)
					*pres /= count;
				++pres;
			}
		} else if (di == 'N') {
			for (j = 0; j < n; ++j) {
				if (pp[j-1] < pp[j]) {
					*pvi = j + 1;
					*pres = pp[j] - pp[j-1];
					if (doMean)
						*pres /= count;
					++pvi;
					++pres;
				}
			}
		} else {
			for (j = 0; j < n; ++j) {
				*pvi = j + 1;
				*pres = pp[j] - pp[j-1] + 1;
				if (doMean)
					*pres /= count;
				++pvi;
				++pres;
			}
		}
	} else {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));

#define CR_COLSUMS_LOOP \
		do { \
			k = 0; \
			if (!doSparse) { \
				if (di == 'N') { \
					for (j = 0; j < n; ++j) { \
						kend = pp[j]; \
						DO_INIT(ZERO); \
						while (k < kend) { DO_INCR; ++k; } \
						DO_SCALE; \
						++pres; \
					} \
				} else { \
					for (j = 0; j < n; ++j) { \
						kend = pp[j]; \
						DO_INIT(ONE); \
						while (k < kend) { DO_INCR; ++k; } \
						DO_SCALE; \
						++pres; \
					} \
				} \
			} else { \
				if (di == 'N') { \
					for (j = 0; j < n; ++j) { \
						kend = pp[j]; \
						if (k < kend) { \
							*pvi = j + 1; \
							DO_INIT(ZERO); \
							while (k < kend) { DO_INCR; ++k; } \
							DO_SCALE; \
							++pvi; \
							++pres; \
						} \
					} \
				} else { \
					for (j = 0; j < n; ++j) { \
						kend = pp[j]; \
						*pvi = j + 1; \
						DO_INIT(ONE); \
						while (k < kend) { DO_INCR; ++k; } \
						DO_SCALE; \
						++pvi; \
						++pres; \
					} \
				} \
			} \
		} while (0)

#define CR_COLSUMS(_CTYPE1_, _PTR1_, _CTYPE2_, _PTR2_) \
		do { \
			_CTYPE1_ *pres = (doSparse) ? _PTR1_(vx) : _PTR1_(res); \
			_CTYPE2_ *px   = _PTR2_(x); \
			CR_COLSUMS_LOOP; \
		} while (0)

		switch (cl[0]) {
		case 'l':

#define ZERO         0.0
#define ONE          1.0
#define DO_INIT(_U_) \
			do { \
				*pres = _U_; \
				if (doCount) \
					count = m; \
			} while (0)
#define DO_INCR \
			do { \
				if (px[k] != NA_LOGICAL) { \
					if (px[k]) *pres += 1.0; \
				} else if (!doNaRm) \
					*pres = NA_REAL; \
				else if (doMean) \
					--count; \
			} while (0)
#define DO_SCALE     if (doMean) *pres /= count

			CR_COLSUMS(double, REAL, int, LOGICAL);
			break;

#undef DO_INCR

		case 'i':

#define DO_INCR \
			do { \
				if (px[k] != NA_INTEGER) \
					*pres += px[k]; \
				else if (!doNaRm) \
					*pres = NA_REAL; \
				else if (doMean) \
					--count; \
			} while (0)

			CR_COLSUMS(double, REAL, int, INTEGER);
			break;

#undef DO_INCR

		case 'd':

#define DO_INCR \
			do { \
				if (!(doNaRm && ISNAN(px[k]))) \
					*pres += px[k]; \
				else if (doMean) \
					--count; \
			} while (0)

			CR_COLSUMS(double, REAL, double, REAL);
			break;

#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_SCALE

		case 'z':

#define ZERO         Matrix_zzero
#define ONE          Matrix_zone
#define DO_INCR \
			do { \
				if (!(doNaRm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) { \
					(*pres).r += px[k].r; \
					(*pres).i += px[k].i; \
				} else if (doMean) \
					--count; \
			} while (0)
#define DO_SCALE \
			do { \
				if (doMean) { \
					(*pres).r /= count; \
					(*pres).i /= count; \
				} \
			} while (0)

			CR_COLSUMS(Rcomplex, COMPLEX, Rcomplex, COMPLEX);
			break;

#undef ZERO
#undef ONE
#undef DO_INIT
#undef DO_INCR
#undef DO_SCALE

		default:
			break;
		}

#undef CR_COLSUMS
#undef CR_COLSUMS_LOOP

		UNPROTECT(1); /* x */
	}

	if (doSparse) {
		if ((cl[0] == 'n' || cl[0] == 'l') && !doMean)
			REPROTECT(vx = coerceVector(vx, INTSXP), pid);

		SET_SLOT(res, Matrix_lengthSym, vl);
		SET_SLOT(res, Matrix_iSym,      vi);
		SET_SLOT(res, Matrix_xSym,      vx);
	} else {
		if ((cl[0] == 'n' || cl[0] == 'l') && !doMean)
			REPROTECT(res = coerceVector(res, INTSXP), pid);

		SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
			nms = VECTOR_ELT(dimnames, margin);
		if (!isNull(nms))
			setAttrib(res, R_NamesSymbol, nms);
		UNPROTECT(1); /* dimnames */
	}

	UNPROTECT(nprotect);
	return res;
}

/* rowSums(<CsparseMatrix>), colSums(<RsparseMatrix>) */
SEXP CRsparse_rowSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse)
{
	static const char *valid[] = { VALID_CSPARSE, VALID_RSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid), nprotect = 0;
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, "CRsparse_rowSums");
	const char *cl = valid[ivalid];

	int doSparse = asLogical(sparse) != 0,
		doNaRm = asLogical(narm) != 0,
		doMean = asLogical(mean) != 0;

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int margin = (cl[2] == 'C') ? 0 : 1,
		*pdim = INTEGER(dim), m = pdim[margin], n = pdim[!margin];
	UNPROTECT(1); /* dim */

	char ul = 'U', di = 'N';
	if (cl[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		UNPROTECT(1); /* uplo */

		if (cl[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			UNPROTECT(1); /* diag */
			if (doSparse && di != 'N')
				warning(_("sparseResult=TRUE inefficient for unit triangular 'x'"));
		}
	}

	SEXP iSym = (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym,
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, iSym));
	nprotect += 2;
	int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k, kend, *pcount = NULL;

	SEXP x = NULL;
	if (cl[0] != 'n') {
		PROTECT(x = GET_SLOT(obj, Matrix_xSym));
		++nprotect;
	}

	PROTECT_INDEX pid;
	SEXP res;
	PROTECT_WITH_INDEX(
		res = allocVector((cl[0] != 'z') ? REALSXP : CPLXSXP, m), &pid);
	++nprotect;

#define CR_ROWSUMS_LOOP \
	do { \
		k = 0; \
		if (cl[1] != 's') { \
			for (j = 0; j < n; ++j) { \
				kend = pp[j]; \
				while (k < kend) { DO_INCR; ++k; } \
			} \
		} else if (ul == ((cl[2] == 'C') ? 'U' : 'L')) { \
			for (j = 0; j < n; ++j) { \
				kend = pp[j]; \
				if (k < kend) { \
					while (kend - k > 1) { DO_INCR_SYMM; ++k; } \
					if (pi[k] == j) \
						DO_INCR; \
					else \
						DO_INCR_SYMM; \
					++k; \
				} \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				kend = pp[j]; \
				if (k < kend) { \
					if (pi[k] == j) \
						DO_INCR; \
					else \
						DO_INCR_SYMM; \
					++k; \
					while (k < kend) { DO_INCR_SYMM; ++k; } \
				} \
			} \
		} \
	} while (0)

#define CR_ROWSUMS_X(_CTYPE1_, _PTR1_, _CTYPE2_, _PTR2_) \
	do { \
		_CTYPE2_ *px = _PTR2_(x); \
		CR_ROWSUMS_N(_CTYPE1_, _PTR1_); \
	} while (0)

#define CR_ROWSUMS_N(_CTYPE1_, _PTR1_) \
	do { \
		_CTYPE1_ *pres = _PTR1_(res), u = (di == 'N') ? ZERO : ONE; \
		if (doNaRm && doMean && cl[0] != 'n') { \
			Matrix_Calloc(pcount, m, int); \
			for (k = 0; k < m; ++k) { \
				pres[k] = u; \
				pcount[k] = n; \
			} \
		} else { \
			for (k = 0; k < m; ++k) \
				pres[k] = u; \
		} \
		CR_ROWSUMS_LOOP; \
	} while (0)

	switch (cl[0]) {
	case 'n':

#define ZERO         0.0
#define ONE          1.0
#define DO_INCR      pres[pi[k]] += 1.0
#define DO_INCR_SYMM \
		do { \
			pres[pi[k]] += 1.0; \
			pres[j]     += 1.0; \
		} while (0)

		CR_ROWSUMS_N(double, REAL);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'l':

#define DO_INCR \
		do { \
			if (px[k] != NA_LOGICAL) { \
				if (px[k]) \
					pres[pi[k]] += 1.0; \
			} else if (!doNaRm) \
				pres[pi[k]] = NA_REAL; \
			else if (doMean) \
				--pcount[pi[k]]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (px[k] != NA_LOGICAL) { \
				if (px[k]) { \
					pres[pi[k]] += 1.0; \
					pres[j]     += 1.0; \
				} \
			} else if (!doNaRm) { \
				pres[pi[k]] = NA_REAL; \
				pres[j]     = NA_REAL; \
			} else if (doMean) { \
				--pcount[pi[k]]; \
				--pcount[j]; \
			} \
		} while (0)

		CR_ROWSUMS_X(double, REAL, int, LOGICAL);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'i':

#define DO_INCR \
		do { \
			if (px[k] != NA_INTEGER) \
				pres[pi[k]] += px[k]; \
			else if (!doNaRm) \
				pres[pi[k]] = NA_REAL; \
			else if (doMean) \
				--pcount[pi[k]]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (px[k] != NA_INTEGER) { \
				pres[pi[k]] += px[k]; \
				pres[j]     += px[k]; \
			} else if (!doNaRm) { \
				pres[pi[k]] = NA_REAL; \
				pres[j]     = NA_REAL; \
			} else if (doMean) { \
				--pcount[pi[k]]; \
				--pcount[j]; \
			} \
		} while (0)

		CR_ROWSUMS_X(double, REAL, int, INTEGER);
		break;

#undef DO_INCR
#undef DO_INCR_SYMM

	case 'd':

#define DO_INCR \
		do { \
			if (!(doNaRm && ISNAN(px[k]))) \
				pres[pi[k]] += px[k]; \
			else if (doMean) \
				--pcount[pi[k]]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (!(doNaRm && ISNAN(px[k]))) { \
				pres[pi[k]] += px[k]; \
				pres[j]     += px[k]; \
			} else if (doMean) { \
				--pcount[pi[k]]; \
				--pcount[j]; \
			} \
		} while (0)

		CR_ROWSUMS_X(double, REAL, double, REAL);
		break;

#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_INCR_SYMM

	case 'z':

#define ZERO         Matrix_zzero
#define ONE          Matrix_zone
#define DO_INCR \
		do { \
			if (!(doNaRm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) { \
				pres[pi[k]].r += px[k].r; \
				pres[pi[k]].i += px[k].i; \
			} else if (doMean) \
				--pcount[pi[k]]; \
		} while (0)
#define DO_INCR_SYMM \
		do { \
			if (!(doNaRm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) { \
				pres[pi[k]].r += px[k].r; \
				pres[pi[k]].i += px[k].i; \
				pres[j].r     += px[k].r; \
				pres[j].i     += px[k].i; \
			} else if (doMean) { \
				--pcount[pi[k]]; \
				--pcount[j]; \
			} \
		} while (0)

		CR_ROWSUMS_X(Rcomplex, COMPLEX, Rcomplex, COMPLEX);
		break;

#undef ZERO
#undef ONE
#undef DO_INCR
#undef DO_INCR_SYMM

	default:
		break;
	}

#undef CR_ROWSUMS
#undef CR_ROWSUMS_LOOP

	if (doMean) {
		if (cl[0] != 'z') {
			double *pres = REAL(res);
			if (doNaRm && cl[0] != 'n') {
				for (k = 0; k < m; ++k)
					pres[k] /= pcount[k];
				Matrix_Free(pcount, m);
			} else {
				for (k = 0; k < m; ++k)
					pres[k] /= n;
			}
		} else {
			Rcomplex *pres = COMPLEX(res);
			if (doNaRm) {
				for (k = 0; k < m; ++k) {
					pres[k].r /= pcount[k];
					pres[k].i /= pcount[k];
				}
				Matrix_Free(pcount, m);
			} else {
				for (k = 0; k < m; ++k) {
					pres[k].r /= n;
					pres[k].i /= n;
				}
			}
		}
	}

	if ((cl[0] == 'n' || cl[0] == 'l') && !doMean)
		REPROTECT(res = coerceVector(res, INTSXP), pid);
	if (doSparse)
		REPROTECT(res = v2spV(res), pid);
	else {
		SEXP dimnames;
		if (cl[1] != 's')
			PROTECT(dimnames = GET_SLOT(obj, Matrix_DimNamesSym));
		else
			PROTECT(dimnames = get_symmetrized_DimNames(obj, -1));
		SEXP nms = VECTOR_ELT(dimnames, margin);
		if (!isNull(nms))
			setAttrib(res, R_NamesSymbol, nms);
		UNPROTECT(1); /* dimnames */
	}

	UNPROTECT(nprotect);
	return res;
}
