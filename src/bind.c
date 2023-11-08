#include "Mdefines.h"
#include "coerce.h"
#include "bind.h"

static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };

static SEXP tagWasVector = NULL;

static
void scanArgs(SEXP args, SEXP exprs, int margin, int level,
              int *rdim, int *rdimnames, char *kind, char *repr)
{
	SEXP a, e, s, tmp;
	int nS4 = 0, nDense = 0,
		anyCsparse = 0, anyRsparse = 0, anyTsparse = 0, anyDiagonal = 0,
		anyN = 0, anyL = 0, anyI = 0, anyD = 0, anyZ = 0,
		i, ivalid, *sdim;
	R_xlen_t slen;
	const char *scl;

	rdim[!margin] = -1;
	rdim[ margin] =  0;
	rdimnames[0] = rdimnames[1] = 0;

	for (a = args; a != R_NilValue; a = CDR(a)) {
		s = CAR(a);
		if (s == R_NilValue)
			continue;
		if (TYPEOF(s) == S4SXP) {
			++nS4;
			ivalid = R_check_class_etc(s, valid);
			if (ivalid < 0) {
				if (margin)
					ERROR_INVALID_CLASS(s, "cbind.Matrix");
				else
					ERROR_INVALID_CLASS(s, "rbind.Matrix");
			}
			scl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

			tmp = GET_SLOT(s, Matrix_DimSym);
			sdim = INTEGER(tmp);
			if (rdim[!margin] < 0)
				rdim[!margin] = sdim[!margin];
			else if (sdim[!margin] != rdim[!margin]) {
				if (margin)
					error(_("number of rows of matrices must match"));
				else
					error(_("number of columns of matrices must match"));
			}
			if (sdim[margin] > INT_MAX - rdim[margin])
				error(_("dimensions cannot exceed %s"), "2^31-1");
			rdim[margin] += sdim[margin];

			if (!rdimnames[0] || !rdimnames[1]) {
				tmp = GET_SLOT(s, Matrix_DimNamesSym);
				if (scl[1] == 's') {
					if (VECTOR_ELT(tmp, 0) != R_NilValue ||
					    VECTOR_ELT(tmp, 1) != R_NilValue)
						rdimnames[0] = rdimnames[1] = 1;
				} else
					for (i = 0; i < 2; ++i)
						if (!rdimnames[i] &&
						    VECTOR_ELT(tmp, i) != R_NilValue)
							rdimnames[i] = 1;
			}

			switch (scl[0]) {
			case 'n':
				anyN = 1;
				break;
			case 'l':
				anyL = 1;
				break;
			case 'i':
				if (scl[2] != 'd')
				anyI = 1;
				break;
			case 'd':
				anyD = 1;
				break;
			case 'z':
				anyZ = 1;
				break;
			default:
				break;
			}

			switch (scl[2]) {
			case 'e':
			case 'y':
			case 'r':
			case 'p':
				++nDense;
				break;
			case 'C':
				anyCsparse = 1;
				break;
			case 'R':
				anyRsparse = 1;
				break;
			case 'T':
			{
				/* defined in ./sparse.c : */
				SEXP Tsparse_aggregate(SEXP);
				SETCAR(a, Tsparse_aggregate(s));
				anyTsparse = 1;
				break;
			}
			case 'i':
				anyDiagonal = 1;
				break;
			case 'd':
				if (INTEGER(GET_SLOT(s, Matrix_marginSym))[0] - 1 != margin) {
					anyN = 1;
					if (margin)
						anyCsparse = 1;
					else
						anyRsparse = 1;
				}
				break;
			default:
				break;
			}
		} else {
			switch (TYPEOF(s)) {
			case LGLSXP:
				anyL = 1;
				break;
			case INTSXP:
				anyI = 1;
				break;
			case REALSXP:
				anyD = 1;
				break;
			case CPLXSXP:
				anyZ = 1;
				break;
			default:
				if (margin)
					ERROR_INVALID_TYPE(s, "cbind.Matrix");
				else
					ERROR_INVALID_TYPE(s, "rbind.Matrix");
				break;
			}

			tmp = getAttrib(s, R_DimSymbol);
			if (TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2) {
				sdim = INTEGER(tmp);
				if (rdim[!margin] < 0)
					rdim[!margin] = sdim[!margin];
				else if (rdim[!margin] != sdim[!margin]) {
					if (margin)
						error(_("number of rows of matrices must match"));
					else
						error(_("number of columns of matrices must match"));
				}
				if (sdim[margin] > INT_MAX - rdim[margin])
					error(_("dimensions cannot exceed %s"), "2^31-1");
				rdim[margin] += sdim[margin];

				if (!rdimnames[0] || !rdimnames[1]) {
					tmp = getAttrib(s, R_DimNamesSymbol);
					if (tmp != R_NilValue)
						for (i = 0; i < 2; ++i)
							if (!rdimnames[i] &&
							    VECTOR_ELT(tmp, i) != R_NilValue)
								rdimnames[i] = 1;
				}
			}
		}
	}

	if (rdim[!margin] < 0) {
		/* Arguments are all vectors or NULL */
		R_xlen_t maxlen = -1;
		for (a = args; a != R_NilValue; a = CDR(a)) {
			s = CAR(a);
			if (s == R_NilValue)
				continue;
			slen = XLENGTH(s);
			if (slen > INT_MAX)
				error(_("dimensions cannot exceed %s"), "2^31-1");
			else if (slen > maxlen)
				maxlen = slen;
		}
		if (maxlen < 0)
			/* Arguments are all NULL */
			return;
		rdim[!margin] = (int) maxlen;
	}

	for (a = args, e = exprs; a != R_NilValue; a = CDR(a), e = CDR(e)) {
		s = CAR(a);
		if ((s == R_NilValue && rdim[!margin] > 0) || TYPEOF(s) == S4SXP)
			continue;
		if (s == R_NilValue)
			rdim[margin] += 1;
		else {
			tmp = getAttrib(s, R_DimSymbol);
			if (TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2)
				continue;
			slen = XLENGTH(s);
			if (slen == 0 && rdim[!margin] > 0)
				continue;
			if (rdim[margin] == INT_MAX)
				error(_("dimensions cannot exceed %s"), "2^31-1");
			rdim[margin] += 1;
			if (slen > rdim[!margin] || rdim[!margin] % (int) slen) {
				if (margin)
					warning(_("number of rows of result is not a multiple of vector length"));
				else
					warning(_("number of columns of result is not a multiple of vector length"));
			}
			if (!rdimnames[!margin] && slen == rdim[!margin]) {
				tmp = getAttrib(s, R_NamesSymbol);
				if (tmp != R_NilValue)
					rdimnames[!margin] = 1;
			}
		}
		if (!rdimnames[margin]) {
			if (TAG(a) != R_NilValue ||
			    level == 2 || (level == 1 && TYPEOF(CAR(e)) == SYMSXP))
				rdimnames[margin] = 1;
		}
	}

	if (anyZ)
		*kind = 'z';
#ifndef MATRIX_ENABLE_IMATRIX
	else if (anyD || anyI)
		*kind = 'd';
#else
	else if (anyD)
		*kind = 'd';
	else if (anyI)
		*kind = 'i';
#endif
	else if (anyL)
		*kind = 'l';
	else if (anyN)
		*kind = 'n';
	else
		*kind = '\0';

	if (nDense == nS4)
		*repr = 'e';
	else if (nDense == 0) {
		if (anyCsparse && anyRsparse)
			*repr = (margin) ? 'C' : 'R';
		else if (anyCsparse)
			*repr = 'C';
		else if (anyRsparse)
			*repr = 'R';
		else if (anyTsparse)
			*repr = 'T';
		else if (anyDiagonal)
			*repr = (margin) ? 'C' : 'R';
		else
			*repr = '\0';
	} else {
		/* The length of the result is at most INT_MAX * INT_MAX,
		   which cannot overflow Matrix_int_fast64_t as long as R
		   builds require sizeof(int) equal to 4
		 */
		Matrix_int_fast64_t nnz = 0, len = 0, snnz = 0, slen = 0;
		for (a = args; a != R_NilValue && nnz < INT_MAX; a = CDR(a)) {
			s = CAR(a);
			if (TYPEOF(s) != S4SXP)
				continue;
			ivalid = R_check_class_etc(s, valid);
			scl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

			PROTECT(tmp = GET_SLOT(s, Matrix_DimSym));
			sdim = INTEGER(tmp);
			slen = (Matrix_int_fast64_t) sdim[0] * sdim[1];

			switch (scl[2]) {
			case 'e':
			case 'y':
			case 'r':
			case 'p':
				snnz = (scl[1] != 't') ? slen : ((slen + sdim[0]) / 2);
				break;
			case 'C':
			case 'R':
			{
				SEXP p = PROTECT(GET_SLOT(s, Matrix_pSym));
				int *pp = INTEGER(p), n = sdim[(scl[2] == 'C') ? 1 : 0];
				snnz = pp[n];
				if (scl[1] == 's') {
					SEXP iSym = (scl[2] == 'C') ? Matrix_iSym : Matrix_jSym,
						i = PROTECT(GET_SLOT(s, iSym));
					int *pi = INTEGER(i), j;
					snnz *= 2;
					if (*CHAR(STRING_ELT(GET_SLOT(s, Matrix_uploSym), 0)) == 'U') {
						for (j = 0; j < n; ++j)
							if (pp[j] < pp[j + 1] && pi[pp[j + 1] - 1] == j)
								--snnz;
					} else {
						for (j = 0; j < n; ++j)
							if (pp[j] < pp[j + 1] && pi[pp[j]] == j)
								--snnz;
					}
					UNPROTECT(1);
				} else if (scl[1] == 't' && *CHAR(STRING_ELT(GET_SLOT(s, Matrix_diagSym), 0)) != 'N')
					snnz += sdim[0];
				UNPROTECT(1);
				break;
			}
			case 'T':
			{
				SEXP i = PROTECT(GET_SLOT(s, Matrix_iSym));
				snnz = XLENGTH(i);
				if (scl[1] == 's') {
					SEXP j = PROTECT(GET_SLOT(s, Matrix_jSym));
					int *pi = INTEGER(i), *pj = INTEGER(j);
					R_xlen_t k = XLENGTH(i);
					snnz *= 2;
					while (k--)
						if (*(pi++) == *(pj++))
							--snnz;
					UNPROTECT(1);
				} else if (scl[1] == 't' && *CHAR(STRING_ELT(GET_SLOT(s, Matrix_diagSym), 0)) != 'N')
					snnz += sdim[0];
				UNPROTECT(1);
				break;
			}
			case 'i':
				snnz = sdim[0];
				break;
			case 'd':
				snnz = XLENGTH(GET_SLOT(s, Matrix_permSym));
				break;
			default:
				break;
			}

			nnz += snnz;
			len += slen;
			UNPROTECT(1);
		}

		if (nnz > INT_MAX || nnz > len / 2)
			*repr = 'e';
		else if (anyCsparse && anyRsparse)
			*repr = (margin) ? 'C' : 'R';
		else if (anyCsparse)
			*repr = 'C';
		else if (anyRsparse)
			*repr = 'R';
		else if (anyTsparse)
			*repr = 'T';
		else
			*repr = (margin) ? 'C' : 'R';
	}

	return;
}

static
void coerceArgs(SEXP args, int margin,
                int *rdim, char kind, char repr)
{
	SEXP a, s, t, tmp;
	int ivalid, isM;
	char scl_[] = "...Matrix";
	const char *scl;

	for (a = args; a != R_NilValue; a = CDR(a)) {
		s = CAR(a);
		t = TAG(a);
		SET_TAG(a, R_NilValue); /* to be replaced only if 's' is a vector */
		if (s == R_NilValue)
			continue;
		PROTECT_INDEX pid;
		PROTECT_WITH_INDEX(s, &pid);
		if (TYPEOF(s) == S4SXP) {
			ivalid = R_check_class_etc(s, valid);
			scl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];
			switch (scl[2]) {
			case 'e':
			case 'y':
			case 'r':
			case 'p':
				switch (repr) {
				case 'e':
					REPROTECT(s = dense_as_kind(s, scl, kind, 0), pid);
					scl_[0] = kind; scl_[1] = scl[1]; scl_[2] = scl[2];
					REPROTECT(s = dense_as_general(
						s, scl_, kindToType(kind) == kindToType(scl[0])), pid);
					break;
				case 'C':
				case 'R':
				case 'T':
					REPROTECT(s = dense_as_sparse(s, scl, repr), pid);
					scl_[0] = scl[0]; scl_[1] = scl[1]; scl_[2] = repr;
					REPROTECT(s = sparse_as_kind(s, scl_, kind), pid);
					scl_[0] = kind;
					REPROTECT(s = sparse_as_general(s, scl_), pid);
					break;
				default:
					break;
				}
				break;
			case 'C':
			case 'R':
			case 'T':
				REPROTECT(s = sparse_as_kind(s, scl, kind), pid);
				scl_[0] = kind; scl_[1] = scl[1]; scl_[2] = scl[2];
				REPROTECT(s = sparse_as_general(s, scl_), pid);
				scl_[1] = 'g';
				switch (repr) {
				case 'e':
					REPROTECT(s = sparse_as_dense(s, scl_, 0), pid);
					break;
				case 'C':
					REPROTECT(s = sparse_as_Csparse(s, scl_), pid);
					break;
				case 'R':
					REPROTECT(s = sparse_as_Rsparse(s, scl_), pid);
					break;
				case 'T':
					REPROTECT(s = sparse_as_Tsparse(s, scl_), pid);
					break;
				default:
					break;
				}
				break;
			case 'i':
				switch (repr) {
				case 'e':
					REPROTECT(s = diagonal_as_dense(s, scl_, kind, 'g', 0, '\0'), pid);
					break;
				case 'C':
				case 'R':
				case 'T':
					REPROTECT(s = diagonal_as_sparse(s, scl_, kind, 'g', repr, '\0'), pid);
					break;
				default:
					break;
				}
				break;
			case 'd':
				switch (repr) {
				case 'e':
					REPROTECT(s = index_as_dense(s, scl, kind), pid);
					break;
				case 'C':
				case 'R':
				case 'T':
					REPROTECT(s = index_as_sparse(s, scl, kind, repr), pid);
					break;
				default:
					break;
				}
				break;
			default:
				break;
			}
		} else {
			tmp = getAttrib(s, R_DimSymbol);
			isM = TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2;
			if (!isM) {
				if (rdim[!margin] > 0 && XLENGTH(s) == 0) {
					UNPROTECT(1);
					continue;
				}
				SET_TAG(a, (t != R_NilValue) ? t : tagWasVector);
			}
			if (TYPEOF(s) != kindToType(kind))
				REPROTECT(s = coerceVector(s, kindToType(kind)), pid);
			if (repr != 'e') {
				if (!isM && XLENGTH(s) != rdim[!margin]) {
					static SEXP replen = NULL;
					if (!replen)
						replen = install("rep_len");
					SEXP lengthout = PROTECT(ScalarInteger(rdim[!margin])),
						call = PROTECT(lang3(replen, s, lengthout));
					REPROTECT(s = eval(call, R_GlobalEnv), pid);
					UNPROTECT(2);
				}
				scl_[1] = 'g';
				scl_[2] = repr;
				REPROTECT(s = matrix_as_sparse(s, scl_, '\0', '\0', !margin), pid);
			}
		}
		SETCAR(a, s);
		UNPROTECT(1);
	}

	return;
}

static
void bindArgs(SEXP args, int margin, SEXP res,
              int *rdim, char kind, char repr)
{
	SEXP a, s;

#define BIND_CASES(_BIND_) \
	do { \
		switch (kind) { \
		case 'l': \
			_BIND_(int, LOGICAL, SHOW); \
			break; \
		case 'i': \
			_BIND_(int, INTEGER, SHOW); \
			break; \
		case 'd': \
			_BIND_(double, REAL, SHOW); \
			break; \
		case 'z': \
			_BIND_(Rcomplex, COMPLEX, SHOW); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	if (repr == 'e') {

		if (rdim[0] == 0 || rdim[1] == 0)
			return;

		int k, m = rdim[0], n = rdim[1];
		R_xlen_t mn = (R_xlen_t) m * n;
		SEXP x = PROTECT(allocVector(kindToType(kind), mn)), tmp;
		SET_SLOT(res, Matrix_xSym, x);

#define BIND_E(_CTYPE_, _PTR_, _MASK_) \
		do { \
			_CTYPE_ *px = _PTR_(x), *ps; \
			for (a = args; a != R_NilValue; a = CDR(a)) { \
				s = CAR(a); \
				if (s == R_NilValue) \
					continue; \
				if (TYPEOF(s) != S4SXP) \
					tmp = getAttrib(s, R_DimSymbol); \
				else { \
					s = GET_SLOT(s, Matrix_xSym); \
					tmp = NULL; \
				} \
				mn = XLENGTH(s); \
				ps = _PTR_(s); \
				if (margin) { \
				if (!tmp || (TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2)) { \
					Matrix_memcpy(px, ps, mn, sizeof(_CTYPE_)); \
					px += mn; \
				} else if (mn >= m) { \
					Matrix_memcpy(px, ps, m , sizeof(_CTYPE_)); \
					px += m; \
				} else if (mn == 1) { \
					_CTYPE_ v = ps[0]; \
					for (k = 0; k < m; ++k) \
						*(px++) = v; \
				} else { \
					int mn_ = (int) mn; \
					for (k = 0; k < rdim[0]; ++k) \
						*(px++) = ps[k % mn_]; \
				} \
				} else { \
				_CTYPE_ *py = px; \
				if (!tmp || (TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2)) { \
					m = (int) (mn / n); \
					for (k = 0; k < n; ++k) { \
						Matrix_memcpy(py, ps, m, sizeof(_CTYPE_)); \
						py += rdim[0]; \
						ps += m; \
					} \
					px += m; \
				} else if (mn >= n) { \
					for (k = 0; k < n; ++k) { \
						*py = *ps; \
						py += rdim[0]; \
						ps += 1; \
					} \
					px += 1; \
				} else if (mn == 1) { \
					_CTYPE_ v = ps[0]; \
					for (k = 0; k < n; ++k) { \
						*py = v; \
						py += rdim[0]; \
					} \
					px += 1; \
				} else { \
					int mn_ = (int) mn; \
					for (k = 0; k < n; ++k) { \
						*py = ps[k % mn_]; \
						py += rdim[0]; \
					} \
					px += 1; \
				} \
				} \
			} \
		} while (0)

		if (kind == 'n')
			BIND_E(int, LOGICAL, SHOW);
		else
			BIND_CASES(BIND_E);
		UNPROTECT(1);

	} else if ((repr == 'C' && margin) || (repr == 'R' && !margin)) {

		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) rdim[margin] + 1));
		int *pp = INTEGER(p);
		SET_SLOT(res, Matrix_pSym, p);

		if (rdim[0] == 0 || rdim[1] == 0) {
			Matrix_memset(pp, 0, (R_xlen_t) rdim[margin] + 1, sizeof(int));
			UNPROTECT(1);
			return;
		}

		SEXP sp;
		int *psp, j, n, nnz = 0;
		*(pp++) = nnz = 0;
		for (a = args; a != R_NilValue; a = CDR(a)) {
			s = CAR(a);
			if (s == R_NilValue)
				continue;
			sp = GET_SLOT(s, Matrix_pSym);
			psp = INTEGER(sp);
			n = (int) (XLENGTH(sp) - 1);
			if (psp[n] > INT_MAX - nnz)
				error(_("%s cannot exceed %s"), "p[length(p)]", "2^31-1");
			for (j = 0; j < n; ++j)
				*(pp++) = nnz = nnz + (psp[j + 1] - psp[j]);
		}

		SEXP i = PROTECT(allocVector(INTSXP, nnz)), si,
			iSym = (repr == 'C') ? Matrix_iSym : Matrix_jSym;
		int *pi = INTEGER(i), *psi;
		SET_SLOT(res, iSym, i);

#define BIND_C1R0(_CTYPE_, _PTR_, _MASK_) \
		do { \
			_MASK_(_CTYPE_ *px = _PTR_(x), *psx); \
			for (a = args; a != R_NilValue; a = CDR(a)) { \
				s = CAR(a); \
				if (s == R_NilValue) \
					continue; \
				PROTECT(sp = GET_SLOT(s, Matrix_pSym)); \
				PROTECT(si = GET_SLOT(s, iSym)); \
				_MASK_(PROTECT(sx = GET_SLOT(s, Matrix_xSym))); \
				psp = INTEGER(sp); \
				psi = INTEGER(si); \
				_MASK_(psx = _PTR_(sx)); \
				n = (int) (XLENGTH(sp) - 1); \
				Matrix_memcpy(pi, psi, psp[n], sizeof(int)); \
				_MASK_(Matrix_memcpy(px, psx, psp[n], sizeof(_CTYPE_))); \
				pi += psp[n]; \
				_MASK_(px += psp[n]); \
				_MASK_(UNPROTECT(1)); \
				UNPROTECT(2); \
			} \
		} while (0)

		if (kind == 'n')
			BIND_C1R0(int, LOGICAL, HIDE);
		else {
			SEXP x = PROTECT(allocVector(kindToType(kind), nnz)), sx;
			SET_SLOT(res, Matrix_xSym, x);
			BIND_CASES(BIND_C1R0);
			UNPROTECT(1);
		}
		UNPROTECT(2);

	} else if ((repr == 'C' && !margin) || (repr == 'R' && margin)) {

		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) rdim[!margin] + 1));
		int *pp = INTEGER(p);
		SET_SLOT(res, Matrix_pSym, p);
		Matrix_memset(pp, 0, (R_xlen_t) rdim[!margin] + 1, sizeof(int));

		if (rdim[0] == 0 || rdim[1] == 0) {
			UNPROTECT(1);
			return;
		}

		SEXP sp;
		int *psp, j, n = rdim[!margin];
		++pp;
		for (a = args; a != R_NilValue; a = CDR(a)) {
			s = CAR(a);
			if (s == R_NilValue)
				continue;
			sp = GET_SLOT(s, Matrix_pSym);
			psp = INTEGER(sp) + 1;
			if (n > 0 && psp[n - 1] > INT_MAX - pp[n - 1])
				error(_("%s cannot exceed %s"), "p[length(p)]", "2^31-1");
			for (j = 0; j < n; ++j)
				pp[j] += psp[j];
		}
		--pp;

		int nnz = pp[n];
		SEXP i = PROTECT(allocVector(INTSXP, nnz)), si,
			iSym = (repr == 'C') ? Matrix_iSym : Matrix_jSym;
		int *pi = INTEGER(i), *psi, *work, k, kend, pos = 0;
		SET_SLOT(res, iSym, i);
		Matrix_Calloc(work, n, int);
		Matrix_memcpy(work, pp, n, sizeof(int));

#define BIND_C0R1(_CTYPE_, _PTR_, _MASK_) \
		do { \
			_MASK_(_CTYPE_ *px = _PTR_(x), *psx); \
			for (a = args; a != R_NilValue; a = CDR(a)) { \
				s = CAR(a); \
				if (s == R_NilValue) \
					continue; \
				PROTECT(sp = GET_SLOT(s, Matrix_pSym)); \
				PROTECT(si = GET_SLOT(s, iSym)); \
				_MASK_(PROTECT(sx = GET_SLOT(s, Matrix_xSym))); \
				psp = INTEGER(sp); \
				psi = INTEGER(si); \
				_MASK_(psx = _PTR_(sx)); \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = psp[j + 1]; \
					while (k < kend) { \
						pi[work[j]] = *(psi++) + pos; \
						_MASK_(px[work[j]] = *(psx++)); \
						work[j]++; \
						++k; \
					} \
				} \
				_MASK_(UNPROTECT(1)); \
				UNPROTECT(2); \
				pos += INTEGER(GET_SLOT(s, Matrix_DimSym))[margin]; \
			} \
		} while (0)

		if (kind == 'n')
			BIND_C0R1(int, LOGICAL, HIDE);
		else {
			SEXP x = PROTECT(allocVector(kindToType(kind), nnz)), sx;
			SET_SLOT(res, Matrix_xSym, x);
			BIND_CASES(BIND_C0R1);
			UNPROTECT(1);
		}
		UNPROTECT(2);
		Matrix_Free(work, n);

	} else if (repr == 'T') {

		if (rdim[0] == 0 || rdim[1] == 0)
			return;

		R_xlen_t k, nnz = 0;
		for (a = args; a != R_NilValue; a = CDR(a)) {
			s = CAR(a);
			if (s == R_NilValue)
				continue;
			k = XLENGTH(GET_SLOT(s, Matrix_iSym));
			if (k > R_XLEN_T_MAX - nnz)
				error(_("attempt to allocate vector of length exceeding %s"),
				      "R_XLEN_T_MAX");
			nnz += k;
		}

		SEXP si, sj,
			i = PROTECT(allocVector(INTSXP, nnz)),
			j = PROTECT(allocVector(INTSXP, nnz));
		int *psi, *psj, *pi = INTEGER(i), *pj = INTEGER(j), pos = 0;
		SET_SLOT(res, Matrix_iSym, i);
		SET_SLOT(res, Matrix_jSym, j);

#define BIND_T(_CTYPE_, _PTR_, _MASK_) \
		do { \
			_MASK_(_CTYPE_ *px = _PTR_(x), *psx); \
			for (a = args; a != R_NilValue; a = CDR(a)) { \
				s = CAR(a); \
				if (s == R_NilValue) \
					continue; \
				PROTECT(si = GET_SLOT(s, Matrix_iSym)); \
				PROTECT(sj = GET_SLOT(s, Matrix_jSym)); \
				_MASK_(PROTECT(sx = GET_SLOT(s, Matrix_xSym))); \
				psi = INTEGER(si); \
				psj = INTEGER(sj); \
				_MASK_(psx = _PTR_(sx)); \
				k = XLENGTH(si); \
				if (margin) { \
					while (k--) { \
						*(pi++) = *(psi++); \
						*(pj++) = *(psj++) + pos; \
						_MASK_(*(px++) = *(psx++)); \
					} \
				} else { \
					while (k--) { \
						*(pi++) = *(psi++) + pos; \
						*(pj++) = *(psj++); \
						_MASK_(*(px++) = *(psx++)); \
					} \
				} \
				_MASK_(UNPROTECT(1)); \
				UNPROTECT(2); \
				pos += INTEGER(GET_SLOT(s, Matrix_DimSym))[margin]; \
			} \
		} while (0)

		if (kind == 'n')
			BIND_T(int, LOGICAL, HIDE);
		else {
			SEXP x = PROTECT(allocVector(kindToType(kind), nnz)), sx;
			SET_SLOT(res, Matrix_xSym, x);
			BIND_CASES(BIND_T);
			UNPROTECT(1);
		}
		UNPROTECT(2);

	} else {

		SEXP p = PROTECT(allocVector(INTSXP, rdim[margin])), sp;
		int *pp = INTEGER(p);
		for (a = args; a != R_NilValue; a = CDR(a)) {
			s = CAR(a);
			if (s == R_NilValue)
				continue;
			sp = GET_SLOT(s, Matrix_permSym);
			Matrix_memcpy(pp, INTEGER(sp), LENGTH(sp), sizeof(int));
			pp += LENGTH(sp);
		}
		SET_SLOT(res, Matrix_permSym, p);
		UNPROTECT(1);
		if (margin)
			INTEGER(GET_SLOT(res, Matrix_marginSym))[0] = 2;

	}

#undef BIND_CASES
#undef BIND_E
#undef BIND_C1R0
#undef BIND_C0R1
#undef BIND_T

	return;
}

static
SEXP bind(SEXP args, SEXP exprs, int margin, int level)
{
	if (!tagWasVector)
		tagWasVector = install(".__WAS_VECTOR__."); /* for now, a hack */

	int rdim[2], rdimnames[2];
	char kind = '\0', repr = '\0';
	scanArgs(args, exprs, margin, level,
	         rdim, rdimnames, &kind, &repr);
	if (rdim[!margin] < 0)
		/* Arguments are all NULL */
		return R_NilValue;
	if (repr == 'e' && (Matrix_int_fast64_t) rdim[0] * rdim[1] > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	char rcl[] = "...Matrix";
	if (kind == '\0' || repr == '\0') {
		if (kind != repr)
			error(_("should never happen ..."));
		rcl[0] = 'i';
		rcl[1] = 'n';
		rcl[2] = 'd';
	} else {
		rcl[0] = kind;
		rcl[1] = 'g';
		rcl[2] = repr;
		coerceArgs(args, margin, rdim, kind, repr);
	}
	SEXP res = PROTECT(newObject(rcl));
	bindArgs(args, margin, res, rdim, kind, repr);

	SEXP dim = PROTECT(GET_SLOT(res, Matrix_DimSym));
	INTEGER(dim)[0] = rdim[0];
	INTEGER(dim)[1] = rdim[1];
	UNPROTECT(1);

	if (rdimnames[0] || rdimnames[1]) {
		SEXP dimnames = PROTECT(GET_SLOT(res, Matrix_DimNamesSym)),
			marnames = R_NilValue, nms[2], nms_, a, e, s, tmp;
		int i, ivalid, r = -1, pos = 0, nprotect = 1;
		const char *scl;
		if (rdimnames[margin]) {
			PROTECT(marnames = allocVector(STRSXP, rdim[margin]));
			++nprotect;
			SET_VECTOR_ELT(dimnames, margin, marnames);
		}
		for (a = args, e = exprs; a != R_NilValue; a = CDR(a), e = CDR(e)) {
			s = CAR(a);
			if (s == R_NilValue && rdim[!margin] > 0)
				continue;
			nms[0] = nms[1] = R_NilValue;
			if (TYPEOF(s) == S4SXP) {
				ivalid = R_check_class_etc(s, valid);
				scl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];
				tmp = GET_SLOT(s, Matrix_DimSym);
				r = INTEGER(tmp)[margin];
				tmp = GET_SLOT(s, Matrix_DimNamesSym);
				if (scl[1] == 's') {
					if ((nms_ = VECTOR_ELT(tmp, 1)) != R_NilValue ||
					    (nms_ = VECTOR_ELT(tmp, 0)) != R_NilValue)
						nms[0] = nms[1] = nms_;
				} else
					for (i = 0; i < 2; ++i)
						nms[i] = VECTOR_ELT(tmp, i);
			} else {
				tmp = getAttrib(s, R_DimSymbol);
				if (TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2) {
					r = INTEGER(tmp)[margin];
					tmp = getAttrib(s, R_DimNamesSymbol);
					if (tmp != R_NilValue)
						for (i = 0; i < 2; ++i)
							nms[i] = VECTOR_ELT(tmp, i);
				} else if (rdim[!margin] == 0 || XLENGTH(s) > 0) {
					r = 1;
					if (rdim[!margin] > 0 && XLENGTH(s) == rdim[!margin])
						nms[!margin] = getAttrib(s, R_NamesSymbol);
				}
			}
			if (TAG(a) != R_NilValue) { /* only if 's' is or was a vector */
				if (TAG(a) != tagWasVector)
					nms[margin] = coerceVector(TAG(a), STRSXP);
				else if (level == 2) {
					PROTECT(nms_ = allocVector(EXPRSXP, 1));
					SET_VECTOR_ELT(nms_, 0, CAR(e));
					nms[margin] = coerceVector(nms_, STRSXP);
					UNPROTECT(1);
				} else if (level == 1 && TYPEOF(CAR(e)) == SYMSXP)
					nms[margin] = coerceVector(CAR(e), STRSXP);
			}
			if (rdimnames[!margin] && nms[!margin] != R_NilValue) {
				SET_VECTOR_ELT(dimnames, !margin, nms[!margin]);
				rdimnames[!margin] = 0;
				if (!rdimnames[margin])
					break;
			}
			if (rdimnames[ margin] && nms[ margin] != R_NilValue)
				for (i = 0; i < r; ++i)
					SET_STRING_ELT(marnames, pos + i,
					               STRING_ELT(nms[margin], i));
			pos += r;
		}
		UNPROTECT(nprotect);
	}

	UNPROTECT(1);
	return res;
}

SEXP R_bind(SEXP args)
{
	SEXP level, margin, exprs;
	args = CDR(args);  level = CAR(args);
	args = CDR(args); margin = CAR(args);
	args = CDR(args);  exprs = CAR(args);
	return bind(CDR(args), CDR(exprs), asInteger(margin), asInteger(level));
}
