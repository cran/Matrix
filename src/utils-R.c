#include "Mdefines.h"
#include "utils-R.h"

SEXP R_Matrix_version(void)
{
	SEXP ans, nms;
	PROTECT(ans = allocVector(INTSXP, 3));
	INTEGER(ans)[0] = MATRIX_PACKAGE_VERSION;
	INTEGER(ans)[1] = MATRIX_ABI_VERSION;
	INTEGER(ans)[2] = MATRIX_SUITESPARSE_VERSION;
	PROTECT(nms = allocVector(STRSXP, 3));
	SET_STRING_ELT(nms, 0, mkChar("package"));
	SET_STRING_ELT(nms, 1, mkChar("abi"));
	SET_STRING_ELT(nms, 2, mkChar("suitesparse"));
	setAttrib(ans, R_NamesSymbol, nms);
	UNPROTECT(2);
	return ans;
}

SEXP R_index_triangle(SEXP n, SEXP packed, SEXP upper, SEXP diag)
{
	SEXP r;
	int i, j, n_ = asInteger(n), packed_ = asLogical(packed),
		upper_ = asLogical(upper), diag_ = asLogical(diag);
	Matrix_int_fast64_t
		nn = (Matrix_int_fast64_t) n_ * n_,
		nx = (packed_) ? n_ + (nn - n_) / 2 : nn,
		nr = (diag_) ? n_ + (nn - n_) / 2 : (nn - n_) / 2;
	if (nx > 0x1.0p+53)
		error(_("indices would exceed %s"), "2^53");
	if (nr > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

#define DO_INDEX \
	do { \
		if (packed_) { \
			if (diag_) { \
				while (k <= nr_) \
					*(pr++) = k++; \
			} else if (upper_) { \
				for (j = 0; j < n_; ++j) { \
					for (i = 0; i < j; ++i) \
						*(pr++) = k++; \
					k++; \
				} \
			} else { \
				for (j = 0; j < n_; ++j) { \
					k++; \
					for (i = j+1; i < n_; ++i) \
						*(pr++) = k++; \
				} \
			} \
		} else if (diag_) { \
			if (upper_) { \
				for (j = 0; j < n_; ++j) { \
					for (i = 0; i <= j; ++i) \
						*(pr++) = k++; \
					k += n_-j-1; \
				} \
			} else { \
				for (j = 0; j < n_; ++j) { \
					k += j; \
					for (i = j; i < n_; ++i) \
						*(pr++) = k++; \
				} \
			} \
		} else { \
			if (upper_) { \
				for (j = 0; j < n_; ++j) { \
					for (i = 0; i < j; ++i) \
						*(pr++) = k++; \
					k += n_-j; \
				} \
			} else { \
				for (j = 0; j < n_; ++j) { \
					k += j+1; \
					for (i = j+1; i < n_; ++i) \
						*(pr++) = k++; \
				} \
			} \
		} \
	} while (0)

	if (nx > INT_MAX) {

		PROTECT(r = allocVector(REALSXP, (R_xlen_t) nr));
		double k = 1.0, nr_ = (double) nr, *pr = REAL(r);

		DO_INDEX;

	} else {

		PROTECT(r = allocVector(INTSXP, (R_xlen_t) nr));
		int k = 1, nr_ = (int) nr, *pr = INTEGER(r);

		DO_INDEX;

	}

#undef DO_INDEX

	UNPROTECT(1);
	return r;
}

SEXP R_index_diagonal(SEXP n, SEXP packed, SEXP upper)
{
	SEXP r;
	int j, n_ = asInteger(n), packed_ = asLogical(packed),
		upper_ = asLogical(upper);
	Matrix_int_fast64_t
		nn = (Matrix_int_fast64_t) n_ * n_,
		nx = (packed_) ? n_ + (nn - n_) / 2 : nn;
	if (nx > 0x1.0p+53)
		error(_("indices would exceed %s"), "2^53");

#define DO_INDEX \
	do { \
		if (!packed_) { \
			for (j = 0; j < n_; ++j) { \
				*(pr++) = k++; \
				k += n_; \
			} \
		} else if (upper_) { \
			for (j = 0; j < n_; ++j) { \
				*(pr++) = k; \
				k += j+2; \
			} \
		} else { \
			for (j = 0; j < n_; ++j) { \
				*(pr++) = k; \
				k += n_-j; \
			} \
		} \
	} while (0)

	if (nx > INT_MAX) {

		PROTECT(r = allocVector(REALSXP, n_));
		double k = 1.0, *pr = REAL(r);

		DO_INDEX;

	} else {

		PROTECT(r = allocVector(INTSXP, n_));
		int k = 1, *pr = INTEGER(r);
		DO_INDEX;

	}

#undef DO_INDEX

	UNPROTECT(1);
	return r;
}

SEXP R_nnz(SEXP x, SEXP countNA, SEXP nnzmax)
{
	int do_countNA = asLogical(countNA);
	R_xlen_t n = XLENGTH(x), nnz = 0;
	double n_ = asReal(nnzmax);
	if (!ISNAN(n_) && n_ >= 0.0 && n_ < (double) n)
		n = (R_xlen_t) n_;

#define DO_NNZ(_CTYPE_, _PTR_, _ISNA_, _ISNZ_, _STRICTLY_ISNZ_) \
	do { \
		_CTYPE_ *px = _PTR_(x); \
		if (do_countNA == NA_LOGICAL) { \
			while (n-- > 0) { \
				if (_ISNA_(*px)) \
					return ScalarInteger(NA_INTEGER); \
				if (_ISNZ_(*px)) \
					++nnz; \
				++px; \
			} \
		} else if (do_countNA != 0) { \
			while (n-- > 0) { \
				if (_ISNZ_(*px)) \
					++nnz; \
				++px; \
			} \
		} else { \
			while (n-- > 0) { \
				if (_STRICTLY_ISNZ_(*px)) \
					++nnz; \
				++px; \
			} \
		} \
	} while (0)

	switch (TYPEOF(x)) {
	case LGLSXP:
		DO_NNZ(int, LOGICAL,
		       ISNA_LOGICAL, ISNZ_LOGICAL, STRICTLY_ISNZ_LOGICAL);
		break;
	case INTSXP:
		DO_NNZ(int, INTEGER,
		       ISNA_INTEGER, ISNZ_INTEGER, STRICTLY_ISNZ_INTEGER);
	break;
	case REALSXP:
		DO_NNZ(double, REAL,
		       ISNA_REAL, ISNZ_REAL, STRICTLY_ISNZ_REAL);
	break;
	case CPLXSXP:
		DO_NNZ(Rcomplex, COMPLEX,
		       ISNA_COMPLEX, ISNZ_COMPLEX, STRICTLY_ISNZ_COMPLEX);
	break;
	default:
		ERROR_INVALID_TYPE(x, __func__);
	}

#undef DO_NNZ

	return (nnz <= INT_MAX)
		? ScalarInteger((int) nnz) : ScalarReal((double) nnz);
}


/* ================================================================== */
/* ================================================================== */


#define TRUE_  ScalarLogical(1)
#define FALSE_ ScalarLogical(0)

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// all0	 <- function(x) !any(is.na(x)) && all(!x) ## ~= allFalse
// allFalse <- function(x) !any(x) && !any(is.na(x)) ## ~= all0
SEXP R_all0(SEXP x) {
	if (!isVectorAtomic(x)) {
		if (length(x) == 0) return TRUE_;
		// Typically S4.  TODO: Call the R code above, instead!
		error(_("Argument must be numeric-like atomic vector"));
	}
	R_xlen_t i, n = XLENGTH(x);
	if (n == 0) return TRUE_;

	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *xx = LOGICAL(x);
		for (i = 0; i < n; i++)
			if (xx[i] == NA_LOGICAL || xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	case INTSXP:
	{
		int *xx = INTEGER(x);
		for (i = 0; i < n; i++)
			if (xx[i] == NA_INTEGER || xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	case REALSXP:
	{
		double *xx = REAL(x);
		for (i = 0; i < n; i++)
			if (ISNAN(xx[i]) || xx[i] != 0.) return FALSE_;
		return TRUE_;
	}
	case RAWSXP:
	{
		unsigned char *xx = RAW(x);
		for (i = 0; i < n; i++)
			if (xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	}
	error(_("Argument must be numeric-like atomic vector"));
	return R_NilValue; // -Wall
}

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// any0 <- function(x) isTRUE(any(x == 0)) ## ~= anyFalse
// anyFalse <- function(x) isTRUE(any(!x)) ## ~= any0
SEXP R_any0(SEXP x) {
	if (!isVectorAtomic(x)) {
		if (length(x) == 0) return FALSE_;
		// Typically S4.  TODO: Call the R code above, instead!
		error(_("Argument must be numeric-like atomic vector"));
	}
	R_xlen_t i, n = XLENGTH(x);
	if (n == 0) return FALSE_;

	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *xx = LOGICAL(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	case INTSXP:
	{
		int *xx = INTEGER(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	case REALSXP:
	{
		double *xx = REAL(x);
		for (i = 0; i < n; i++) if (xx[i] == 0.) return TRUE_;
		return FALSE_;
	}
	case RAWSXP:
	{
		unsigned char *xx = RAW(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	}
	error(_("Argument must be numeric-like atomic vector"));
	return R_NilValue; // -Wall
}

#undef TRUE_
#undef FALSE_

// Almost "Cut n Paste" from ...R../src/main/array.c  do_matrix() :
// used in ../R/Matrix.R as
//
// .External(Mmatrix,
//		 data, nrow, ncol, byrow, dimnames,
//		 missing(nrow), missing(ncol))
SEXP Mmatrix(SEXP args)
{
	SEXP vals, ans, snr, snc, dimnames;
	int nr = 1, nc = 1, byrow, miss_nr, miss_nc;
	R_xlen_t lendat;

	args = CDR(args); /* skip 'name' */
	vals = CAR(args); args = CDR(args);
	/* Supposedly as.vector() gave a vector type, but we check */
	switch (TYPEOF(vals)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
	case STRSXP:
	case RAWSXP:
	case EXPRSXP:
	case VECSXP:
		break;
	default:
		error(_("'data' must be of a vector type"));
	}
	lendat = XLENGTH(vals);
	snr = CAR(args); args = CDR(args);
	snc = CAR(args); args = CDR(args);
	byrow = asLogical(CAR(args)); args = CDR(args);
	if (byrow == NA_INTEGER)
		error(_("invalid '%s' argument"), "byrow");
	dimnames = CAR(args);
	args = CDR(args);
	miss_nr = asLogical(CAR(args)); args = CDR(args);
	miss_nc = asLogical(CAR(args));

	if (!miss_nr) {
		if (!isNumeric(snr)) error(_("non-numeric matrix extent"));
		nr = asInteger(snr);
		if (nr == NA_INTEGER)
			error(_("invalid 'nrow' value (too large or NA)"));
		if (nr < 0)
			error(_("invalid 'nrow' value (< 0)"));
	}
	if (!miss_nc) {
		if (!isNumeric(snc)) error(_("non-numeric matrix extent"));
		nc = asInteger(snc);
		if (nc == NA_INTEGER)
			error(_("invalid 'ncol' value (too large or NA)"));
		if (nc < 0)
			error(_("invalid 'ncol' value (< 0)"));
	}
	if (miss_nr && miss_nc) {
		if (lendat > INT_MAX) error("data is too long");
		nr = (int) lendat;
	} else if (miss_nr) {
		if (lendat > (double) nc * INT_MAX) error("data is too long");
		nr = (int) ceil((double) lendat / (double) nc);
	} else if (miss_nc) {
		if (lendat > (double) nr * INT_MAX) error("data is too long");
		nc = (int) ceil((double) lendat / (double) nr);
	}

	if (lendat > 0) {
		R_xlen_t nrc = (R_xlen_t) nr * nc;
		if (lendat > 1 && nrc % lendat != 0) {
			if ((lendat > nr && (lendat / nr) * nr != lendat) ||
			    (lendat < nr && (nr / lendat) * lendat != nr))
				warning(_("data length [%lld] is not a sub-multiple "
				          "or multiple of the number of rows [%d]"),
				        (long long)lendat, nr);
			else if ((lendat > nc && (lendat / nc) * nc != lendat) ||
				 (lendat < nc && (nc / lendat) * lendat != nc))
				warning(_("data length [%lld] is not a sub-multiple "
				          "or multiple of the number of columns [%d]"),
					(long long)lendat, nc);
		} else if (lendat > 1 && nrc == 0)
			warning(_("data length exceeds size of matrix"));
	}

#ifndef LONG_VECTOR_SUPPORT
	if ((double) nr * (double) nc > INT_MAX)
		error(_("too many elements specified"));
#endif

	PROTECT(ans = allocMatrix(TYPEOF(vals), nr, nc));
	if (isVector(vals)) {
	    if(lendat)
		copyMatrix(ans, vals, (Rboolean) byrow);
	    else { /* fill with NAs */
		R_xlen_t N = (R_xlen_t) nr * nc, i;
		switch (TYPEOF(vals)) {
		case STRSXP:
			for (i = 0; i < N; i++)
				SET_STRING_ELT(ans, i, NA_STRING);
			break;
		case LGLSXP:
			for (i = 0; i < N; i++)
				LOGICAL(ans)[i] = NA_LOGICAL;
			break;
		case INTSXP:
			for (i = 0; i < N; i++)
				INTEGER(ans)[i] = NA_INTEGER;
			break;
		case REALSXP:
			for (i = 0; i < N; i++)
				REAL(ans)[i] = NA_REAL;
			break;
		case CPLXSXP:
		{
			/* Initialization must work whether Rcomplex is typedef-ed
			   to a struct { R < 4.3.0 } or to a union { R >= 4.3.0 }
			*/
			Rcomplex zna = { .r = NA_REAL, .i = 0.0 };
			for (i = 0; i < N; i++)
				COMPLEX(ans)[i] = zna;
			break;
		}
		case RAWSXP:
			// FIXME:  N may overflow size_t !!
			memset(RAW(ans), 0, N);
			break;
		default:
			/* don't fill with anything */
			;
		}
	    }
	}
	if (!isNull(dimnames)&& length(dimnames) > 0)
		ans = dimnamesgets(ans, dimnames);
	UNPROTECT(1);
	return ans;
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
static
int *expand_cmprPt(int ncol, const int mp[], int mj[])
{
	int j;
	for (j = 0; j < ncol; j++) {
		int j2 = mp[j+1], jj;
		for (jj = mp[j]; jj < j2; jj++)
			mj[jj] = j;
	}
	return mj;
}

/** Return a 2 column matrix  '' cbind(i, j) ''  of 0-origin index vectors (i,j)
 *  which entirely correspond to the (i,j) slots of
 *  as(x, "TsparseMatrix") :
 */
SEXP compressed_non_0_ij(SEXP x, SEXP colP)
{
    int col = asLogical(colP); /* 1 if "C"olumn compressed;  0 if "R"ow */
    SEXP ans, indSym = col ? Matrix_iSym : Matrix_jSym;
    SEXP indP = PROTECT(GET_SLOT(x, indSym)),
	 pP   = PROTECT(GET_SLOT(x, Matrix_pSym));
    int i, *ij;
    int nouter = INTEGER(GET_SLOT(x, Matrix_DimSym))[col ? 1 : 0],
	n_el   = INTEGER(pP)[nouter]; /* is only == length(indP), if the
				     inner slot is not over-allocated */

    ij = INTEGER(ans = PROTECT(allocMatrix(INTSXP, n_el, 2)));
    /* expand the compressed margin to 'i' or 'j' : */
    expand_cmprPt(nouter, INTEGER(pP), &ij[col ? n_el : 0]);
    /* and copy the other one: */
    if (col)
	for(i = 0; i < n_el; i++)
	    ij[i] = INTEGER(indP)[i];
    else /* row compressed */
	for(i = 0; i < n_el; i++)
	    ij[i + n_el] = INTEGER(indP)[i];

    UNPROTECT(3);
    return ans;
}

SEXP Matrix_expand_pointers(SEXP pP)
{
	int n = length(pP) - 1;
	int *p = INTEGER(pP);
	SEXP ans = PROTECT(allocVector(INTSXP, p[n]));

	expand_cmprPt(n, p, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param ij: 2-column integer matrix
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd(SEXP ij, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
	SEXP ans;
	int *ij_di = NULL, n, nprot=1;
	int check_bounds = Rf_asLogical(chk_bnds), one_ind = Rf_asLogical(orig_1);

	if (TYPEOF(di) != INTSXP) {
		di = PROTECT(coerceVector(di, INTSXP));
		nprot++;
	}
	if (TYPEOF(ij) != INTSXP) {
		ij = PROTECT(coerceVector(ij, INTSXP));
		nprot++;
	}
	if (!isMatrix(ij) ||
	    (ij_di = INTEGER(getAttrib(ij, R_DimSymbol)))[1] != 2)
		error(_("Argument ij must be 2-column integer matrix"));
	n = ij_di[0];
	int *Di = INTEGER(di), *IJ = INTEGER(ij),
		*j_ = IJ+n;/* pointer offset! */

	if ((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
		ans = PROTECT(allocVector(REALSXP, n));
		double *ii = REAL(ans), nr = (double) Di[0];

#define do_ii_FILL(_i_, _j_) \
		int i; \
		if (check_bounds) { \
			for (i = 0; i < n; i++) { \
				if (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER) \
					ii[i] = NA_INTEGER; \
				else { \
					register int i_i, j_i; \
					if (one_ind) { \
						i_i = _i_[i]-1; \
						j_i = _j_[i]-1; \
					} else { \
						i_i = _i_[i]; \
						j_i = _j_[i]; \
					} \
					if (i_i < 0 || i_i >= Di[0]) \
						error(_("subscript 'i' out of bounds in M[ij]")); \
					if (j_i < 0 || j_i >= Di[1]) \
						error(_("subscript 'j' out of bounds in M[ij]")); \
					ii[i] = i_i + j_i * nr; \
				} \
			} \
		} else { \
			for (i = 0; i < n; i++) \
				ii[i] = (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER) \
					? NA_INTEGER \
					: ((one_ind) \
					   ? ((_i_[i]-1) + (_j_[i]-1) * nr) \
					   :   _i_[i] + _j_[i] * nr); \
		}

		do_ii_FILL(IJ, j_);
	} else {
	ans = PROTECT(allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];

	do_ii_FILL(IJ, j_);
	}
	UNPROTECT(nprot);
	return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param i: integer vector
 * @param j: integer vector of same length as 'i'
 * @param orig_1: logical: if TRUE, "1-origin" otherwise "0-origin"
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
	SEXP ans;
	int n = LENGTH(i), nprot = 1;
	int check_bounds = Rf_asLogical(chk_bnds), one_ind = Rf_asLogical(orig_1);

	if (TYPEOF(di)!= INTSXP) {
		di = PROTECT(coerceVector(di,INTSXP));
		nprot++;
	}
	if (TYPEOF(i) != INTSXP) {
		i = PROTECT(coerceVector(i, INTSXP));
		nprot++;
	}
	if (TYPEOF(j) != INTSXP) {
		j = PROTECT(coerceVector(j, INTSXP));
		nprot++;
	}
	if (LENGTH(j) != n)
		error(_("i and j must be integer vectors of the same length"));

	int *Di = INTEGER(di), *i_ = INTEGER(i), *j_ = INTEGER(j);

	if ((Di[0] * (double) Di[1]) >= 1 + (double) INT_MAX) { /* need double */
		ans = PROTECT(allocVector(REALSXP, n));
		double *ii = REAL(ans), nr = (double) Di[0];

		do_ii_FILL(i_, j_);
	} else {
		ans = PROTECT(allocVector(INTSXP, n));
		int *ii = INTEGER(ans), nr = Di[0];

		do_ii_FILL(i_, j_);
	}
	UNPROTECT(nprot);
	return ans;
}
#undef do_ii_FILL

#define _rle_d_
#include "t_rle.c"
#undef _rle_d_

#define _rle_i_
#include "t_rle.c"
#undef _rle_i_
