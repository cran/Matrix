#include "Lapack-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "idz.h"
#include "coerce.h"
#include "dense.h"
#include "sparse.h"
#include "matmult.h"

static
void matmultDim(SEXP x, SEXP y, int *xtrans, int *ytrans, int *ztrans,
                int *m, int *n, int *v)
{
	*xtrans = (*xtrans) ? 1 : 0;
	*ytrans = (*ytrans) ? 1 : 0;
	*ztrans = (*ztrans) ? 1 : 0;
	if (y == R_NilValue) {
		SEXP
			xdim = (TYPEOF(x) == S4SXP)
			? GET_SLOT(x, Matrix_DimSym) : getAttrib(x, R_DimSymbol);
		if (TYPEOF(xdim) == INTSXP && LENGTH(xdim) == 2) {
			*v = 0;
			*m = *n = INTEGER(xdim)[(*xtrans) ? 1 : 0];
		} else if (XLENGTH(x) <= INT_MAX) {
			*v = 1;
			*m = *n = (*xtrans) ? 1 : LENGTH(x);
		} else
			error(_("dimensions cannot exceed %s"), "2^31-1");
		*ytrans = (*xtrans) ? 0 : 1;
	} else {
		/* MJ: So that I don't lose my mind ... : */
		if (*ztrans) {
			int tmp = !(*xtrans); *xtrans = !(*ytrans); *ytrans = tmp;
			SEXP s = x; x = y; y = s;
		}
		SEXP
			xdim = (TYPEOF(x) == S4SXP)
			? GET_SLOT(x, Matrix_DimSym) : getAttrib(x, R_DimSymbol),
			ydim = (TYPEOF(y) == S4SXP)
			? GET_SLOT(y, Matrix_DimSym) : getAttrib(y, R_DimSymbol);
		int xm, xn, ym, yn, x2, y2;
		xm = xn = ym = yn = -1;
		x2 = TYPEOF(xdim) == INTSXP && LENGTH(xdim) == 2;
		y2 = TYPEOF(ydim) == INTSXP && LENGTH(ydim) == 2;
		if (x2) {
			int *pxdim = INTEGER(xdim);
			xm = pxdim[0];
			xn = pxdim[1];
		} else if (XLENGTH(x) > INT_MAX)
			error(_("dimensions cannot exceed %s"), "2^31-1");
		if (y2) {
			int *pydim = INTEGER(ydim);
			ym = pydim[0];
			yn = pydim[1];
		} else if (XLENGTH(y) > INT_MAX)
			error(_("dimensions cannot exceed %s"), "2^31-1");
		/* MJ: R's do_matprod behaves quite asymmetrically ... what a pain */
		if (x2 && y2)
			*v = 0;
		else if (y2) {
			*v = (*ztrans) ? 2 : 1;
			int k = (*ytrans) ? yn : ym, xl = LENGTH(x);
			if (k == xl || (k == 1 && !(*xtrans))) {
				xm = (int) xl;
				xn = 1;
				*xtrans = (k == xl) ? 1 : 0;
			}
		} else if (x2) {
			*v = (*ztrans) ? 1 : 2;
			int k = (*xtrans) ? xm : xn, yl = LENGTH(y);
			if (*ytrans) {
				if (xm == 1 || xn == 1) {
					ym = (int) yl;
					yn = 1;
					*ytrans = (((*xtrans) ? xn : xm) == 1) ? 0 : 1;
				}
			} else {
				if (k == yl || k == 1) {
					ym = (int) yl;
					yn = 1;
					*ytrans = (k == yl) ? 0 : 1;
				}
			}
		} else {
			*v = 3;
			int xl = LENGTH(x), yl = LENGTH(y);
			if (*xtrans) {
				xm = xl;
				xn = 1;
				ym = yl;
				yn = 1;
				*ytrans = xl == 1;
			} else if (*ytrans) {
				xm = xl;
				xn = 1;
				ym = yl;
				yn = 1;
				/* *xtrans = 0; */
			} else {
				xm = 1;
				xn = xl;
				ym = (xl == 1) ? 1 : yl;
				yn = (xl == 1) ? yl : 1;
			}
		}
		if (((*xtrans) ? xm : xn) != ((*ytrans) ? yn : ym))
			error(_("non-conformable arguments"));
		*m = (*xtrans) ? xn : xm;
		*n = (*ytrans) ? ym : yn;
		if (*ztrans) {
			int tmp = !(*xtrans); *xtrans = !(*ytrans); *ytrans = tmp;
			tmp = *m; *m = *n; *n = tmp;
		}
	}
	return;
}

static
void matmultDN(SEXP dest, SEXP asrc, int ai, SEXP bsrc, int bi) {
	SEXP s;
	if (!isNull(s = VECTOR_ELT(asrc, ai)))
		SET_VECTOR_ELT(dest, 0, s);
	if (!isNull(s = VECTOR_ELT(bsrc, bi)))
		SET_VECTOR_ELT(dest, 1, s);
	PROTECT(asrc = getAttrib(asrc, R_NamesSymbol));
	PROTECT(bsrc = getAttrib(bsrc, R_NamesSymbol));
	if (!isNull(asrc) || !isNull(bsrc)) {
		SEXP destnms = PROTECT(allocVector(STRSXP, 2));
		if (!isNull(asrc) && *CHAR(s = STRING_ELT(asrc, ai)) != '\0')
			SET_STRING_ELT(destnms, 0, s);
		if (!isNull(bsrc) && *CHAR(s = STRING_ELT(bsrc, bi)) != '\0')
			SET_STRING_ELT(destnms, 1, s);
		setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return;
}

/* op(<dge>) * op(<dge>) */
static
SEXP dgeMatrix_matmult(SEXP a, SEXP b, int atrans, int btrans)
{
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int *padim = INTEGER(adim), am = padim[0], an = padim[1],
		rm = (atrans) ? an : am, rk = (atrans) ? am : an;

	if (b == R_NilValue) {

		if ((Matrix_int_fast64_t) rm * rm > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");

		SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

		char rcl[] = ".poMatrix";
		rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
		SEXP r = PROTECT(newObject(rcl));

		SEXP rdim = GET_SLOT(r, Matrix_DimSym);
		int *prdim = INTEGER(rdim);
		prdim[0] = prdim[1] = rm;

		SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
			rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
		symDN(rdimnames, adimnames, (atrans) ? 1 : 0);
		UNPROTECT(2); /* rdimnames, adimnames */

		if (rm > 0) {
		SEXP rx = PROTECT(allocVector(TYPEOF(ax), (R_xlen_t) rm * rm));
#ifdef MATRIX_ENABLE_ZMATRIX
		if (TYPEOF(ax) == CPLXSXP) {
		Rcomplex *prx = COMPLEX(rx);
		Matrix_memset(prx, 0, (R_xlen_t) rm * rm, sizeof(Rcomplex));
		if (rk > 0) {
			Rcomplex *pax = COMPLEX(ax),
				zero = Matrix_zzero, one = Matrix_zone;
			F77_CALL(zsyrk)(
				"U", (atrans) ? "T" : "N", &rm, &rk,
				&one, pax, &am, &zero, prx, &rm FCONE FCONE);
		}
		} else {
#endif
		double *prx = REAL(rx);
		Matrix_memset(prx, 0, (R_xlen_t) rm * rm, sizeof(double));
		if (rk > 0) {
			double *pax = REAL(ax),
				zero = 0.0, one = 1.0;
			F77_CALL(dsyrk)(
				"U", (atrans) ? "T" : "N", &rm, &rk,
				&one, pax, &am, &zero, prx, &rm FCONE FCONE);
		}
#ifdef MATRIX_ENABLE_ZMATRIX
		}
#endif
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
		}

		UNPROTECT(2); /* r, ax */
		return r;

	} else {

		SEXP bdim = GET_SLOT(b, Matrix_DimSym);
		int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
			rn = (btrans) ? bm : bn;

		if (rk != ((btrans) ? bn : bm))
			error(_("non-conformable arguments"));
		if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");

		SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

		char rcl[] = ".geMatrix";
		rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
		SEXP r = PROTECT(newObject(rcl));

		SEXP rdim = GET_SLOT(r, Matrix_DimSym);
		int *prdim = INTEGER(rdim);
		prdim[0] = rm;
		prdim[1] = rn;

		SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
			bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
			rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
		matmultDN(rdimnames,
		          adimnames, (atrans) ? 1 : 0,
		          bdimnames, (btrans) ? 0 : 1);
		UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

		if (rm > 0 && rn > 0) {
		SEXP rx = PROTECT(allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
#ifdef MATRIX_ENABLE_ZMATRIX
		if (TYPEOF(ax) == CPLXSXP) {
		Rcomplex *prx = COMPLEX(rx);
		if (rk == 0)
			Matrix_memset(prx, 0, (R_xlen_t) rm * rn, sizeof(Rcomplex));
		else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx),
				zero = Matrix_zzero, one = Matrix_zone;
			F77_CALL(zgemm)(
				(atrans) ? "T" : "N", (btrans) ? "T" : "N", &rm, &rn, &rk,
				&one, pax, &am, pbx, &bm, &zero, prx, &rm FCONE FCONE);
			UNPROTECT(1); /* bx */
		}
		} else {
#endif
		double *prx = REAL(rx);
		if (rk == 0)
			Matrix_memset(prx, 0, (R_xlen_t) rm * rn, sizeof(double));
		else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			double *pax = REAL(ax), *pbx = REAL(bx),
				zero = 0.0, one = 1.0;
			F77_CALL(dgemm)(
				(atrans) ? "T" : "N", (btrans) ? "T" : "N", &rm, &rn, &rk,
				&one, pax, &am, pbx, &bm, &zero, prx, &rm FCONE FCONE);
			UNPROTECT(1); /* bx */
		}
#ifdef MATRIX_ENABLE_ZMATRIX
		}
#endif
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
		}

		UNPROTECT(2); /* r, ax */
		return r;

	}
}

/* <dsy> * op(<dge>)  or  op(<dge>) * <dsy> */
static
SEXP dsyMatrix_matmult(SEXP a, SEXP b, int aleft, int btrans)
{
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int rk = INTEGER(adim)[0];

	SEXP bdim = GET_SLOT(b, Matrix_DimSym);
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans) ? bn : bm, rn = (btrans) ? bm : bn;

	if (rk != ((aleft == btrans) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = ".geMatrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;

	SEXP adimnames = PROTECT(get_symmetrized_DimNames(a, -1)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aleft)
	matmultDN(rdimnames, adimnames,      0, bdimnames, !btrans);
	else
	matmultDN(rdimnames, bdimnames, btrans, adimnames,       1);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	if (rm > 0 && rn > 0) {
	SEXP auplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
		bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
	char aul = *CHAR(STRING_ELT(auplo, 0));
	int i, d = (aleft) ? rn : rm,
		binc = (aleft) ? bm : 1, bincp = (aleft) ? 1 : bm,
		rinc = (aleft) ? 1 : rm, rincp = (aleft) ? rm : 1;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		zero = Matrix_zero, one = Matrix_zone;
	if (!btrans)
		F77_CALL(zsymm)(
			(aleft) ? "L" : "R", &aul, &rm, &rn,
			&one, pax, &rk, pbx, &bm, &zero, prx, &rm FCONE FCONE);
	else {
		for (i = 0; i < d; ++i) {
			F77_CALL(zsymv)(
				&aul, &rk, &one, pax, &rk, pbx, &binc, &zero, prx, &rinc FCONE);
			pbx += bincp;
			prx += rincp;
		}
	}
	} else {
#endif
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		zero = 0.0, one = 1.0;
	if (!btrans)
		F77_CALL(dsymm)(
			(aleft) ? "L" : "R", &aul, &rm, &rn,
			&one, pax, &rk, pbx, &bm, &zero, prx, &rm FCONE FCONE);
	else {
		for (i = 0; i < d; ++i) {
			F77_CALL(dsymv)(
				&aul, &rk, &one, pax, &rk, pbx, &binc, &zero, prx, &rinc FCONE);
			pbx += bincp;
			prx += rincp;
		}
	}
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(3); /* rx, bx, auplo */
	}

	UNPROTECT(2); /* r, ax */
	return r;
}

/* <dsp> * op(<dge>)  or  op(<dge>) * <dsp> */
static
SEXP dspMatrix_matmult(SEXP a, SEXP b, int aleft, int btrans)
{
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int rk = INTEGER(adim)[0];

	SEXP bdim = GET_SLOT(b, Matrix_DimSym);
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans) ? bn : bm, rn = (btrans) ? bm : bn;

	if (rk != ((aleft == btrans) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = ".geMatrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;

	SEXP adimnames = PROTECT(get_symmetrized_DimNames(a, -1)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aleft)
	matmultDN(rdimnames, adimnames,      0, bdimnames, !btrans);
	else
	matmultDN(rdimnames, bdimnames, btrans, adimnames,       1);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	if (rm > 0 && rn > 0) {
	SEXP auplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
		bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
	char aul = *CHAR(STRING_ELT(auplo, 0));
	int i, d = (aleft) ? rn : rm,
		binc = (aleft == btrans) ? bm : 1, bincp = (aleft == btrans) ? 1 : bm,
		rinc = (aleft          ) ? 1 : rm, rincp = (aleft          ) ? rm : 1;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		zero = Matrix_zzero, one = Matrix_zone;
	for (i = 0; i < d; ++i) {
		F77_CALL(zspmv)(
			&aul, &rk, &one, pax, pbx, &binc, &zero, prx, &rinc FCONE);
		pbx += bincp;
		prx += rincp;
	}
	} else {
#endif
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		zero = 0.0, one = 1.0;
	for (i = 0; i < d; ++i) {
		F77_CALL(dspmv)(
			&aul, &rk, &one, pax, pbx, &binc, &zero, prx, &rinc FCONE);
		pbx += bincp;
		prx += rincp;
	}
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(3); /* rx, bx, auplo */
	}

	UNPROTECT(2); /* r, ax */
	return r;
}

/* op(<dtr>) * op(<dge>)  or  op(<dge>) * op(<dtr>) */
static
SEXP dtrMatrix_matmult(SEXP a, SEXP b, int aleft, int atrans, int btrans,
                       int triangular)
{
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int rk = INTEGER(adim)[0];

	SEXP bdim = GET_SLOT(b, Matrix_DimSym);
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans) ? bn : bm, rn = (btrans) ? bm : bn;

	if (rk != ((aleft == btrans) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	rcl[1] = (triangular > 0) ? 't' : 'g';
	rcl[2] = (triangular > 0) ? 'r' : 'e';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;

	SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aleft)
	matmultDN(rdimnames, adimnames, atrans, bdimnames, !btrans);
	else
	matmultDN(rdimnames, bdimnames, btrans, adimnames, !atrans);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	SEXP auplo = GET_SLOT(a, Matrix_uploSym);
	char aul = *CHAR(STRING_ELT(auplo, 0));
	if (triangular > 0 && ((atrans) ? aul == 'U' : aul != 'U')) {
		if (atrans)
			auplo = mkString("L");
		PROTECT(auplo);
		SET_SLOT(r, Matrix_uploSym, auplo);
		UNPROTECT(1); /* auplo */
	}

	SEXP adiag = GET_SLOT(a, Matrix_diagSym);
	char adi = *CHAR(STRING_ELT(adiag, 0));
	if (triangular > 1 && adi != 'N') {
		PROTECT(adiag);
		SET_SLOT(r, Matrix_diagSym, adiag);
		UNPROTECT(1); /* adiag */
	}

	if (rm > 0 && rn > 0) {
	SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		one = Matrix_zone;
	if (!btrans)
		Matrix_memcpy(prx, pbx, (R_xlen_t) bm * bn, sizeof(Rcomplex));
	else
		ztranspose2(prx, pbx, bm, bn);
	F77_CALL(ztrmm)(
		(aleft) ? "L" : "R", &aul, (atrans) ? "T" : "N", &adi, &rm, &rn,
		&one, pax, &rk, prx, &rm FCONE FCONE FCONE FCONE);
	} else {
#endif
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		one = 1.0;
	if (!btrans)
		Matrix_memcpy(prx, pbx, (R_xlen_t) bm * bn, sizeof(double));
	else
		dtranspose2(prx, pbx, bm, bn);
	F77_CALL(dtrmm)(
		(aleft) ? "L" : "R", &aul, (atrans) ? "T" : "N", &adi, &rm, &rn,
		&one, pax, &rk, prx, &rm FCONE FCONE FCONE FCONE);
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(2); /* rx, bx */
	}

	UNPROTECT(2); /* r, ax */
	return r;
}

/* op(<dtp>) * op(<dge>)  or  op(<dge>) * op(<dtp>) */
static
SEXP dtpMatrix_matmult(SEXP a, SEXP b, int aleft, int atrans, int btrans,
                       int triangular)
{
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int rk = INTEGER(adim)[0];

	SEXP bdim = GET_SLOT(b, Matrix_DimSym);
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans) ? bn : bm, rn = (btrans) ? bm : bn;

	if (rk != ((aleft == btrans) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	rcl[1] = (triangular > 0) ? 't' : 'g';
	rcl[2] = (triangular > 0) ? 'r' : 'e';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;

	SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aleft)
	matmultDN(rdimnames, adimnames, atrans, bdimnames, !btrans);
	else
	matmultDN(rdimnames, bdimnames, btrans, adimnames, !atrans);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	SEXP auplo = GET_SLOT(a, Matrix_uploSym);
	char aul = *CHAR(STRING_ELT(auplo, 0));
	if (triangular > 0 && ((atrans) ? aul == 'U' : aul != 'U')) {
		if (atrans)
			auplo = mkString("L");
		PROTECT(auplo);
		SET_SLOT(r, Matrix_uploSym, auplo);
		UNPROTECT(1); /* auplo */
	}

	SEXP adiag = GET_SLOT(a, Matrix_diagSym);
	char adi = *CHAR(STRING_ELT(adiag, 0));
	if (triangular > 1 && adi != 'N') {
		PROTECT(adiag);
		SET_SLOT(r, Matrix_diagSym, adiag);
		UNPROTECT(1); /* adiag */
	}

	if (rm > 0 && rn > 0) {
	SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
	int i, rinc = (aleft) ? 1 : rm, rincp = (aleft) ? rm : 1;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx);
	if (!btrans)
		Matrix_memcpy(prx, pbx, (R_xlen_t) bm * bn, sizeof(Rcomplex));
	else
		ztranspose2(prx, pbx, bm, bn);
	for (i = 0; i < rn; ++i) {
		F77_CALL(ztpmv)(
			&aul, (aleft == atrans) ? "T" : "N", &adi, &rk,
			pax, prx, &rinc FCONE FCONE FCONE);
		prx += rincp;
	}
	} else {
#endif
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx);
	if (!btrans)
		Matrix_memcpy(prx, pbx, (R_xlen_t) bm * bn, sizeof(double));
	else
		dtranspose2(prx, pbx, bm, bn);
	for (i = 0; i < rn; ++i) {
		F77_CALL(dtpmv)(
			&aul, (aleft == atrans) ? "T" : "N", &adi, &rk,
			pax, prx, &rinc FCONE FCONE FCONE);
		prx += rincp;
	}
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(2); /* rx, bx */
	}

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP R_dense_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans)
{
	int xtrans_ = LOGICAL(xtrans)[0], ytrans_ = LOGICAL(ytrans)[0],
		ztrans_ = 0, m, n, v;
	matmultDim(x, y, &xtrans_, &ytrans_, &ztrans_, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(x, &xpid);
	PROTECT_WITH_INDEX(y, &ypid);

	if (TYPEOF(x) != S4SXP) {
		REPROTECT(x = matrix_as_dense(x, ",ge", '\0', '\0', xtrans_, 0), xpid);
		if (v == 1) {
			/* Vector: discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
			               (xtrans_) ? 1 : 0, R_NilValue);
			xtrans_ = 0;
		}
	}
	if (TYPEOF(y) != S4SXP && y != R_NilValue) {
		REPROTECT(y = matrix_as_dense(y, ",ge", '\0', '\0', ytrans_, 0), ypid);
		if (v == 2) {
			/* Vector: discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym),
			               (ytrans_) ? 1 : 0, R_NilValue);
			ytrans_ = 0;
		}
	}

	static const char *valid[] = { VALID_DENSE, "" };
	const char *xcl = NULL, *ycl = NULL;
	int ivalid;
	ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	xcl = valid[ivalid];
	if (y != R_NilValue) {
	ivalid = R_check_class_etc(y, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(y, __func__);
	ycl = valid[ivalid];
	}

	char kind = (xcl[0] == 'z' || (y != R_NilValue && ycl[0] == 'z'))
		? 'z' : 'd';
	if (xcl[0] != kind) {
		REPROTECT(x = dense_as_kind(x, xcl, kind, 0), xpid);
		xcl = valid[R_check_class_etc(x, valid)];
	}
	if (y != R_NilValue) {
	if (ycl[0] != kind) {
		REPROTECT(y = dense_as_kind(y, ycl, kind, 0), ypid);
		ycl = valid[R_check_class_etc(y, valid)];
	}
	}

	if (y == R_NilValue) {
		REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
		x = dgeMatrix_matmult(x, y, xtrans_, !xtrans_);
	} else if (xcl[1] == 'g' && ycl[1] == 'g') {
		x = dgeMatrix_matmult(x, y, xtrans_, ytrans_);
	} else if (xcl[1] == 'g' || ycl[1] == 'g') {
		x = (xcl[1] == 'g')
			? ((ycl[1] == 's')
			   ? ((ycl[2] != 'p')
			      ? dsyMatrix_matmult(y, x, 0, xtrans_)
			      : dspMatrix_matmult(y, x, 0, xtrans_))
			   : ((ycl[2] != 'p')
			      ? dtrMatrix_matmult(y, x, 0, ytrans_, xtrans_, 0)
			      : dtpMatrix_matmult(y, x, 0, ytrans_, xtrans_, 0)))
			: ((xcl[1] == 's')
			   ? ((xcl[2] != 'p')
			      ? dsyMatrix_matmult(x, y, 1, ytrans_)
			      : dspMatrix_matmult(x, y, 1, ytrans_))
			   : ((xcl[2] != 'p')
			      ? dtrMatrix_matmult(x, y, 1, xtrans_, ytrans_, 0)
			      : dtpMatrix_matmult(x, y, 1, xtrans_, ytrans_, 0)));
	} else if (xcl[1] == 's' && ycl[1] == 's') {
		if (xcl[2] == 'p' && ycl[2] == 'p') {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = dspMatrix_matmult(x, y, 1, ytrans_);
		} else if (xcl[2] == 'p') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = dsyMatrix_matmult(y, x, 0, xtrans_);
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = dsyMatrix_matmult(x, y, 1, ytrans_);
		}
	} else if (xcl[1] == 's' || ycl[1] == 's') {
		if (xcl[1] == 's') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = (ycl[2] != 'p')
				? dtrMatrix_matmult(y, x, 0, ytrans_, 0, 0)
				: dtpMatrix_matmult(y, x, 0, ytrans_, 0, 0);
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = (xcl[2] != 'p')
				? dtrMatrix_matmult(x, y, 1, xtrans_, 0, 0)
				: dtpMatrix_matmult(x, y, 1, xtrans_, 0, 0);
		}
	} else {
		SEXP
			xuplo = PROTECT(GET_SLOT(x, Matrix_uploSym)),
			yuplo = PROTECT(GET_SLOT(y, Matrix_uploSym)),
			xdiag = PROTECT(GET_SLOT(x, Matrix_diagSym)),
			ydiag = PROTECT(GET_SLOT(y, Matrix_diagSym));
		char
			xul = *CHAR(STRING_ELT(xuplo, 0)),
			yul = *CHAR(STRING_ELT(yuplo, 0)),
			xdi = *CHAR(STRING_ELT(xdiag, 0)),
			ydi = *CHAR(STRING_ELT(ydiag, 0));
		if (xtrans_)
			xul = (xul == 'U') ? 'L' : 'U';
		if (ytrans_)
			yul = (yul == 'U') ? 'L' : 'U';
		int triangular = (xul != yul) ? 0 : ((xdi != ydi || xdi == 'N') ? 1 : 2);
		UNPROTECT(4); /* ydiag, xdiag, yuplo, xuplo */

		if (xcl[2] == 'p' && ycl[2] == 'p') {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = dtpMatrix_matmult(x, y, 1, xtrans_, ytrans_, triangular);
		} else if (xcl[2] == 'p') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = dtrMatrix_matmult(y, x, 0, ytrans_, xtrans_, triangular);
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = dtrMatrix_matmult(x, y, 1, xtrans_, ytrans_, triangular);
		}
	}

	UNPROTECT(2); /* y, x */
	return x;
}

/* boolean: op(op(<.gC>) & op(<.gC>)) */
/* numeric: op(op(<dgC>) * op(<dgC>)) */
static
SEXP dgCMatrix_dgCMatrix_matmult(SEXP x, SEXP y, int xtrans, int ytrans,
                                 int ztrans, int triangular, int boolean)
{
	PROTECT_INDEX zpid;
	SEXP z;
	char zcl[] = "..CMatrix";
	zcl[0] = (boolean) ? 'n' : 'd';
	if (y == R_NilValue) {
		zcl[1] = 's';
		cholmod_sparse *X = M2CHS(x, !boolean);
		if (X->xtype == CHOLMOD_COMPLEX)
			error(_("'%s' does not support complex matrices"), "cholmod_aat");
		if (xtrans)
			X = cholmod_transpose(X, !boolean, &c);
		cholmod_sparse *Z = cholmod_aat(X, (int *) NULL, 0, !boolean, &c);
		if (xtrans)
			cholmod_free_sparse(&X, &c);
		if (!Z->sorted)
			cholmod_sort(Z, &c);
		X = cholmod_copy(Z, (ztrans) ? -1 : 1, 1, &c);
		cholmod_free_sparse(&Z, &c);
		Z = X;
		PROTECT_WITH_INDEX(z = CHS2M(Z, !boolean, zcl[1]), &zpid);
		cholmod_free_sparse(&Z, &c);
		SEXP xdimnames = PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
			zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
		symDN(zdimnames, xdimnames, (xtrans) ? 1 : 0);
		UNPROTECT(2); /* zdimnames, xdimnames */
		if (ztrans) {
			SEXP uplo = PROTECT(mkString("L"));
			SET_SLOT(z, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
	} else {
		zcl[1] = (triangular != 0) ? 't' : 'g';
		cholmod_sparse
			*X = M2CHS(x, !boolean),
			*Y = M2CHS(y, !boolean);
		if (X->xtype == CHOLMOD_COMPLEX || Y->xtype == CHOLMOD_COMPLEX)
			error(_("'%s' does not support complex matrices"), "cholmod_ssmult");
		if (((xtrans) ? X->nrow : X->ncol) != ((ytrans) ? Y->ncol : Y->nrow))
			error(_("non-conformable arguments"));
		if (xtrans)
			X = cholmod_transpose(X, !boolean, &c);
		if (ytrans)
			Y = cholmod_transpose(Y, !boolean, &c);
		cholmod_sparse *Z = cholmod_ssmult(X, Y, 0, !boolean, 1, &c);
		if (xtrans)
			cholmod_free_sparse(&X, &c);
		if (ytrans)
			cholmod_free_sparse(&Y, &c);
		PROTECT_WITH_INDEX(z = CHS2M(Z, !boolean, zcl[1]), &zpid);
		cholmod_free_sparse(&Z, &c);
		SEXP xdimnames = PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
			ydimnames = PROTECT(GET_SLOT(y, Matrix_DimNamesSym)),
			zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
		matmultDN(zdimnames,
		          xdimnames, (xtrans) ? 1 : 0,
		          ydimnames, (ytrans) ? 0 : 1);
		UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */
		if (triangular < 0) {
			SEXP uplo = PROTECT(mkString("L"));
			SET_SLOT(z, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
		if (triangular < -1 || triangular > 1)
			REPROTECT(z = sparse_diag_N2U(z, zcl), zpid);
	}
	if (ztrans)
		REPROTECT(z = sparse_transpose(z, zcl, 1), zpid);
	UNPROTECT(1); /* z */
	return z;
}

/* op(op(<d[gs]C>) * op(<dge>)) */
static
SEXP dgCMatrix_dgeMatrix_matmult(SEXP x, SEXP y, int xtrans, int ytrans,
                                 int ztrans, int triangular, int symmetric)
{
	SEXP z;
	char zcl[] = "...Matrix";
	cholmod_sparse *X = M2CHS(x, 1);
	cholmod_dense  *Y = M2CHD(y, ytrans);
	zcl[0] = (X->xtype == CHOLMOD_COMPLEX || Y->xtype == CHOLMOD_COMPLEX)
		? 'z' : 'd';
	zcl[1] = (triangular) ? 't' : 'g';
	zcl[2] = (triangular) ? 'r' : 'e';
	X->stype = symmetric;
	if (((xtrans) ? X->nrow : X->ncol) != Y->nrow) {
		if (ytrans)
			R_Free(Y->x);
		error(_("non-conformable arguments"));
	}
	int m = (int) ((xtrans) ? X->ncol : X->nrow), n = (int) Y->ncol;
	if ((Matrix_int_fast64_t) m * n > R_XLEN_T_MAX) {
		if (ytrans)
			R_Free(Y->x);
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	}
	cholmod_dense *Z = (cholmod_dense *) R_alloc(1, sizeof(cholmod_dense));
	memset(Z, 0, sizeof(cholmod_dense));
	Z->nrow = (size_t) m;
	Z->ncol = (size_t) n;
	Z->d = Z->nrow;
	Z->nzmax = Z->nrow * Z->ncol;
	Z->xtype = Y->xtype;
	Z->dtype = Y->dtype;
	double alpha[2] = { 1.0, 0.0 }, beta[2] = { 0.0, 0.0 };
	if (ztrans) {
		if (Z->xtype == CHOLMOD_COMPLEX)
			Z->x = R_Calloc(Z->nzmax, Rcomplex);
		else
			Z->x = R_Calloc(Z->nzmax, double);
		cholmod_sdmult(X, xtrans, alpha, beta, Y, Z, &c);
		PROTECT(z = CHD2M(Z, ztrans, zcl[1]));
		R_Free(Z->x);
	} else {
		PROTECT(z = newObject(zcl));
		SEXP zdim = GET_SLOT(z, Matrix_DimSym);
		INTEGER(zdim)[0] = m;
		INTEGER(zdim)[1] = n;
		SEXP zx;
		if (Z->xtype == CHOLMOD_COMPLEX) {
			PROTECT(zx = allocVector(CPLXSXP, (R_xlen_t) m * n));
			Z->x = COMPLEX(zx);
		} else {
			PROTECT(zx = allocVector(REALSXP, (R_xlen_t) m * n));
			Z->x = REAL(zx);
		}
		cholmod_sdmult(X, xtrans, alpha, beta, Y, Z, &c);
		SET_SLOT(z, Matrix_xSym, zx);
		UNPROTECT(1); /* zx */
	}
	if (ytrans)
		R_Free(Y->x);
	SEXP xdimnames = (symmetric)
		? PROTECT(get_symmetrized_DimNames(x, -1))
		: PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
		ydimnames = PROTECT(GET_SLOT(y, Matrix_DimNamesSym)),
		zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
	if (ztrans)
	matmultDN(zdimnames,
	          ydimnames, (ytrans) ? 0 : 1,
	          xdimnames, (xtrans) ? 1 : 0);
	else
	matmultDN(zdimnames,
	          xdimnames, (xtrans) ? 1 : 0,
	          ydimnames, (ytrans) ? 0 : 1);
	UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */
	if (triangular != 0 && ztrans == (triangular > 0)) {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(z, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (triangular < -1 || triangular > 1) {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(z, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}
	UNPROTECT(1); /* z */
	return z;
}

SEXP R_sparse_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans, SEXP ztrans,
                      SEXP boolean)
{
	if (TYPEOF(boolean) != LGLSXP || LENGTH(boolean) < 1)
		error(_("invalid '%s' to '%s'"), "boolean", __func__);
	int boolean_ = LOGICAL(boolean)[0];

	int xtrans_ = LOGICAL(xtrans)[0], ytrans_ = LOGICAL(ytrans)[0],
		ztrans_ = LOGICAL(ztrans)[0], m, n, v;
	matmultDim(x, y, &xtrans_, &ytrans_, &ztrans_, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(x, &xpid);
	PROTECT_WITH_INDEX(y, &ypid);

	if (TYPEOF(x) != S4SXP) {
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(x = matrix_as_dense( x, ",ge", '\0', '\0', xtrans_, 0), xpid);
		else if (!xtrans_)
		REPROTECT(x = matrix_as_sparse(x, "ngC", '\0', '\0', xtrans_   ), xpid);
		else
		REPROTECT(x = matrix_as_sparse(x, "ngR", '\0', '\0', xtrans_   ), xpid);
		if (v == 1) {
			/* Discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
			               (xtrans_) ? 1 : 0, R_NilValue);
			xtrans_ = 0;
		}
	}
	if (TYPEOF(y) != S4SXP && y != R_NilValue) {
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(y = matrix_as_dense( y, ",ge", '\0', '\0', ytrans_, 0), ypid);
		else if (!ytrans_)
		REPROTECT(y = matrix_as_sparse(y, "ngC", '\0', '\0', ytrans_   ), ypid);
		else
		REPROTECT(y = matrix_as_sparse(y, "ngR", '\0', '\0', ytrans_   ), ypid);
		if (v == 2) {
			/* Discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym),
			               (ytrans_) ? 1 : 0, R_NilValue);
			ytrans_ = 0;
		}
	}

	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, VALID_DENSE, "" };
	const char *xcl = NULL, *ycl = NULL;
	int ivalid;
	ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	xcl = valid[ivalid];
	if (y != R_NilValue) {
	ivalid = R_check_class_etc(y, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(y, __func__);
	ycl = valid[ivalid];
	}
	if (boolean_ == NA_LOGICAL)
		boolean_ = xcl[0] == 'n' && (y == R_NilValue || ycl[0] == 'n');
	char kind = (boolean_) ? 'n' :
		((xcl[0] == 'z' || (y != R_NilValue && ycl[0] == 'z')) ? 'z' : 'd');

	if (xcl[2] != 'C' && xtrans_) {
		if (xcl[2] != 'R' && xcl[2] != 'T') {
			REPROTECT(x = dense_as_sparse(x, xcl, 'R'), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (xcl[1] != 's' || xcl[2] != 'T') {
			REPROTECT(x = sparse_transpose(x, xcl, 1), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		xtrans_ = 0;
	}
	if (xcl[2] != 'C') {
		if (xcl[2] != 'R' && xcl[2] != 'T')
			REPROTECT(x = dense_as_sparse(x, xcl, 'C'), xpid);
		else
			REPROTECT(x = sparse_as_Csparse(x, xcl), xpid);
		xcl = valid[R_check_class_etc(x, valid)];
	}
	if (xcl[1] == 's')
		xtrans_ = 0;
	if (xcl[0] != kind) {
		if (boolean_)
			REPROTECT(x = sparse_drop0(x, xcl, 0.0), xpid);
		else {
			REPROTECT(x = sparse_as_kind(x, xcl, kind), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
	}

	if (y == R_NilValue) {
		REPROTECT(x = sparse_as_general(x, xcl), xpid);
		x = dgCMatrix_dgCMatrix_matmult(
			x, y, xtrans_, !xtrans_, ztrans_, 0, boolean_);
		UNPROTECT(2); /* y, x */
		return x;
	}

	int triangular = 0;
	if (xcl[1] == 't' && ycl[1] == 't') {
		SEXP
			xuplo = PROTECT(GET_SLOT(x, Matrix_uploSym)),
			yuplo = PROTECT(GET_SLOT(y, Matrix_uploSym)),
			xdiag = PROTECT(GET_SLOT(x, Matrix_diagSym)),
			ydiag = PROTECT(GET_SLOT(y, Matrix_diagSym));
		char
			xul = *CHAR(STRING_ELT(xuplo, 0)),
			yul = *CHAR(STRING_ELT(yuplo, 0)),
			xdi = *CHAR(STRING_ELT(xdiag, 0)),
			ydi = *CHAR(STRING_ELT(ydiag, 0));
		if (xtrans_)
			xul = (xul == 'U') ? 'L' : 'U';
		if (ytrans_)
			yul = (yul == 'U') ? 'L' : 'U';
		triangular = (xul != yul) ? 0 : ((xdi != ydi || xdi == 'N') ? 1 : 2);
		if (xul != 'U')
			triangular = -triangular;
		UNPROTECT(4); /* ydiag, xdiag, yuplo, xuplo */
	}

	if (!boolean_ && ycl[2] != 'C' && ycl[2] != 'R' && ycl[2] != 'T') {
		int symmetric = xcl[1] == 's';
		if (symmetric) {
			SEXP xuplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
			char xul = *CHAR(STRING_ELT(xuplo, 0));
			if (xul != 'U')
				symmetric = -1;
			UNPROTECT(1); /* xuplo */
		}
		if (ycl[0] != kind) {
			REPROTECT(y = dense_as_kind(y, ycl, kind, 0), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
		if (xcl[1] == 't')
			REPROTECT(x = sparse_diag_U2N(x, xcl), xpid);
		x = dgCMatrix_dgeMatrix_matmult(
			x, y, xtrans_, ytrans_, ztrans_, triangular, symmetric);
		UNPROTECT(2); /* y, x */
		return x;
	}

	if (ycl[2] != 'C' && ytrans_) {
		if (ycl[2] != 'R' && ycl[2] != 'T') {
			REPROTECT(y = dense_as_sparse(y, ycl, 'R'), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (ycl[1] != 's' || ycl[2] != 'T') {
			REPROTECT(y = sparse_transpose(y, ycl, 1), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		ytrans_ = 0;
	}
	if (ycl[2] != 'C') {
		if (ycl[2] != 'R' && ycl[2] != 'T')
			REPROTECT(y = dense_as_sparse(y, ycl, 'C'), ypid);
		else
			REPROTECT(y = sparse_as_Csparse(y, ycl), ypid);
		ycl = valid[R_check_class_etc(y, valid)];
	}
	if (ycl[1] == 's')
		ytrans_ = 0;
	if (ycl[0] != kind) {
		if (boolean_)
			REPROTECT(y = sparse_drop0(y, ycl, 0.0), ypid);
		else {
			REPROTECT(y = sparse_as_kind(y, ycl, kind), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
	}

	REPROTECT(x = sparse_as_general(x, xcl), xpid);
	REPROTECT(y = sparse_as_general(y, ycl), ypid);
	x = dgCMatrix_dgCMatrix_matmult(
		x, y, xtrans_, ytrans_, ztrans_, triangular, boolean_);
	UNPROTECT(2); /* y, x */
	return x;
}

#define MULTIPLY_COMPLEX(_X_, _D_) \
	do { \
		tmp = (_X_); \
		(_X_).r = tmp.r * (_D_).r - tmp.i * (_D_).i; \
		(_X_).i = tmp.r * (_D_).i + tmp.i * (_D_).r; \
	} while (0)
#define MULTIPLY_REAL(_X_, _D_) \
	(_X_) = (_X_)  * (_D_)
#define MULTIPLY_LOGICAL(_X_, _D_) \
	(_X_) = (_X_) && (_D_)

#define SCALE_CASES(_J_) \
	do { \
		switch (TYPEOF(d)) { \
		case CPLXSXP: \
		{ \
			Rcomplex tmp; \
			SCALE(Rcomplex, COMPLEX, MULTIPLY_COMPLEX, _J_); \
			break; \
		} \
		case REALSXP: \
			SCALE(double, REAL, MULTIPLY_REAL, _J_); \
			break; \
		default: \
			SCALE(int, LOGICAL, MULTIPLY_LOGICAL, _J_); \
			break; \
		} \
	} while (0)

static
void dense_colscale(SEXP obj, SEXP d, int m, int n, char uplo, char diag)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = XLENGTH(x) < (R_xlen_t) m * n;

#define SCALE(_CTYPE_, _PTR_, _OP_, _J_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pd = _PTR_(d); \
		if (uplo == '\0') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < m; ++i) { \
					_OP_(*px, pd[_J_]); \
					++px; \
				} \
			} \
		} else if (uplo == 'U') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i <= j; ++i) { \
					_OP_(*px, pd[_J_]); \
					++px; \
				} \
				if (!packed) \
					px += m - j - 1; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				for (i = j; i < m; ++i) { \
					_OP_(*px, pd[_J_]); \
					++px; \
				} \
			} \
		} \
		if (diag != '\0' && diag != 'N') { \
			px = _PTR_(x); \
			if (!packed) { \
				R_xlen_t m1a = (R_xlen_t) m + 1; \
				for (j = 0; j < n; ++j, px += m1a, pd += 1) \
					*px = *pd; \
			} else if (uplo == 'U') { \
				for (j = 0; j < n; px += (++j)+1, pd += 1) \
					*px = *pd; \
			} else { \
				for (j = 0; j < n; px += m-(j++), pd += 1) \
					*px = *pd; \
			} \
		} \
	} while (0)

	SCALE_CASES(j);
	return;
}

static
void dense_rowscale(SEXP obj, SEXP d, int m, int n, char uplo, char diag)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = XLENGTH(x) < (R_xlen_t) m * n;
	SCALE_CASES(i);

#undef SCALE

	return;
}

/* boolean: <lgC> & <ldi>  or  <ldi> & <dgR> */
/* numeric: <dgC> * <ddi>  or  <ddi> * <dgR> */
static
void Csparse_colscale(SEXP obj, SEXP d)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp = INTEGER(p) + 1, n = (int) (XLENGTH(p) - 1), j, k = 0, kend;
	UNPROTECT(2); /* p, x */

#define SCALE(_CTYPE_, _PTR_, _OP_, _J_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pd = _PTR_(d); \
		for (j = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				_OP_(*px, *pd); \
				++px; \
				++k; \
			} \
			++pd; \
		} \
	} while (0)

	SCALE_CASES();

#undef SCALE

	return;
}

/* boolean: <ldi> & <lgC>  or  <lgR> & <ldi> */
/* numeric: <ddi> * <dgC>  or  <dgR> * <ddi> */
static
void Csparse_rowscale(SEXP obj, SEXP d, SEXP iSym)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, iSym));
	int *pi = INTEGER(i), k, nnz = INTEGER(p)[XLENGTH(p) - 1];
	UNPROTECT(3); /* i, p, x */

#define SCALE(_CTYPE_, _PTR_, _OP_, _J_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pd = _PTR_(d); \
		for (k = 0; k < nnz; ++k) { \
			_OP_(*px, pd[*pi]); \
			++px; \
			++pi; \
		} \
	} while (0)

	SCALE_CASES();
	return;
}

/* boolean: <ldi> & <lgT>  or  <lgT> & <ldi> */
/* numeric: <ddi> * <dgT>  or  <dgT> * <ddi> */
static
void Tsparse_rowscale(SEXP obj, SEXP d, SEXP iSym)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		i = PROTECT(GET_SLOT(obj, iSym));
	int *pi = INTEGER(i);
	R_xlen_t k, nnz = XLENGTH(i);
	UNPROTECT(2); /* i, x */
	SCALE_CASES();

#undef SCALE

	return;
}

SEXP R_diagonal_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans,
                        SEXP boolean)
{
	SEXP x_ = x, y_ = y; /* for later pointer comparison */

	if (TYPEOF(boolean) != LGLSXP || LENGTH(boolean) < 1)
		error(_("invalid '%s' to '%s'"), "boolean", __func__);
	int boolean_ = LOGICAL(boolean)[0];

	int xtrans_ = LOGICAL(xtrans)[0], ytrans_ = LOGICAL(ytrans)[0],
		ztrans_ = 0, m, n, v;
	matmultDim(x, y, &xtrans_, &ytrans_, &ztrans_, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(x, &xpid);
	PROTECT_WITH_INDEX(y, &ypid);

	if (TYPEOF(x) != S4SXP) {
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(x = matrix_as_dense(x, ",ge", '\0', '\0', xtrans_, 2), xpid);
		else
		REPROTECT(x = matrix_as_dense(x, "nge", '\0', '\0', xtrans_, 2), xpid);
		if (v == 1) {
			/* Vector: discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
			               (xtrans_) ? 1 : 0, R_NilValue);
			xtrans_ = 0;
		}
	}
	if (TYPEOF(y) != S4SXP) {
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(y = matrix_as_dense(y, ",ge", '\0', '\0', ytrans_, 2), ypid);
		else
		REPROTECT(y = matrix_as_dense(y, "nge", '\0', '\0', ytrans_, 2), ypid);
		if (v == 2) {
			/* Vector: discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym),
			               (ytrans_) ? 1 : 0, R_NilValue);
			ytrans_ = 0;
		}
	}

	static const char *valid[] = {
		VALID_DIAGONAL,
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, VALID_DENSE, "" };
	const char *xcl = NULL, *ycl = NULL;
	int ivalid;
	ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	xcl = valid[ivalid];
	if (xcl[1] == 's')
		xtrans_ = 0;
	ivalid = R_check_class_etc(y, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(y, __func__);
	ycl = valid[ivalid];
	if (ycl[1] == 's')
		ytrans_ = 0;
	if (boolean_ == NA_LOGICAL)
		boolean_ = xcl[0] == 'n' && ycl[0] == 'n';
	char kind = (boolean_) ? 'n' :
		((xcl[0] == 'z' || ycl[0] == 'z') ? 'z' : 'd');

	int margin = -1, unit = -1;
	if (xcl[2] == 'i') {
		margin = 0;
		unit = *CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0)) != 'N';
	} else if (ycl[2] == 'i') {
		margin = 1;
		unit = *CHAR(STRING_ELT(GET_SLOT(y, Matrix_diagSym), 0)) != 'N';
	} else
		error(_("should never happen ..."));

	char ks = (boolean_) ? 'l' : kind, kd = kind;
	switch (xcl[2]) {
	case 'i':
		if (!unit && xcl[0] != ks) {
			REPROTECT(x = diagonal_as_kind(x, xcl, ks), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		break;
	case 'C':
	case 'R':
	case 'T':
		if (xcl[0] != ks) {
			REPROTECT(x = sparse_as_kind(x, xcl, ks), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (!unit && xcl[1] == 's') {
			REPROTECT(x = sparse_as_general(x, xcl), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		} else if (!unit && xcl[1] == 't')
			REPROTECT(x = sparse_diag_U2N(x, xcl), xpid);
		if (xtrans_) {
			REPROTECT(x = sparse_transpose(x, xcl, 0), xpid);
			xtrans_ = 0;
		}
		break;
	default:
		if (xcl[0] != kd) {
			REPROTECT(x = dense_as_kind(x, xcl, kd, 1), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (!unit && xcl[1] == 's') {
			REPROTECT(x = dense_as_general(x, xcl, x == x_), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (xtrans_) {
			REPROTECT(x = dense_transpose(x, xcl), xpid);
			xtrans_ = 0;
		}
		break;
	}
	switch (ycl[2]) {
	case 'i':
		if (!unit && ycl[0] != ks) {
			REPROTECT(y = diagonal_as_kind(y, ycl, ks), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		break;
	case 'C':
	case 'R':
	case 'T':
		if (ycl[0] != ks) {
			REPROTECT(y = sparse_as_kind(y, ycl, ks), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (!unit && ycl[1] == 's') {
			REPROTECT(y = sparse_as_general(y, ycl), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		} else if (!unit && ycl[1] == 't')
			REPROTECT(y = sparse_diag_U2N(y, ycl), ypid);
		if (ytrans_) {
			REPROTECT(y = sparse_transpose(y, ycl, 0), ypid);
			ytrans_ = 0;
		}
		break;
	default:
		if (ycl[0] != kd) {
			REPROTECT(y = dense_as_kind(y, ycl, kd, 1), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (!unit && ycl[1] == 's') {
			REPROTECT(y = dense_as_general(y, ycl, y == y_), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (ytrans_) {
			REPROTECT(y = dense_transpose(y, ycl), ypid);
			ytrans_ = 0;
		}
		break;
	}

	SEXP z;
	PROTECT_INDEX zpid;
	const char *zcl = (margin == 0) ? ycl : xcl;
	PROTECT_WITH_INDEX(z = newObject(zcl), &zpid);

	SEXP zdim = PROTECT(GET_SLOT(z, Matrix_DimSym));
	int *pzdim = INTEGER(zdim);
	pzdim[0] = m;
	pzdim[1] = n;
	UNPROTECT(1); /* zdim */

	SEXP xdimnames = PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
		ydimnames = PROTECT(GET_SLOT(y, Matrix_DimNamesSym)),
		zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
	matmultDN(zdimnames,
	          xdimnames, (xtrans_) ? 1 : 0,
	          ydimnames, (ytrans_) ? 0 : 1);
	UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */

	char ul = '\0', di = '\0';
	if (zcl[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT((margin == 0) ? y : x, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(z, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		if (zcl[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT((margin == 0) ? y : x, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N' && unit)
			SET_SLOT(z, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
		}
	}

	if (zcl[2] == 'C' || zcl[2] == 'R' || zcl[2] == 'T') {
		if (zcl[2] != 'T') {
			SEXP p = PROTECT(GET_SLOT((margin == 0) ? y : x, Matrix_pSym));
			SET_SLOT(z, Matrix_pSym, p);
			UNPROTECT(1); /* p */
		}
		if (zcl[2] != 'R') {
			SEXP i = PROTECT(GET_SLOT((margin == 0) ? y : x, Matrix_iSym));
			SET_SLOT(z, Matrix_iSym, i);
			UNPROTECT(1); /* i */
		}
		if (zcl[2] != 'C') {
			SEXP j = PROTECT(GET_SLOT((margin == 0) ? y : x, Matrix_jSym));
			SET_SLOT(z, Matrix_jSym, j);
			UNPROTECT(1); /* j */
		}
	}

	SEXP x0 = PROTECT(GET_SLOT((margin == 0) ? y : x, Matrix_xSym));
	if (unit || ((margin == 0) ? y != y_ : x != x_))
		SET_SLOT(z, Matrix_xSym, x0);
	else {
		SEXP x1 = PROTECT(allocVector(TYPEOF(x0), XLENGTH(x0)));
		switch (kind) {
		case 'z':
			Matrix_memcpy(COMPLEX(x1), COMPLEX(x0), XLENGTH(x0), sizeof(Rcomplex));
			break;
		case 'd':
			Matrix_memcpy(   REAL(x1),    REAL(x0), XLENGTH(x0), sizeof(  double));
			break;
		default:
			Matrix_memcpy(LOGICAL(x1), LOGICAL(x0), XLENGTH(x0), sizeof(     int));
			break;
		}
		SET_SLOT(z, Matrix_xSym, x1);
		UNPROTECT(1); /* x1 */
	}
	UNPROTECT(1); /* x0 */

	if (!unit) {
		SEXP d = PROTECT(GET_SLOT((margin == 0) ? x : y, Matrix_xSym));
		switch (zcl[2]) {
		case 'C':
			if (margin == 0)
				Csparse_rowscale(z, d, Matrix_iSym);
			else
				Csparse_colscale(z, d);
			break;
		case 'R':
			if (margin == 0)
				Csparse_colscale(z, d);
			else
				Csparse_rowscale(z, d, Matrix_jSym);
			break;
		case 'T':
			if (margin == 0)
				Tsparse_rowscale(z, d, Matrix_iSym);
			else
				Tsparse_rowscale(z, d, Matrix_jSym);
			break;
		default:
			if (margin == 0)
				  dense_rowscale(z, d, m, n, ul, di);
			else
				  dense_colscale(z, d, m, n, ul, di);
			break;
		}
		UNPROTECT(1); /* d */
	}

	if (boolean_ && (zcl[2] == 'C' || zcl[2] == 'R' || zcl[2] == 'T')) {
		REPROTECT(z = sparse_drop0(z, zcl, 0.0), zpid);
		REPROTECT(z = sparse_as_kind(z, zcl, 'n'), zpid);
	}

	UNPROTECT(3); /* z, y, x */
	return z;
}
