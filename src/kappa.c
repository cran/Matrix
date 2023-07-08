#include "kappa.h"

static char La_norm_type(SEXP s)
{
#define ARGNAME "type"
	if (TYPEOF(s) != STRSXP)
		error(_("argument '%s' is not of type \"character\""), ARGNAME);
	if (LENGTH(s) == 0)
		error(_("argument '%s' has length 0"), ARGNAME);
	const char *type = CHAR(STRING_ELT(s, 0));
	if (type[0] == '\0' || type[1] != '\0')
		error(_("argument '%s' (\"%s\") does not have string length 1"),
		      ARGNAME, type);
	char type_ = '\0';
	switch (type[0]) {
	case 'M':
	case 'm':
		type_ = 'M';
		break;
	case 'O':
	case 'o':
	case '1':
		type_ = 'O';
		break;
	case 'I':
	case 'i':
		type_ = 'I';
		break;
	case 'F':
	case 'f':
	case 'E':
	case 'e':
		type_ = 'F';
		break;
	default:
		error(_("argument '%s' (\"%s\") is not \"M\", \"O\", \"1\", \"I\", \"F\", or \"E\""),
		      ARGNAME, type);
		break;
	}
	return type_;
#undef ARGNAME
}

static char La_rcond_type(SEXP s)
{
#define ARGNAME "norm"
	if (TYPEOF(s) != STRSXP)
		error(_("argument '%s' is not of type \"character\""), ARGNAME);
	if (LENGTH(s) == 0)
		error(_("argument '%s' has length 0"), ARGNAME);
	const char *type = CHAR(STRING_ELT(s, 0));
	if (type[0] == '\0' || type[1] != '\0')
		error(_("argument '%s' (\"%s\") does not have string length 1"),
		      ARGNAME, type);
	char type_ = '\0';
	switch (type[0]) {
	case 'O':
	case 'o':
	case '1':
		type_ = 'O';
		break;
	case 'I':
	case 'i':
		type_ = 'I';
		break;
	default:
		error(_("argument '%s' (\"%s\") is not \"O\", \"1\", or \"I\""),
		      ARGNAME, type);
		break;
	}
	return type_;
#undef ARGNAME
}

SEXP dgeMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */
	if (m == 0 || n == 0)
		return ScalarReal(0.0);

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'I')
		work = (double *) R_alloc((size_t) m, sizeof(double));
	norm =
	F77_CALL(dlange)(&type_, &m, &n, REAL(x), &m,
	                 work FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dgeMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	UNPROTECT(1); /* dim */
	if (m != n)
		error(_("rcond(x) is undefined: 'x' is not square"));
	if (n == 0)
		error(_("rcond(x) is undefined: 'x' has length 0"));

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond,
		*work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	int info, *iwork = (int *) R_alloc((size_t) n, sizeof(int));
	norm =
	F77_CALL(dlange)(&type_, &n, &n, REAL(x), &n,
	                 work FCONE);
	F77_CALL(dgecon)(&type_, &n,     REAL(y), &n, &norm, &rcond,
	                 work, iwork, &info FCONE);
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}

SEXP dtrMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
		diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0],
		diag_ = CHAR(STRING_ELT(diag, 0))[0];
	UNPROTECT(2); /* uplo, diag */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	norm =
	F77_CALL(dlantr)(&type_, &uplo_, &diag_, &n, &n, REAL(x), &n,
	                 work FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dtrMatrix_rcond(SEXP obj, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		error(_("rcond(x) is undefined: 'x' has length 0"));

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
		diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0],
		diag_ = CHAR(STRING_ELT(diag, 0))[0];
	UNPROTECT(2); /* uplo, diag */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double rcond,
		*work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int info, *iwork = (int *) R_alloc((size_t) n, sizeof(int));
	F77_CALL(dtrcon)(&type_, &uplo_, &diag_, &n, REAL(x), &n, &rcond,
	                 work, iwork, &info FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(rcond);
}

SEXP dtpMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
		diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0],
		diag_ = CHAR(STRING_ELT(diag, 0))[0];
	UNPROTECT(2); /* uplo, diag */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	norm =
	F77_CALL(dlantp)(&type_, &uplo_, &diag_, &n, REAL(x),
	                 work FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dtpMatrix_rcond(SEXP obj, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		error(_("rcond(x) is undefined: 'x' has length 0"));

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
		diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0],
		diag_ = CHAR(STRING_ELT(diag, 0))[0];
	UNPROTECT(2); /* uplo, diag */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double rcond,
		*work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int info, *iwork = (int *) R_alloc((size_t) n, sizeof(int));
	F77_CALL(dtpcon)(&type_, &uplo_, &diag_, &n, REAL(x), &rcond,
	                 work, iwork, &info FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(rcond);
}

SEXP dsyMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];
	UNPROTECT(1); /* uplo */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'O' || type_ == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	norm =
	F77_CALL(dlansy)(&type_, &uplo_, &n, REAL(x), &n,
	                 work FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dsyMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		error(_("rcond(x) is undefined: 'x' has length 0"));

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];
	UNPROTECT(1); /* uplo */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym)),
		pivot = PROTECT(GET_SLOT(trf, Matrix_permSym));
	double norm, rcond,
		*work = (double *) R_alloc((size_t) 2 * n, sizeof(double));
	int info, *iwork = (int *) R_alloc((size_t) n, sizeof(int));
	norm =
	F77_CALL(dlansy)(&type_, &uplo_, &n, REAL(x), &n,
	                 work FCONE FCONE);
	F77_CALL(dsycon)(        &uplo_, &n, REAL(y), &n,
	                 INTEGER(pivot), &norm, &rcond,
	                 work, iwork, &info FCONE);
	UNPROTECT(3); /* x, y, pivot */

	return ScalarReal(rcond);
}

SEXP dspMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];
	UNPROTECT(1); /* uplo */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'O' || type_ == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	norm =
	F77_CALL(dlansp)(&type_, &uplo_, &n, REAL(x),
	                 work FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dspMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		error(_("rcond(x) is undefined: 'x' has length 0"));

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];
	UNPROTECT(1); /* uplo */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym)),
		pivot = PROTECT(GET_SLOT(trf, Matrix_permSym));
	double norm, rcond,
		*work = (double *) R_alloc((size_t) 2 * n, sizeof(double));
	int info, *iwork = (int *) R_alloc((size_t) n, sizeof(int));
	norm =
	F77_CALL(dlansp)(&type_, &uplo_, &n, REAL(x),
	                 work FCONE FCONE);
	F77_CALL(dspcon)(        &uplo_, &n, REAL(y),
	                 INTEGER(pivot), &norm, &rcond,
	                 work, iwork, &info FCONE);
	UNPROTECT(3); /* x, y, pivot */

	return ScalarReal(rcond);
}

SEXP dpoMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		error(_("rcond(x) is undefined: 'x' has length 0"));

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];
	UNPROTECT(1); /* uplo */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond,
		*work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int info, *iwork = (int *) R_alloc((size_t) n, sizeof(int));
	norm =
	F77_CALL(dlansy)(&type_, &uplo_, &n, REAL(x), &n,
	                 work FCONE FCONE);
	F77_CALL(dpocon)(        &uplo_, &n, REAL(y), &n, &norm, &rcond,
	                 work, iwork, &info FCONE);
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}

SEXP dppMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	if (n == 0)
		error(_("rcond(x) is undefined: 'x' has length 0"));

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];
	UNPROTECT(1); /* uplo */

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond,
		*work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int info, *iwork = (int *) R_alloc((size_t) n, sizeof(int));
	norm =
	F77_CALL(dlansp)(&type_, &uplo_, &n, REAL(x),
	                 work FCONE FCONE);
	F77_CALL(dppcon)(        &uplo_, &n, REAL(y), &norm, &rcond,
	                 work, iwork, &info FCONE);
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}
