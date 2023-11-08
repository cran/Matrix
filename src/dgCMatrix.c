#include "Mdefines.h"
#include "cs-etc.h"
#include "cholmod-etc.h"
#include "dgCMatrix.h"

/* TODO: support NCOL(b) > 1                   */
SEXP dgCMatrix_lusol(SEXP a, SEXP b)
{
	Matrix_cs *A = M2CXS(a, 1);
	MCS_XTYPE_SET(MCS_REAL);
	PROTECT(b = (TYPEOF(b) == REALSXP) ?
		duplicate(b) : coerceVector(b, REALSXP));
	if (A->m != A->n || A->m <= 0)
		error(_("'%s' is empty or not square"), "a");
	if (LENGTH(b) != A->m)
		error(_("dimensions of '%s' and '%s' are inconsistent"), "a", "b");
	if (!Matrix_cs_lusol(1, A, REAL(b), 1e-07))
		error(_("'%s' failed"), "cs_lusol");
	UNPROTECT(1);
	return b;
}

/* called from package MatrixModels's R code : */
/* TODO: support NCOL(b) > 1                   */
/* TODO: give result list(L, coef, Xty, resid) */
SEXP dgCMatrix_qrsol(SEXP a, SEXP b, SEXP order)
{
	/* FIXME? 'cs_qrsol' supports underdetermined systems.  */
	/*        We should only require LENGTH(b) = max(m, n). */
	int order_ = asInteger(order);
	if (order_ < 0 || order_ > 3)
		order_ = 0;
	Matrix_cs *A = M2CXS(a, 1);
	MCS_XTYPE_SET(MCS_REAL);
	PROTECT(b = (TYPEOF(b) == REALSXP)
		? duplicate(b) : coerceVector(b, REALSXP));
	if (LENGTH(b) != A->m)
		error(_("dimensions of '%s' and '%s' are inconsistent"), "a", "b");
	if (A->n <= 0 || A->n > A->m)
		error(_("%s(%s, %s) requires m-by-n '%s' with m >= n > 0"),
		      "dgCMatrix_qrsol", "a", "b", "a");
	if (!Matrix_cs_qrsol(order_, A, REAL(b)))
		error(_("'%s' failed"), "cs_qrsol");
	if (A->n < A->m) {
		SEXP tmp = allocVector(REALSXP, A->n);
		Matrix_memcpy(REAL(tmp), REAL(b), A->n, sizeof(double));
		b = tmp;
	}
	UNPROTECT(1);
	return b;
}

/* called from package MatrixModels's R code : */
/* TODO: support NCOL(b) > 1                   */
SEXP dgCMatrix_cholsol(SEXP at, SEXP b)
{
	/* Find least squares solution of A * X = B, given A' and B : */
	cholmod_sparse *At = M2CHS(at, 1);
	PROTECT(b = coerceVector(b, REALSXP));
	if (LENGTH(b) != At->ncol)
		error(_("dimensions of '%s' and '%s' are inconsistent"), "at", "b");
	if (At->ncol <= 0 || At->ncol < At->nrow)
		error(_("%s(%s, %s) requires m-by-n '%s' with n >= m > 0"),
		      "dgCMatrix_cholsol", "at", "b", "at");
	double zero[] = { 0.0, 0.0 }, one[] = {1.0, 0.0}, mone[] = { -1.0, 0.0 };

	/* L * L' = A' * A */
	cholmod_factor *L = cholmod_analyze(At, &c);
	if (!cholmod_factorize(At, L, &c))
		error(_("'%s' failed"), "cholmod_factorize");

	cholmod_dense *B = (cholmod_dense *) R_alloc(1, sizeof(cholmod_dense));
	memset(B, 0, sizeof(cholmod_dense));
	B->nrow = B->d = B->nzmax = LENGTH(b);
	B->ncol = 1;
	B->x = REAL(b);
	B->dtype = CHOLMOD_DOUBLE;
	B->xtype = CHOLMOD_REAL;

	/* A' * B = 1 * A' * B + 0 * <in/out> */
	cholmod_dense *AtB = cholmod_allocate_dense(
		At->nrow, 1, At->nrow, CHOLMOD_REAL, &c);
	if (!cholmod_sdmult(At, 0, one, zero, B, AtB, &c))
		error(_("'%s' failed"), "cholmod_sdmult");

	/* C := solve(A' * A, A' * B) = solve(L', solve(L, A' * B)) */
	cholmod_dense *C = cholmod_solve(CHOLMOD_A, L, AtB, &c);
	if (!C)
		error(_("'%s' failed"), "cholmod_solve");

	/* R := A * A' * C - B = 1 * (A')' * A' * X + (-1) * B */
	cholmod_dense *R = cholmod_copy_dense(B, &c);
	if (!cholmod_sdmult(At, 1, mone, one, C, R, &c))
		error(_("'%s' failed"), "cholmod_sdmult");

	const char *nms[] = {"L", "coef", "Xty", "resid", ""};
	SEXP ans = PROTECT(Rf_mkNamed(VECSXP, nms)), tmp;
	/* L : */
	PROTECT(tmp = CHF2M(L, 1));
	SET_VECTOR_ELT(ans, 0, tmp);
	/* coef : */
	PROTECT(tmp = allocVector(REALSXP, At->nrow));
	Matrix_memcpy(REAL(tmp),   C->x, At->nrow, sizeof(double));
	SET_VECTOR_ELT(ans, 1, tmp);
	/* Xty : */
	PROTECT(tmp = allocVector(REALSXP, At->nrow));
	Matrix_memcpy(REAL(tmp), AtB->x, At->nrow, sizeof(double));
	SET_VECTOR_ELT(ans, 2, tmp);
	/* resid : */
	PROTECT(tmp = allocVector(REALSXP, At->ncol));
	Matrix_memcpy(REAL(tmp),   R->x, At->ncol, sizeof(double));
	SET_VECTOR_ELT(ans, 3, tmp);

	cholmod_free_factor(&  L, &c);
	cholmod_free_dense (&AtB, &c);
	cholmod_free_dense (&  C, &c);
	cholmod_free_dense (&  R, &c);

	UNPROTECT(6);
	return ans;
}
