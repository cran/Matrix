#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#ifndef R_MATRIX_INLINE
# define R_MATRIX_INLINE
#endif

/* ==== cholmod.h =================================================== */

#include "cholmod.h"

#ifdef __cplusplus
extern "C" {
#endif

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(aat)(CHM_SP A, int *fset, size_t fsize, int mode,
                      CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, int *, size_t, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, int *, size_t, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_aat");
	return fun(A, fset, fsize, mode, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(add)(CHM_SP A, CHM_SP B, double alpha[2], double beta[2],
                      int values, int sorted, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, CHM_SP, double[2], double[2],
	                    int, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, CHM_SP, double[2], double[2],
		                 int, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_add");
	return fun(A, B, alpha, beta, values, sorted, Common);
}

R_MATRIX_INLINE CHM_DN attribute_hidden
R_MATRIX_CHOLMOD(allocate_dense)(size_t nrow, size_t ncol, size_t d, int xtype,
                                 CHM_CM Common)
{
	static CHM_DN(*fun)(size_t, size_t, size_t, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(size_t, size_t, size_t, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_allocate_dense");
	return fun(nrow, ncol, d, xtype, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(allocate_sparse)(size_t nrow, size_t ncol, size_t nzmax,
                                  int sorted, int packed, int stype, int xtype,
                                  CHM_CM Common)
{
	static CHM_SP(*fun)(size_t, size_t, size_t, int, int, int, int,
	                    CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(size_t, size_t, size_t, int, int, int, int,
		                 CHM_CM))
			R_GetCCallable("Matrix", "cholmod_allocate_sparse");
	return fun(nrow, ncol, nzmax, sorted, packed, stype, xtype, Common);
}

R_MATRIX_INLINE CHM_TR attribute_hidden
R_MATRIX_CHOLMOD(allocate_triplet)(size_t nrow, size_t ncol, size_t nzmax,
                                   int stype, int xtype, CHM_CM Common)
{
	static CHM_TR(*fun)(size_t, size_t, size_t, int, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_TR(*)(size_t, size_t, size_t, int, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_allocate_triplet");
	return fun(nrow, ncol, nzmax, stype, xtype, Common);
}

R_MATRIX_INLINE CHM_FR attribute_hidden
R_MATRIX_CHOLMOD(analyze)(CHM_SP A, CHM_CM Common)
{
	static CHM_FR(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_SP,CHM_CM))
			R_GetCCallable("Matrix", "cholmod_analyze");
	return fun(A, Common);
}

R_MATRIX_INLINE CHM_FR attribute_hidden
R_MATRIX_CHOLMOD(analyze_p)(CHM_SP A, int *Perm, int *fset, size_t fsize,
                            CHM_CM Common)
{
	static CHM_FR(*fun)(CHM_SP, int *, int *, size_t, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_SP, int *, int *, size_t, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_analyze_p");
	return fun(A, Perm, fset, fsize, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(band_inplace)(int k1, int k2, int mode, CHM_SP A,
                               CHM_CM Common)
{
	static int(*fun)(int, int, int, CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(int, int, int, CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_band_inplace");
	return fun(k1, k2, mode, A, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(change_factor)(int to_xtype, int to_ll, int to_super,
                                int to_packed, int to_monotonic,
                                CHM_FR L, CHM_CM Common)
{
	static int(*fun)(int, int, int, int, int, CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
	fun = (int(*)(int, int, int, int, int, CHM_FR, CHM_CM))
		R_GetCCallable("Matrix", "cholmod_change_factor");
	return fun(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(copy)(CHM_SP A, int stype, int mode, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, int, int, CHM_CM) = NULL;
	if (fun == NULL)
	fun = (CHM_SP(*)(CHM_SP, int, int, CHM_CM))
		R_GetCCallable("Matrix", "cholmod_copy");
	return fun(A, stype, mode, Common);
}

R_MATRIX_INLINE CHM_DN attribute_hidden
R_MATRIX_CHOLMOD(copy_dense)(CHM_DN  A, CHM_CM Common)
{
	static CHM_DN(*fun)(CHM_DN, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(CHM_DN, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_copy_dense");
	return fun(A, Common);
}

R_MATRIX_INLINE CHM_FR attribute_hidden
R_MATRIX_CHOLMOD(copy_factor)(CHM_FR L, CHM_CM Common)
{
	static CHM_FR(*fun)(CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_copy_factor");
	return fun(L, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(copy_sparse)(CHM_SP A, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_copy_sparse");
	return fun(A, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(defaults)(CHM_CM Common)
{
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_CM))
			R_GetCCallable("Matrix", "cholmod_defaults");
	return fun(Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(dense_to_sparse)(CHM_DN X, int values, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_DN, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_DN, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_dense_to_sparse");
	return fun(X, values, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(factor_to_sparse)(CHM_FR L, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_factor_to_sparse");
	return fun(L, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(factorize)(CHM_SP A, CHM_FR L, CHM_CM Common)
{
	static int(*fun)(CHM_SP, CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_factorize");
	return fun(A, L, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(factorize_p)(CHM_SP A, double beta[2], int *fset,
                      size_t fsize, CHM_FR L, CHM_CM Common)
{
	static int(*fun)(CHM_SP, double[2], int *, size_t, CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, double[2], int *, size_t, CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_factorize_p");
	return fun(A, beta, fset, fsize, L, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(finish)(CHM_CM Common)
{
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_CM))
			R_GetCCallable("Matrix", "cholmod_finish");
	return fun(Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(free_dense)(CHM_DN *A, CHM_CM Common)
{
	static int(*fun)(CHM_DN *, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_DN *,CHM_CM))
			R_GetCCallable("Matrix", "cholmod_free_dense");
	return fun(A, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(free_factor)(CHM_FR *L, CHM_CM Common)
{
	static int(*fun)(CHM_FR *,CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_FR *, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_free_factor");
	return fun(L, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(free_sparse)(CHM_SP *A, CHM_CM Common)
{
	static int(*fun)(CHM_SP *, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP *, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_free_sparse");
	return fun(A, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(free_triplet)(CHM_TR *T, CHM_CM Common)
{
	static int(*fun)(CHM_TR *, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_TR *,CHM_CM))
			R_GetCCallable("Matrix", "cholmod_free_triplet");
	return fun(T, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(nnz)(CHM_SP A, CHM_CM Common)
{
	static int(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_nnz");
	return fun(A, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(scale)(CHM_DN S, int scale, CHM_SP A, CHM_CM Common)
{
	static int(*fun)(CHM_DN, int, CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_DN, int, CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_scale");
	return fun(S, scale, A, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(sdmult)(CHM_SP A, int transpose,
                         double alpha[2], double beta[2],
                         CHM_DN X, CHM_DN Y, CHM_CM Common)
{
	static int(*fun)(CHM_SP, int, double[2], double[2],
	                 CHM_DN, CHM_DN, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, int, double[2], double[2],
		              CHM_DN, CHM_DN, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_sdmult");
	return fun(A, transpose, alpha, beta, X, Y, Common);
}

R_MATRIX_INLINE CHM_DN attribute_hidden
R_MATRIX_CHOLMOD(solve)(int sys, CHM_FR L, CHM_DN B, CHM_CM Common)
{
	static CHM_DN(*fun)(int, CHM_FR, CHM_DN, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(int, CHM_FR, CHM_DN, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_solve");
	return fun(sys, L, B, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(solve2)(int sys, CHM_FR L, CHM_DN B,
                         CHM_DN *X_Handle, CHM_DN *Y_Handle, CHM_DN *E_Handle,
                         CHM_CM Common)
{
	static int(*fun)(int, CHM_FR, CHM_DN, CHM_SP,
	                 CHM_DN *, CHM_SP *, CHM_DN *, CHM_DN *, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(int, CHM_FR, CHM_DN, CHM_SP,
		              CHM_DN *, CHM_SP *, CHM_DN *, CHM_DN *, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_solve2");
	return fun(sys, L, B, NULL, X_Handle, NULL, Y_Handle, E_Handle, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(sort)(CHM_SP A, CHM_CM Common)
{
	static int(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_sort");
	return fun(A, Common);
}

R_MATRIX_INLINE CHM_DN attribute_hidden
R_MATRIX_CHOLMOD(sparse_to_dense)(CHM_SP A, CHM_CM Common)
{
	static CHM_DN(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_sparse_to_dense");
	return fun(A, Common);
}

R_MATRIX_INLINE CHM_TR attribute_hidden
R_MATRIX_CHOLMOD(sparse_to_triplet)(CHM_SP A, CHM_CM Common)
{
	static CHM_TR(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_TR(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_sparse_to_triplet");
	return fun(A, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(speye)(size_t nrow, size_t ncol, int xtype, CHM_CM Common)
{
	static CHM_SP(*fun)(size_t, size_t, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(size_t, size_t, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_speye");
	return fun(nrow, ncol, xtype, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(spsolve)(int sys, CHM_FR L, CHM_SP B, CHM_CM Common)
{
	static CHM_SP(*fun)(int, CHM_FR, CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(int,CHM_FR, CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_spsolve");
	return fun(sys, L, B, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(ssmult)(CHM_SP A, CHM_SP B,
                 int stype, int values, int sorted, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, CHM_SP, int, int, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, CHM_SP, int, int, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_ssmult");
	return fun(A, B, stype, values, sorted, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(submatrix)(CHM_SP A, int *rset, int rsize, int *cset,
                            int csize, int values, int sorted, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, int *, int, int *,
	                    int, int, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, int *, int, int *,
		                 int, int, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_submatrix");
	return fun(A, rset, rsize, cset, csize, values, sorted, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(transpose)(CHM_SP A, int values, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_transpose");
	return fun(A, values, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(triplet_to_sparse)(CHM_TR T, int nzmax, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_TR, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_TR, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_triplet_to_sparse");
	return fun(T, nzmax, Common);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(updown)(int update, CHM_SP C, CHM_FR L, CHM_CM Common)
{
	static int(*fun)(int, CHM_SP, CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(int, CHM_SP, CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_updown");
	return fun(update, C, L, Common);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
R_MATRIX_CHOLMOD(vertcat)(CHM_SP A, CHM_SP B, int values, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, CHM_SP, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, CHM_SP, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_vertcat");
	return fun(A, B, values, Common);
}


/* ---- cholmod_start ----------------------------------------------- */
/* NB: keep synchronized with analogues in ../../src/chm_common.c     */

#if 0
static int attribute_hidden
R_MATRIX_CHOLMOD(print_function)(const char *fmt, ...)
{
	va_list(ap);
	va_start(ap, fmt);
	Rprintf((char *) fmt, ap);
	va_end(ap);
	return 0;
}
#endif

R_MATRIX_INLINE void attribute_hidden
R_MATRIX_CHOLMOD(error_handler)(int status, const char *file, int line,
                                const char *message)
{
	/* NB: Matrix itself uses cholmod_common_env(ini|set|get) to preserve
	   settings through error calls.  Consider defining *your* own error
	   handler and restoring the instance of cholmod_common that *you* use.
	*/

	if (status < 0)
		Rf_error("CHOLMOD error '%s' at file '%s', line %d",
		         message, file, line);
	else
		Rf_warning("CHOLMOD warning '%s' at file '%s', line %d",
		           message, file, line);
}

R_MATRIX_INLINE int attribute_hidden
R_MATRIX_CHOLMOD(start)(CHM_CM Common)
{
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_CM))
			R_GetCCallable("Matrix", "cholmod_start");
	int ans = fun(Common);
#if 0
	/* No longer, with SuiteSparse 5.7.1 : */
	Common->print_function =
# if 0
		R_MATRIX_CHOLMOD(print_function);
# else
		NULL;
# endif
#endif
	Common->error_handler = R_MATRIX_CHOLMOD(error_handler);
	return ans;
}

#ifdef	__cplusplus
}
#endif


/* ==== cholmod-utils.h ============================================= */

#ifndef R_MATRIX_NO_CHOLMOD_UTILS

#include "cholmod-utils.h"

#ifdef __cplusplus
extern "C" {
#endif

R_MATRIX_INLINE CHM_FR attribute_hidden
M_sexp_as_cholmod_factor(CHM_FR L, SEXP from)
{
	static CHM_FR(*fun)(CHM_FR, SEXP) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_FR, SEXP))
			R_GetCCallable("Matrix", "sexp_as_cholmod_factor");
	return fun(L, from);
}

R_MATRIX_INLINE CHM_SP attribute_hidden
M_sexp_as_cholmod_sparse(CHM_SP A, SEXP from,
                         Rboolean checkUnit, Rboolean sortInPlace)
{
	static CHM_SP(*fun)(CHM_SP, SEXP, Rboolean, Rboolean) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, SEXP, Rboolean, Rboolean))
			R_GetCCallable("Matrix", "sexp_as_cholmod_sparse");
	return fun(A, from, checkUnit, sortInPlace);
}

R_MATRIX_INLINE CHM_TR attribute_hidden
M_sexp_as_cholmod_triplet(CHM_TR A, SEXP from,
                          Rboolean checkUnit)
{
	static CHM_TR(*fun)(CHM_TR, SEXP, Rboolean) = NULL;
	if (fun == NULL)
		fun = (CHM_TR(*)(CHM_TR, SEXP, Rboolean))
			R_GetCCallable("Matrix", "sexp_as_cholmod_triplet");
	return fun(A, from, checkUnit);
}

R_MATRIX_INLINE CHM_DN attribute_hidden
M_sexp_as_cholmod_dense(CHM_DN A, SEXP from)
{
	static CHM_DN(*fun)(CHM_DN, SEXP) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(CHM_DN, SEXP))
			R_GetCCallable("Matrix", "sexp_as_cholmod_dense");
	return fun(A, from);
}

R_MATRIX_INLINE CHM_DN attribute_hidden
M_numeric_as_cholmod_dense(CHM_DN A, double *data, int nrow, int ncol)
{
	static CHM_DN(*fun)(CHM_DN, double *, int, int) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(CHM_DN, double *, int, int))
			R_GetCCallable("Matrix", "numeric_as_cholmod_dense");
	return fun(A, data, nrow, ncol);
}

R_MATRIX_INLINE SEXP attribute_hidden
M_cholmod_factor_as_sexp(CHM_FR L, int doFree)
{
	static SEXP(*fun)(CHM_FR, int) = NULL;
	if (fun == NULL)
		fun = (SEXP(*)(CHM_FR, int))
			R_GetCCallable("Matrix", "cholmod_factor_as_sexp");
	return fun(L, doFree);
}

R_MATRIX_INLINE SEXP attribute_hidden
M_cholmod_sparse_as_sexp(CHM_SP A, int doFree,
                         int ttype, int doLogic, const char *diagString,
                         SEXP dimnames)
{
	static SEXP(*fun)(CHM_SP, int, int, int, const char *, SEXP) = NULL;
	if (fun == NULL)
		fun = (SEXP(*)(CHM_SP, int, int, int, const char *, SEXP))
			R_GetCCallable("Matrix", "cholmod_sparse_as_sexp");
	return fun(A, doFree, ttype, doLogic, diagString, dimnames);
}

R_MATRIX_INLINE SEXP attribute_hidden
M_cholmod_triplet_as_sexp(CHM_TR A, int doFree,
                          int ttype, int doLogic, const char *diagString,
                          SEXP dimnames)
{
	static SEXP(*fun)(CHM_TR, int, int, int, const char *, SEXP) = NULL;
	if (fun == NULL)
		fun = (SEXP(*)(CHM_TR, int, int, int, const char *, SEXP))
			R_GetCCallable("Matrix", "cholmod_triplet_as_sexp");
	return fun(A, doFree, ttype, doLogic, diagString, dimnames);
}

R_MATRIX_INLINE SEXP attribute_hidden
M_cholmod_dense_as_sexp(CHM_DN A, int doFree)
{
	static SEXP(*fun)(CHM_DN, int) = NULL;
	if (fun == NULL)
		fun = (SEXP(*)(CHM_DN, int))
			R_GetCCallable("Matrix", "cholmod_dense_as_sexp");
	return fun(A, doFree);
}

R_MATRIX_INLINE double attribute_hidden
M_cholmod_factor_ldetA(CHM_FR L)
{
	static double(*fun)(CHM_FR) = NULL;
	if (fun == NULL)
		fun = (double(*)(CHM_FR))
			R_GetCCallable("Matrix", "cholmod_factor_ldetA");
	return fun(L);
}

R_MATRIX_INLINE CHM_FR attribute_hidden
M_cholmod_factor_update(CHM_FR L, CHM_SP A, double beta)
{
	static CHM_FR(*fun)(CHM_FR, CHM_SP, double) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_FR, CHM_SP, double))
			R_GetCCallable("Matrix", "cholmod_factor_update");
	return fun(L, A, beta);
}

#ifdef __cplusplus
}
#endif

#endif /* !defined(R_MATRIX_NO_CHOLMOD_UTILS) */
