#include <Rconfig.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "Matrix.h"

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

cholmod_dense attribute_hidden
*M_as_cholmod_dense(SEXP x)
{
    static cholmod_dense*(*fun)(SEXP) = NULL;
    if(fun == NULL)
	fun = (cholmod_dense*(*)(SEXP))
	    R_GetCCallable("Matrix", "as_cholmod_dense");
    return fun(x);
}

cholmod_factor attribute_hidden
*M_as_cholmod_factor(SEXP x)
{
    static cholmod_factor*(*fun)(SEXP) = NULL;
    if(fun == NULL)
	fun = (cholmod_factor*(*)(SEXP))
	    R_GetCCallable("Matrix", "as_cholmod_factor");
    return fun(x);
}

cholmod_sparse attribute_hidden
*M_as_cholmod_sparse(SEXP x)
{
    static cholmod_sparse*(*fun)(SEXP)= NULL;
    if(fun == NULL)
	fun = (cholmod_sparse*(*)(SEXP))
	    R_GetCCallable("Matrix", "as_cholmod_sparse");
    return fun(x);
}

SEXP attribute_hidden
M_chm_factor_to_SEXP(cholmod_factor *f, int dofree)
{
    static SEXP(*fun)(cholmod_factor*,int) = NULL;
    if(fun == NULL)
	fun = (SEXP(*)(cholmod_factor*,int))
	    R_GetCCallable("Matrix", "chm_factor_to_SEXP");
    return fun(f, dofree);
}

SEXP attribute_hidden
M_chm_sparse_to_SEXP(cholmod_sparse *a, int dofree,
		     int uploT, int Rkind, char *diag, SEXP dn)
{
    static SEXP(*fun)(cholmod_sparse*, int,
		     int, int, char*, SEXP) = NULL;
    if(fun == NULL)
	fun = (SEXP(*)(cholmod_sparse*, int, int, int, char*, SEXP))
	    R_GetCCallable("Matrix", "chm_sparse_to_SEXP");
    return fun(a, dofree, uploT, Rkind, diag, dn);
}

cholmod_sparse attribute_hidden
*M_cholmod_aat(cholmod_sparse *A, int *fset, size_t fsize,
	       int mode, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(cholmod_sparse*,int*,size_t,
				 int,cholmod_common*) = NULL;
    if(fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_sparse*,int*,size_t,
				  int,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_aat");
    return fun(A, fset, fsize, mode, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_add(cholmod_sparse *A, cholmod_sparse *B,
	       double alpha[2], double beta[2], int values,
	       int sorted, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(cholmod_sparse*,cholmod_sparse*,
				 double*,double*,int,int,
				 cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_sparse*,cholmod_sparse*,
				  double*,double*,int,int,
				  cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_add");
    return fun(A, B, alpha, beta, values, sorted, Common);
}

cholmod_dense attribute_hidden
*M_cholmod_allocate_dense(size_t nrow, size_t ncol, size_t d,
			  int xtype, cholmod_common *Common)
{
    static cholmod_dense*(*fun)(size_t,size_t,size_t,
				int,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_dense*(*)(size_t,size_t,size_t,
				 int,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_allocate_dense");
    return fun(nrow, ncol, d, xtype, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_allocate_sparse(size_t nrow, size_t ncol, size_t nzmax,
			  int sorted, int packed, int stype,
			  int xtype, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(size_t,size_t,size_t,int,int,
				 int,int,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)
	       (size_t,size_t,size_t,int,int,int,int,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_allocate_sparse");
    return fun(nrow,ncol,nzmax,sorted,packed,stype,xtype,Common);
}

cholmod_triplet attribute_hidden
*M_cholmod_allocate_triplet(size_t nrow, size_t ncol, size_t nzmax,
			    int stype, int xtype, cholmod_common *Common)
{
    static cholmod_triplet*(*fun)(size_t,size_t,size_t,
				 int,int,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_triplet*(*)
	       (size_t,size_t,size_t,int,int,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_allocate_triplet");
    return fun(nrow,ncol,nzmax,stype,xtype,Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_triplet_to_sparse(cholmod_triplet *T, int nzmax,
			     cholmod_common *Common)
{
    static cholmod_sparse*(*fun)
	(cholmod_triplet*,int,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_triplet*,int,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_triplet_to_sparse");
    return fun(T, nzmax, Common);
}

cholmod_triplet attribute_hidden
*M_cholmod_sparse_to_triplet(cholmod_sparse *A, cholmod_common *Common)
{
    static cholmod_triplet*(*fun)
	(cholmod_sparse*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_triplet*(*)(cholmod_sparse*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_sparse_to_triplet");
    return fun(A, Common);
}

cholmod_dense attribute_hidden
*M_cholmod_sparse_to_dense(cholmod_sparse *A, cholmod_common *Common)
{
    static cholmod_dense*(*fun)
	(cholmod_sparse*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_dense*(*)(cholmod_sparse*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_sparse_to_dense");
    return fun(A, Common);
}

cholmod_factor attribute_hidden
*M_cholmod_analyze(cholmod_sparse *A, cholmod_common *Common)
{
    static cholmod_factor*(*fun)(cholmod_sparse*,cholmod_common*) = NULL;
    if (fun == NULL)
    fun = (cholmod_factor*(*)(cholmod_sparse*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_analyze");
    return fun(A, Common);
}

cholmod_factor attribute_hidden
*M_cholmod_analyze_p(cholmod_sparse *A, int *Perm, int *fset,
		     size_t fsize, cholmod_common *Common)
{
    static cholmod_factor*(*fun)(cholmod_sparse*,int*,int*,size_t,
				 cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_factor*(*)(cholmod_sparse*,int*,int*,
				  size_t,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_analyze_p");
    return fun(A, Perm, fset, fsize, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_copy(cholmod_sparse *A, int stype,
			       int mode, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)
	(cholmod_sparse*,int,int,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_sparse*,int,int,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_copy");
    return fun(A, stype, mode, Common);
}

attribute_hidden
cholmod_dense* M_cholmod_copy_dense(cholmod_dense *A, cholmod_common *Common)
{
    static cholmod_dense*(*fun)(cholmod_dense*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_dense*(*)(cholmod_dense*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_copy_dense");
    return fun(A, Common);
}

cholmod_factor attribute_hidden
*M_cholmod_copy_factor(cholmod_factor *L, cholmod_common *Common)
{
    static cholmod_factor*(*fun)(cholmod_factor*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_factor*(*)(cholmod_factor*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_copy_factor");
    return fun(L, Common);
}

int attribute_hidden
M_cholmod_change_factor(int to_xtype, int to_ll, int to_super, int to_packed,
			int to_monotonic, cholmod_factor *L, cholmod_common *Common)
{
    static int(*fun)(int,int,int,int,int,cholmod_factor*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(int,int,int,int,int,cholmod_factor*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_change_factor");
    return fun(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_copy_sparse(cholmod_sparse *A, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(cholmod_sparse*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_sparse*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_copy_sparse");
    return fun(A, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_factor_to_sparse(cholmod_factor *L, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(cholmod_factor*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_factor*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_factor_to_sparse");
    return fun(L, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_submatrix(cholmod_sparse *A, int *rset, int rsize, int *cset,
		     int csize, int values, int sorted, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(cholmod_sparse*,int*,int,int*,int,
				 int,int,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_sparse*,int*,int,int*,
				  int,int,int,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_submatrix");
    return fun(A, rset, rsize, cset, csize, values, sorted, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_dense_to_sparse(cholmod_dense *X, int values, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(cholmod_dense*,int,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_dense*,int,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_dense_to_sparse");
    return fun(X, values, Common);
}

int attribute_hidden
M_cholmod_factorize(cholmod_sparse *A, cholmod_factor *L,
		    cholmod_common *Common)
{
    static int(*fun)(cholmod_sparse*,cholmod_factor*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_sparse*,cholmod_factor*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_factorize");
    return fun(A, L, Common);
}

int attribute_hidden
M_cholmod_finish(cholmod_common *Common)
{

    static int(*fun)(cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_finish");
    return fun(Common);
}

int attribute_hidden
M_cholmod_sort(cholmod_sparse *A, cholmod_common *Common)
{
    static int(*fun)(cholmod_sparse*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_sparse*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_sort");
    return fun(A, Common);
}

int attribute_hidden
M_cholmod_free_dense(cholmod_dense **A, cholmod_common *Common)
{
    static int(*fun)(cholmod_dense**,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_dense**,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_free_dense");
    return fun(A, Common);
}

int attribute_hidden
M_cholmod_free_factor(cholmod_factor **L, cholmod_common *Common)
{
    static int(*fun)(cholmod_factor**,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_factor**,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_free_factor");
    return fun(L, Common);
}

int attribute_hidden
M_cholmod_free_sparse(cholmod_sparse **A, cholmod_common *Common)
{
    static int(*fun)(cholmod_sparse**,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_sparse**,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_free_sparse");
    return fun(A, Common);
}

int attribute_hidden
M_cholmod_free_triplet(cholmod_triplet **T, cholmod_common *Common)
{
    static int(*fun)(cholmod_triplet**,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_triplet**,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_free_triplet");
    return fun(T, Common);
}

long attribute_hidden
M_cholmod_nnz(cholmod_sparse *A, cholmod_common *Common)
{
    static long(*fun)(cholmod_sparse*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (long(*)(cholmod_sparse*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_nnz");
    return fun(A, Common);
}

int attribute_hidden
M_cholmod_sdmult(cholmod_sparse *A, int transpose,
		 double alpha [2], double beta [2],
		 cholmod_dense *X, cholmod_dense *Y,
		 cholmod_common *Common)
{
    static int(*fun)(cholmod_sparse*,int,double*,double*,
		     cholmod_dense*,cholmod_dense*,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_sparse*,int,double*,double*,
		      cholmod_dense*,cholmod_dense*,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_sdmult");
    return fun(A, transpose, alpha, beta, X, Y, Common);
}

cholmod_dense attribute_hidden
*M_cholmod_solve(int sys, cholmod_factor *L,
			       cholmod_dense *B, cholmod_common *Common)
{
    static cholmod_dense*(*fun)(int,cholmod_factor*,cholmod_dense*,
				cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_dense*(*)(int,cholmod_factor*,cholmod_dense*,
				 cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_solve");
    return fun(sys, L, B, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_speye(size_t nrow, size_t ncol,
		 int xtype, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(size_t,size_t,int,cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(size_t,size_t,int,cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_speye");
    return fun(nrow, ncol, xtype, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_spsolve(int sys, cholmod_factor *L,
		   cholmod_sparse *B, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(int,cholmod_factor*,
				 cholmod_sparse*, cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(int,cholmod_factor*,
				  cholmod_sparse*, cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_spsolve");
    return fun(sys, L, B, Common);
}

void attribute_hidden
M_R_cholmod_error(int status, char *file, int line, char *message)
{
        error("Cholmod error `%s' at file:%s, line %d", message, file, line);
}
    
int attribute_hidden
M_R_cholmod_start(cholmod_common *Common)
{
    int val;
    static int(*fun)(cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_start");
    val = fun(Common);
    Common->print_function = Rprintf;
    Common->error_handler = M_R_cholmod_error;
    return val;
}

cholmod_sparse attribute_hidden
*M_cholmod_transpose(cholmod_sparse *A, int values, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(cholmod_sparse*,int,
				 cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_sparse*,int,
				  cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_transpose");
    return fun(A, values, Common);
}

cholmod_sparse attribute_hidden
*M_cholmod_vertcat(cholmod_sparse *A, cholmod_sparse *B,
		   int values, cholmod_common *Common)
{
    static cholmod_sparse*(*fun)(cholmod_sparse*, cholmod_sparse*,
				 int, cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (cholmod_sparse*(*)(cholmod_sparse*,cholmod_sparse*,
				  int, cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_vertcat");
    return fun(A, B, values, Common);
}

SEXP attribute_hidden
M_dpoMatrix_chol(SEXP x)
{
    static SEXP(*fun)(SEXP) = NULL;
    if (fun == NULL)
	fun = (SEXP(*)(SEXP))
	    R_GetCCallable("Matrix", "dpoMatrix_chol");
    return fun(x);
}

cholmod_dense attribute_hidden
*M_numeric_as_chm_dense(double *v, int n)
{
    static cholmod_dense*(*fun)(double*,int) = NULL;
    if (fun == NULL)
	fun = (cholmod_dense*(*)(double*,int))
	    R_GetCCallable("Matrix", "numeric_as_chm_dense");
    return fun(v, n);
}

int attribute_hidden
M_cholmod_scale(cholmod_dense *S, int scale, cholmod_sparse *A,
		cholmod_common *Common)
{
    static int(*fun)(cholmod_dense*,int,cholmod_sparse*,
		     cholmod_common*) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_dense*,int,cholmod_sparse*,
				  cholmod_common*))
	    R_GetCCallable("Matrix", "cholmod_scale");
    return fun(S, scale, A, Common);
}
