#include "taucs_utils.h"

/** 
 * Create a pointer to a taucs_ccs_matrix from an R object that
 * inherits from class cscMatrix according to the flags.
 * 
 * @param A Pointer to an object that inherits from cscMatrix
 * @param flags taucs flags describing the matrix
 * 
 * @return A taucs_ccs_matrix pointer to the existing storage (no copying).
 */
taucs_ccs_matrix* csc_taucs_ptr(SEXP A, int flags)
{
    taucs_ccs_matrix *ans =
	(taucs_ccs_matrix *) R_alloc(1, sizeof(taucs_ccs_matrix));
    int *dims = INTEGER(GET_SLOT(A, Matrix_DimSym));

    ans->flags = flags;
    ans->m = dims[0];
    ans->n = dims[1];
    ans->colptr = INTEGER(GET_SLOT(A, Matrix_pSym));
    ans->rowind = INTEGER(GET_SLOT(A, Matrix_iSym));
    if (flags & TAUCS_DOUBLE)
	ans->values.d = REAL(GET_SLOT(A, Matrix_xSym));
    if (flags & TAUCS_DCOMPLEX)
	ans->values.z =
	    (taucs_dcomplex *) COMPLEX(GET_SLOT(A, Matrix_zSym));
    return ans;
}

/** 
 * Copy a taucs_ccs_matrix to an R object of the appropriate class and
 * free the storage used by the taucs_ccs_matrix.
 * 
 * @param tm A pointer to a taucs_ccs_matrix
 * 
 * @return An R object of class "cscMatrix" or "sscMatrix" or "tscMatrix"
 */

SEXP mat_from_taucs(taucs_ccs_matrix *tm)
{
    SEXP ans;
    char *cls;
    int nnz = tm->colptr[tm->n];

    cls = "cscMatrix";
    if (tm->flags & TAUCS_SYMMETRIC) cls = "sscMatrix";
    if (tm->flags & TAUCS_TRIANGULAR) cls = "tscMatrix";
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cls)));
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, tm->n + 1));
    Memcpy(INTEGER(GET_SLOT(ans, Matrix_pSym)), tm->colptr, tm->n + 1);
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    Memcpy(INTEGER(GET_SLOT(ans, Matrix_iSym)), tm->rowind, nnz);
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    Memcpy(REAL(GET_SLOT(ans, Matrix_xSym)), tm->values.d, nnz);
    cscMatrix_set_Dim(ans, tm->m);
    taucs_dccs_free(tm);
    UNPROTECT(1);
    return ans;
}

taucs_ccs_matrix* copy_csc_to_taucs(SEXP A, int typ)
{
    SEXP pslot = GET_SLOT(A, Matrix_pSym),
	islot = GET_SLOT(A, Matrix_iSym);
    int *dims = INTEGER(GET_SLOT(A, Matrix_DimSym));
    taucs_ccs_matrix *ans =
	taucs_ccs_create(dims[0], dims[1], length(islot), typ);

    Memcpy(ans->colptr, INTEGER(pslot), length(pslot));
    Memcpy(ans->rowind, INTEGER(islot), length(islot));
    if (typ & TAUCS_DOUBLE)
	Memcpy(ans->values.d, REAL(GET_SLOT(A, Matrix_xSym)),
	       length(islot));
    if (typ & TAUCS_DCOMPLEX)
	Memcpy(ans->values.d,
	       (taucs_dcomplex *) COMPLEX(GET_SLOT(A, Matrix_zSym)),
	       length(islot));
    return ans;
}
    


/* Utilities for the TAUCS library */
				/* timers */
double taucs_wtime() { return 0.0; }
double taucs_ctime() { return 0.0; }
				/* memory allocation */
#undef malloc
#undef calloc
#undef realloc
#undef free

void* taucs_malloc_stub (size_t size)               { return malloc(size); }
void* taucs_calloc_stub (size_t nmemb, size_t size) { return calloc(nmemb,size); }
void* taucs_realloc_stub(void* ptr, size_t size)    { return realloc(ptr,size); }
void  taucs_free_stub   (void* ptr)                 { free(ptr); }

double taucs_allocation_amount()   { return 0.0; }
int    taucs_allocation_count()    { return 0; }
int    taucs_allocation_attempts() { return 0; }
void   taucs_allocation_assert_clean() {}
void   taucs_allocation_mark_clean() {}
void   taucs_allocation_induce_failure(int i) {}
				/* logging */
int
taucs_printf(char *fmt, ...)
{
    return 0;
}
				/* arithmetic constants */
double taucs_get_nan() { return R_NaN; }
double taucs_dzero_const     =  0.0;
double taucs_done_const      =  1.0;
double taucs_dminusone_const = -1.0;
