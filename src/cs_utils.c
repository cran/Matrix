#include "cs_utils.h"

/* Borrowed from one of Tim Davis' examples in the CSparse Demo directory */
/* 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
static int is_sym (cs *A)
{
    int is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i ;
    if (m != n) return (0) ;
    is_upper = 1 ;
    is_lower = 1 ;
    for (j = 0 ; j < n ; j++)
    {
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	    if (Ai [p] > j) is_upper = 0 ;
	    if (Ai [p] < j) is_lower = 0 ;
	}
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
}

static R_INLINE int
check_class(char *class, char **valid)
{
    int ans;
    for (ans = 0; ; ans++) {
	if (!strlen(valid[ans])) return -1;
	if (!strcmp(class, valid[ans])) return ans;
    }
}

/**
 * Create a cs object with the contents of x.  Note that
 * the result should *not* be freed with cs_spfree.  Use
 * Free on the result.
 *
 * @param x pointer to an object that inherits from CsparseMatrix
 *
 * @return pointer to a cs object that contains pointers
 * to the slots of x.
 */
cs *Matrix_as_cs(SEXP x)
{
    cs *ans = Calloc(1, cs);
    char *valid[] = {"dgCMatrix", "dsCMatrix", "dtCMatrix", ""};
    int *dims,
	ctype = check_class(CHAR(asChar(getAttrib(x, R_ClassSymbol))), valid);
    SEXP islot;

    if (ctype < 0) error("invalid class of object to Matrix_as_cs");
				/* dimensions and nzmax */
    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    ans->m = dims[0]; ans->n = dims[1];
    islot = GET_SLOT(x, Matrix_iSym);
    ans->nz = -1;		/* indicates compressed column storage */
    ans->nzmax = LENGTH(islot);
    ans->i = INTEGER(islot);
    ans->p = INTEGER(GET_SLOT(x, Matrix_pSym));
    ans->x = REAL(GET_SLOT(x, Matrix_xSym));

    return ans;
}

/**
 * Copy the contents of a to an appropriate CsparseMatrix object and,
 * optionally, free a or free both a and the pointers to its contents.
 *
 * @param a matrix to be converted
 * @param cl the name of the S4 class of the object to be generated
 * @param dofree 0 - don't free a; > 0 cs_free a; < 0 Free a
 *
 * @return SEXP containing a copy of a
 */
SEXP Matrix_cs_to_SEXP(cs *a, const char *cl, int dofree)
{
    SEXP ans;
    int *dims, sym = !strcmp(cl, "dsCMatrix"), tr = !strcmp(cl, "dtCMatrix");
    
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cl)));
				/* allocate and copy common slots */
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = a->m; dims[1] = a->n;
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, a->n + 1)),
	   a->p, a->n + 1);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, a->nz)),
	   a->i, a->nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, a->nz)),
	   a->x, a->nz);
    if (sym || tr) {
	int uplo = is_sym(a);
	if (!uplo) error("cs matrix not compatible with class");
	SET_SLOT(ans, Matrix_diagSym, mkString("N"));
	SET_SLOT(ans, Matrix_uploSym, mkString(uplo < 0 ? "L" : "U"));
    }
    if (dofree > 0) cs_free((void *) a);
    if (dofree < 0) Free(a);
    UNPROTECT(1);
    return ans;
}

