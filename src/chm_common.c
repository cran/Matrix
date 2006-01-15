#include "chm_common.h"
#include "Mutils.h"

static R_INLINE int
check_class(char *class, char **valid)
{
    int ans;
    for (ans = 0; ; ans++) {
	if (!strlen(valid[ans])) return -1;
	if (!strcmp(class, valid[ans])) return ans;
    }
}

cholmod_sparse *as_cholmod_sparse(SEXP x)
{
    cholmod_sparse *ans = Calloc(1, cholmod_sparse);
    char *valid[] = {"dgCMatrix", "dsCMatrix", "dtCMatrix",
		     "lgCMatrix", "lsCMatrix", "ltCMatrix",
		     "zgCMatrix", "zsCMatrix", "ztCMatrix",
		     ""};
    int *dims, ctype = check_class(CHAR(asChar(getAttrib(x, R_ClassSymbol))),
				   valid);
    SEXP islot;

    if (ctype < 0) error("invalid class of object to as_cholmod_sparse");
				/* characteristics of the system */
    ans->itype = CHOLMOD_INT;
    ans->dtype = CHOLMOD_DOUBLE;
    ans->packed = TRUE;
    ans->sorted = TRUE;
    ans->x = ans->z = ans->nz = (void *) NULL;
				/* dimensions and nzmax */
    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    ans->nrow = dims[0];
    ans->ncol = dims[1];
    islot = GET_SLOT(x, Matrix_iSym);
    ans->nzmax = LENGTH(islot);
				/* slots always present */
    ans->i = (void *) INTEGER(islot);
    ans->p = (void *) INTEGER(GET_SLOT(x, Matrix_pSym));
				/* set the xtype and any elements */
    switch(ctype) {
    case 0:
    case 1:
    case 2:
	ans->xtype = CHOLMOD_REAL;
	ans->x = (void *) REAL(GET_SLOT(x, Matrix_xSym));
	break;
    case 3:
    case 4:
    case 5:
	ans->xtype = CHOLMOD_PATTERN;
	break;
    case 6:
    case 7:
    case 8:
	ans->xtype = CHOLMOD_COMPLEX;
	ans->x = (void *) COMPLEX(GET_SLOT(x, Matrix_xSym));
	break;
    }
				/* set the stype */
    switch(ctype % 3) {
    case 0: ans->stype = 0; break;
    case 1:
	ans->stype =
	    (!strcmp(CHAR(asChar(getAttrib(x, Matrix_uploSym))), "U")) ?
	    1 : -1;
	break;
    case 2: error("triangular matrices not yet mapped to CHOLMOD");
    }

    return ans;
}

/**
 * Copy the contents of a to an appropriate TsparseMatrix object and,
 * optionally, free a or free both a and its the pointers to its contents.
 *
 * @param a matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_sparse_to_SEXP(cholmod_sparse *a, int dofree)
{
    SEXP ans;
    char *cl = "";		/* -Wall */
    int *dims, nnz = cholmod_nnz(a, &c);
				/* ensure a is sorted and packed */
    if (!a->sorted || !a->packed) cholmod_sort(a, &c);
				/* determine the class of the result */
    switch(a->xtype){
    case CHOLMOD_PATTERN:
	cl = (a->stype) ? "lsCMatrix" : "lgCMatrix";
	break;
    case CHOLMOD_REAL:
	cl = (a->stype) ? "dsCMatrix" : "dgCMatrix";
	break;
    case CHOLMOD_COMPLEX:
	cl = (a->stype) ? "zsCMatrix" : "zgCMatrix";
	break;
    default: error("unknown xtype in cholmod_sparse object");
    }
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cl)));
				/* allocate and copy common slots */
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = a->nrow; dims[1] = a->ncol;
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, a->ncol + 1)),
	   (int *) a->p, a->ncol + 1);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nnz)),
	   (int *) a->i, nnz);
				/* copy data slot if present */
    if (a->xtype == CHOLMOD_REAL)
	Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nnz)),
	       (double *) a->x, nnz);
    if (a->xtype == CHOLMOD_COMPLEX)
	error("complex sparse matrix code not yet written");
/* 	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, nnz)), */
/* 	       (complex *) a->x, nnz); */
				/* set symmetry attributes */
    if (a->stype)
	SET_SLOT(ans, Matrix_uploSym,
		 mkString((a->stype > 0) ? "U" : "L"));
    if (dofree > 0) cholmod_free_sparse(&a, &c);
    if (dofree < 0) Free(a);
    UNPROTECT(1);
    return ans;
}

/**
 * Create a cholmod_triplet object with the contents of x.  Note that
 * the result should *not* be freed with cholmod_triplet_free.  Use
 * free or Free on the result.
 *
 * @param x pointer to an object that inherits from TsparseMatrix
 *
 * @return pointer to a cholmod_triplet object that contains pointers
 * to the slots of x.
 */
cholmod_triplet *as_cholmod_triplet(SEXP x)
{
    cholmod_triplet *ans = Calloc(1, cholmod_triplet);
    char *valid[] = {"dgTMatrix", "dsTMatrix", "dtTMatrix",
		     "lgTMatrix", "lsTMatrix", "ltTMatrix",
		     "zgTMatrix", "zsTMatrix", "ztTMatrix",
		     ""};
    int *dims, ctype = check_class(CHAR(asChar(getAttrib(x, R_ClassSymbol))),
				   valid);
    SEXP islot;

    if (ctype < 0) error("invalid class of object to as_cholmod_triplet");
				/* characteristics of the system */
    ans->itype = CHOLMOD_INT;
    ans->dtype = CHOLMOD_DOUBLE;
    ans->x = ans->z = (void *) NULL;
				/* dimensions and nzmax */
    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    ans->nrow = dims[0];
    ans->ncol = dims[1];
    islot = GET_SLOT(x, Matrix_iSym);
    ans->nnz = ans->nzmax = LENGTH(islot);
				/* slots always present */
    ans->i = (void *) INTEGER(islot);
    ans->j = (void *) INTEGER(GET_SLOT(x, Matrix_jSym));
				/* set the xtype and any elements */
    switch(ctype) {
    case 0:
    case 1:
    case 2:
	ans->xtype = CHOLMOD_REAL;
	ans->x = (void *) REAL(GET_SLOT(x, Matrix_xSym));
	break;
    case 3:
    case 4:
    case 5:
	ans->xtype = CHOLMOD_PATTERN;
	break;
    case 6:
    case 7:
    case 8:
	ans->xtype = CHOLMOD_COMPLEX;
	ans->x = (void *) COMPLEX(GET_SLOT(x, Matrix_xSym));
	break;
    }
				/* set the stype */
    switch(ctype % 3) {
    case 0: ans->stype = 0; break;
    case 1:
	ans->stype =
	    (!strcmp(CHAR(asChar(getAttrib(x, Matrix_uploSym))), "U")) ?
	    1 : -1;
	break;
    case 2: error("triangular matrices not yet mapped to CHOLMOD");
    }

    return ans;
}

/**
 * Copy the contents of a to an appropriate TsparseMatrix object and,
 * optionally, free a or free both a and its the pointers to its contents.
 *
 * @param a matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_triplet_to_SEXP(cholmod_triplet *a, int dofree)
{
    SEXP ans;
    char *cl = "";		/* -Wall */
    int *dims;
				/* determine the class of the result */
    switch(a->xtype){
    case CHOLMOD_PATTERN:
	cl = (a->stype) ? "lsTMatrix" : "lgTMatrix";
	break;
    case CHOLMOD_REAL:
	cl = (a->stype) ? "dsTMatrix" : "dgTMatrix";
	break;
    case CHOLMOD_COMPLEX:
	cl = (a->stype) ? "zsTMatrix" : "zgTMatrix";
	break;
    default: error("unknown xtype in cholmod_triplet object");
    }
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cl)));
				/* allocate and copy common slots */
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = a->nrow; dims[1] = a->ncol;
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, a->nnz)),
	   (int *) a->i, a->nnz);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_jSym, INTSXP, a->nnz)),
	   (int *) a->j, a->nnz);
				/* copy data slot if present */
    if (a->xtype == CHOLMOD_REAL)
	Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, a->nnz)),
	       (double *) a->x, a->nnz);
    if (a->xtype == CHOLMOD_COMPLEX)
	error("complex sparse matrix code not yet written");
/* 	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, a->nnz)), */
/* 	       (complex *) a->x, a->nz); */
				/* set symmetry attributes */
    if (a->stype)
	SET_SLOT(ans, Matrix_uploSym,
		 mkString((a->stype > 0) ? "U" : "L"));
    if (dofree > 0) cholmod_free_triplet(&a, &c);
    if (dofree < 0) Free(a);
    UNPROTECT(1);
    return ans;
}

/**
 * Create a cholmod_dense object with the contents of x.  Note that
 * the result should *not* be freed with cholmod_dense_free.  Use
 * free or Free on the result.
 *
 * @param x pointer to an object that inherits from (denseMatrix ^ generalMatrix)
 *
 * @return pointer to a cholmod_dense object that contains a pointer
 * to the contents of x.
 */
cholmod_dense *as_cholmod_dense(SEXP x)
{
    cholmod_dense *ans = Calloc(1, cholmod_dense);
    char *valid[] = {"dmatrix", "dgeMatrix",
		     "lmatrix", "lgeMatrix",
		     "zmatrix", "zgeMatrix", ""},
	*cl = CHAR(asChar(getAttrib(x, R_ClassSymbol)));
    int dims[2], ctype = check_class(cl, valid);

    if (ctype < 0) {
	ctype = (isReal(x) ? 0 :
		 (isLogical(x) ? 2 :
		  (isComplex(x) ? 4 : -1)));
	if (isMatrix(x)) Memcpy(dims, INTEGER(getAttrib(x, R_DimSymbol)), 2);
	else {dims[0] = LENGTH(x); dims[1] = 1;}
    } else Memcpy(dims, INTEGER(GET_SLOT(x, Matrix_DimSym)), 2);
    if (ctype < 0) error("invalid class of object to as_cholmod_dense");
				/* characteristics of the system */
    ans->dtype = CHOLMOD_DOUBLE;
    ans->x = ans->z = (void *) NULL;
				/* dimensions and nzmax */
    ans->d = ans->nrow = dims[0];
    ans->ncol = dims[1];
    ans->nzmax = dims[0] * dims[1];
				/* set the xtype and any elements */
    switch(ctype / 2) {
    case 0:
	ans->xtype = CHOLMOD_REAL;
	ans->x = (void *) REAL((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x);
	break;
    case 1:
	ans->xtype = CHOLMOD_PATTERN;
	ans->x = (void *) LOGICAL((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x);
	break;
    case 3:
	ans->xtype = CHOLMOD_COMPLEX;
	ans->x = (void *) COMPLEX((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x);
	break;
    }

    return ans;
}

/**
 * Copy the contents of a to an appropriate denseMatrix object and,
 * optionally, free a or free both a and its pointer to its contents.
 *
 * @param a matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_dense_to_SEXP(cholmod_dense *a, int dofree)
{
    SEXP ans;
    char *cl;
    int *dims, ntot;
				/* determine the class of the result */
    cl = (a->xtype == CHOLMOD_PATTERN) ? "lgeMatrix" :
	((a->xtype == CHOLMOD_REAL) ? "dgeMatrix" :
	 ((a->xtype == CHOLMOD_COMPLEX) ? "zgeMatrix" : ""));
    if (!strlen(cl)) error("unknown xtype");

    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cl)));
				/* allocate and copy common slots */
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = a->nrow; dims[1] = a->ncol;
    ntot = dims[0] * dims[1];
    if (a->d == a->nrow) {	/* copy data slot if present */
	if (a->xtype == CHOLMOD_REAL)
	    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, ntot)),
		   (double *) a->x, ntot);
	if (a->xtype == CHOLMOD_COMPLEX)
	    error("complex sparse matrix code not yet written");
/* 	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, a->nnz)), */
/* 	       (complex *) a->x, a->nz); */
    } else error("code for cholmod_dense with holes not yet written");

    if (dofree > 0) cholmod_free_dense(&a, &c);
    if (dofree < 0) Free(a);
    UNPROTECT(1);
    return ans;
}

cholmod_dense *numeric_as_chm_dense(double *v, int n)
{
    cholmod_dense *ans = Calloc(1, cholmod_dense);

    ans->d = ans->nzmax = ans->nrow = n;
    ans->ncol = 1;
    ans->x = (void *) v;
    ans->xtype = CHOLMOD_REAL;
    ans->dtype = CHOLMOD_DOUBLE;
    return ans;
}

/**
 * Create a cholmod_factor object from the contents of x.  Note that
 * the result should *not* be freed with cholmod_free_factor.  Use
 * Free on the result.
 *
 * @param x pointer to an object that inherits from ddenseMatrix
 *
 * @return pointer to a cholmod_dense object that contains a pointer
 * to the contents of x.
 */
cholmod_factor *as_cholmod_factor(SEXP x)
{
    cholmod_factor *ans = Calloc(1, cholmod_factor);
    char *valid[] = {"dCHMsuper", "dCHMsimpl", "lCHMsuper", "lCHMsimpl", ""};
    int *type = INTEGER(GET_SLOT(x, install("type"))),
	ctype = check_class(CHAR(asChar(getAttrib(x, R_ClassSymbol))), valid);
    SEXP tmp;

    if (ctype < 0) error("invalid class of object to as_cholmod_factor");
				/* characteristics of the system */
    ans->itype = CHOLMOD_INT;
    ans->dtype = CHOLMOD_DOUBLE;
    ans->z = (void *) NULL;
    ans->xtype = (ctype < 2) ? CHOLMOD_REAL : CHOLMOD_PATTERN;

				/* unravel the type */
    ans->ordering = type[0];
    ans->is_ll = (type[1] ? 1 : 0);
    ans->is_super = (type[2] ? 1 : 0);
    ans->is_monotonic = (type[3] ? 1 : 0);
				/* check for consistency */
    if ((!(ans->is_ll)) && ans->is_super)
	error(_("Supernodal LDL' decomposition not available"));
    if ((!type[2]) ^ (ctype % 2))
	error(_("Supernodal/simplicial class inconsistent with type flags"));
				/* slots always present */
    tmp = GET_SLOT(x, Matrix_permSym);
    ans->minor = ans->n = LENGTH(tmp); ans->Perm = INTEGER(tmp);
    if (ctype < 2) {
	tmp = GET_SLOT(x, Matrix_xSym);
	ans->x = REAL(tmp);
    } else ans->x = (void*)NULL;
    if (ans->is_super) {	/* supernodal factorization */
	ans->xsize = LENGTH(tmp);
	ans->maxcsize = type[4]; ans->maxesize = type[5];
	ans->i = (int*)NULL;
	tmp = GET_SLOT(x, install("super"));
	ans->nsuper = LENGTH(tmp) - 1; ans->super = INTEGER(tmp);
	/* Move these checks to the dCHMfactor_validate function */
	if (ans->nsuper < 1)
	    error(_("Number of supernodes must be positive when is_super is TRUE"));
	tmp = GET_SLOT(x, install("pi"));
	if (LENGTH(tmp) != ans->nsuper + 1)
	    error(_("Lengths of super and pi must be equal"));
	ans->pi = INTEGER(tmp);
	tmp = GET_SLOT(x, install("px"));
	if (LENGTH(tmp) != ans->nsuper + 1)
	    error(_("Lengths of super and px must be equal"));
	ans->px = INTEGER(tmp);
	tmp = GET_SLOT(x, install("s"));
	ans->ssize = LENGTH(tmp); ans->s = INTEGER(tmp);
    } else {
	ans->nzmax = LENGTH(tmp);
	ans->p = INTEGER(GET_SLOT(x, Matrix_pSym));
	ans->i = INTEGER(GET_SLOT(x, Matrix_iSym));
	ans->ColCount = INTEGER(GET_SLOT(x, install("colcount")));
	ans->nz = INTEGER(GET_SLOT(x, install("nz")));
	ans->next = INTEGER(GET_SLOT(x, install("nxt")));
	ans->prev = INTEGER(GET_SLOT(x, install("prv")));
    }
    return ans;
}

/**
 * Copy the contents of f to an appropriate dCHMfactor object and,
 * optionally, free f or free both f and its pointer to its contents.
 *
 * @param a matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_factor_to_SEXP(cholmod_factor *f, int dofree)
{
    SEXP ans;
    int *type;
    char *class = (char*) NULL;	/* -Wall */

    switch(f->xtype) {
    case CHOLMOD_REAL:
	class = f->is_super ? "dCHMsuper" : "dCHMsimpl";
	break;
    case CHOLMOD_PATTERN:
	class = f->is_super ? "dCHMsuper" : "dCHMsimpl";
	break;
    default:
	error(_("f->xtype of %d not recognized"), f->xtype);
    }
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(class)));
    if (f->minor < f->n)
	error(_("CHOLMOD factorization was unsuccessful"));
				/* copy component of known length */
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_permSym, INTSXP, f->n)),
	   (int*)f->Perm, f->n);
    type = INTEGER(ALLOC_SLOT(ans, install("type"), INTSXP, f->is_super ? 6 : 4));
    type[0] = f->ordering; type[1] = f->is_ll;
    type[2] = f->is_super; type[3] = f->is_monotonic;
    if (f->is_super) {
	type[4] = f->maxcsize; type[5] = f->maxesize;
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("super"), INTSXP, f->nsuper + 1)),
	       (int*)f->super, f->nsuper+1);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("pi"), INTSXP, f->nsuper + 1)),
	       (int*)f->pi, f->nsuper + 1);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("px"), INTSXP, f->nsuper + 1)),
	       (int*)f->px, f->nsuper + 1);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("s"), INTSXP, f->ssize)),
	       (int*)f->s, f->ssize);
	Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, f->xsize)),
	       (double*)f->x, f->xsize);
    } else {
	Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, f->nzmax)),
	   (int*)f->i, f->nzmax);
	Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, f->n + 1)),
	   (int*)f->p, f->n + 1);
	Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, f->nzmax)),
	       (double*)f->x, f->nzmax);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("colcount"), INTSXP, f->n)),
	       (int*)f->ColCount, f->n);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("nxt"), INTSXP, f->n + 2)),
	       (int*)f->next, f->n + 2);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("prv"), INTSXP, f->n + 2)),
	       (int*)f->prev, f->n + 2);
    }
    if (dofree > 0) cholmod_free_factor(&f, &c);
    if (dofree < 0) Free(f);
    UNPROTECT(1);
    return ans;
}

SEXP CHMfactor_validate(SEXP obj) /* placeholder */
{
    return ScalarLogical(1);
}

SEXP CHMsimpl_validate(SEXP obj) /* placeholder */
{
    return ScalarLogical(1);
}

SEXP CHMsuper_validate(SEXP obj) /* placeholder */
{
    return ScalarLogical(1);
}

