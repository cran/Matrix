#include "chm_common.h"
#include "Mutils.h"

cholmod_common c;

/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * In most cases this function is called through the macro AS_CHM_SP.
 * It is unusual to call it directly.
 *
 * @param ans a CHM_SP pointer
 * @param x pointer to an object that inherits from CsparseMatrix
 *
 * @return ans containing pointers to the slots of x.
 */
CHM_SP as_cholmod_sparse(CHM_SP ans, SEXP x)
{
    char *valid[] = {"dgCMatrix", "dsCMatrix", "dtCMatrix",
		     "lgCMatrix", "lsCMatrix", "ltCMatrix",
		     "ngCMatrix", "nsCMatrix", "ntCMatrix",
		     "zgCMatrix", "zsCMatrix", "ztCMatrix",
		     ""};
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	ctype = Matrix_check_class(class_P(x), valid);
    SEXP islot = GET_SLOT(x, Matrix_iSym);

    if (ctype < 0) error("invalid class of object to as_cholmod_sparse");
    memset(ans, 0, sizeof(cholmod_sparse)); /* zero the struct */

    ans->itype = CHOLMOD_INT;	/* characteristics of the system */
    ans->dtype = CHOLMOD_DOUBLE;
    ans->packed = TRUE;
    ans->sorted = TRUE;
    ans->x = ans->z = ans->nz = (void *) NULL;
				/* dimensions and nzmax */
    ans->nrow = dims[0];
    ans->ncol = dims[1];

    ans->nzmax = LENGTH(islot);
				/* slots always present */
    ans->i = (void *) INTEGER(islot);
    ans->p = (void *) INTEGER(GET_SLOT(x, Matrix_pSym));

#define AS_CHM_FINISH								\
    				/* set the xtype and any elements */		\
    switch(ctype / 3) {								\
    case 0: /* "d" */								\
	ans->xtype = CHOLMOD_REAL;						\
	ans->x = (void *) REAL(GET_SLOT(x, Matrix_xSym));			\
	break;									\
    case 1: /* "l" */								\
	ans->xtype = CHOLMOD_REAL;						\
	ans->x = (void *) REAL(coerceVector(GET_SLOT(x, Matrix_xSym), REALSXP)); \
	break;									\
    case 2: /* "n" */								\
	ans->xtype = CHOLMOD_PATTERN;						\
	break;									\
    case 3: /* "z" */								\
	ans->xtype = CHOLMOD_COMPLEX;						\
	ans->x = (void *) COMPLEX(GET_SLOT(x, Matrix_xSym));			\
	break;									\
    }										\
				/* set the stype */				\
    switch(ctype % 3) {								\
    case 0: /* g(eneral) */							\
	ans->stype = 0; break;							\
    case 1: /* s(ymmetric) */							\
	ans->stype =								\
	    (!strcmp(CHAR(asChar(getAttrib(x, Matrix_uploSym))), "U")) ?	\
	    1 : -1;								\
	break;									\
    case 2: /* t(riangular) */							\
	ans->stype = 0; break; /* Note that triangularity property is lost */	\
/*	error("triangular matrices not yet mapped to CHOLMOD"); */		\
    }										\
    return ans

    AS_CHM_FINISH;
}

/**
 * Copy the contents of a to an appropriate CsparseMatrix object and,
 * optionally, free a or free both a and its the pointers to its contents.
 *
 * @param a matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 * @param uploT 0 - not triangular; > 0 upper triangular; < 0 lower
 * @param Rkind - vector type to store for a->xtype == CHOLMOD_REAL,
 *                0 - REAL; 1 - LOGICAL
 * @param diag character string suitable for the diag slot of a
 *          triangular matrix (not accessed if uploT == 0).
 * @param dn either R_NilValue or an SEXP suitable for the Dimnames slot.
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_sparse_to_SEXP(CHM_SP a, int dofree, int uploT, int Rkind,
			const char* diag, SEXP dn)
{
    SEXP ans;
    char *cl = "";		/* -Wall */
    int *dims, nnz;

    PROTECT(dn);  /* dn is usually UNPROTECTed before the call */
				/* ensure a is sorted and packed */
    if (!a->sorted || !a->packed) cholmod_sort(a, &c);
    nnz = cholmod_nnz(a, &c);
				/* determine the class of the result */
    switch(a->xtype){
    case CHOLMOD_PATTERN:
	cl = uploT ? "ntCMatrix": ((a->stype) ? "nsCMatrix" : "ngCMatrix");
	break;
    case CHOLMOD_REAL:
	switch(Rkind) {
	case 0:
	    cl = uploT ? "dtCMatrix": ((a->stype) ? "dsCMatrix" : "dgCMatrix");
	    break;
	case 1:
	    cl = uploT ? "ltCMatrix": ((a->stype) ? "lsCMatrix" : "lgCMatrix");
	    break;
	}
	break;
    case CHOLMOD_COMPLEX:
	cl = uploT ? "ztCMatrix": ((a->stype) ? "zsCMatrix" : "zgCMatrix");
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
    if (a->xtype == CHOLMOD_REAL) {
	int i, *m_x;
	switch(Rkind) {
	case 0:
	    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nnz)),
		   (double *) a->x, nnz);
	    break;
	case 1:
	    m_x = LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, nnz));
	    for (i=0; i < nnz; i++)
		m_x[i] = (int) ((double *) a->x)[i];
	    break;
	}
    }
    else if (a->xtype == CHOLMOD_COMPLEX)
	error("complex sparse matrix code not yet written");
/* 	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, nnz)), */
/* 	       (complex *) a->x, nnz); */
    if (uploT) {		/* slots for triangularMatrix */
	if (a->stype) error("Symmetric and triangular both set");
	SET_SLOT(ans, Matrix_uploSym, mkString((uploT > 0) ? "U" : "L"));
	SET_SLOT(ans, Matrix_diagSym, mkString(diag));
    }
    if (a->stype)		/* slot for symmetricMatrix */
	SET_SLOT(ans, Matrix_uploSym,
		 mkString((a->stype > 0) ? "U" : "L"));
    if (dofree > 0) cholmod_free_sparse(&a, &c);
    if (dofree < 0) Free(a);
    if (dn != R_NilValue)
	SET_SLOT(ans, Matrix_DimNamesSym, duplicate(dn));

    UNPROTECT(2);
    return ans;
}

/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * In most cases this function is called through the macro AS_CHM_SP.
 * It is unusual to call it directly.
 *
 * @param ans a CHM_TR pointer
 * @param x pointer to an object that inherits from TsparseMatrix
 *
 * @return ans containing pointers to the slots of x.
 */
CHM_TR as_cholmod_triplet(CHM_TR ans, SEXP x)
{
    char *valid[] = {"dgTMatrix", "dsTMatrix", "dtTMatrix",
		     "lgTMatrix", "lsTMatrix", "ltTMatrix",
		     "ngTMatrix", "nsTMatrix", "ntTMatrix",
		     "zgTMatrix", "zsTMatrix", "ztTMatrix",
		     ""};
    int *dims, ctype = Matrix_check_class(class_P(x), valid);
    SEXP islot;

    if (ctype < 0) error("invalid class of object to as_cholmod_triplet");
    memset(ans, 0, sizeof(cholmod_triplet)); /* zero the struct */

    ans->itype = CHOLMOD_INT;	/* characteristics of the system */
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

    AS_CHM_FINISH;
}

/**
 * Copy the contents of a to an appropriate TsparseMatrix object and,
 * optionally, free a or free both a and its the pointers to its contents.
 *
 * @param a matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 * @param uploT 0 - not triangular; > 0 upper triangular; < 0 lower
 * @param Rkind - vector type to store for a->xtype == CHOLMOD_REAL,
 *                0 - REAL; 1 - LOGICAL
 * @param diag character string suitable for the diag slot of a
 *          triangular matrix (not accessed if uploT == 0).
 * @param dn either R_NilValue or an SEXP suitable for the Dimnames slot.
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_triplet_to_SEXP(CHM_TR a, int dofree, int uploT, int Rkind,
			 const char* diag, SEXP dn)
{
    SEXP ans;
    char *cl = "";		/* -Wall */
    int *dims;

    PROTECT(dn);  /* dn is usually UNPROTECTed before the call */
				/* determine the class of the result */
    switch(a->xtype) {
    case CHOLMOD_PATTERN:
	cl = uploT ? "ntTMatrix" :
	    ((a->stype) ? "nsTMatrix" : "ngTMatrix");
	break;
    case CHOLMOD_REAL:
	switch(Rkind) {
	case 0:
	    cl = uploT ? "dtTMatrix" :
		((a->stype) ? "dsTMatrix" : "dgTMatrix");
	    break;
	case 1:
	    cl = uploT ? "ltTMatrix" :
		((a->stype) ? "lsTMatrix" : "lgTMatrix");
	    break;
	}
	break;
    case CHOLMOD_COMPLEX:
	cl = uploT ? "ztTMatrix" :
	    ((a->stype) ? "zsTMatrix" : "zgTMatrix");
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
    if (a->xtype == CHOLMOD_REAL) {
	int i, *m_x;
	switch(Rkind) {
	case 0:
	    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, a->nnz)),
		   (double *) a->x, a->nnz);
	    break;
	case 1:
	    m_x= LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, a->nnz));
	    for (i=0; i < a->nnz; i++)
		m_x[i] = (int) ((double *) a->x)[i];
/* 	    Memcpy(LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, a->nnz)), */
/* 		   (int *) a->x, a->nnz); */
	    break;
	}
    }
    else if (a->xtype == CHOLMOD_COMPLEX)
	error("complex sparse matrix code not yet written");
/* 	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, a->nnz)), */
/* 	       (complex *) a->x, a->nz); */
    if (uploT) {		/* slots for triangularMatrix */
	if (a->stype) error("Symmetric and triangular both set");
	SET_SLOT(ans, Matrix_uploSym, mkString((uploT > 0) ? "U" : "L"));
	SET_SLOT(ans, Matrix_diagSym, mkString(diag));
    }
				/* set symmetry attributes */
    if (a->stype)
	SET_SLOT(ans, Matrix_uploSym,
		 mkString((a->stype > 0) ? "U" : "L"));
    if (dofree > 0) cholmod_free_triplet(&a, &c);
    if (dofree < 0) Free(a);
    if (dn != R_NilValue)
	SET_SLOT(ans, Matrix_DimNamesSym, duplicate(dn));
    UNPROTECT(2);
    return ans;
}

/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * In most cases this function is called through the macro AS_CHM_DN.
 * It is unusual to call it directly.
 *
 * @param ans a CHM_DN pointer.
 * @param x pointer to an object that inherits from (denseMatrix ^ generalMatrix)
 *
 * @return ans containing pointers to the slots of x.
 */
CHM_DN as_cholmod_dense(CHM_DN ans, SEXP x)
{
#define _AS_cholmod_dense_1						\
    char *valid[] = {"dmatrix", "dgeMatrix",				\
		     "lmatrix", "lgeMatrix",				\
		     "nmatrix", "ngeMatrix",				\
		     "zmatrix", "zgeMatrix", ""};			\
    int dims[2], ctype = Matrix_check_class(class_P(x), valid), nprot = 0; \
									\
    if (ctype < 0) {		/* not a classed matrix */		\
	if (isMatrix(x)) Memcpy(dims, INTEGER(getAttrib(x, R_DimSymbol)), 2); \
	else {dims[0] = LENGTH(x); dims[1] = 1;}			\
	if (isInteger(x)) {						\
	    x = PROTECT(coerceVector(x, REALSXP));			\
	    nprot++;							\
	}								\
	ctype = (isReal(x) ? 0 :					\
		 (isLogical(x) ? 2 : /* logical -> default to "l", not "n" */ \
		  (isComplex(x) ? 6 : -1)));				\
    } else Memcpy(dims, INTEGER(GET_SLOT(x, Matrix_DimSym)), 2);	\
    if (ctype < 0) error("invalid class of object to as_cholmod_dense");\
    memset(ans, 0, sizeof(cholmod_dense)); /* zero the struct */        \
                                                                        \
    ans->dtype = CHOLMOD_DOUBLE; /* characteristics of the system */	\
    ans->x = ans->z = (void *) NULL;					\
				/* dimensions and nzmax */		\
    ans->d = ans->nrow = dims[0];					\
    ans->ncol = dims[1];						\
    ans->nzmax = dims[0] * dims[1];					\
				/* set the xtype and any elements */	\
    switch(ctype / 2) {							\
    case 0: /* "d" */							\
	ans->xtype = CHOLMOD_REAL;					\
	ans->x = (void *) REAL((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x); \
	break

    _AS_cholmod_dense_1;

    case 1: /* "l" */
	ans->xtype = CHOLMOD_REAL;
	ans->x =
	    (void *) REAL(coerceVector((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x,
				       REALSXP));
	break;
    case 2: /* "n" */
	ans->xtype = CHOLMOD_PATTERN;
	ans->x = (void *) LOGICAL((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x);
	break;

#define _AS_cholmod_dense_2						\
    case 3: /* "z" */							\
	ans->xtype = CHOLMOD_COMPLEX;					\
	ans->x = (void *) COMPLEX((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x); \
	break;								\
    }									\
    UNPROTECT(nprot);							\
    return ans

    _AS_cholmod_dense_2;
}

/* version of as_cholmod_dense() that produces a cholmod_dense matrix
 * with REAL 'x' slot -- i.e. treats "nMatrix" as "lMatrix" -- as only difference;
 * Not just via a flag in as_cholmod_dense() since that has fixed API */
CHM_DN as_cholmod_x_dense(cholmod_dense *ans, SEXP x)
{
    _AS_cholmod_dense_1;

    case 1: /* "l" */
    case 2: /* "n" (no NA in 'x', but *has* 'x' slot => treat as "l" */
	ans->xtype = CHOLMOD_REAL;
	ans->x =
	    (void *) REAL(coerceVector((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x,
				       REALSXP));
	break;

    _AS_cholmod_dense_2;
}

#undef _AS_cholmod_dense_1
#undef _AS_cholmod_dense_2

void R_cholmod_error(int status, char *file, int line, char *message)
{
    error(_("Cholmod error '%s' at file:%s, line %d"), message, file, line);
}

/* just to get 'int' instead of 'void' as required by CHOLMOD's print_function */
static
int R_cholmod_printf(const char* fmt, ...)
{
    va_list(ap);

    va_start(ap, fmt);
    Rprintf((char *)fmt, ap);
    va_end(ap);
    return 0;
}

/**
 * Initialize the CHOLMOD library and replace the print and error functions
 * by R-specific versions.
 *
 * @param c pointer to a cholmod_common structure to be initialized
 *
 * @return CHOLMOD_OK if successful
 */
int R_cholmod_start(CHM_CM c)
{
    int res;
    if (!(res = cholmod_start(c)))
	error(_("Unable to initialize cholmod: error code %d"), res);
    c->print_function = R_cholmod_printf; /* Rprintf gives warning */
    c->error_handler = R_cholmod_error;
    return TRUE;
}

/**
 * Copy the contents of a to an appropriate denseMatrix object and,
 * optionally, free a or free both a and its pointer to its contents.
 *
 * @param a matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 * @param Rkind type of R matrix to be generated (special to this function)
 * @param dn   -- dimnames [list(.,.) or NULL]
 *
 * @return SEXP containing a copy of a
 */

/* FIXME: should also have args  (int uploST, char *diag) */

SEXP chm_dense_to_SEXP(CHM_DN a, int dofree, int Rkind, SEXP dn)
{
    SEXP ans;
    char *cl = ""; /* -Wall */
    int *dims, ntot;

    PROTECT(dn); /* << (why? -- just cut&paste from chm_dense_to_mat.. below*/

    switch(a->xtype) {		/* determine the class of the result */
/* CHOLMOD_PATTERN never happens because cholmod_dense can't :
 *     case CHOLMOD_PATTERN:
 * 	cl = "ngeMatrix"; break;
 */
    case CHOLMOD_REAL:
	switch(Rkind) { /* -1: special for this function! */
	case -1: cl = "ngeMatrix"; break;
	case 0:	 cl = "dgeMatrix"; break;
	case 1:	 cl = "lgeMatrix"; break;
	default: error("unknown 'Rkind'");
	}
	break;
    case CHOLMOD_COMPLEX:
	cl = "zgeMatrix"; break;
    default:
	error("unknown xtype");
    }

    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cl)));
				/* allocate and copy common slots */
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = a->nrow; dims[1] = a->ncol;
    ntot = dims[0] * dims[1];
    if (a->d == a->nrow) {	/* copy data slot -- always present in dense(!) */
	if (a->xtype == CHOLMOD_REAL) {
	    int i, *m_x;
	    switch(Rkind) {
	    case 0:
		Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, ntot)),
		       (double *) a->x, ntot);
		break;
	    case -1: /* nge*/
	    case 1:  /* lge*/
		m_x = LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, ntot));
		for (i=0; i < ntot; i++)
		    m_x[i] = (int) ((double *) a->x)[i];
/* 		Memcpy(LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, ntot)), */
/* 		       (int *) a->x, ntot); */
		break;
	    }
	}
	else if (a->xtype == CHOLMOD_COMPLEX)
	    error("complex sparse matrix code not yet written");
/*	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, ntot)), */
/*	       (complex *) a->x, ntot); */
    } else error("code for cholmod_dense with holes not yet written");

    if (dofree > 0) cholmod_free_dense(&a, &c);
    if (dofree < 0) Free(a);
    if (dn != R_NilValue)
	SET_SLOT(ans, Matrix_DimNamesSym, duplicate(dn));
    UNPROTECT(2);
    return ans;
}

/**
 * Copy the contents of a to a matrix object and, optionally, free a
 * or free both a and its pointer to its contents.
 *
 * @param a cholmod_dense structure to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 * @param dn either R_NilValue or an SEXP suitable for the Dimnames slot.
 *
 * @return SEXP containing a copy of a as a matrix object
 */
SEXP chm_dense_to_matrix(CHM_DN a, int dofree, SEXP dn)
{
    SEXP ans;
    SEXPTYPE typ;

    PROTECT(dn);
				/* determine the class of the result */
    typ = (a->xtype == CHOLMOD_PATTERN) ? LGLSXP :
	((a->xtype == CHOLMOD_REAL) ? REALSXP :
	 ((a->xtype == CHOLMOD_COMPLEX) ? CPLXSXP : NILSXP));
    if (typ == NILSXP) error("unknown xtype");

    ans = PROTECT(allocMatrix(typ, a->nrow, a->ncol));
    if (a->d == a->nrow) {	/* copy data slot if present */
	if (a->xtype == CHOLMOD_REAL)
	    Memcpy(REAL(ans), (double *) a->x, a->nrow * a->ncol);
	else if (a->xtype == CHOLMOD_COMPLEX)
	    error("complex sparse matrix code not yet written");
	else if (a->xtype == CHOLMOD_PATTERN)
	    error("don't know if a dense pattern matrix makes sense");
/* 	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, a->nnz)), */
/* 	       (complex *) a->x, a->nz); */
    } else error("code for cholmod_dense with holes not yet written");

    if (dofree > 0) cholmod_free_dense(&a, &c);
    if (dofree < 0) Free(a);
    if (dn != R_NilValue)
        setAttrib(ans, R_DimNamesSymbol, duplicate(dn));
    UNPROTECT(2);
    return ans;
}

CHM_DN numeric_as_chm_dense(CHM_DN ans, double *v, int nr, int nc)
{
    ans->d = ans->nrow = nr;
    ans->ncol = nc;
    ans->nzmax = nr * nc;
    ans->x = (void *) v;
    ans->xtype = CHOLMOD_REAL;
    ans->dtype = CHOLMOD_DOUBLE;
    return ans;
}

/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * In most cases this function is called through the macro AS_CHM_FR.
 * It is unusual to call it directly.
 *
 * @param ans an CHM_FR object
 * @param x pointer to an object that inherits from CHMfactor
 *
 * @return ans containing pointers to the slots of x.
 */
CHM_FR as_cholmod_factor(CHM_FR ans, SEXP x)
{
    char *valid[] = {"dCHMsuper", "dCHMsimpl", "nCHMsuper", "nCHMsimpl", ""};
    int *type = INTEGER(GET_SLOT(x, install("type"))),
	ctype = Matrix_check_class(class_P(x), valid);
    SEXP tmp;

    if (ctype < 0) error("invalid class of object to as_cholmod_factor");
    memset(ans, 0, sizeof(cholmod_factor)); /* zero the struct */

    ans->itype = CHOLMOD_INT;	/* characteristics of the system */
    ans->dtype = CHOLMOD_DOUBLE;
    ans->z = (void *) NULL;
    ans->xtype = (ctype < 2) ? CHOLMOD_REAL : CHOLMOD_PATTERN;

    ans->ordering = type[0];	/* unravel the type */
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
    ans->ColCount = INTEGER(GET_SLOT(x, install("colcount")));
    ans->z = ans->x = (void *) NULL;
    if (ctype < 2) {
	tmp = GET_SLOT(x, Matrix_xSym);
	ans->x = REAL(tmp);
    }
    if (ans->is_super) {	/* supernodal factorization */
	ans->xsize = LENGTH(tmp);
	ans->maxcsize = type[4]; ans->maxesize = type[5];
	ans->i = (int*)NULL;
	tmp = GET_SLOT(x, install("super"));
	ans->nsuper = LENGTH(tmp) - 1; ans->super = INTEGER(tmp);
	/* Move these checks to the CHMfactor_validate function */
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
	ans->nz = INTEGER(GET_SLOT(x, install("nz")));
	ans->next = INTEGER(GET_SLOT(x, install("nxt")));
	ans->prev = INTEGER(GET_SLOT(x, install("prv")));
    }
    if (!cholmod_check_factor(ans, &c))
	error(_("failure in as_cholmod_factor"));
    return ans;
}

/**
 * Copy the contents of f to an appropriate dCHMfactor object and,
 * optionally, free f or free both f and its pointer to its contents.
 *
 * @param f cholmod_factor object to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 Free a
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_factor_to_SEXP(CHM_FR f, int dofree)
{
    SEXP ans;
    int *dims, *type;
    char *class = (char*) NULL;	/* -Wall */

    switch(f->xtype) {
    case CHOLMOD_REAL:
	class = f->is_super ? "dCHMsuper" : "dCHMsimpl";
	break;
    case CHOLMOD_PATTERN:
	class = f->is_super ? "nCHMsuper" : "nCHMsimpl";
	break;
    default:
	error(_("f->xtype of %d not recognized"), f->xtype);
    }
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(class)));
    if (f->minor < f->n)
	error(_("CHOLMOD factorization was unsuccessful"));
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = dims[1] = f->n;
				/* copy component of known length */
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_permSym, INTSXP, f->n)),
	   (int*)f->Perm, f->n);
    Memcpy(INTEGER(ALLOC_SLOT(ans, install("colcount"), INTSXP, f->n)),
	   (int*)f->ColCount, f->n);
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
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("nz"), INTSXP, f->n)),
	       (int*)f->nz, f->n);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("nxt"), INTSXP, f->n + 2)),
	       (int*)f->next, f->n + 2);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("prv"), INTSXP, f->n + 2)),
	       (int*)f->prev, f->n + 2);

    }
    if(dofree) {
	if (dofree > 0) cholmod_free_factor(&f, &c);
	else /* dofree < 0 */ Free(f);
    }
    UNPROTECT(1);
    return ans;
}

/* Placeholders; TODO: use checks above (search "CHMfactor_validate"): */

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

