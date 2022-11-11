#ifndef MATRIX_INLINES_H
#define MATRIX_INLINES_H

/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.  The validity of this function
 * depends on SET_SLOT not duplicating val when NAMED(val) == 0.  If
 * this behavior changes then ALLOC_SLOT must use SET_SLOT followed by
 * GET_SLOT to ensure that the value returned is indeed the SEXP in
 * the slot.
 * NOTE:  GET_SLOT(x, what)        :== R_do_slot       (x, what)
 * ----   SET_SLOT(x, what, value) :== R_do_slot_assign(x, what, value)
 * and the R_do_slot* are in src/main/attrib.c
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, R_xlen_t length)
{
    SEXP val = allocVector(type, length);

    SET_SLOT(obj, nm, val);
    return val;
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
static R_INLINE
int* expand_cmprPt(int ncol, const int mp[], int mj[])
{
    int j;
    for (j = 0; j < ncol; j++) {
	int j2 = mp[j+1], jj;
	for (jj = mp[j]; jj < j2; jj++) mj[jj] = j;
    }
    return mj;
}

/** Inverse Permutation
 * C version of   .inv.perm.R <- function(p) { p[p] <- seq_along(p) ; p }
 */
static R_INLINE
SEXP inv_permutation(SEXP p_, SEXP zero_p, SEXP zero_res)
{
    int np = 1;
    if(!isInteger(p_)) {p_ = PROTECT(coerceVector(p_, INTSXP)); np++; }
    int *p = INTEGER(p_), n = LENGTH(p_);
    SEXP val = PROTECT(allocVector(INTSXP, n));
    int *v = INTEGER(val), p_0 = asLogical(zero_p), r_0 = asLogical(zero_res);
    if(!p_0) v--; // ==> use 1-based indices
    // shorter (but not 100% sure if ok: is LHS always eval'ed *before* RHS ?) :
    // for(int i=0; i < n; ) v[p[i]] = ++i;
    for(int i=0; i < n; ) {
	int j = p[i]; v[j] = (r_0) ? i++ : ++i;
    }
    UNPROTECT(np);
    return val;
}

/**
 * Return the 0-based index of a string match in a vector of strings
 * terminated by an empty string.  Returns -1 for no match.
 * Is  __cheap__ :  __not__ looking at superclasses --> better use  R_check_class_etc(obj, *)
 *
 * @param x string to match
 * @param valid vector of possible matches terminated by an empty string
 *
 * @return index of match or -1 for no match
 */
static R_INLINE
int strmatch(char *x, const char **valid)
{
    int ans = 0;
    while (strlen(valid[ans]) > 0) {
	if (strcmp(x, valid[ans]) == 0)
	    return ans;
	++ans;
    }
    return -1;
}

static R_INLINE
int strmatch2(const char *x, SEXP valid)
{
    int i, n = LENGTH(valid);
    for (i = 0; i < n; ++i)
	if (strcmp(x, CHAR(STRING_ELT(valid, i))) == 0)
	    return i;
    return -1;
}

#endif /* MATRIX_INLINES_H */
