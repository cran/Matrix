/* C-level Methods for the ``abstract Index'' class
 *
 * Note: this heavily builds on ideas and code from  Jens Oehlschlaegel,
 * ----  as implemented (in the GPL'ed part of) package 'ff'.
 */

#include "abIndex.h"

/**
 * Integer RLE (Run Length Encoding) -- only when it's worth
 *
 * @param x  R vector   which can be coerced to "integer"
 *
 * @return NULL or a valid R object of class "rle"
 */
SEXP Matrix_int_rle(SEXP x_)
{
    int n = LENGTH(PROTECT(x_ = coerceVector(x_, INTSXP)));
    if (n < 3)
	return R_NilValue;
    else {
	register int lv,ln, i, c = 0;
	int n2 = n / 3;
	/* upper bound: ==> max RAM requirement 2 x n2, (= 2/3 n);
	 * using 2 instead of 3 would need 50% more time, have max
	 * RAM requirement 2.5x for savings of any size */
	int *x = INTEGER(x_);
	int *val, *len;
	const char *res_nms[] = {"lengths", "values", ""};
	SEXP ans;

	len = Calloc(n2, int);
	val = Calloc(n2, int);

	lv = x[0];
	ln = 1;
	for(i = 1; i < n; i++) {
	    if (x[i] == lv) {
		ln++;
	    } else {
		val[c] = lv;
		len[c] = ln;
		c++;
		if (c == n2) { /* reached the "efficiency bound" */
		    Free(len);
		    Free(val);
		    UNPROTECT(1);
		    return R_NilValue;
		}
		lv = x[i];
		ln = 1;
	    }
	}
	val[c] = lv;
	len[c] = ln;
	c++;

	ans = PROTECT(Matrix_make_named(VECSXP, res_nms));
	SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, c)); /* lengths */
	SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, c)); /* values */
	Memcpy(INTEGER(VECTOR_ELT(ans, 0)), len, c);
	Memcpy(INTEGER(VECTOR_ELT(ans, 1)), val, c);

	setAttrib(ans, R_ClassSymbol, mkString("rle"));

	Free(len);
	Free(val);
	UNPROTECT(2);
	return ans;
    }
}
