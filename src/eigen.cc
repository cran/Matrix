//  Copyright (C) 2000-2000 the R Development Core Team

#include "eigen.h"

LaSymmEigenDouble::LaSymmEigenDouble(const LaMatDouble& a, bool findVecs = true) :
    vals(a.size(0))
{
    char jobz = (findVecs) ? 'V' : 'N';
    int n = a.size(0);
    if (a.size(1) != n)
	throw(LaException("LaSymmEigenDouble : only square matrices allowed"));

    LaGenMatDouble aa(a);
    int lwork = 5 * n, info;
    VectorDouble work(lwork);
    F77_CALL(dsyev)(jobz, 'L', n, &aa(0,0), aa.gdim(0), &vals(0),
		    &work(0), lwork, info);
    if (info != 0)
	throw(LaException("LaSymmEigenDouble : non-zero info returned by dsyev"));
    if (findVecs) vecs.copy(aa);
}

SEXP LaSymmEigenDouble::asSEXP() const
{
    SEXP ret = PROTECT(allocVector(VECSXP, 2));
    SEXP nm = PROTECT(allocVector(STRSXP, 2));
    STRING(nm)[0] = mkChar("values");
    STRING(nm)[1] = mkChar("vectors");
    setAttrib(ret, R_NamesSymbol, nm);
    SEXP classes = PROTECT(allocVector(STRSXP, 1));
    STRING(classes)[0] = mkChar("eigen.Matrix");
    setAttrib(ret, R_ClassSymbol, classes);
    VECTOR(ret)[0] = vals.asSEXP();
    VECTOR(ret)[1] = vecs.asSEXP();
    UNPROTECT(3);
    return ret;
}

    
