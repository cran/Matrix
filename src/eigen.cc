//  Copyright (C) 2000-2000 the R Development Core Team

#include "eigen.h"

LaSymmEigenDouble::LaSymmEigenDouble(const LaMatDouble& a,
				     const char uplo,
				     const bool findVecs = true) :
    vals(a.size(0))
{
    char jobz = (findVecs) ? 'V' : 'N';
    int n = a.size(0);
    if (a.size(1) != n)
	throw(LaException("LaSymmEigenDouble : only square matrices allowed"));

    LaGenMatDouble aa(a);
    int lwork = 5 * n, info;
    VectorDouble work(lwork);
    F77_CALL(dsyev)(jobz, uplo, n, &aa(0,0), aa.gdim(0), &vals(0),
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
    SEXP classes = PROTECT(allocVector(STRSXP, 2));
    STRING(classes)[0] = mkChar("eigen.Matrix");
    STRING(classes)[1] = mkChar("decomp");
    setAttrib(ret, R_ClassSymbol, classes);
    VECTOR(ret)[0] = vals.asSEXP();
    VECTOR(ret)[1] = vecs.asSEXP();
    UNPROTECT(3);
    return ret;
}

    
LaGenEigenDouble::LaGenEigenDouble(const LaMatDouble& a,
				   bool leftEV = true,
				   bool rightEV = true) :
    wR(a.size(0)), wI(a.size(0))
{
    char jobVL = (leftEV) ? 'V' : 'N', jobVR = (rightEV) ? 'V' : 'N';
    int n = a.size(0);
    if (a.size(1) != n)
	throw(LaException("LaGenEigenDouble : only square matrices allowed"));

    LaGenMatDouble aa(a);
    int lwork = 16 * n, info;
    VectorDouble work(lwork);
    if (leftEV) left.resize(n, n);
    if (rightEV) right.resize(n, n);
    F77_CALL(dgeev)(jobVL, jobVR, n, &aa(0,0), aa.gdim(0), &wR(0), &wI(0),
		    &left(0,0), n, &right(0,0), n, &work(0), lwork, info);
    if (info != 0)
	throw(LaException("LaGenEigenDouble : non-zero info returned by dgeev"));
    complexVectors_ = false;
    for (int i = 0; i < n; i++)
	if (wI(i) != 0) { complexVectors_ = true; break; }
}

static SEXP unscramble(const LaVectorDouble& imaginary,
		       const LaGenMatDouble& vecs)
{
    int n = vecs.size(1);
    SEXP s = allocMatrix(CPLXSXP, n, n);

    for (int j = 0; j < n; j++) {
	if (imaginary(j) != 0) {
	    int j1 = j + 1;
	    for (int i = 0; i < n; i++) {
		COMPLEX(s)[i+n*j].r = COMPLEX(s)[i+n*j1].r = vecs(i, j);
		COMPLEX(s)[i+n*j1].i = -(COMPLEX(s)[i+n*j].i = vecs(i, j1));
	    }
	    j = j1;
	} else {
	    for (int i = 0; i < n; i++) {
		COMPLEX(s)[i+n*j].r = vecs(i, j);
		COMPLEX(s)[i+n*j].i = 0.0;
	    }
	}
    }
    return s;
}

SEXP LaGenEigenDouble::asSEXP() const
{
    SEXP ret = PROTECT(allocVector(VECSXP, 2));
    SEXP nm = PROTECT(allocVector(STRSXP, 2));
    STRING(nm)[0] = mkChar("values");
    STRING(nm)[1] = mkChar("vectors");
    setAttrib(ret, R_NamesSymbol, nm);
    SEXP classes = PROTECT(allocVector(STRSXP, 2));
    STRING(classes)[0] = mkChar("eigen.Matrix");
    STRING(classes)[1] = mkChar("decomp");
    setAttrib(ret, R_ClassSymbol, classes);
    SEXP vecs = PROTECT(allocVector(VECSXP, 2));
    STRING(nm)[0] = mkChar("left");
    STRING(nm)[1] = mkChar("right");
    setAttrib(vecs, R_NamesSymbol, nm);

    if (complexVectors()) {
	int n = wR.size();
	SEXP val = allocVector(CPLXSXP, n);
	for (int i = 0; i < n; i++) {
	    COMPLEX(val)[i].r = wR(i);
	    COMPLEX(val)[i].i = wI(i);
	}
	VECTOR(ret)[0] = val;
	if (left.size(0) == n)
	    VECTOR(vecs)[0] = unscramble(wI, left);
	else VECTOR(vecs)[0] = R_NilValue;

	if (right.size(0) == n)
	    VECTOR(vecs)[1] = unscramble(wI, right);
	else VECTOR(vecs)[1] = R_NilValue;

    } else {
	VECTOR(ret)[0] = wR.asSEXP();
	VECTOR(vecs)[0] = left.asSEXP();
	VECTOR(vecs)[1] = right.asSEXP();
    }
    VECTOR(ret)[1] = vecs;
    UNPROTECT(4);
    return ret;
}

    
