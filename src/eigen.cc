//  Copyright (C) 2000-2000, 2002 the R Development Core Team

#include "eigen.h"

LaSymmEigenDouble::LaSymmEigenDouble(const LaMatDouble& a,
				     const char uplo,
				     const bool findVecs)
    : vals(a.size(0)), vecs()
{
    char jobz = (findVecs) ? 'V' : 'N';
    int n = a.size(0);
    if (a.size(1) != n)
	throw(LaException("LaSymmEigenDouble : only square matrices allowed"));

    LaGenMatDouble temp(a);
    int lwork = 5 * n, info;
    VectorDouble work(lwork);
    F77_CALL(dsyev)(jobz, uplo, n, &temp(0, 0), temp.gdim(0), &vals(0),
		    &work(0), lwork, info);
    if (info != 0)
	throw(LaException("LaSymmEigenDouble : non-zero info returned by dsyev"));
    if (findVecs) vecs.copy(temp);
}

SEXP LaSymmEigenDouble::asSEXP() const
{
    SEXP ret = PROTECT(allocVector(VECSXP, 2));
    SEXP nm = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nm, 0, mkChar("values"));
    SET_STRING_ELT(nm, 1, mkChar("vectors"));
    setAttrib(ret, R_NamesSymbol, nm);
    SEXP classes = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(classes, 0, mkChar("eigen.Matrix"));
    SET_STRING_ELT(classes, 1, mkChar("decomp"));
    setAttrib(ret, R_ClassSymbol, classes);
    SET_VECTOR_ELT(ret, 0, vals.asSEXP());
    SET_VECTOR_ELT(ret, 1, vecs.asSEXP());
    UNPROTECT(3);
    return ret;
}
    
LaGenEigenDouble::LaGenEigenDouble(const LaMatDouble& a,
				   bool leftEV,
				   bool rightEV,
				   char balanc,
				   char sense)
    : wR(a.size(0)), wI(a.size(0)), scale(a.size(0)), rcondE(a.size(0)),
      rcondV(a.size(0))
{
    char jobVL = (leftEV) ? 'V' : 'N', jobVR = (rightEV) ? 'V' : 'N';
    int n = a.size(0);
    if (a.size(1) != n)
	throw(LaException("LaGenEigenDouble : only square matrices allowed"));

    LaGenMatDouble temp(a);
    int lwork = 5*(n*n + 2*n), info;
    VectorDouble work(lwork);
    VectorInt iwork(2*n);
    if (leftEV) left.resize(n, n);
    if (rightEV) right.resize(n, n);
    if (sense != 'E' && sense != 'B') rcondE.resize(0);
    if (sense != 'V' && sense != 'B') rcondV.resize(0);
    F77_CALL(dgeevx)(balanc, jobVL, jobVR, sense, n, &temp(0,0), temp.gdim(0),
		     &wR(0), &wI(0), &left(0,0), n, &right(0,0), n,
		     &ilo, &ihi, &scale(0), &abnrm, &rcondE(0), &rcondV(0),
		     &work(0), lwork, &iwork(0), info);
    if (info != 0)
	throw(LaException("LaGenEigenDouble : non-zero info returned by dgeevx"));
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
    SEXP ret = PROTECT(allocVector(VECSXP, 3));
    SEXP nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("values"));
    SET_STRING_ELT(nm, 1, mkChar("vectors"));
    SET_STRING_ELT(nm, 2, mkChar("rcond"));
    setAttrib(ret, R_NamesSymbol, nm);
    SEXP classes = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(classes, 0, mkChar("eigen.Matrix"));
    SET_STRING_ELT(classes, 1, mkChar("decomp"));
    setAttrib(ret, R_ClassSymbol, classes);
    SEXP vecs = PROTECT(allocVector(VECSXP, 2));
    SEXP vnms = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(vnms, 0, mkChar("left"));
    SET_STRING_ELT(vnms, 1, mkChar("right"));
    setAttrib(vecs, R_NamesSymbol, vnms);
    if (complexVectors()) {
	int n = wR.size();
	SEXP val = allocVector(CPLXSXP, n);
	for (int i = 0; i < n; i++) {
	    COMPLEX(val)[i].r = wR(i);
	    COMPLEX(val)[i].i = wI(i);
	}
	SET_VECTOR_ELT(ret, 0, val);

	if (left.size(0) == n)
	    SET_VECTOR_ELT(vecs, 0, unscramble(wI, left));
	else SET_VECTOR_ELT(vecs, 0, R_NilValue);

	if (right.size(0) == n)
	    SET_VECTOR_ELT(vecs, 1, unscramble(wI, right));
	else SET_VECTOR_ELT(vecs, 1, R_NilValue);

    } else {
	SET_VECTOR_ELT(ret, 0, wR.asSEXP());
	SET_VECTOR_ELT(vecs, 0, left.asSEXP());
	SET_VECTOR_ELT(vecs, 1, right.asSEXP());
    }
    SET_VECTOR_ELT(ret, 1, vecs);
    SEXP rcond = PROTECT(allocVector(VECSXP, 2)),
	rnms = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(rnms, 0, mkChar("values"));
    SET_STRING_ELT(rnms, 1, mkChar("vectors"));
    setAttrib(rcond, R_NamesSymbol, rnms);
    SET_VECTOR_ELT(rcond, 0, rcondE.asSEXP());
    SET_VECTOR_ELT(rcond, 1, rcondV.asSEXP());
    SET_VECTOR_ELT(ret, 2, rcond);
    UNPROTECT(7);
    return ret;
}

    
