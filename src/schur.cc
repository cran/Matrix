//  Copyright (C) 2000-2000 the R Development Core Team

#include "schur.h"

static int dummy(const double& wr, const double& wi) { return 0; }

LaGenSchurDouble::LaGenSchurDouble(const LaMatDouble& a,
				   bool jobV = true) :
    a(a), wR(a.size(0)), wI(a.size(0)), vecs(a)
{
    char jobVS = (jobV) ? 'V' : 'N';
    int n = a.size(0);
    if (a.size(1) != n)
	throw(LaException("LaGenSchurDouble : only square matrices allowed"));

    int lwork = 5*n, info, sdim;
    double rconde, rcondv;	// never used
    VectorDouble work(lwork);
    VectorInt iwork(0), bwork(0);
    F77_CALL(dgeesx)(jobVS, 'N', &dummy, 'N', n, &a(0,0), a.gdim(0),
		     sdim, &wR(0), &wI(0), &vecs(0,0), n, rconde, rcondv,
		     &work(0), lwork, &iwork(0), 1, bwork, info);
    if (info != 0)
	throw(LaException("LaGenSchurDouble : non-zero info returned by dgeesx"));
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

SEXP LaGenSchurDouble::asSEXP() const
{
    SEXP ret = PROTECT(allocVector(VECSXP, 3));
    SEXP nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("values"));
    SET_STRING_ELT(nm, 1, mkChar("schur"));
    SET_STRING_ELT(nm, 2, mkChar("vectors"));
    setAttrib(ret, R_NamesSymbol, nm);
    SEXP classes = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(classes, 0, mkChar("schur.Matrix"));
    SET_STRING_ELT(classes, 1, mkChar("decomp"));
    setAttrib(ret, R_ClassSymbol, classes);
    if (complexVectors()) {
	int n = wR.size();
	SEXP val = allocVector(CPLXSXP, n);
	for (int i = 0; i < n; i++) {
	    COMPLEX(val)[i].r = wR(i);
	    COMPLEX(val)[i].i = wI(i);
	}
	SET_VECTOR_ELT(ret, 0, val);
    } else {
	SET_VECTOR_ELT(ret, 0, wR.asSEXP());
    }
    SET_VECTOR_ELT(ret, 1, a.asSEXP());
    SET_VECTOR_ELT(ret, 2, vecs.asSEXP());
    UNPROTECT(3);
    return ret;
}

    
