//
// Copyright (C) 2000-2000 the R Development Core Team

#include "lapack++.h"
#include "eigen.h"
#include "lavi.h"
 
#include <iostream>
#include <cstdlib>
#include <new>
using namespace std;

static int count = 0;
static bool setNewHandler = false;

void out_of_memory() {
    cerr << "memory exhausted after " << count 
	 << " allocations!" << endl;
    exit(1);
}

static LaVectorInt* piv2perm(const LaVectorInt& piv)
{				// transform a pivot vector to the permutation
    LaVectorInt *perm;
    int n = piv.size();

    perm = new LaVectorInt(n);
    for (int i = 0; i < n; i++) (*perm)(i) = i + 1;
    for (int i = 0; i < n; i++) {
	int tmp = (*perm)(i);
	(*perm)(i) = (*perm)(piv(i) - 1);
	(*perm)(piv(i) - 1) = tmp;
    }
    return perm;
}

static int isMMatrix(SEXP s)
{
    if (isObject(s)) {
	SEXP classes = getAttrib(s, R_ClassSymbol);
	for (int i = 0; i < Rf_length(classes); i++)
	    if (!strcmp(CHAR(STRING_ELT(classes, i)), "Matrix")) return 1;
    }
    return 0;
}

static int checkClass(SEXP classes, char *cname)
{
    for (int i = 0; i < Rf_length(classes); i++)
	if (!strcmp(CHAR(STRING_ELT(classes, i)), cname)) return 1;
    return 0;
}
    
// copy an R matrix into a suitable class inheriting from LaMatrix
static LaMatDouble* asLaMatrix(SEXP x)
{
    if (isComplex(x)) error("Complex Matrix classes not yet implemented");
    if (isMMatrix(x)) {
	SEXP classes = getAttrib(x, R_ClassSymbol);
	if (checkClass(classes, "UnitLowerTriangular"))
	    return new LaUnitLowerTriangMatDouble(x);
	if (checkClass(classes, "UnitUpperTriangular"))
	    return new LaUnitUpperTriangMatDouble(x);
	if (checkClass(classes, "LowerTriangular"))
	    return new LaLowerTriangMatDouble(x);
	if (checkClass(classes, "UpperTriangular"))
	    return new LaUpperTriangMatDouble(x);
	if (checkClass(classes, "Hermitian"))
	    return new LaSymmMatDouble(x);
	return new LaGenMatDouble(x);
    }
    if (isMatrix(x)) return new LaGenMatDouble(x);
    if (isNumeric(x)) return new LaVectorDouble(x);
    error("object must be a matrix or a numeric vector a Classed Matrix");
    return new LaGenMatDouble(); // to keep -Wall happy
}

// create a reference to an R matrix in a suitable class
// inheriting from LaMatrix
LaMatDouble* asLaRef(SEXP x)
{
    if (isComplex(x)) error("Complex Matrix classes not yet implemented");
    LaMatDouble *val = 0;
    // We can't coerce because we will lose track of the pointer.
    if (!isReal(x)) error("object must be a numeric matrix or vector");
    if (isMMatrix(x)) {
	SEXP classes = getAttrib(x, R_ClassSymbol);
	if (checkClass(classes, "UnitLowerTriangular")) {
	    val = new LaUnitLowerTriangMatDouble();
	    val->ref(x);
	    return val;
	}
	if (checkClass(classes, "UnitUpperTriangular")) {
	    val = new LaUnitUpperTriangMatDouble();
	    val->ref(x);
	    return val;
	}
	if (checkClass(classes, "LowerTriangular")) {
	    val = new LaLowerTriangMatDouble();
	    val->ref(x);
	    return val;
	}
	if (checkClass(classes, "UpperTriangular")) {
	    val = new LaUpperTriangMatDouble();
	    val->ref(x);
	    return val;
	}
	if (checkClass(classes, "Hermitian")) {
	    val = new LaSymmMatDouble();
	    val->ref(x);
	    return val;
	}
	val = new LaGenMatDouble();
	val->ref(x);
	return val;
    }
    if (isMatrix(x)) { 
	val = new LaGenMatDouble();
	val->ref(x);
	return val;
    }
    val = new LaVectorDouble();
    val->ref(x);
    return val;
}
	
//inline double abs(x) { (x < 0) ? -x : x; }
extern "C" {
    SEXP R_LapackPP_det(SEXP x, SEXP logarithm)
    {
	if (isComplex(x)) error("Complex Matrix classes not yet implemented");
	LaGenMatDouble xx(x);
	int n = xx.size(0);
	if (xx.size(1) != n) {error("x must be a square matrix");}
	int lwork = 5 * n, info;
	LaVectorDouble work(lwork), tau(n);
	F77_CALL(dgeqrf)(n, n, &xx(0,0), xx.gdim(0), &tau(0), &work(0),
			 lwork, info);
	int sign = (n % 2) ? 1 : -1,
	    useLog = LOGICAL(coerceVector(logarithm, LGLSXP))[0];
	double modulus;
	if (useLog) {
	    modulus = 0.0;
	    for (int i = 0; i < n; i++) {
		modulus += log(abs(xx(i,i)));
		if (xx(i,i) < 0) sign = -sign;
	    }
	} else {
	    modulus = 1.0;
	    for (int i = 0; i < n; i++) { modulus *= xx(i, i); }
	    if (modulus < 0) {
		modulus = - modulus;
		sign = -sign;
	    }
	}
	SEXP val = PROTECT(allocVector(VECSXP, 2));
	SEXP nm = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(nm, 0, mkChar("modulus"));
	SET_STRING_ELT(nm, 1, mkChar("sign"));
	setAttrib(val, R_NamesSymbol, nm);
	SET_VECTOR_ELT(val, 0, ScalarReal(modulus));
	setAttrib(VECTOR_ELT(val, 0), install("logarithm"), ScalarLogical(useLog));
	SET_VECTOR_ELT(val, 1, ScalarInteger(sign));
	setAttrib(val, R_ClassSymbol, ScalarString(mkChar("det")));
	UNPROTECT(2);
	return val;
    }
	
    SEXP R_LapackPP_eigen(SEXP x, SEXP vectors)
    {
	LaMatDouble *aa = 0;
	LaEigenDouble *eig;
	try {
	    aa = asLaMatrix(x);
	    bool vecs = (LOGICAL(coerceVector(vectors, LGLSXP))[0] != 0) ? true
		: false;
	    if (aa->size(0) != aa->size(1)) {
		delete aa;
		error("Eigenvalue/eigenvector decompositions are defined only for square matrices");
	    }
	    eig = aa->eigen(vecs, vecs);
	    SEXP val = eig->asSEXP();
	    delete aa;
	    delete eig;
	    return val;
	} catch(LaException xcp) {
	    delete aa;
	    delete eig;
	    error(xcp.what());
	    return R_NilValue;	// to keep -Wall happy
	}
    }    

    SEXP R_LapackPP_lu(SEXP x, SEXP normComp)
    {
	LaGenMatDouble *xx = 0;
	LaLUFactorDouble *fact = 0;
	LaVectorInt *perm = 0;
	try {
	    xx = new LaGenMatDouble(x);
	    fact = new LaLUFactorDouble(*xx);
	    SEXP val = PROTECT(allocVector(VECSXP, 3));
	    SET_VECTOR_ELT(val, 0, fact->L().asSEXP());
	    SET_VECTOR_ELT(val, 1, fact->U().asSEXP());
	    perm = piv2perm(fact->pivot());
	    SET_VECTOR_ELT(val, 2, perm->asSEXP());
	    SEXP nm = PROTECT(allocVector(STRSXP, 3));
	    SET_STRING_ELT(nm, 0, mkChar("l"));
	    SET_STRING_ELT(nm, 1, mkChar("u"));
	    SET_STRING_ELT(nm, 2, mkChar("permutation"));
	    setAttrib(val, R_NamesSymbol, nm);
	    UNPROTECT(1);
	    setAttrib(val, R_ClassSymbol, ScalarString(mkChar("lu.Matrix")));
	    SEXP norms = PROTECT(allocVector(VECSXP, 2));
	    nm = PROTECT(allocVector(STRSXP, 2));
	    SET_STRING_ELT(nm, 0, mkChar("One"));
	    SET_STRING_ELT(nm, 1, mkChar("Infinity"));
	    setAttrib(norms, R_NamesSymbol, nm);
	    UNPROTECT(1);
	    PROTECT(normComp = coerceVector(normComp, LGLSXP));
	    if (Rf_length(normComp) > 0 && LOGICAL(normComp)[0] != 0) {
		SET_VECTOR_ELT(norms, 0, ScalarReal(xx->norm('O')));
	    } else {
		SET_VECTOR_ELT(norms, 0, R_NilValue);
	    }
	    if (Rf_length(normComp) > 1 && LOGICAL(normComp)[1] != 0) {
		SET_VECTOR_ELT(norms, 1, ScalarReal(xx->norm('I')));
	    } else {
		SET_VECTOR_ELT(norms, 1, R_NilValue);
	    }
	    setAttrib(val, install("norms"), norms);
	    delete xx; delete fact; delete perm; UNPROTECT(3);
	    return val;
	} catch(LaException xcp) {
	    delete xx; delete fact;
	    error(xcp.what());
	    return R_NilValue;
	}
    }
	    
    SEXP R_LapackPP_luH(SEXP x, SEXP lower, SEXP normComp)
    {
	LaSymmMatDouble *xx = 0;
	LaBunchKaufmanFactorDouble *fact = 0;
	try {
	    xx = new LaSymmMatDouble(x);
	    fact = new LaBunchKaufmanFactorDouble(*xx);
	    SEXP val = PROTECT(allocVector(VECSXP, 2));
	    SET_VECTOR_ELT(val, 0, fact->decomp().asSEXP());
	    SET_VECTOR_ELT(val, 1, fact->pivot().asSEXP());
	    SEXP nm = PROTECT(allocVector(STRSXP, 2));
	    SET_STRING_ELT(nm, 0, mkChar("decomp"));
	    SET_STRING_ELT(nm, 2, mkChar("permutation"));
	    setAttrib(val, R_NamesSymbol, nm);
	    setAttrib(val, R_ClassSymbol, ScalarString(mkChar("lu.Hermitian")));
	    delete xx; delete fact; UNPROTECT(2); return val;
	} catch(LaException xcp) {
	    delete xx; delete fact;
	    error(xcp.what());
	    return R_NilValue;
	}
    }
		
    SEXP R_LapackPP_norm(SEXP a, SEXP which)
    {
	LaMatDouble *aa = 0;
	try {
	    if (!isString(which))
		error("R_LapackPP_norm : which should be of mode character");
	    aa = asLaMatrix(a);
	    SEXP val =
		PROTECT(ScalarReal(aa->norm(CHAR(STRING_ELT(which, 0))[0])));
	    setAttrib(val, R_ClassSymbol, ScalarString(mkChar("norm")));
	    delete aa;
	    UNPROTECT(1);
	    return val;
	} catch(LaException xcp) {
	    delete aa;
	    error(xcp.what());
	    return R_NilValue;	// to keep -Wall happy
	}
    }    

    SEXP R_LapackPP_rcond(SEXP a, SEXP which)
    {
	LaMatDouble *aa = 0;
	try {
	    if (!isString(which))
		error("R_LapackPP_rcond : which should be of mode character");
	    aa = asLaMatrix(a);
	    SEXP val =
		PROTECT(ScalarReal(aa->rcond(CHAR(STRING_ELT(which, 0))[0])));
	    setAttrib(val, R_ClassSymbol, ScalarString(mkChar("rcond")));
	    delete aa;
	    UNPROTECT(1);
	    return val;
	} catch(LaException xcp) {
	    delete aa;
	    error(xcp.what());
	    return R_NilValue;
	}
    }

    SEXP R_LapackPP_solve(SEXP a, SEXP b)
    {
	LaMatDouble *aa = 0, *bb = 0;
	try {
	    if (b == R_NilValue) { // calculate inverse
		aa = asLaMatrix(a);
		if (aa->size(0) != aa->size(1))
		    error("Only square matrices can be inverted");
		SEXP val = aa->solve()->asSEXP();
		delete aa;
		return val;
	    }
	    bb = asLaMatrix(b); // a pointer to a copy of b
	    aa = asLaMatrix(a); // a pointer to a copy of a
	    SEXP val = aa->solve(*bb).asSEXP();
	    delete aa;
	    delete bb;
	    return val;
	} catch (LaException xcp) {
	    delete aa;
	    delete bb;
	    error(xcp.what());
	    return R_NilValue;
	}
    }
   
    SEXP R_LapackPP_svd(SEXP x, SEXP nu, SEXP nv)
    {
	if (isComplex(x)) error("Complex Matrix classes not yet implemented");
	try {
	    LaGenMatDouble a;
	    int nnu = INTEGER(coerceVector(nu, INTSXP))[0],
		nnv = INTEGER(coerceVector(nv, INTSXP))[0];
	    a.ref(x);		// it gets copied in the SVD constructor
	    SVD sv(a, nnu, nnv);
	    SEXP val = PROTECT(allocVector(VECSXP, 3));
	    SEXP nm = PROTECT(allocVector(STRSXP, 3));
	    SET_STRING_ELT(nm, 0, mkChar("d"));
	    SET_STRING_ELT(nm, 1, mkChar("u"));
	    SET_STRING_ELT(nm, 2, mkChar("vt"));
	    setAttrib(val, R_NamesSymbol, nm);
	    SET_VECTOR_ELT(val, 0, sv.getS().asSEXP());
	    if (nnu > 0) 
		SET_VECTOR_ELT(val, 1, sv.getU().asSEXP());
	    else
		SET_VECTOR_ELT(val, 1, R_NilValue);
	    if (nnv > 0)
		SET_VECTOR_ELT(val, 2, sv.getVT().asSEXP());
	    else
		SET_VECTOR_ELT(val, 2, R_NilValue);
	    UNPROTECT(2);
	    return val;
	} catch (LaException xcp) {
	    error(xcp.what());
	    return R_NilValue;
	}
    }
}
