//
// Copyright (C) 2000-2000 the R Development Core Team

#include "lapack++.h"
 
static int isMMatrix(SEXP s)
{
    if (isObject(s)) {
	SEXP classes = getAttrib(s, R_ClassSymbol);
	for (int i = 0; i < Rf_length(classes); i++)
	    if (!strcmp(CHAR(STRING(classes)[i]), "Matrix")) return 1;
    }
    return 0;
}

static int checkClass(SEXP classes, char *cname)
{
    for (int i = 0; i < Rf_length(classes); i++)
	if (!strcmp(CHAR(STRING(classes)[i]), cname)) return 1;
    return 0;
}
    
// copy an R matrix into a suitable class inheriting from LaMatrix
static LaMatDouble* asLaMatrix(SEXP x)
{
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
	
extern "C" {
    SEXP R_LapackPP_norm(SEXP a, SEXP which)
    {
	try {
	    if (!isString(which))
		error("R_LapackPP_norm : which should be of mode character");
	    LaMatDouble *aa = asLaMatrix(a);
	    SEXP val = ScalarReal(aa->norm(CHAR(STRING(which)[0])[0]));
	    delete aa;
	    return val;
	} catch(LaException xcp) {
	    error(xcp.what());
	    return R_NilValue;
	}
    }    

    SEXP R_LapackPP_rcond(SEXP a, SEXP which)
    {
	try {
	    if (!isString(which))
		error("R_LapackPP_rcond : which should be of mode character");
	    LaMatDouble *aa = asLaMatrix(a);
	    SEXP val = ScalarReal(aa->rcond(CHAR(STRING(which)[0])[0]));
	    delete aa;
	    return val;
	} catch(LaException xcp) {
	    error(xcp.what());
	    return R_NilValue;
	}
    }

    SEXP R_LapackPP_solve(SEXP a, SEXP b)
    {
	try {
	    LaMatDouble *bb = asLaMatrix(b); // a pointer to a copy of b
	    LaMatDouble *aa = asLaMatrix(a); // a pointer to a copy of a
	    LaGenMatDouble x;
	    SEXP val = PROTECT(allocMatrix(REALSXP, aa->size(1), bb->size(1)));
	    x.ref(val);
	    aa->solve(x, *bb);
	    UNPROTECT(1);
	    delete aa;
	    delete bb;
	    return val;
	} catch (LaException xcp) {
	    error(xcp.what());
	    return R_NilValue;
	}
    }
   
    SEXP R_LapackPP_svd(SEXP x, SEXP nu, SEXP nv)
    {
	try {
	    LaGenMatDouble a;
	    a.ref(x);		// it gets copied in the SVD constructor
	    SVD sv(a, INTEGER(coerceVector(nu, INTSXP))[0], 
		   INTEGER(coerceVector(nv, INTSXP))[0]);
	    SEXP val = PROTECT(allocVector(VECSXP, 3));
	    SEXP nm = PROTECT(allocVector(STRSXP, 3));
	    STRING(nm)[0] = mkChar("d");
	    STRING(nm)[1] = mkChar("u");
	    STRING(nm)[2] = mkChar("vt");
	    setAttrib(val, R_NamesSymbol, nm);
	    VECTOR(val)[0] = sv.getS().asSEXP();
	    VECTOR(val)[1] = sv.getU().asSEXP();
	    VECTOR(val)[2] = sv.getVT().asSEXP();
	    UNPROTECT(2);
	    return val;
	} catch (LaException xcp) {
	    error(xcp.what());
	    return R_NilValue;
	}
    }
}
