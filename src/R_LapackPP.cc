//
// Copyright (C) 2000-2000 the R Development Core Team

#include "lapack++.h"
#include "SVD.h"
 
int isMMatrix(SEXP s)
{
    if (isObject(s)) {
	SEXP classes = getAttrib(s, R_ClassSymbol);
	for (int i = 0; i < Rf_length(classes); i++)
	    if (!strcmp(CHAR(STRING(classes)[i]), "Matrix")) return 1;
    }
    return 0;
}

int checkClass(SEXP classes, char *cname)
{
    for (int i = 0; i < Rf_length(classes); i++)
	if (!strcmp(CHAR(STRING(classes)[i]), cname)) return 1;
    return 0;
}
    
LaMatrix* asLaMat(SEXP x)
{   // copy an R matrix into a suitable class inheriting from LaMatrix
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

LaMatrix* asLaRef(SEXP x)
{   // create a reference to an R matrix into a suitable class
    // inheriting from LaMatrix
    LaMatrix *val = 0;
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
//    SEXP R_LapackPP_solve2(SEXP a, SEXP b)
    SEXP R_LapackPP_solve(SEXP a, SEXP b)
    {
	LaMatrix *bb = asLaMat(b); // a pointer to a copy of b
	LaMatrix *aa = asLaMat(a); // a pointer to a copy of a
	LaGenMatDouble x;
	SEXP val = PROTECT(allocMatrix(REALSXP, aa->size(1), bb->size(1)));
	x.ref(val);
	aa->solve(x, *bb);
	UNPROTECT(1);
	delete aa;
	delete bb;
	return val;
    }
    
//      SEXP R_LapackPP_solve1(SEXP a)
//      {
//  	LaMatrix *aa = asLaMat(a);
//  	int m = aa->size(0), n = aa->size(1);
//  	if (m != n) error("only square matrices can be inverted");
//  	LaGenMatDouble ainv = aa->solve(); // fix this
//  	SEXP val = PROTECT(allocMatrix(REALSXP, m, m));
//  	Memcpy(REAL(val), &(ainv)(0,0), m * m);
//  	delete aa;
//  	UNPROTECT(1);
//  	return val;
//      }
   
    SEXP R_LapackPP_norm(SEXP a, SEXP which)
    {
	if (!isString(which))
	    error("R_LapackPP_norm : which should be of mode character");
	LaMatrix *aa = asLaMat(a);
	SEXP val = ScalarReal(aa->norm(CHAR(STRING(which)[0])[0]));
	delete aa;
	return val;
    }    

    SEXP R_LapackPP_svd(SEXP x, SEXP nu, SEXP nv)
    {
	LaGenMatDouble a(x);	// create a copy of x
	int ncu = INTEGER(coerceVector(nu, INTSXP))[0];
	if (ncu != 0 && ncu != a.size(0) && ncu != a.size(1))
	    error("nu must be one of 0, %d, or %d", a.size(0), a.size(1));
	int ncvt = INTEGER(coerceVector(nv, INTSXP))[0];
	if (ncvt != 0 && ncvt != a.size(1))
	    error("nv must be 0 or %d", a.size(1));
				// establish the return value
	SEXP val = PROTECT(allocVector(VECSXP, 3));
	SEXP nm = PROTECT(allocVector(STRSXP, 3));
	STRING(nm)[0] = mkChar("d");
	STRING(nm)[1] = mkChar("u");
	STRING(nm)[2] = mkChar("vt");
	setAttrib(val, R_NamesSymbol, nm);
	VECTOR(val)[0] = allocVector(REALSXP, min(a.size(0), a.size(1)));
	VectorDouble vv(REAL(VECTOR(val)[0]), min(a.size(0), a.size(1)));
				// create lapack++ objects and decompose
	LaGenMatDouble u, v;
	u.ref(VECTOR(val)[1] = allocMatrix(REALSXP, a.size(0), ncu));
	v.ref(VECTOR(val)[2] = allocMatrix(REALSXP, a.size(1), ncvt));
	SVD sv(a, vv, u, v);
	if (sv.getInfo() != 0)
	    error("dgesvd returned an info of %d", sv.getInfo());
	UNPROTECT(2);
	return val;
    }
}
