//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
//
// Modifications Copyright (C) 2000-2000 the R Development Core Team
//

#include "lafnames.h"
#include LA_LOWER_TRIANG_MAT_DOUBLE_H
#include "blas3.h"
#include "vi.h"

double LaLowerTriangMatDouble::outofbounds_ = 0; // initialize outofbounds.

LaLowerTriangMatDouble& LaLowerTriangMatDouble::copy(const LaMatDouble &ob)
{
    int M = ob.size(0);
    
    resize(ob);
    for (int i = 0; i < M; i++)
	for (int j = 0; j <= i; j++)
	    (*this)(i,j) = ob(i,j);

    return *this;
}

LaMatDouble& LaLowerTriangMatDouble::operator=(double s)
{
    int N = size(1);

    for (int j = 0; j < N; j++)
	for (int i = 0; i >= j; i++)
	    (*this)(i,j) = s;
    return *this;
}

ostream& LaLowerTriangMatDouble::printMatrix(ostream& s) const
{
    if (*info_) {   // print out only matrix info, not actual values
	*info_ = 0; // reset the flag
	s << "(" << size(0) << "x" << size(1) << ") " ;
	s << "Indices: " << index(0) << " " << index(1);
	s << " #ref: " << ref_count();
	s << " sa:" << shallow();
    } else {
	int M = size(0);
	int N = size(1);

	for (int i = 0; i < M; i++) {
	    for (int j = 0; j < N; j++) {
		if (i >= j)
		    s << (*this)(i,j) << "  ";
	    }
	    s << endl;
	}
    }
    return s;
}

LaLowerTriangMatDouble& LaLowerTriangMatDouble::solve() const
{				// inverse
    LaLowerTriangMatDouble *inv; //create a copy to return
    inv = new LaLowerTriangMatDouble(*this); 

    int info;
    F77_CALL(dtrtri)('L', 'N', inv->size(0), &(*inv)(0,0), inv->gdim(0), info);
    if (info != 0)
	throw(LaException("LaLowerTriangMatDouble::solve()",
			  "Non-zero return code from dtrtri"));
    return *inv;
}

LaMatDouble& LaLowerTriangMatDouble::solve(LaMatDouble& B) const
{				// in-place solution
    F77_CALL(dtrsm)('L', 'L', 'N', 'N', size(0), B.size(1), 1.0,
		    &data_(0,0), gdim(0), &B(0,0), B.gdim(0));
    return B;
}

LaMatDouble& LaLowerTriangMatDouble::solve(LaMatDouble& X, const LaMatDouble& B) const
{
    X.inject(B);
    return solve(X);
}

double LaLowerTriangMatDouble::norm(char which) const
{
    VectorDouble work(size(0));	// only needed for Infinity norm
    return F77_CALL(dlantr)(which, 'L', 'N', size(0), size(1),
			    &(*this)(0,0), gdim(0), &work(0));
}

double LaLowerTriangMatDouble::rcond(char which) const
{
    double val;
    VectorDouble work(3 * size(0));
    int info;
    VectorInt iwork(size(0));
    F77_CALL(dtrcon)(which, 'L', 'N', size(0), &(*this)(0,0),
		     gdim(0), val, &work(0), &iwork(0), info);
    return val;
}

SEXP LaLowerTriangMatDouble::asSEXP() const
{
    SEXP val = allocMatrix(REALSXP, size(0), size(1));
    F77_CALL(dlacpy)('L', size(0), size(1), &(*this)(0,0), gdim(0),
		     REAL(val), size(0));
    SEXP classes = allocVector(STRSXP, 2);
    STRING(classes)[0] = mkChar("LowerTriangular");
    STRING(classes)[1] = mkChar("Matrix");
    setAttrib(val, R_ClassSymbol, classes);
    return val;
}
