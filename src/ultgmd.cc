//
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.
//
// Modifications Copyright (C) 2000-2000 the R Development Core Team
//

#include "lafnames.h"
#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H
#include "lapackd.h"
#include "blas3.h"
#include "vi.h"

double LaUnitLowerTriangMatDouble::outofbounds_ = 0; // initialize outofbounds.

double& LaUnitLowerTriangMatDouble::operator()(int i, int j)
{

#ifdef UNIT_LOWER_INDEX_CHK
    if (i <= j) {
	warning("index to Lower Triular matrix out of range!\n i = %d  j = %d",
		i, j);
    }
#endif
    if ((j==0)&&(i==0))  // this special case allows us to pass unit matrices
	return data_(0,0); // by accessing the first element of the matrix
    else if (i<=j)       // which under normal circumstances would return 
	return outofbounds_; //outofbounds_.
    else
	return data_(i,j);
}

double& LaUnitLowerTriangMatDouble::operator()(int i, int j) const
{

#ifdef UNIT_LOWER_INDEX_CHK
    if (i<=j) {
	warning("index to Lower Triular matrix out of range!\n i = %d  j = %d",
		i, j);
    }
#endif

    if ((j==0)&&(i==0))  // this special case allows us to pass unit matrices
	return data_(0,0); // by accessing the first element of the matrix
    else if (i<=j)       // which under normal circumstances would return
	return outofbounds_; //outofbounds_.
    else
	return data_(i,j);
}

LaUnitLowerTriangMatDouble& LaUnitLowerTriangMatDouble::copy(const LaMatDouble &ob)
{
    int M = ob.size(0);

    resize(ob);
    for (int i = 0; i < M; i++)
	for (int j = 0; j < i; j++)
	    (*this)(i,j) = ob(i,j);

    return *this;
}

LaMatDouble& LaUnitLowerTriangMatDouble::operator=(double s)
{
    int N = size(1);
  
    for (int j = 0; j < N; j++)
	for (int i = 0; i > j; i++)
	    (*this)(i,j) = s;

    return *this;
}

ostream& LaUnitLowerTriangMatDouble::printMatrix(ostream& s) const
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
		if (i > j)
		    s << (*this)(i,j) << "  ";
	    }
	    s << endl;
	}
    }
    return s;
}

LaUnitLowerTriangMatDouble* LaUnitLowerTriangMatDouble::solve() const
{				// inverse
    LaUnitLowerTriangMatDouble *inv; //create a copy to return
    inv = new LaUnitLowerTriangMatDouble(*this); 

    int info;
    F77_CALL(dtrtri)('L', 'U', inv->size(0), &(*inv)(0,0), inv->gdim(0), info);
    if (info != 0)		// can't happen - "this" is unit triangular
	throw(LaException("LaUnitLowerTriangMatDouble::solve()",
			  "Non-zero return code from dtrtri"));
    return inv;
}

LaMatDouble& LaUnitLowerTriangMatDouble::solve(LaMatDouble& B) const
{				// in-place solution
    F77_CALL(dtrsm)('L', 'L', 'N', 'U', size(0), B.size(1), 1.0,
		    &data_(0,0), gdim(0), &B(0,0), B.gdim(0));
    return B;
}

LaMatDouble& LaUnitLowerTriangMatDouble::solve(LaMatDouble& X, const LaMatDouble& B) const
{
    X.inject(B);
    return solve(X);
}
    
double LaUnitLowerTriangMatDouble::norm(char which) const
{
    VectorDouble work(size(0)); // only needed for Infinity norm
    double val = F77_CALL(dlantr)(which, 'L', 'U', size(0), size(1),
				  &(*this)(0,0), gdim(0), &work(0));
    return val;
}

double LaUnitLowerTriangMatDouble::rcond(char which) const
{
    double val;
    VectorDouble work(3 * size(0));
    int info;
    VectorInt iwork(size(0));
    F77_CALL(dtrcon)(which, 'L', 'U', size(0), &(*this)(0,0),
		     gdim(0), val, &work(0), &iwork(0), info);
    return val;
}

SEXP LaUnitLowerTriangMatDouble::asSEXP() const
{
    int m = size(0), n = size(1);
    SEXP val = PROTECT(allocMatrix(REALSXP, m, n));
    F77_CALL(dlacpy)('L', m, n, &(*this)(0,0), gdim(0),
		     REAL(val), m);
    int ldiag = (m < n) ? m : n;
    for (int i = 0; i < ldiag; i++) // ensure the diagonal entries are 1.0
	REAL(val)[i * (m + 1)] = 1.0;
    SEXP classes = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(classes, 0, mkChar("UnitLowerTriangular"));
    SET_STRING_ELT(classes, 1, mkChar("LowerTriangular"));
    SET_STRING_ELT(classes, 2, mkChar("Matrix"));
    setAttrib(val, R_ClassSymbol, classes);
    UNPROTECT(2);
    return val;
}
