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
#include LA_SYMM_MAT_DOUBLE_H 
#include "vi.h"

LaSymmMatDouble& LaSymmMatDouble::copy(const LaMatDouble &ob)
{
    int M = ob.size(0);

    // current scheme in copy() is to detach the left-hand-side
    // from whatever it was pointing to.
    //

    resize(ob);
    for (int i = 0; i < M; i++)
	for (int j = 0; j <= i; j++)
	    (*this)(i,j) = ob(i,j);
    return *this;
}

ostream& LaSymmMatDouble::printMatrix(ostream& s) const
{
    if (*info_) {   // print out only matrix info, not actual values
	*info_ = 0; // reset the flag
	s << "(" << size(0) << "x" << size(1) << ") " ;
	s << "Indices: " << index(0) << " " << index(1);
	s << " #ref: " << ref_count() ;
	s << " sa:" << shallow();
    } else {
	int M = size(0);
	int N = size(1);
	for (int i = 0; i < M; i++) {
	    for (int j = 0; j < N; j++)
		s << (*this)(i,j) << " ";
	    s << endl;
	}
    }
    return s;
}

LaSymmMatDouble::operator LaGenMatDouble()
{
  int M = size(0);
  int N = size(1);

  LaGenMatDouble G(M,N);

  for (int i = 0; i< M; i++)
    for (int j = 0; j < N; j++)
        G(i,j) = (*this)(i,j);

  return G;
}

LaSymmMatDouble::operator LaLowerTriangMatDouble()
{
  int M = size(0);
  int N = size(1);

  LaLowerTriangMatDouble Lower(M,N);

  Lower.copy(lower_data_);

  return Lower;
}

LaSymmMatDouble& LaSymmMatDouble::solve() const
{				// inverse
    LaSymmMatDouble *inv; //create a copy to return
    inv = new LaSymmMatDouble(*this); 

    int info;
    F77_CALL(dtrtri)('U', 'N', inv->size(0), &(*inv)(0,0), inv->gdim(0), info);
    if (info != 0)
	throw(LaException("LaSymmMatDouble::solve()",
			  "Non-zero return code from dtrtri"));
    return *inv;
}

LaMatDouble& LaSymmMatDouble::solve(LaMatDouble& B) const
{				// in-place solution
    VectorInt ipiv(size(0));
    int lwork = 5 * size(0), info;
    VectorDouble work(lwork);

    F77_CALL(dsysv)('L', size(0), B.size(1),
		    &lower_data_(0,0), gdim(0), &ipiv(0),
		    &B(0,0), B.gdim(0), &work(0), lwork, info);
    return B;
}

LaMatDouble& LaSymmMatDouble::solve(LaMatDouble& X, const LaMatDouble& B) const
{
    X.inject(B);
    return solve(X);
}

double LaSymmMatDouble::norm(char which) const
{
    VectorDouble work(size(0));	// only needed for Infinity norm
    return F77_CALL(dlansy)(which, 'L', size(0),
			    &(*this)(0,0), gdim(0), &work(0));
}

double LaSymmMatDouble::rcond(char which) const
{
    double val;
    VectorDouble work(5 * size(0));
    int info;
    VectorInt ipiv(size(0)), iwork(size(0));
    LaSymmMatDouble th(*this);	// create a copy to pass

    F77_CALL(dsytrf)('L', th.size(0), &th(0,0), th.gdim(0), &ipiv(0),
		     &work(0), 5*size(0), info);
    F77_CALL(dsycon)('L', th.size(0), &th(0,0), th.gdim(0), &ipiv(0),
		     norm('O'), val, &work(0), &iwork(0), info);
    return val;
}

SEXP LaSymmMatDouble::asSEXP() const
{
    int n = size(0);
    SEXP val = allocMatrix(REALSXP, n, n);
    F77_CALL(dlacpy)('L', n, n, &(*this)(0,0), gdim(0),
		     REAL(val), n);
    for (int i = 1; i < n; i++) // symmetrize the result
	for (int j = 0; j < i; j++)
	    REAL(val)[i * n + j] = REAL(val)[j * n + i];
    SEXP classes = allocVector(STRSXP, 2);
    STRING(classes)[0] = mkChar("Hermitian");
    STRING(classes)[1] = mkChar("Matrix");
    setAttrib(val, R_ClassSymbol, classes);
    return val;
}
