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

ostream& LaSymmMatDouble::printMatrix(ostream& s) const
{
    if (*info_) {   // print out only matrix info, not actual values
	*info_ = 0; // reset the flag
	s << "(" << size(0) << "x" << size(1) << ") " ;
	s << "Indices: " << index(0) << " " << index(1);
	s << " #ref: " << ref_count() ;
	s << " sa:" << shallow();
	s << " uplo:" << uplo();
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

//  LaSymmMatDouble::operator LaLowerTriangMatDouble()
//  {
//    LaLowerTriangMatDouble Lower;

//    if (uplo_ == 'U') {
//        Lower.copy(*this);
//    } else Lower.ref(*data_.lower);
//    return Lower;
//  }

//  LaSymmMatDouble::operator LaUpperTriangMatDouble()
//  {
//    LaUpperTriangMatDouble Upper;

//    if (uplo_ == 'L') {
//        Upper.copy(*this);
//    } else Upper.ref(*data_.upper);
//    return Upper;
//  }

double LaSymmMatDouble::norm(char which) const
{
    VectorDouble work(size(0));	// only needed for Infinity norm
    return F77_CALL(dlansy)(which, uplo(), size(0),
			    &(*this)(0,0), gdim(0), &work(0));
}

SEXP LaSymmMatDouble::asSEXP() const
{
    int n = size(0);
    SEXP val = PROTECT(allocMatrix(REALSXP, n, n));
    LaGenMatDouble tmp(REAL(val), n, n);
    tmp.inject(*this);
    SEXP classes = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(classes, 0, mkChar("Hermitian"));
    SET_STRING_ELT(classes, 1, mkChar("Matrix"));
    setAttrib(val, R_ClassSymbol, classes);
    UNPROTECT(2);
    return val;
}

