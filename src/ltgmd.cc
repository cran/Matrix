//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
//
// Modifications Copyright (C) 2000-2000 the R Development Core Team
//

#include "lafnames.h"
#include LA_LOWER_TRIANG_MAT_DOUBLE_H
#include "blas3.h"

double LaLowerTriangMatDouble::outofbounds_ = 0; // initialize outofbounds.

LaMatrix& LaLowerTriangMatDouble::copy(const LaMatrix &ob)
{
    if (debug()) {
	cout << " ob: " << ob.info() << endl;
    }

    int M = ob.size(0);
    
    resize(ob);
    for (int i = 0; i < M; i++)
	for (int j = 0; j <= i; j++)
	    (*this)(i,j) = ob(i,j);

    if (debug()) {
	cout << " *this: " << this->info() << endl;
    }

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

LaMatrix& LaLowerTriangMatDouble::solve() const
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

LaMatrix& LaLowerTriangMatDouble::solve(LaMatrix& B) const
{				// in-place solution
    F77_CALL(dtrsm)('L', 'L', 'N', 'N', size(0), B.size(1), 1.0,
		    &data_(0,0), gdim(0), &B(0,0), B.gdim(0));
    return B;
}

LaMatrix& LaLowerTriangMatDouble::solve(LaMatrix& X, const LaMatrix& B) const
{
    X.inject(B);
    return solve(X);
}

double LaLowerTriangMatDouble::norm(char which) const
{
    assert(which == 'M' || which == 'm' ||
	   which == '1' || which == 'O' || which == 'o' ||
	   which == 'I' || which == 'i' ||
	   which == 'F' || which == 'f' ||
	   which == 'E' || which == 'e');
    double *work = new double[size(0)]; // only needed for Infinity norm
    double val = F77_CALL(dlantr)(which, 'L', 'N', size(0), size(1),
				  &(*this)(0,0), gdim(0), work);
    delete[] work;
    return val;
}
