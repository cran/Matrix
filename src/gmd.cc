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

#include "lafnames.h"
#include LA_GEN_MAT_DOUBLE_H
#include "gfd.h"
#include "qr.h"

LaGenMatDouble::~LaGenMatDouble()
{
    delete solver;
}

LaGenMatDouble::LaGenMatDouble() : v(0)
{
    dim[0] = dim[1] = 0;
    sz[0] = sz[1] = 0;
    *info_ = 0;
    shallow_= 0;
    solver = 0;
}

LaGenMatDouble::LaGenMatDouble(int m, int n) :v(m*n)
{
    ii[0](0,m-1);
    ii[1](0,n-1);
    dim[0] = m;
    dim[1] = n;
    sz[0] = m;
    sz[1] = n;
    *info_ = 0;
    shallow_= 0;
    solver = 0;
}

LaGenMatDouble::LaGenMatDouble(double *d, int m, int n) :v(d, m*n)
{
    if (debug())
    {
        cout << ">>> LaGenMatDouble::LaGenMatDouble(double *d, int m, int n) :v(d, m*n)\n";
    }

    ii[0](0,m-1);
    ii[1](0,n-1);
    dim[0] = m;
    dim[1] = n;
    sz[0] = m;
    sz[1] = n;
    *info_ = 0;
    shallow_= 0;  
    solver = 0;
    
    if (debug())
    {
        cout << "<<< LaGenMatDouble::LaGenMatDouble(double *d, int m, int n) :v(d, m*n)\n";
    }
}

LaMatDouble& LaGenMatDouble::operator=(double s)
{
    if (debug())
    {
        cout << ">>> LaGenMatDouble& LaGenMatDouble::operator=(double s)\n";
        cout << "       s = " << s << endl;
    }

    for (int j = 0; j < size(1); j++)
    {
        for (int i = 0; i < size(0); i++)
        {
            (*this)(i,j) = s;
        }
    }

    if (solver != 0)
	clearDecomposition();

    if (debug())
    {
	cout << "  *this: " << this->info() << endl;
	cout << " >>> LaGenMatDouble& LaGenMatDouble::operator=(double s)\n";
    }
    return *this;
}

LaMatrix& LaGenMatDouble::ref(const LaGenMatDouble& s)
{
    if (this == &s) return *this; // handle trivial M.ref(M) case
    
    ii[0] = s.ii[0];
    ii[1] = s.ii[1];
    dim[0] = s.dim[0];
    dim[1] = s.dim[1];
    sz[0] = s.sz[0];
    sz[1] = s.sz[1];
    shallow_ = 0;
    v.ref(s.v);

    if (solver != 0)
	clearDecomposition();

    return *this;
}

LaMatrix& LaGenMatDouble::ref(SEXP x)
{				// create a reference to the data
    if (!isMatrix(x)) error("x must be a matrix");
    int *dims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    LaGenMatDouble tmp(REAL(coerceVector(x, REALSXP)), dims[0], dims[1]);
    return ref(tmp);
}

LaMatrix&  LaGenMatDouble::resize(int m, int n)
{
    if (debug())
    {
        cout << ">>> resize("<< m << "," << n << ")\n";
    }

    // first, reference 0x0 matrix, potentially freeing memory
    // this allows one to resize a matrix > 1/2 of the available
    // memory

    LaGenMatDouble tmp1(0,0);
    ref(tmp1);
    
    // now, reference an MxN matrix
    LaGenMatDouble tmp(m,n);
    ref(tmp);
    
    if (debug())
    {
	cout << "<<< resize: *this: " << this->info() << endl;
    }

    if (solver != 0)
	clearDecomposition();

    return *this;
}

LaGenMatDouble::LaGenMatDouble(const LaGenMatDouble& X) : v(0)
{
    if (X.debug())
    {
        cout << ">>> LaGenMatDouble::LaGenMatDouble(const LaGenMatDouble& X)\n";
        cout << "X: " << X.info() << endl;
	
    }
    solver = 0;
    debug_ = X.debug_;
    shallow_ = 0;		// do not perpetuate shallow copies, otherwise
				//  B = A(I,J) does not work properly...
    if (X.shallow_) {
	v.ref(X.v);
	dim[0] = X.dim[0];
	dim[1] = X.dim[1];
	sz[0] = X.sz[0];
	sz[1] = X.sz[1];
	ii[0] = X.ii[0];
	ii[1] = X.ii[1];
    } else {
	if (X.debug())
	    cout << endl << "Data is being copied!\n" << endl;
	v.resize(X.size(0)*X.size(1)); 
	ii[0](0,X.size(0)-1);
	ii[1](0,X.size(1)-1);
	dim[0] = sz[0] = X.size(0);
	dim[1] = sz[1] = X.size(1);
	int M = X.size(0), N = X.size(1);
	for (int j = 0; j < N; j++)
	    for (int i = 0; i < M; i++)
		(*this)(i,j) = X(i,j);
    }

    if (debug()) {
        cout << "*this: " << info() << endl;
        cout << "<<< LaGenMatDouble::LaGenMatDouble(const LaGenMatDouble& X)\n";
    }
}

LaGenMatDouble::LaGenMatDouble(SEXP x) : v(0)
{				// constructor performs a copy
    if (!isMatrix(x)) error("x must be a matrix");
    int *dims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    int m = dims[0];
    int n = dims[1];
    solver = 0;
    LaGenMatDouble tmp(REAL(coerceVector(x, REALSXP)), m, n);
    copy(tmp);
}

LaMatrix& LaGenMatDouble::copy(const LaMatrix& X) 
{
    if (debug())
    {
        cout << ">>> LaGenMatDouble& LaGenMatDouble::copy(const LaGenMatDouble& X)\n";
        cout << " X: " << X.info() << endl;
    }
    
				// current scheme in copy() is to
				// detach the left-hand-side 
				// from whatever it was pointing to.
    resize(X);
    
    int M = X.size(0), N = X.size(1);
    for (int i = 0; i < M; i++)
	for (int j = 0; j < N; j++)
	    (*this)(i,j) = X(i,j);
    
    if ((debug()))
    {
        cout << " *this: " << this->info() << endl;
        cout << " <<< LaGenMatDouble& LaGenMatDouble::copy(const LaGenMatDouble& s)\n";
    }
    return *this;
}

LaMatrix& LaGenMatDouble::inject(const LaMatrix& s)
{
    assert(s.size(0) == size(0));
    assert(s.size(1) == size(1));
    
    int M=size(0), N=size(1);
    for (int j = 0; j < N; j++)
        for (int i = 0; i < M; i++)
            (*this)(i,j) = s(i,j);


    if (solver != 0)
	clearDecomposition();

    return *this;
}

LaGenMatDouble LaGenMatDouble::operator()(const LaIndex& II, const LaIndex& JJ) const
{
    if (debug())
    {
	cout << ">>> LaGenMatDouble::operator(const LaIndex& const LaIndex&)\n";
    }
    LaIndex I, J;

    if (II.null()) {
        I(0,size(0)-1);
    } else {
        I = II;
    }
    if (JJ.null()) {
        J(0,size(1)-1);
    } else {
        J = JJ;
    }

    assert(I.inc() != 0);
    assert(J.inc() != 0);

    if (I.inc() > 0) {
        assert(I.start() >= 0);
        assert(I.start() <= I.end());
        assert(I.end() < size(0));
    } else {			// I.inc() < 0
        assert(I.start() < size(0));
        assert(I.start() >= I.end());
        assert(I.end() >= 0);
    }

    if (J.inc() > 0) {
        assert(J.start() >= 0);
        assert(J.start() <= J.end());
        assert(J.end() < size(1));
    } else {			// J.inc() < 0
        assert(J.start() < size(1));
        assert(J.start() >= J.end());
        assert(J.end() >= 0);
    }

    LaGenMatDouble tmp;
    
    tmp.dim[0] = dim[0];
    tmp.dim[1] = dim[1];
    tmp.sz[0] = (I.end() - I.start())/I.inc() + 1;
    tmp.sz[1] = (J.end() - J.start())/J.inc() + 1;
    
    tmp.ii[0].start() =  ii[0].start() + I.start()*ii[0].inc();
    tmp.ii[0].inc() = ii[0].inc() * I.inc();
    tmp.ii[0].end() = (I.end() - I.start())/ I.inc() * tmp.ii[0].inc() 
	+ tmp.ii[0].start();
    
    tmp.ii[1].start() =  ii[1].start() + J.start()*ii[1].inc();
    tmp.ii[1].inc() = ii[1].inc() * J.inc();
    tmp.ii[1].end() = (J.end() - J.start())/ J.inc() * tmp.ii[1].inc() 
	+ tmp.ii[1].start();
    
    tmp.v.ref(v);
    tmp.shallow_assign();
    
    if (debug())
    {
	cout << "    return value: " << tmp.info() << endl;
	cout << "<<< LaGenMatDouble::operator(const LaIndex& const LaIndex&)\n";
    }
    return tmp;
}

LaGenMatDouble LaGenMatDouble::operator()(const LaIndex& II, const LaIndex& JJ) 
{
    if (debug())
    {
	cout << ">>> LaGenMatDouble::operator(const LaIndex& const LaIndex&)\n";
    }
    LaIndex I, J;

    if (II.null()) {
        I(0,size(0)-1);
    } else {
        I = II;
    }
    if (JJ.null()) {
        J(0,size(1)-1);
    } else {
        J = JJ;
    }

    assert(I.inc() != 0);
    assert(J.inc() != 0);
    
    if (I.inc() > 0) {
        assert(I.start() >= 0);
        assert(I.start() <= I.end());
        assert(I.end() < size(0));
    } else {			// I.inc() < 0
        assert(I.start() < size(0));
        assert(I.start() >= I.end());
        assert(I.end() >= 0);
    }

    if (J.inc() > 0) {
        assert(J.start() >= 0);
        assert(J.start() <= J.end());
        assert(J.end() < size(1));
    } else {			// J.inc() < 0
        assert(J.start() < size(1));
        assert(J.start() >= J.end());
        assert(J.end() >= 0);
    }
    
    LaGenMatDouble tmp;
    
    tmp.dim[0] = dim[0];
    tmp.dim[1] = dim[1];
    tmp.sz[0] = (I.end() - I.start())/I.inc() + 1;
    tmp.sz[1] = (J.end() - J.start())/J.inc() + 1;
    
    tmp.ii[0].start() =  ii[0].start() + I.start()*ii[0].inc();
    tmp.ii[0].inc() = ii[0].inc() * I.inc();
    tmp.ii[0].end() = (I.end() - I.start())/ I.inc() * tmp.ii[0].inc() 
	+ tmp.ii[0].start();
    
    tmp.ii[1].start() =  ii[1].start() + J.start()*ii[1].inc();
    tmp.ii[1].inc() = ii[1].inc() * J.inc();
    tmp.ii[1].end() = (J.end() - J.start())/ J.inc() * tmp.ii[1].inc() 
	+ tmp.ii[1].start();
    
    tmp.v.ref(v);
    tmp.shallow_assign();
    
    if (debug())
    {
	cout << "    return value: " << tmp.info() << endl;
	cout << "<<< LaGenMatDouble::operator(const LaIndex& const LaIndex&)\n";
    }
    return tmp;
}

double LaGenMatDouble::norm(char which) const
{
    assert(which == 'M' || which == 'm' ||
	   which == '1' || which == 'O' || which == 'o' ||
	   which == 'I' || which == 'i' ||
	   which == 'F' || which == 'f' ||
	   which == 'E' || which == 'e');
    double *work = new double[size(0)]; // only for Frobenius norm
    double val = F77_CALL(dlange)(which, size(0), size(1),
				  &(*this)(0,0), gdim(0), work);
    delete[] work;
    return val;
}

inline void LaGenMatDouble::doDecomposition() const
{
    if (solver != 0)
	delete solver;
    if (size(0) == size(1))
	solver = new LaLUFactorDouble(*this);
    else solver = new LaQRFactorDouble(*this);
}

ostream& LaGenMatDouble::printMatrix(ostream& s) const
{
    if (*info_)     // print out only matrix info, not actual values
    {
        *info_ = 0; // reset the flag
        s << "(" << size(0) << "x" << size(1) << ") " ;
        s << "Indices: " << index(0) << " " << index(1);
        s << " #ref: " << ref_count();
        s << " shallow:" << shallow_  ;
    } else {
        for (int i = 0; i < size(0); i++)
        {
            for (int j = 0; j < size(1); j++) { s << (*this)(i,j) << "  "; }
            s << "\n";
        }
    }
    return s;
}

ostream& LaGenMatDouble::Info(ostream& s)
{
    LaMatDouble::Info(s);
    s << "#ref: " << ref_count() << endl;
    return s;
};

