// -*- c++ -*-
//      LAPACK++ (V. 1.0a Beta)
//      (C) 1992-1996 All Rights Reserved.
//
//  Modifications Copyright (C) 2000-2000 the R Development Core Team
//
//      Lapack++ Rectangular Matrix Class
//
//      Dense (nonsingular) matrix, assumes no special structure or properties.
//
//      ) allows 2-d indexing
//      ) non-unit strides
//      ) inject assignment
//      ) cout << A.info()  prints out internal states of A
//      ) indexing via A(i,j) where i,j are either integers or
//              LaIndex         

#ifndef _LA_GEN_MAT_DOUBLE_H_
#define _LA_GEN_MAT_DOUBLE_H_

#include "lamatrix.h"
#include "lapack.h"
#include VECTOR_DOUBLE_H
#include "solvable.h"

class LaGenMatDouble : public LaMatDouble
{
private:
    VectorDouble    v;
    LaIndex         ii[2];
    int             dim[2];	// size of original matrix, not submatrix
    int             sz[2];	// size of this submatrix
    mutable solvable* solver;   // LU or QR decomposition

public:
    // constructors and destructors
    LaGenMatDouble();
    LaGenMatDouble(int, int);
    LaGenMatDouble(double*, int, int);
    LaGenMatDouble(const LaGenMatDouble&);
    LaGenMatDouble(SEXP);
    virtual ~LaGenMatDouble();

    //  Indices and access operations 
    int size(int d) const       // submatrix size
	{ return sz[d]; }
    int gdim(int d) const       // global dimensions
	{ return dim[d]; }
    LaIndex index(int d) const
	{ return ii[d]; }
    int ref_count() const
	{ return v.ref_count(); }
    double* addr() const        // begining addr of data space
	{ return v.addr(); }
    inline double& operator()(int i, int j)
	{
#ifdef LA_BOUNDS_CHECK
	    assert(i>=0);
	    assert(i<size(0));
	    assert(j>=0);
	    assert(j<size(1));
#endif
	    return v( dim[0]*(ii[1].start() + j*ii[1].inc()) + 
		      ii[0].start() + i*ii[0].inc());
	}
    inline double& operator()(int i, int j) const
	{
#ifdef LA_BOUNDS_CHECK
	    assert(i>=0);
	    assert(i<size(0));
	    assert(j>=0);
	    assert(j<size(1));
#endif
	    return v( dim[0]*(ii[1].start() + j*ii[1].inc()) + 
		      ii[0].start() + i*ii[0].inc());
	}
    LaGenMatDouble operator()(const LaIndex& I, const LaIndex& J) ;
    LaGenMatDouble operator()(const LaIndex& I, const LaIndex& J) const;

    LaMatDouble& operator=(double s);

    LaMatrix& resize(const LaMatrix& m)
	{ return resize(m.size(0), m.size(1)); }
    LaMatrix& resize(int m, int n);
    LaMatrix& ref(const LaGenMatDouble& s);
    LaMatrix& ref(SEXP);
    LaMatrix& inject(const LaMatrix& s);
    LaMatrix& copy(const LaMatrix& s);

    double norm(char which) const;
    inline void doDecomposition() const;
    void clearDecomposition() const
	{ delete solver; solver = 0; }
    LaMatrix& solve() const
	{ if (solver == 0) doDecomposition(); return solver->solve(); };
    LaMatrix& solve(LaMatrix& B) const
	{ if (solver == 0) doDecomposition(); return solver->solve(B); };
    LaMatrix& solve(LaMatrix& X, const LaMatrix& B) const
	{ if (solver == 0) doDecomposition(); return solver->solve(X, B); };

    //* I/O *//
    ostream& printMatrix(ostream&) const;
    ostream& Info(ostream& s);
};				// End of LaGenMatDouble Class

#endif 

// _LA_GEN_MAT_H_
