// -*- c++ -*-
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
// Modifications Copyright (C) 2000-2002 the R Development Core Team

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
#include "factor.h"

class LaGenMatDouble : public LaMatDouble
{
private:
    VectorDouble    v;
    LaIndex         ii[2];
    int             dim[2];	// size of original matrix, not submatrix
    int             sz[2];	// size of this submatrix
    mutable Factor* solver;	// LU or QR decomposition

public:
				// constructors and destructors
    LaGenMatDouble();
    LaGenMatDouble(int, int);
    LaGenMatDouble(double*, int, int);
    explicit LaGenMatDouble(const LaMatDouble&);
    LaGenMatDouble(const LaGenMatDouble&);
    LaGenMatDouble& operator=(const LaGenMatDouble&);
    explicit LaGenMatDouble(SEXP);
    ~LaGenMatDouble();

				//  Indices and access operations 
    int size(int d) const       // submatrix size
	{ return sz[d]; }
    int gdim(int d) const       // global dimensions
	{ return dim[d]; }
    LaIndex index(int d) const
	{ return ii[d]; }
    int ref_count() const
	{ return v.ref_count(); }
    double* addr()                    // begining addr of data space
	{ return v.addr() + dim[0]*ii[1].start() + ii[0].start(); }
    const double* addr() const        // begining addr of data space
	{ return v.addr() + dim[0]*ii[1].start() + ii[0].start(); }
    inline double& operator()(int i, int j)
	{
#ifdef LA_BOUNDS_CHECK
	    if (!(i>=0)) throw(LaException("assert failed : i>=0"));
	    if (!(i<size(0))) throw(LaException("assert failed : i<size(0)"));
	    if (!(j>=0)) throw(LaException("assert failed : j>=0"));
	    if (!(j<size(1))) throw(LaException("assert failed : j<size(1)"));
#endif
	    return v( dim[0]*(ii[1].start() + j*ii[1].inc()) + 
		      ii[0].start() + i*ii[0].inc());
	}
    inline double operator()(int i, int j) const
	{
#ifdef LA_BOUNDS_CHECK
	    if (!(i>=0)) throw(LaException("assert failed : i>=0"));
	    if (!(i<size(0))) throw(LaException("assert failed : i<size(0)"));
	    if (!(j>=0)) throw(LaException("assert failed : j>=0"));
	    if (!(j<size(1))) throw(LaException("assert failed : j<size(1)"));
#endif
	    return v( dim[0]*(ii[1].start() + j*ii[1].inc()) + 
		      ii[0].start() + i*ii[0].inc());
	}
    LaGenMatDouble operator()(const LaIndex& I, const LaIndex& J) ;
    LaGenMatDouble operator()(const LaIndex& I, const LaIndex& J) const;

    LaMatDouble& operator=(double s);

    LaGenMatDouble& resize(const LaMatDouble& m)
	{ return resize(m.size(0), m.size(1)); }
    LaGenMatDouble& resize(int m, int n);
    LaGenMatDouble& ref(const LaGenMatDouble& s);
    LaGenMatDouble& ref(SEXP);
    LaGenMatDouble& inject(const LaMatDouble& s);
    LaGenMatDouble& copy(const LaMatDouble& s);
    LaGenMatDouble* clone() const;
    SEXP asSEXP() const;

    double norm(char which) const;
    double rcond(char which) const;
    void doDecomposition() const;
    Factor& clearDecomposition() const
	{ delete solver; solver = 0; return *solver; }
    LaGenMatDouble* solve() const
	{ if (solver == 0) doDecomposition();
	  return dynamic_cast<LaGenMatDouble*>(solver->solve()); };
    LaMatDouble& solve(LaMatDouble& B) const
	{ if (solver == 0) doDecomposition(); return solver->solve(B); };
    LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const
	{ if (solver == 0) doDecomposition(); return solver->solve(X, B); };


    //* I/O *//
    std::ostream& printMatrix(std::ostream&) const;
    std::ostream& Info(std::ostream& s);
};				// End of LaGenMatDouble Class

#endif 

// _LA_GEN_MAT_H_
