//
//  Copyright (C) 2000-2000 the R Development Core Team
//

#ifndef _LA_MATRIX_H_
#define _LA_MATRIX_H_

inline int min(int a, int b) { return (a > b) ? b : a; }
inline int max(int a, int b) { return (a > b) ? a : b; }

#include "lafnames.h"
#include LA_INDEX_H
#include LA_EXCEPTION_H

extern "C" {
#include <R.h>
#include <Rinternals.h>
}

#include <iostream.h>

class LaMatrix
{
protected:
    static int      debug_;	// trace all entry and exits into methods and 
				// operators of this class.  This variable is
				// explicitly initalized in gmd.cc

    static int      *info_;	// print matrix info only, not values
				// originally 0, set to 1, and then
				// reset to 0 after use.
				// use as in
				//
				//    cout << B.info() << endl;
				//
		// this *info_ member is unique in that it really isn't
                // part of the matrix info, just a flag as to how
                // to print it.   We've included in this beta release
                // as part of our testing, but we do not expect it 
                // to be user accessable.
                // It has to be declared as global static
                // so that we may monitor expresssions like
                // X::(const &X) and still utilize without violating
                // the "const" condition.
                // Because this *info_ is used at most one at a time,
                // there is no harm in keeping only one copy of it,
                // also, we do not need to malloc free space every time
                // we call a matrix constructor.

    
    int             shallow_;	// set flag to '0' in order to return matrices
				// by value from functions without unecessary
				// copying.
    inline LaMatrix& shallow_assign()
	{ shallow_ = 1; return *this; }

public:
    virtual ~LaMatrix() {}
				//  Indices and access operations 
    virtual int size(int d) const = 0; // submatrix size
    virtual int inc(int d) const       // explicit increment
	{ return index(d).inc(); }
    virtual int gdim(int d) const = 0; // global dimensions
    virtual int start(int d) const
	{ return index(d).start(); }
    virtual int end(int d) const
	{ return index(d).end(); }
    virtual LaIndex index(int d) const = 0; // index
    
    virtual double& operator()(int i, int j) = 0;
    virtual double& operator()(int i, int j) const = 0;
    virtual LaMatrix& operator=(const LaMatrix& s)
	{ return inject(s); }

    virtual LaMatrix& resize(int m, int n) = 0;
    virtual LaMatrix& resize(const LaMatrix &A) = 0;

    virtual LaMatrix& ref(SEXP) = 0;
    virtual LaMatrix& inject(const LaMatrix& s) = 0;
    virtual LaMatrix& copy(const LaMatrix& s) = 0;
    
    virtual double norm(char which) const = 0;
    virtual void doDecomposition(){};
    virtual LaMatrix& solve() const = 0;
    virtual LaMatrix& solve(LaMatrix& B) const = 0;
    virtual LaMatrix& solve(LaMatrix& X, const LaMatrix& B) const = 0;

    virtual int shallow() const          // read global shallow flag
        { return shallow_;}
    virtual int debug() const            // read global debug flag
	{ return debug_; }
    virtual int debug(int d)             // set global debug flag
	{ return debug_ = d; }
    virtual const LaMatrix& info() const
	{ int *t = info_; *t = 1; return *this; }
    
    //* I/O *//
    virtual ostream& printMatrix(ostream&) const = 0;
    friend ostream& operator<<(ostream& s, const LaMatrix& mat)
	{ return mat.printMatrix(s); }
    virtual ostream& Info(ostream& s);
};

class LaMatDouble : public LaMatrix
{
public:
    virtual ~LaMatDouble() { }
				//  Indices and access operations 
    virtual double* addr() const = 0;// begining addr of data space
    
    virtual LaMatDouble& operator=(double s) = 0;

    ostream& Info(ostream& s);
};

class LaMatInt : public LaMatrix
{
public:
    virtual ~LaMatInt() { }
				//  Indices and access operations 
    virtual int* addr() const = 0;// begining addr of data space
    
    virtual LaMatInt& operator=(int s) = 0;

    ostream& Info(ostream& s);
};

#endif
