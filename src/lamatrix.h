// -*- c++ -*-
//
//  Copyright (C) 2000-2001 the R Development Core Team

#ifndef _LA_MATRIX_H_
#define _LA_MATRIX_H_

inline int min(int a, int b) { return (a > b) ? b : a; }
inline int max(int a, int b) { return (a > b) ? a : b; }

#include "lafnames.h"
#include LA_INDEX_H
#include LA_EXCEPTION_H
class LaEigenDouble;

extern "C" {
#include <R.h>
#include <Rinternals.h>
}

#include <iostream.h>

class LaMatrix
{
protected:
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
    
    virtual LaMatrix& resize(int m, int n) = 0;
    virtual LaMatrix* clone() const = 0;

    virtual double norm(char which) const = 0;
    virtual double rcond(char which) const = 0;
    virtual void doDecomposition() const { };
    virtual SEXP asSEXP() const = 0; // copy the matrix to an SEXP

    virtual int shallow() const	// read global shallow flag
        { return shallow_;}
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
    
    virtual double& operator()(int i, int j) = 0;
    virtual const double& operator()(int i, int j) const = 0;
    virtual LaMatDouble& operator=(const LaMatDouble& s)
	{ return inject(s); }
    virtual LaMatDouble& operator=(double s) = 0;

    LaMatDouble& resize(int m, int n) = 0;
    virtual LaMatDouble& resize(const LaMatDouble &A) = 0;

    virtual LaMatDouble& ref(SEXP) = 0;
    virtual LaMatDouble& inject(const LaMatDouble& s) = 0;
    virtual LaMatDouble& copy(const LaMatDouble& s) = 0;
    virtual LaMatDouble* clone() const = 0;

    virtual LaMatDouble* solve() const = 0;
    virtual LaMatDouble& solve(LaMatDouble& B) const = 0;
    virtual LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const = 0;

    virtual LaEigenDouble* eigen(bool leftEV = true, bool rightEV = true,
				 char balanc = 'B', char rcond = 'N');

    ostream& Info(ostream& s);
};

class LaMatInt : public LaMatrix
{
public:
    virtual ~LaMatInt() { }
				//  Indices and access operations 
    virtual int* addr() const = 0;// begining addr of data space
    
    virtual int& operator()(int i, int j) = 0;
    virtual const int& operator()(int i, int j) const = 0;
    virtual LaMatInt& operator=(const LaMatInt& s)
	{ return inject(s); }
    virtual LaMatInt& operator=(int s) = 0;
 
    LaMatInt& resize(int m, int n) = 0;
    virtual LaMatInt& resize(const LaMatInt &A) = 0;

    virtual LaMatInt& ref(SEXP) = 0;
    virtual LaMatInt& inject(const LaMatInt& s) = 0;
    virtual LaMatInt& copy(const LaMatInt& s) = 0;

    ostream& Info(ostream& s);
};

#endif
