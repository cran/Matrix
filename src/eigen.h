// -*- c++ -*-
//
//  Copyright (C) 2000-2000 the R Development Core Team

#ifndef _LA_EIGEN_H
#define _LA_EIGEN_H

#include "lamatrix.h"
#include LA_VECTOR_DOUBLE_H
#include "orthonormal.h"

class LaEigen			// eigenvalue/eigenvector decompositions
{
public:
    virtual ~LaEigen() {}
    virtual SEXP asSEXP() const = 0; // copy the decomposition to an SEXP
    virtual LaMatrix& values() const = 0; // will always be a LaVectorDouble
    virtual LaMatrix& vectors() const = 0;
    virtual LaMatrix& vectors(char side = 'L') const = 0;
};

class LaEigenDouble : public LaEigen
{
public:
    virtual ~LaEigenDouble() {}
};

class LaSymmEigenDouble : public LaEigenDouble
{
protected:
    LaVectorDouble vals;
    LaOrthogonalMatDouble vecs;
public:
    LaSymmEigenDouble(const LaMatDouble& a, bool findVecs = true);
				// accessor methods
    LaMatrix& values() { return vals; }
    LaMatrix& vectors() { return vecs; }
    LaMatrix& vectors(char side = 'L') { return vecs; }

    SEXP asSEXP() const;
};

#endif // _LA_EIGEN_H
