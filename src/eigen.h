// -*- c++ -*-
//
//  Copyright (C) 2000-2000, 2002 the R Development Core Team

#ifndef _LA_EIGEN_H_
#define _LA_EIGEN_H_

#include "lamatrix.h"
#include LA_VECTOR_DOUBLE_H
#include LA_VECTOR_INT_H
#include "orthonormal.h"

class LaEigen			// eigenvalue/eigenvector decompositions
{
public:
    virtual ~LaEigen() {}
    virtual SEXP asSEXP() const = 0; // copy the decomposition to an SEXP
    virtual LaMatrix& values() = 0; // will always be a LaVectorDouble
    virtual LaMatrix& vectors() = 0;
    virtual LaMatrix& vectors(char side = 'L') = 0;
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
    LaSymmEigenDouble(const LaMatDouble& a, const char uplo,
                      const bool findVecs = true);
				// accessor methods
    LaMatrix& values() { return vals; }
    LaMatrix& vectors() { return vecs; }
    LaMatrix& vectors(char side = 'L') { return vecs; }

    SEXP asSEXP() const;
};

class LaGenEigenDouble : public LaEigenDouble
{
protected:
    LaVectorDouble wR;
    LaVectorDouble wI;
    LaVectorDouble scale;
    LaVectorDouble rcondE;
    LaVectorDouble rcondV;
    LaGenMatDouble left;
    LaGenMatDouble right;
    double abnrm;
    int ilo, ihi;
    bool complexVectors_;
public:
    LaGenEigenDouble(const LaMatDouble& a, bool leftEV = true,
		     bool rightEV = true, char balanc = 'B', char sense = 'N');
				// accessor methods
    bool complexVectors() const { return complexVectors_; }
    LaMatrix& valuesR() { return wR; }
    LaMatrix& valuesI() { return wI; }
    LaMatrix& values()
	{
	    if (complexVectors())
		throw(LaException("Can not return complex values"));
	    return wR;
	}
    LaMatrix& vectors() { if (left.size(0) > 0) return left; return right; }
    LaMatrix& vectors(char side = 'L') {
	if (side == 'L') return left;
	return right;
    }

    SEXP asSEXP() const;
};
    
#endif // _LA_EIGEN_H_
