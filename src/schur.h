// -*- c++ -*-
//
//  Copyright (C) 2000-2000 the R Development Core Team

#ifndef _LA_SCHUR_H_
#define _LA_SCHUR_H_

#include "lamatrix.h"
#include LA_VECTOR_DOUBLE_H
#include LA_VECTOR_INT_H
#include "orthonormal.h"

class LaSchur			// eigenvalue/eigenvector decompositions
{
public:
    virtual ~LaSchur() {}
    virtual SEXP asSEXP() const = 0; // copy the decomposition to an SEXP
    virtual LaMatrix& values() = 0; // will always be a LaVectorDouble
    virtual LaMatrix& vectors() = 0;
};

class LaSchurDouble : public LaSchur
{
public:
    virtual ~LaSchurDouble() {}
};

class LaGenSchurDouble : public LaSchurDouble
{
protected:
    LaGenMatDouble a;
    LaVectorDouble wR;
    LaVectorDouble wI;
    LaGenMatDouble vecs;
    bool complexVectors_;
public:
    LaGenSchurDouble(LaMatDouble& a, bool jobV = true);
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
    LaMatrix& vectors() { return vecs; }

    SEXP asSEXP() const;
};
    
#endif // _LA_SCHUR_H_
