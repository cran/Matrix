// -*- c++ -*-
//
//  Copyright (C) 2000-2000, 2002 the R Development Core Team

#ifndef _SYMMETRIC_FACTOR_H
#define _SYMMETRIC_FACTOR_H

#include "factor.h"
#include "lavi.h"
class LaSymmMatDouble;

class LaSymmFactor : public Factor {
 public:
    virtual ~LaSymmFactor() { };

    virtual LaMatDouble* solve() const = 0;
    virtual LaMatDouble& solve(LaMatDouble&) const = 0;
    virtual LaMatDouble& solve(LaMatDouble&, const LaMatDouble&) const = 0;
    virtual double rcond(double infnorm) const = 0;

    virtual const LaSymmMatDouble& decomp() const = 0;
    virtual char uplo() const = 0;
    virtual bool singular() const = 0;
    virtual const LaVectorInt& pivot() const = 0;
    virtual LaSymmFactor& ref(LaSymmMatDouble&) = 0;
    virtual LaSymmFactor& ref(const LaSymmFactor&) = 0;
};

#endif
