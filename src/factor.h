// -*- c++ -*-
//
//  Copyright (C) 2000-2000 the R Development Core Team

#ifndef _FACTOR_H
#define _FACTOR_H

class Factor {
 public:
    virtual LaMatDouble& solve() const = 0;
    virtual LaMatDouble& solve(LaMatDouble&) const = 0;
    virtual LaMatDouble& solve(LaMatDouble&, const LaMatDouble&) const = 0;
};

#endif
