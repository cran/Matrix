#ifndef _SOLVABLE_H
#define _SOLVABLE_H

class solvable {
 public:
    virtual LaMatrix& solve() const = 0;
    virtual LaMatrix& solve(LaMatrix&) const = 0;
    virtual LaMatrix& solve(LaMatrix&, const LaMatrix&) const = 0;
};

#endif
