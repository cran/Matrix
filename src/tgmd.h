// -*- c++ -*-

#ifndef _LA_TRIANG_MAT_DOUBLE_H_
#define _LA_TRIANG_MAT_DOUBLE_H_

#include "lafnames.h"
#include LA_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_UPPER_TRIANG_MAT_DOUBLE_H

class LaTriangMatDouble : public LaMatDouble
{
    char uplo_;

    union {
        LaLowerTriangMatDouble* lower_;
        LaUpperTriangMatDouble* upper_;
    };

public:
				// constructors
    LaTriangMatDouble(char uplo = 'U')
        : uplo_(uplo),
          upper_(uplo_ == 'U'?
                 (new LaUpperTriangMatDouble()):
                 reinterpret_cast<LaUpperTriangMatDouble*>(
                     new LaLowerTriangMatDouble())) {}
    LaTriangMatDouble(int i, int j, char uplo = 'U')
        : uplo_(uplo),
          upper_(uplo_ == 'U'?
                 (new LaUpperTriangMatDouble(i, j)):
                 reinterpret_cast<LaUpperTriangMatDouble*>(
                     new LaLowerTriangMatDouble(i, j))) {}
    LaTriangMatDouble(double* d, int i, int j, char uplo = 'U')
        : uplo_(uplo),
          upper_(uplo_ == 'U'?
                 (new LaUpperTriangMatDouble(d, i, j)):
                 reinterpret_cast<LaUpperTriangMatDouble*>(
                     new LaLowerTriangMatDouble(d, i, j))) {}
    LaTriangMatDouble(const LaTriangMatDouble& A)
        : uplo_(A.uplo_),
          upper_(uplo_ == 'U'?
                 (new LaUpperTriangMatDouble(*A.upper_)):
                 reinterpret_cast<LaUpperTriangMatDouble*>(
                     new LaLowerTriangMatDouble(*A.lower_))) {}
    LaTriangMatDouble(const LaUpperTriangMatDouble& A)
        : uplo_('U'), upper_(new LaUpperTriangMatDouble())
	{
	    upper_->ref(A);
	}
    LaTriangMatDouble(const LaLowerTriangMatDouble& A)
        : uplo_('L'), lower_(new LaLowerTriangMatDouble())
	{
	    lower_->ref(A);
	}
    explicit LaTriangMatDouble(SEXP s, char uplo = 'U')
        : uplo_(uplo),
          upper_(uplo_ == 'U'?
                 (new LaUpperTriangMatDouble(s)):
                 reinterpret_cast<LaUpperTriangMatDouble*>(
                     new LaLowerTriangMatDouble(s))) {}

				// destructor
    ~LaTriangMatDouble()
	{
	    if (uplo_ == 'U')
		delete upper_;
	    else delete lower_;
	}

    char uplo() const { return uplo_; }
    int size(int d) const	// submatrix size
	{ return (uplo_ == 'U')?upper_->size(d):lower_->size(d); }
    int gdim(int d) const	// global dimensions
	{ return (uplo_ == 'U')?upper_->gdim(d):lower_->gdim(d); }
    LaIndex index(int d) const	// return indices of matrix.
	{ return (uplo_ == 'U')?upper_->index(d):lower_->index(d); }
    int ref_count() const	// return ref_count of matrix.
	{ return (uplo_ == 'U')?upper_->ref_count()
	      :lower_->ref_count(); };
    const double* addr() const	// return address of matrix.
	{ return (uplo_ == 'U')?upper_->addr():lower_->addr(); };
    double* addr()	// return address of matrix.
	{ return (uplo_ == 'U')?upper_->addr():lower_->addr(); };

				// operators
    double& operator()(int i,int j) 
	{
	    if (uplo_ == 'U')
            return (*upper_)(i, j);
	    else return (*lower_)(j, i);
	}
    double operator()(int i, int j) const
	{
	    if (uplo_ == 'U')
		return (*upper_)(i, j);
	    else return (*lower_)(j, i);
	}
    LaMatDouble& operator=(double s)
	{ 
	    if (uplo_ == 'U')
		*upper_ = s;
	    else *lower_ = s;
	    return *this;
	}

    LaTriangMatDouble& inject(const LaMatDouble& A)
	{ 
	    if (uplo_ == 'U')
            upper_->inject(A);
	    else lower_->inject(A);
	    return *this;
	}
    LaTriangMatDouble& resize(const LaMatDouble& A)
	{ return resize(A.size(0), A.size(1)); }
    LaTriangMatDouble& resize(int m, int n)
	{ 
	    if (uplo_ == 'U')
            upper_->resize(m, n);
	    else lower_->resize(m, n);
	    return *this;
	}
    LaTriangMatDouble& ref(const LaTriangMatDouble& A)
	{ 
	    if (uplo_ != A.uplo_) {
            if (uplo_ == 'U') {
                delete upper_;
                lower_ = new LaLowerTriangMatDouble();
                lower_->ref(*A.lower_);
            } else {
                delete lower_;
                upper_ = new LaUpperTriangMatDouble();
                upper_->ref(*A.upper_);
            }
            uplo_ = A.uplo_;
	    } else if (uplo_ == 'U')
            upper_->ref(*A.upper_);
	    else lower_->ref(*A.lower_);

	    return *this;
	}
    LaTriangMatDouble& ref(SEXP s)
	{
	    if (uplo_ == 'U')
            upper_->ref(s);
	    else lower_->ref(s);
	    return *this;
	}
    LaTriangMatDouble& copy(const LaMatDouble& A)
	{
	    if (uplo_ == 'U')
            upper_->copy(A);
	    else lower_->copy(A);
	    return *this;
	}
    LaTriangMatDouble* clone() const
	{
	    LaTriangMatDouble* ans = new LaTriangMatDouble(uplo());
	    if (uplo() == 'U')
            ans->upper_ = upper_->clone();
	    else ans->lower_ = lower_->clone();
	    return ans;
	}

    std::ostream& printMatrix(std::ostream& os) const
	{
	    if (uplo_ == 'U')
            upper_->printMatrix(os);
	    else lower_->printMatrix(os);
	    return os;
	}

//      operator LaGenMatDouble()
//  	{
//  	    if (uplo_ == 'U')
//  		return *data_.upper;
//  	    else return *data_.lower;
//  	}

//    operator LaLowerTriangMatDouble();
//    operator LaUpperTriangMatDouble();

				// linear equation solvers
    LaTriangMatDouble* solve() const
	{
	    LaTriangMatDouble* ans = new LaTriangMatDouble('U');
	    delete ans->upper_;
	    if (uplo() == 'U') {
            ans->uplo_ = 'U';
            ans->upper_ = upper_->solve();
	    } else {
            ans->uplo_ = 'L';
            ans->lower_ = lower_->solve();
	    }
	    return ans;
	}
    LaMatDouble& solve(LaMatDouble& B) const
	{
	    if (uplo() == 'U') {
            return upper_->solve(B);
	    } else {
            return lower_->solve(B);
	    }
	}
    LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const
	{
	    if (uplo() == 'U') {
            return upper_->solve(X, B);
	    } else {
            return lower_->solve(X, B);
	    }
	}
    
				// matrix norms, etc.
    double norm(char which) const
	{
	    if (uplo() == 'U') {
            return upper_->norm(which);
	    } else {
            return lower_->norm(which);
	    }
	}
    double rcond(char which) const
	{
	    if (uplo() == 'U') {
            return upper_->rcond(which);
	    } else {
            return lower_->rcond(which);
	    }
	}
    SEXP asSEXP() const
	{
	    if (uplo() == 'U') {
            return upper_->asSEXP();
	    } else {
            return lower_->asSEXP();
	    }
	}
};

#endif
