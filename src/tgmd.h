// -*- c++ -*-

#ifndef _LA_TRIANG_MAT_DOUBLE_H_
#define _LA_TRIANG_MAT_DOUBLE_H_

#include "lafnames.h"
#include LA_LOWER_TRIANG_MAT_DOUBLE_H
#include LA_UPPER_TRIANG_MAT_DOUBLE_H

class LaTriangMatDouble : public LaMatDouble
{
    char uplo_;

    union LowerUpper {
	LaLowerTriangMatDouble *lower;
	LaUpperTriangMatDouble *upper;
    } data_;

public:
				// constructors
    LaTriangMatDouble(char uplo = 'U')
	: uplo_(uplo)
	{
	    if (uplo_ == 'U')
		data_.upper = new LaUpperTriangMatDouble();
	    else data_.lower = new LaLowerTriangMatDouble();
	}
    LaTriangMatDouble(int i, int j, char uplo = 'U')
	: uplo_(uplo)
	{
	    if (uplo_ == 'U')
		data_.upper = new LaUpperTriangMatDouble(i, j);
	    else data_.lower = new LaLowerTriangMatDouble(i, j);
	}
    LaTriangMatDouble(double* d, int i, int j, char uplo = 'U')
	: uplo_(uplo)
	{
	    if (uplo_ == 'U')
		data_.upper = new LaUpperTriangMatDouble(d, i, j);
	    else data_.lower = new LaLowerTriangMatDouble(d, i, j);
	}
    LaTriangMatDouble(const LaTriangMatDouble& A)
	: uplo_(A.uplo_)
	{
	    if (uplo_ == 'U')
		data_.upper = new LaUpperTriangMatDouble(*A.data_.upper);
	    else data_.lower = new LaLowerTriangMatDouble(*A.data_.lower);
	}
    LaTriangMatDouble(const LaUpperTriangMatDouble& A)
	: uplo_('U')
	{
	    data_.upper = new LaUpperTriangMatDouble();
	    data_.upper->ref(A);
	}
    LaTriangMatDouble(const LaLowerTriangMatDouble& A)
	: uplo_('L')
	{
	    data_.lower = new LaLowerTriangMatDouble();
	    data_.lower->ref(A);
	}
    explicit LaTriangMatDouble(SEXP s, char uplo = 'U')
	: uplo_(uplo)
	{
	    if (uplo_ == 'U')
		data_.upper = new LaUpperTriangMatDouble(s);
	    else data_.lower = new LaLowerTriangMatDouble(s);
	}

				// destructor
    ~LaTriangMatDouble()
	{
	    if (uplo_ == 'U')
		delete data_.upper;
	    else delete data_.lower;
	}

    char uplo() const { return uplo_; }
    int size(int d) const	// submatrix size
	{ return (uplo_ == 'U')?data_.upper->size(d):data_.lower->size(d); }
    int gdim(int d) const	// global dimensions
	{ return (uplo_ == 'U')?data_.upper->gdim(d):data_.lower->gdim(d); }
    LaIndex index(int d) const	// return indices of matrix.
	{ return (uplo_ == 'U')?data_.upper->index(d):data_.lower->index(d); }
    int ref_count() const	// return ref_count of matrix.
	{ return (uplo_ == 'U')?data_.upper->ref_count()
	      :data_.lower->ref_count(); };
    double* addr() const	// return address of matrix.
	{ return (uplo_ == 'U')?data_.upper->addr():data_.lower->addr(); };

				// operators
    double& operator()(int i,int j) 
	{
	    if (uplo_ == 'U')
		return (*data_.upper)(i, j);
	    else return (*data_.lower)(j, i);
	}
    const double& operator()(int i, int j) const
	{
	    if (uplo_ == 'U')
		return (*data_.upper)(i, j);
	    else return (*data_.lower)(j, i);
	}
    LaMatDouble& operator=(double s)
	{ 
	    if (uplo_ == 'U')
		*data_.upper = s;
	    else *data_.lower = s;
	    return *this;
	}

    LaTriangMatDouble& inject(const LaMatDouble& A)
	{ 
	    if (uplo_ == 'U')
		data_.upper->inject(A);
	    else data_.lower->inject(A);
	    return *this;
	}
    LaTriangMatDouble& resize(const LaMatDouble& A)
	{ return resize(A.size(0), A.size(1)); }
    LaTriangMatDouble& resize(int m, int n)
	{ 
	    if (uplo_ == 'U')
		data_.lower->resize(m, n);
	    else data_.lower->resize(m, n);
	    return *this;
	}
    LaTriangMatDouble& ref(const LaTriangMatDouble& A)
	{ 
	    if (uplo_ != A.uplo_) {
		if (uplo_ == 'U') {
		    delete data_.upper;
		    data_.lower = new LaLowerTriangMatDouble();
		    data_.lower->ref(*A.data_.lower);
		} else {
		    delete data_.lower;
		    data_.upper = new LaUpperTriangMatDouble();
		    data_.upper->ref(*A.data_.upper);
		}
		uplo_ = A.uplo_;
	    } else if (uplo_ == 'U')
		data_.upper->ref(*A.data_.upper);
	    else data_.lower->ref(*A.data_.lower);

	    return *this;
	}
    LaTriangMatDouble& ref(SEXP s)
	{
	    if (uplo_ == 'U')
		data_.upper->ref(s);
	    else data_.lower->ref(s);
	    return *this;
	}
    LaTriangMatDouble& copy(const LaMatDouble& A)
	{
	    if (uplo_ == 'U')
		data_.upper->copy(A);
	    else data_.lower->copy(A);
	    return *this;
	}
    LaTriangMatDouble* clone() const
	{
	    LaTriangMatDouble* ans = new LaTriangMatDouble(uplo());
	    if (uplo() == 'U')
		ans->data_.upper = data_.upper->clone();
	    else ans->data_.lower = data_.lower->clone();
	    return ans;
	}

    ostream &printMatrix(ostream& os) const
	{
	    if (uplo_ == 'U')
		data_.upper->printMatrix(os);
	    else data_.lower->printMatrix(os);
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
	    delete ans->data_.upper;
	    if (uplo() == 'U') {
		ans->uplo_ = 'U';
		ans->data_.upper = data_.upper->solve();
	    } else {
		ans->uplo_ = 'L';
		ans->data_.lower = data_.lower->solve();
	    }
	    return ans;
	}
    LaMatDouble& solve(LaMatDouble& B) const
	{
	    if (uplo() == 'U') {
		return data_.upper->solve(B);
	    } else {
		return data_.lower->solve(B);
	    }
	}
    LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const
	{
	    if (uplo() == 'U') {
		return data_.upper->solve(X, B);
	    } else {
		return data_.lower->solve(X, B);
	    }
	}
    
				// matrix norms, etc.
    double norm(char which) const
	{
	    if (uplo() == 'U') {
		return data_.upper->norm(which);
	    } else {
		return data_.lower->norm(which);
	    }
	}
    double rcond(char which) const
	{
	    if (uplo() == 'U') {
		return data_.upper->rcond(which);
	    } else {
		return data_.lower->rcond(which);
	    }
	}
    SEXP asSEXP() const
	{
	    if (uplo() == 'U') {
		return data_.upper->asSEXP();
	    } else {
		return data_.lower->asSEXP();
	    }
	}
};

#endif
