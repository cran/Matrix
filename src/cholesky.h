// -*- c++ -*-

#ifndef _LA_CHOLESKY_H_
#define _LA_CHOLESKY_H_

#include "lafnames.h"
#include "symm_factor.h"
#include "tgmd.h"
#include "symd.h"

class LaSpdMatDouble;

class LaCholeskyFactorDouble : public LaSymmFactor
{
    LaSpdMatDouble* decomp_;

public:
				// constructor
    inline LaCholeskyFactorDouble();
    inline explicit LaCholeskyFactorDouble(const LaSpdMatDouble&);
    inline LaCholeskyFactorDouble(const LaCholeskyFactorDouble& F);
    inline LaCholeskyFactorDouble& operator=(const LaCholeskyFactorDouble& F);

    inline ~LaCholeskyFactorDouble();

				// extractor methods for components
    inline const LaSymmMatDouble& decomp() const;
    inline const LaSpdMatDouble& localDecomp() const;
    inline char uplo() const;
    bool singular() const
	{ throw(LaException("singularity of decomposition is unknown")); };
    const LaVectorInt& pivot() const
	{ throw(LaException("pivoting not used in LaCholeskyFactorization")); };

    inline double rcond(double infnorm) const;

				// linear equation solvers
    inline LaMatDouble* solve() const;// inverse
    inline LaMatDouble& solve(LaMatDouble& B) const; // in-place solution
    inline LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const;

				// operators
    inline LaCholeskyFactorDouble& ref(const LaCholeskyFactorDouble& F);
    LaCholeskyFactorDouble& ref(const LaSymmFactor& F)
    {
        return ref(dynamic_cast<const LaCholeskyFactorDouble&>(F));
    }
    inline LaCholeskyFactorDouble& ref(LaSpdMatDouble&);
    inline LaCholeskyFactorDouble& ref(LaSymmMatDouble& A);
};

#include "spdmd.h"

// constructor/destructor functions

inline LaCholeskyFactorDouble::LaCholeskyFactorDouble()
    : decomp_(0)
{
}

inline LaCholeskyFactorDouble::LaCholeskyFactorDouble(const LaCholeskyFactorDouble& F)
    : decomp_(new LaSpdMatDouble)
{
    ref(F);
}

inline LaCholeskyFactorDouble&
LaCholeskyFactorDouble::operator=(const LaCholeskyFactorDouble& F)
{
    ref(F);
    return *this;
}

inline LaCholeskyFactorDouble::LaCholeskyFactorDouble(const LaSpdMatDouble& A)
    : decomp_(new LaSpdMatDouble)
{
    LaSpdMatDouble A1;
    A1.copy(A);
    ref(A1);
}

inline LaCholeskyFactorDouble::~LaCholeskyFactorDouble()
{
    delete decomp_;
}

inline const LaSymmMatDouble& LaCholeskyFactorDouble::decomp() const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));
    return *decomp_;
}

inline const LaSpdMatDouble& LaCholeskyFactorDouble::localDecomp() const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));
    return *decomp_;
}

inline char LaCholeskyFactorDouble::uplo() const
{
    if (decomp_ == 0)
        throw(LaException("No decomposition present"));
    return decomp_->uplo();
}

// operators
inline LaCholeskyFactorDouble& LaCholeskyFactorDouble::ref(const LaCholeskyFactorDouble& F)
{
    if (decomp_ == 0)
        decomp_ = new LaSpdMatDouble();
    decomp_->ref(F.localDecomp());
    return *this;
}

inline LaCholeskyFactorDouble& LaCholeskyFactorDouble::ref(LaSpdMatDouble& A)
{
    if(A.inc(0) != 1 || A.inc(1) != 1)	
	throw(LaException("LaCholeskyFactorDouble::ref(const LaSpdMatDouble&)",
			  "input matrix has non unit increment"));
    if (decomp_ == 0)
	decomp_ = new LaSpdMatDouble();
    decomp_->ref(A);
    int info;
    F77_CALL(dpotrf)(uplo(), A.size(0), &A(0, 0), A.gdim(0), info);
    if (info < 0)
	throw(LaException("LaCholeskyFactorDouble::ref(const LaSpdMatDouble&)",
			  "illegal input"));
    if (info > 0)
	throw(LaException("LaCholeskyFactorDouble::ref(const LaSpdMatDouble&)",
			  "non positive definite matrix"));
    return *this;
}

inline LaCholeskyFactorDouble& LaCholeskyFactorDouble::ref(LaSymmMatDouble& A)
{
    return ref(dynamic_cast<LaSpdMatDouble&>(A));
}

inline double LaCholeskyFactorDouble::rcond(double infnorm) const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));
    double ans;
    VectorDouble work(3*localDecomp().size(0));
    VectorInt iwork(localDecomp().size(0));
    int info;
    F77_CALL(dpocon)(uplo(), localDecomp().size(0), localDecomp().addr(),
                     localDecomp().gdim(0), infnorm, ans, work.addr(),
                     iwork.addr(), info);
    if (info < 0)
	throw(LaException("LaCholeskyFactorDouble::rcond(double)",
			  "illegal input"));
    return ans;
}

inline LaMatDouble* LaCholeskyFactorDouble::solve() const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));
    int info;
    LaSpdMatDouble* ans = localDecomp().clone();
    F77_CALL(dpotri)(uplo(), localDecomp().size(0), &(*ans)(0,0),
		     localDecomp().gdim(0), info);
    if (info < 0)
	throw(LaException("LaCholeskyFactorDouble::solve(const LaSpdMatDouble&)",
			  "illegal input"));
    if (info > 0)
	throw(LaException("LaCholeskyFactorDouble::solve(const LaSpdMatDouble&)",
			  "non positive definite matrix"));
    return ans;
}

inline LaMatDouble& LaCholeskyFactorDouble::solve(LaMatDouble& B) const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));
    int info;
    F77_CALL(dpotrs)(uplo(), localDecomp().size(0), B.size(1),
                     localDecomp().addr(), localDecomp().gdim(0),
                     B.addr(), B.gdim(0), info);
    return B;
}

inline LaMatDouble& LaCholeskyFactorDouble::solve(LaMatDouble& X,
                                                  const LaMatDouble& B ) const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));
    X.inject(B);
    return solve(X);
}

#endif
