// -*- c++ -*-

#ifndef _LA_BUNCH_KAUFMAN_H_
#define _LA_BUNCH_KAUFMAN_H_

#include "lafnames.h"
#include "symm_factor.h"
#include "tgmd.h"
#include "lavi.h"
#include "lautil.h"

class LaSymmMatDouble;

class LaBunchKaufmanFactorDouble : public LaSymmFactor
{
    LaSymmMatDouble* decomp_;
    bool           uplo_;
    LaVectorInt    pivot_;
    bool           singular_;

public:
				// constructor
    LaBunchKaufmanFactorDouble()
	{ decomp_ = 0; };
    inline explicit LaBunchKaufmanFactorDouble(const LaSymmMatDouble&);
    LaBunchKaufmanFactorDouble(const LaBunchKaufmanFactorDouble& F);

    inline ~LaBunchKaufmanFactorDouble();

				// extractor methods for components
    const LaSymmMatDouble& decomp() const
	{
	    if (decomp_ == 0)
		throw(LaException("No decomposition present"));
	    return *decomp_;
	}
    const LaVectorInt& pivot() const
	{ return pivot_; };
    bool singular() const
	{ return singular_; };
    inline char uplo() const;

    inline double rcond(double infnorm) const;

				// linear equation solvers
    inline LaMatDouble* solve() const;// inverse
    inline LaMatDouble& solve(LaMatDouble& B) const; // in-place solution
    inline LaMatDouble& solve(LaMatDouble& X, const LaMatDouble& B) const;

				// operators
    inline LaBunchKaufmanFactorDouble& ref(const LaBunchKaufmanFactorDouble& F);
    inline LaBunchKaufmanFactorDouble& ref(LaSymmMatDouble&);
};

#include "symd.h"

inline LaBunchKaufmanFactorDouble::LaBunchKaufmanFactorDouble(const LaSymmMatDouble& A)
{
    decomp_ = new LaSymmMatDouble;
    LaSymmMatDouble A1;
    A1.copy(A);
    ref(A1);
}

inline LaBunchKaufmanFactorDouble::LaBunchKaufmanFactorDouble(const LaBunchKaufmanFactorDouble& F)
{
    decomp_ = new LaSymmMatDouble;
    ref(F);
}

inline LaBunchKaufmanFactorDouble::~LaBunchKaufmanFactorDouble()
{
    delete decomp_;
}

inline char LaBunchKaufmanFactorDouble::uplo() const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));
    return decomp_->uplo();
}

// operators
inline LaBunchKaufmanFactorDouble& LaBunchKaufmanFactorDouble::ref(const LaBunchKaufmanFactorDouble& F)
{
    if (decomp_ == 0)
	decomp_ = new LaSymmMatDouble();
    decomp_->ref(F.decomp());
    pivot_.ref(F.pivot());
    singular_ = F.singular();
    return *this;
}

inline LaBunchKaufmanFactorDouble& LaBunchKaufmanFactorDouble::ref(LaSymmMatDouble& A)
{
    if(A.inc(0) != 1 || A.inc(1) != 1)	
	throw(LaException("LaBunchKaufmanFactorDouble::ref(const LaSymmMatDouble&)",
			  "input matrix has non unit increment"));
    if (decomp_ == 0)
	decomp_ = new LaSymmMatDouble();
    decomp_->ref(A);
    pivot_.resize(A.size(0));
    int info;
    char uplo_str[2];
    uplo_str[0] = uplo();
    uplo_str[1] = '\0';
    int lwork = decomp().size(0)*F77_NAME(ilaenv)(1, "DSYTRF", uplo_str,
					     decomp().size(0), -1, -1, -1);
    VectorDouble work(lwork);
    F77_CALL(dsytrf)(uplo(), decomp().size(0), &decomp()(0, 0), decomp().gdim(0), &pivot_(0),
		     &work(0), lwork, info);
    if (info < 0)
	throw(LaException("LaBunchKaufmanFactorDouble::ref(const LaSymmMatDouble&)",
			  "illegal input"));
    singular_ = info > 0;
    return *this;
}

inline double LaBunchKaufmanFactorDouble::rcond(double infnorm) const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));

    double ans;

    VectorDouble work(5*decomp().size(0));
    VectorInt iwork(decomp().size(0));
    int info;
    F77_CALL(dsycon)(uplo(), decomp().size(0), &(decomp()(0,0)),
		     decomp().gdim(0), &(pivot()(0)), infnorm, ans,
		     &work(0), &iwork(0), info);
    if (info < 0)
	throw(LaException("LaSymmMatDouble::rcond(char which)",
			  "illegal input"));
    return ans;
}

inline LaMatDouble* LaBunchKaufmanFactorDouble::solve() const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));

    if (singular())
	throw(LaException("LaBunchKaufmanFactorDouble::solve()",
			  "singular matrix"));
    int info;
    LaSymmMatDouble *ans = decomp().clone();
    VectorDouble work(decomp().size(0));
    F77_CALL(dsytri)(uplo(), decomp().size(0), &(*ans)(0,0), decomp().gdim(0),
		     &pivot_(0), &work(0), info);
    if (info < 0)
	throw(LaException("LaBunchKaufmanFactorDouble::solve()",
			  "illegal input"));
    if (info > 0)
	throw(LaException("LaBunchKaufmanFactorDouble::solve()",
			  "singular matrix"));
    return ans;
}

inline LaMatDouble& LaBunchKaufmanFactorDouble::solve(LaMatDouble& B) const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));

    int info;
    F77_CALL(dsytrs)(uplo(), decomp().size(0), B.size(1), &decomp()(0,0),
		     decomp().gdim(0), &pivot_(0), &B(0,0), B.gdim(0), info);
    return B;
}

inline LaMatDouble& LaBunchKaufmanFactorDouble::solve(LaMatDouble& X,
						  const LaMatDouble& B ) const
{
    if (decomp_ == 0)
	throw(LaException("No decomposition present"));

    X.inject(B);
    return solve(X);
}

#endif
