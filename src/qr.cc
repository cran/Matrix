#include "qr.h"
#include "lapackd.h"
#include "vd.h"

#ifdef length
#undef length
#endif

#ifdef append
#undef append
#endif

#include <valarray>

LaQRFactorDouble& LaQRFactorDouble::ref(const LaGenMatDouble& A)
{
    if(A.inc(0) != 1 || A.inc(1) != 1)
	throw(LaException("LaQRFactorDouble::ref(const LaGenMatDouble&)",
			  "input matrix has non unit increment"));
    qr_.ref(A);
    R_.ref(qr_);
    LaVectorInt pivot(A.size(1));
    pivot_.ref(pivot);
    LaVectorDouble qraux(min(A.size(0), A.size(1)));
    qraux_.ref(qraux);
    rank_ = -1;

    int lwork = 2*qr_.size(1)+(qr_.size(1)+1)*F77_NAME(ilaenv)(1, "DGEQP3",
							       "", qr_.size(0),
							       qr_.size(1),
							       -1, -1);
    VectorDouble work(lwork);
    int info;
    F77_CALL(dgeqp3)(qr_.size(0), qr_.size(1), &qr_(0, 0), qr_.gdim(0),
		     &pivot_(0), &qraux_(0), &work(0), lwork, info);
    if (info < 0)
	throw(LaException("LaQRFactorDouble::ref(const LaGenMatDouble&)",
			  "illegal input"));
    return *this;
}

LaMatDouble& LaQRFactorDouble::applyQ(LaMatDouble& y, bool left,
				   bool transpose) const
{
    int info;
    char opts[] = "LT";
    if (!left) opts[0] = 'R';
    if (!transpose) opts[1] = 'N';
    int lwork = y.size(1)*F77_NAME(ilaenv)(1, "DORMQR",
					   opts, y.size(0), y.size(1),
					   qr_.size(0), -1);
    std::valarray<double> work(lwork);
    F77_CALL(dormqr)(left?'L':'R', transpose?'T':'N', y.size(0), y.size(1),
                     qr_.size(0), qr_.addr(), qr_.gdim(0),
                     qraux_.addr(), y.addr(), y.gdim(0), &work[0], lwork,
                     info);
    if (info < 0)
	throw(LaException("LaQRFactorDouble::applyQ",
			  "illegal input "));
    return y;
}

LaGenMatDouble* LaQRFactorDouble::solve() const
{
    if (qr_.size(0) != qr_.size(1))
	throw(LaException("singular matrix"));
    LaUpperTriangMatDouble* Rinv = R_.solve();
    LaGenMatDouble* inv = new LaGenMatDouble(&(*Rinv)(0, 0), Rinv->size(0),
					     Rinv->size(1));
    delete Rinv;
    applyQ(*inv, false, true);
    return inv;
}

LaMatDouble& LaQRFactorDouble::solve(LaMatDouble& B) const
{
//    dynamic_cast<LaGenMatDouble&>(B);
    applyQ(B);
    R_.solve(B);
    F77_CALL(dlaswp)(B.size(1), B.addr(), B.gdim(0), pivot_.start(),
                     pivot_.end(), pivot_.addr(), pivot_.inc());
    return B;
}

LaMatDouble& LaQRFactorDouble::solve(LaMatDouble& X, const LaMatDouble& B ) const
{
//    dynamic_cast<LaGenMatDouble&>(X);
    if (!(X.size(1) == B.size(1))) throw(LaException("assert failed : X.size(1) == B.size(1)"));

    if (X.size(0) == B.size(0)) {
	X.inject(B);
	return solve(X);
    }

    LaGenMatDouble BB;
    BB.copy(B);

    applyQ(BB);
    LaGenMatDouble XX = BB(LaIndex(0, X.size(0) - 1),
			   LaIndex(0, X.size(1) - 1));
    X.inject(XX);
    R_.solve(X);
    F77_CALL(dlaswp)(X.size(1), X.addr(), X.gdim(0), pivot_.start(),
                     pivot_.end(), pivot_.addr(), pivot_.inc());
    return X;
}

//
//     Determine RANK using incremental condition estimation
//     This code is borrowed from the linpack-3.0 dgelsy FORTRAN routine
//
int LaQRFactorDouble::rank(double rcond) {
    if (rank_ < 0) {

	int mn = min(qr_.size(0), qr_.size(1));

	VectorDouble work1(mn);
	VectorDouble work2(mn);
	work1(0) = work2(0) = 1.0;

	double smax = qr_(0, 0);
	if (smax < 0) smax = -smax;
	double smin = smax;
	rank_ = (smin == 0.0) ? 0 : 1;

	for (int i = 1; i < mn; i++) {
	    double sminpr, smaxpr, s1, s2, c1, c2;

	    F77_CALL(dlaic1)(2, i, &work1(0), smin, &qr_(1, i),
			     qr_(i, i), sminpr, s1, c1);
	    F77_CALL(dlaic1)(1, i, &work2(0), smax, &qr_(1, i),
			     qr_(i, i), smaxpr, s2, c2);

	    for (int j = 0; j < i; j++) {
		work1(j) = s1*work1(j);
		work2(j) = s2*work2(j);
	    }
	    work1(i) = c1;
	    work2(i) = c2;
	    smin = sminpr;
	    smax = smaxpr;
	    if (smaxpr*rcond < sminpr)
		rank_++;
	}
    }
    return rank_;
}
