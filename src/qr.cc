#include "qr.h"
#include "lapackd.h"
#include "vd.h"

LaQRFactorDouble& LaQRFactorDouble::ref(const LaGenMatDouble& A)
{
    assert(A.inc(0) == 1 && A.inc(1) == 1);
    qr_.ref(A);
    R_.ref(qr_);
    LaVectorInt pivot(A.size(1));
    pivot_.ref(pivot);
    LaVectorDouble qraux(min(A.size(0), A.size(1)));
    qraux_.ref(qraux);
    rank_ = 0;

    double* work = new double[3*qr_.size(1)];
    int info;
    F77_CALL(dgeqpf)(qr_.size(0), qr_.size(1), &qr_(0, 0), qr_.gdim(0),
		     &pivot_(0), &qraux_(0), work, info);
    delete[] work;
    return *this;
}

LaMatrix& LaQRFactorDouble::solve(LaMatrix& B) const
{
//    dynamic_cast<LaGenMatDouble&>(B);
    int info;
    int nb = F77_NAME(ilaenv)(1, "DORMQR",
			      "LT", B.size(0), B.size(1),
			      qr_.size(0), -1);
    double* work = new double[B.size(1)*nb];
    F77_CALL(dormqr)('L', 'T', B.size(0), B.size(1),
		     qr_.size(0), &qr_(0, 0), qr_.gdim(0),
		     &qraux_(0), &B(0, 0), B.gdim(0), work, B.size(1)*nb,
		     info);
    delete[] work;
    R_.solve(B);
    F77_CALL(dlaswp)(B.size(1), &B(0, 0), B.gdim(0), pivot_.start(),
		     pivot_.end(), &pivot_(0), pivot_.inc());
    return B;
}

LaMatrix& LaQRFactorDouble::solve(LaMatrix& X, const LaMatrix& B ) const
{
//    dynamic_cast<LaGenMatDouble&>(X);
    assert(X.size(1) == B.size(1));

    if (X.size(0) == B.size(0)) {
	X.inject(B);
	return solve(X);
    }

    LaGenMatDouble BB;
    BB.copy(B);

    int info;
    int nb = F77_NAME(ilaenv)(1, "DORMQR",
			      "LT", BB.size(0), BB.size(1),
			      qr_.size(0), -1);
    double* work = new double[BB.size(1)*nb];
    F77_CALL(dormqr)('L', 'T', BB.size(0), BB.size(1),
		     qr_.size(0), &qr_(0, 0), qr_.gdim(0),
		     &qraux_(0), &BB(0, 0), BB.gdim(0), work, BB.size(1)*nb,
		     info);
    delete[] work;
    LaGenMatDouble XX = BB(LaIndex(0, X.size(0) - 1),
			   LaIndex(0, X.size(1) - 1));
    X.inject(XX);
    R_.solve(X);
    F77_CALL(dlaswp)(X.size(1), &X(0, 0), X.gdim(0), pivot_.start(),
		     pivot_.end(), &pivot_(0), pivot_.inc());
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
