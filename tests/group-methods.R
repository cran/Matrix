### Testing the group methods

library(Matrix)
set.seed(2001)

mm <- Matrix(rnorm(50 * 7), nc = 7)
xpx <- crossprod(mm)# -> "factors" in mm !
round(xpx, 3) # works via "Math2"

y <- rnorm(nrow(mm))
xpy <- crossprod(mm, y)
res <- solve(xpx, xpy)
signif(res, 4) # 7 x 1 Matrix

## exp(): component wise
signif(dd <- (expm(xpx) - exp(xpx)) / 1e34, 3)# 7 x 7

stopifnot(validObject(xpx),
          validObject(xpy),
          validObject(dd))

## "Math" also, for log() and [l]gamma() which need special treatment
stopifnot(identical(exp(res)@x, exp(res@x)),
          identical(log(abs(res))@x, log(abs((res@x)))),
          identical(lgamma(res)@x, lgamma(res@x)))


###--- sparse matrices ---------

m <- Matrix(c(0,0,2:0), 3,5)
(mC <- as(m, "dgCMatrix"))
sm <- sin(mC)
stopifnot(class(sm) == class(mC), class(mC) == class(mC^2),
          dim(sm) == dim(mC),
          class(0 + 100*mC) == class(mC),
          all.equal(0.1 * ((0 + 100*mC)/10), mC),
          all.equal(sqrt(mC ^ 2), mC),
          all.equal(m^m, mC^mC),
          identical(mC^2, mC * mC),
          identical(mC*2, mC + mC)
          )

x <- Matrix(rbind(0,cbind(0, 0:3,0,0,-1:2,0),0))
x # sparse
stopifnot(is(show(x + 10*t(x)), "sparseMatrix"))
(px <- Matrix(x^x - 1))#-> sparse again
stopifnot(px@i == c(3,4,1,4),
          px@x == c(3,26,-2,3))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
