library(Matrix)

set.seed(123)
mm <- Matrix(rnorm(500 * 150), nc = 150)
stopifnot(validObject(mm))
xpx <- crossprod(mm)# alters mm !
stopifnot(validObject(mm))
stopifnot(validObject(xpx))
str(mm)
str(xpx)
xpy <- crossprod(mm, rnorm(500))
res <- solve(xpx, xpy)
str(xpx)#
stopifnot(validObject(xpx))

proc.time()
