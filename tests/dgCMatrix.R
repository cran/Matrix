library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))

data(mm)
stopifnot(##is(mm) == c("dgCMatrix", "dMatrix", "Matrix"),
          dim(mm) == (dm <- c(1850, 712)),
          identical(dimnames(mm), list(NULL,NULL)))
str(mm)
tmm <- t(mm)
str(tmm)

str(mTm <- crossprod(mm))
mmT  <- crossprod(tmm) ## should be the same as
mmT. <- tcrossprod(mm)
## but not quite: even length( * @ x ) differs!
str(mmT, max=2)# much larger than mTm (i.e less sparse)
str(mmT., max=2)# x slot is currently slightly larger --> improve tcrossprod()?
system.time(ae <- all.equal(as(mmT.,"matrix"), as(mmT,"matrix"), tol = 1e-14))
## 4-5 seconds on a 850 MHz, P III
stopifnot(ae)

stopifnot(validObject(tmm), dim(tmm) == dm[2:1],
          validObject(mTm), dim(mTm) == dm[c(2,2)],
          validObject(mmT), dim(mmT) == dm[c(1,1)],
          identical(as(tmm, "matrix"), t(as(mm, "matrix"))))

## from a bug report by Guissepe Ragusa <gragusa@ucsd.edu>
set.seed(101)
for(i in 1:10) {
    A <- matrix(rnorm(400), nrow = 100, ncol = 4)
    A[A < +1] <- 0
    Acsc <- as(A, "dgCMatrix")
    A <- as(A, "dgeMatrix")
    b <- matrix(rnorm(400), nrow = 4, ncol = 100)
    B <- as(b, "dgeMatrix")
    stopifnot(is.all.equal4(A %*% B, Acsc %*% B,
                            A %*% b, Acsc %*% b),
              is.all.equal4(b %*% A, b %*% Acsc,
                            B %*% A, B %*% Acsc))
}

proc.time() # for ``statistical reasons''
