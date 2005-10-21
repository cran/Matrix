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
    A[A < +1] <- 0 ; Am <- A
    Acsc <- as(Am, "dgCMatrix")
    A    <- as(Am, "dgeMatrix")
    b <- matrix(rnorm(400), nrow = 4, ncol = 100)
    B <- as(b, "dgeMatrix")
    assert.EQ.mat(A %*% B, Am %*%  b, tol=0)
    assert.EQ.mat(B %*% A,  b %*% Am, tol=0)
    stopifnot(identical(A, as(Acsc, "dgeMatrix")),
              identical(Acsc, as(A, "dgCMatrix")),
              is.all.equal4(A %*% B, Acsc %*% B,
                            A %*% b, Acsc %*% b),
              is.all.equal4(b %*% A, b %*% Acsc,
                            B %*% A, B %*% Acsc))
}

###--- dgTMatrix {was ./dgTMatrix.R } -------

### Use ``non-unique'' versions of dgTMatrix objects

N <- 200
set.seed(1)
i <- as.integer(round(runif (N, 0, 100)))
j <- as.integer(3* rpois (N, lam=15))
x <- round(rnorm(N), 2)
which(duplicated(cbind(i,j))) # 8 index pairs are duplicated

m1 <- new("dgTMatrix", Dim = c(max(i)+1:1, max(j)+1:1), i = i, j = j, x = x)
mc <- as(m1, "dgCMatrix")
m2 <- as(mc, "dgTMatrix")## the same as 'm1' but without duplicates

stopifnot(!isTRUE(all.equal(m1, m2)),
          all.equal(as(m1,"matrix"), as(m2,"matrix"), tol=1e-15),
          all.equal(crossprod(m1), crossprod(m2), tol=1e-15),
          identical(mc, as(m2, "dgCMatrix")))

### -> uniq* functions now in ../R/Auxiliaries.R
(t2 <- system.time(um2 <- Matrix:::uniq(m1)))



proc.time() # for ``statistical reasons''
