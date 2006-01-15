library(Matrix)

### Do all kinds of object creation and coercion

source(system.file("test-tools.R", package = "Matrix"))

## the empty ones:
chk.matrix(new("dgeMatrix"))

## "dge"
assertError( new("dgeMatrix", Dim = c(2,2), x= 1:4) )# double 'Dim'
assertError( new("dgeMatrix", Dim = as.integer(c(2,2)), x= 1:4) )# int 'x'
assertError( new("dgeMatrix", Dim = 2:2, x=as.double(1:4)) )# length(Dim) !=2
assertError( new("dgeMatrix", Dim = as.integer(c(2,2)), x= as.double(1:5)))

chk.matrix(m1 <- Matrix(1:6, ncol=2))
chk.matrix(m2 <- Matrix(1:7, ncol=3)) # a (desired) warning
stopifnot(unique(is(m1)) == c("dgeMatrix", "ddenseMatrix", "generalMatrix",
	    "dMatrix", "denseMatrix", "Matrix", "compMatrix"),
	  dim(t(m1)) == 2:3,
	  identical(m1, t(t(m1))))
c.nam <- paste("C",1:2, sep='')
dimnames(m1) <- list(NULL, c.nam)
stopifnot(colnames(m1) == c.nam,
	  identical(dimnames(t(m1)), list(c.nam, NULL)),
	  identical(m1, t(t(m1))))

## an example of *named* dimnames
(t34N <- as(unclass(table(x = gl(3,4), y=gl(4,3))), "dgeMatrix"))
stopifnot(identical(dimnames(t34N),
		    dimnames(as(t34N, "matrix"))))

## "dpo"
chk.matrix(cm <- crossprod(m1))
chk.matrix(cp <- as(cm, "dppMatrix"))# 'dpp'
chk.matrix(as(cm, "dsyMatrix"))
chk.matrix(dcm <- as(cm, "dgeMatrix"))
chk.matrix(mcm <- as(cm, "dMatrix"))
##BUG - FIXME: stopifnot(identical(dcm, mcm))
##---	------	mcm[2,1] is garbage
chk.matrix(as(cm, "Matrix"))

## Cholesky
chk.matrix(ch <- chol(cm))
chk.matrix(ch2 <- chol(as(cm, "dsyMatrix")))
#not yet{FIXME}: chk.matrix(ch3 <- chol(as(cm, "dgeMatrix")))
stopifnot(all.equal(as(ch, "matrix"), as(ch2, "matrix")))

### Very basic	triangular matrix stuff

assertError( new("dtrMatrix", Dim = c(2,2), x= 1:4) )# double 'Dim'
if(paste(R.version$major, R.version$minor, sep=".") >= "2.0.1")
assertError( new("dtrMatrix", Dim = as.integer(c(2,2)), x= 1:4) )# int 'x'
## This caused a segfault (before revision r1172 in ../src/dtrMatrix.c):
assertError( new("dtrMatrix", Dim = 2:2, x=as.double(1:4)) )# length(Dim) !=2
assertError( new("dtrMatrix", Dim = as.integer(c(2,2)), x= as.double(1:5)))

tr22 <- new("dtrMatrix", Dim = as.integer(c(2,2)), x=as.double(1:4))
tt22 <- t(tr22)
(tPt <- tr22 + tt22)
stopifnot(identical(10 * tPt, tPt * 10),
	  (t.22 <- (tr22 / .5)* .5)@x == c(1,0,3,4),
	  TRUE) ## not yet: class(t.22) == "dtrMatrix")

## non-square triagonal Matrices --- are forbidden ---
assertError(new("dtrMatrix", Dim = 2:3,
                x=as.double(1:6), uplo="L", diag="U"))
