library(Matrix)

### Do all kinds of object creation and coercion

chk.matrix <- function(M) {
    ## check object; including coercion to "matrix" :
    cl <- class(M)
    cat("class ", dQuote(cl), " [",nrow(M)," x ",ncol(M),"]; slots (",
        paste(slotNames(M), collapse=","), ")\n", sep='')
    stopifnot(validObject(M),
              dim(M) == c(nrow(M), ncol(M)),
              identical(dim(m <- as(M, "matrix")), dim(M))
              )
}

## Make sure errors are signaled
assertError <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- try(expr, silent = TRUE)
    if(!inherits(t.res, "try-error"))
        stop(d.expr, "\n\t did not give an error", call. = FALSE)
    invisible(t.res)
}

## the empty ones:
chk.matrix(new("dgeMatrix"))

## "dge"
assertError( new("dgeMatrix", Dim = c(2,2), x= 1:4) )# double 'Dim'
assertError( new("dgeMatrix", Dim = as.integer(c(2,2)), x= 1:4) )# int 'x'
assertError( new("dgeMatrix", Dim = 2:2, x=as.double(1:4)) )# length(Dim) !=2
assertError( new("dgeMatrix", Dim = as.integer(c(2,2)), x= as.double(1:5)))

chk.matrix(m1 <- Matrix(1:6, ncol=2))
chk.matrix(m2 <- Matrix(1:7, ncol=3)) # a warning
stopifnot(is(m1) == c("dgeMatrix", "ddenseMatrix", "dMatrix", "Matrix"),
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
chk.matrix(as(cm, "dsyMatrix"))
chk.matrix(dcm <- as(cm, "dgeMatrix"))
chk.matrix(mcm <- as(cm, "dMatrix"))
#BUG - FIXME: stopifnot(identical(dcm, mcm))
try( chk.matrix(as(cm, "Matrix")) )# gives an error: "Matrix" has NULL 'dim()'

## Cholesky
chk.matrix(ch <- chol(cm))
#if(FALSE)# fails for Doug in R-devel (2005-06-06) :
chk.matrix(ch2 <- chol(as(cm, "dsyMatrix")))
#not yet{FIXME}: chk.matrix(ch3 <- chol(as(cm, "dgeMatrix")))
#if(FALSE)# ...R-devel
stopifnot(all.equal(as(ch, "matrix"), as(ch2, "matrix")))

### Very basic  triangular matrix stuff

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


## non-square  triagonal Matrices --- should this be forbidden anyway? ---
tru <- new("dtrMatrix", Dim = 2:3, x=as.double(1:6), uplo="L", diag="U")
trn <- new("dtrMatrix", Dim = 2:3, x=as.double(1:6), uplo="L", diag="N")
tru + trn  # a 'dgeMatrix'

as(t(tru),"dgeMatrix")
as(t(trn),"dgeMatrix")
as(t(t(tru)), "dgeMatrix")# pretty non sense
