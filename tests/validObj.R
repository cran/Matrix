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

## "dMatrix"
str(new("Matrix"))

## "dge"
assertError( new("dgeMatrix", Dim = c(2,2), x= 1:4) )# double 'Dim'
if(FALSE)## FIXME: this creates an integer '@ x' !
assertError( new("dgeMatrix", Dim = as.integer(c(2,2)), x= 1:4) )# int 'x'
assertError( new("dgeMatrix", Dim = 2:2, x=as.double(1:4)) )# length(Dim) !=2
assertError( new("dgeMatrix", Dim = as.integer(c(2,2)), x= as.double(1:5)))

chk.matrix(m1 <- Matrix(1:6, ncol=2))
chk.matrix(m2 <- Matrix(1:7, ncol=3)) # a warning

## "dpo"
chk.matrix(cm <- crossprod(m1))
chk.matrix(as(cm, "dsyMatrix"))
chk.matrix(as(cm, "dgeMatrix"))
try( chk.matrix(as(cm, "Matrix")) ) # gives an error

## Cholesky
chk.matrix(ch <- chol(cm))
## FIXME:
try( chk.matrix(ch2 <- chol(as(cm, "dsyMatrix"))) ) # should not give an error
try( chk.matrix(ch3 <- chol(as(cm, "dgeMatrix"))) ) # nor that one

### Very basic  triangular matrix stuff

assertError( new("dtrMatrix", Dim = c(2,2), x= 1:4) )# double 'Dim'
if(FALSE)## FIXME: this creates an integer '@ x' !
assertError( new("dtrMatrix", Dim = as.integer(c(2,2)), x= 1:4) )# int 'x'
if(FALSE)## FIXME: this causes a segfault
assertError( new("dtrMatrix", Dim = 2:2, x=as.double(1:4)) )# length(Dim) !=2
assertError( new("dtrMatrix", Dim = as.integer(c(2,2)), x= as.double(1:5)))

tr22 <- new("dtrMatrix", Dim = as.integer(c(2,2)), x=as.double(1:4))
try( t(tr22) ) # fails -- FIXME

## non-square
tru <- new("dtrMatrix", Dim = 2:3, x=as.double(1:6), uplo="L", diag="U")
trn <- new("dtrMatrix", Dim = 2:3, x=as.double(1:6), uplo="L", diag="N")
try( tru + trn ) # not yet

try( t(tru) ) ## FIXME !
try( t(trn) ) ## FIXME

