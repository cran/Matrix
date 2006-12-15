#### Will be sourced by several R scripts in ../tests/

paste0 <- function(...) paste(..., sep = '')

identical3 <- function(x,y,z)	identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d) identical(a,b) && identical3(b,c,d)

as.mat <- function(m) {
    ## as(., "matrix")	but with no extraneous empty dimnames
    m <- as(m, "matrix")
    if(identical(dimnames(m), list(NULL,NULL)))
	dimnames(m) <- NULL
    m
}

asD <- function(m) { ## as "Dense"
    if(canCoerce(m, "denseMatrix")) as(m, "denseMatrix")
    else if(canCoerce(m, (cl <- paste(.M.kind(m), "denseMatrix", sep=''))))
        as(m, cl)
    else if(canCoerce(m, "dgeMatrix")) as(m, "dgeMatrix")
    else stop("cannot coerce to a typical dense Matrix")
}

assert.EQ.mat <- function(M, m, tol = if(show) 0 else 1e-15, show=FALSE) {
    ## Purpose: check equality of  'Matrix' M with  'matrix' m
    ## ----------------------------------------------------------------------
    ## Arguments: M: is(., "Matrix")
    ##            m: is(., "matrix")
    ##            show: if TRUE, return (and hence typically print) all.equal(...)
    MM <- as.mat(M)                     # as(M, "matrix")
    if(is.logical(MM) && is.numeric(m))
	storage.mode(MM) <- "integer"
    attr(MM, "dimnames") <- attr(m, "dimnames") <- NULL
    if(show) all.equal(MM, m, tol = tol)
    else stopifnot(all.equal(MM, m, tol = tol))
}

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

is.all.equal3 <- function(x,y,z, tol = .Machine$double.eps^0.5)
    isTRUE(all.equal(x,y, tol=tol)) && isTRUE(all.equal(y,z, tol=tol))

is.all.equal4 <- function(x,y,z,u, tol = .Machine$double.eps^0.5)
    is.all.equal3(x,y,z, tol=tol) && isTRUE(all.equal(z,u, tol=tol))


Qidentical <- function(x,y) {
    ## quasi-identical:
    if(class(x) != class(y)) return(FALSE)
    slts <- slotNames(x)
    if("factors" %in% slts) { ## allow one empty and one non-empty 'factors'
        slts <- slts[slts != "factors"]
        ## if both are not empty, they must be the same:
        if(length(xf <- x@factors) && length(yf <- y@factors))
            if(!identical(xf, yf)) return(FALSE)
    }
    for(sl in slts)
        if(!identical(slot(x,sl), slot(y,sl)))
            return(FALSE)
    TRUE
}

## The relative error typically returned by all.equal:
relErr <- function(target, current)
    mean(abs(target - current)) / mean(abs(target))

is.R22 <- (paste(R.version$major, R.version$minor, sep=".") >= "2.2")
