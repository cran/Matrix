## METHODS FOR GENERIC: c
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c.Matrix <- function(...) {
    if(nargs() == 0L)
        return(NULL)
    args <- lapply(list(...), as.vector)
    unlist(args, FALSE, TRUE)
}

c.sparseVector <- function(...) {
    N <- nargs()
    if(N == 0L)
        return(NULL)
    args        <- lapply(list(...), as, "sparseVector")
    args.length <- vapply(args, slot, 0, "length")
    args.i      <- lapply(args, slot,         "i")
    args.nnz    <- lengths(args.i, FALSE)

    s <- c("n", "l", "i", "d", "z")
    i <- match(vapply(args, .M.kind, ""), s)
    k <- range(i)
    n <- sum(args.length)
    a <- if(n - 1 <= .Machine$integer.max) as.integer else as.double

    r <- new(paste0(s[k[2L]], "sparseVector"))
    r@length <- a(n)
    r@i <- a(unlist(args.i, FALSE, FALSE)) +
        rep.int(cumsum(c(0L, a(args.length)[-N])), args.nnz)
    if(k[2L] > 1L) {
        if(k[1L] > 1L)
            args.x <- lapply(args, slot, "x")
        else {
            pattern <- i == 1L
            args.x <- vector("list", N)
            args.x[!pattern] <- lapply(args    [!pattern],    slot,      "x")
            args.x[ pattern] <- lapply(args.nnz[ pattern], rep.int, x = TRUE)
        }
        r@x <- unlist(args.x, FALSE, FALSE)
    }
    r
}

## These are insufficient as dispatch only consides the first argument,
## which need not be a Matrix or sparseVector:
if(FALSE) {
setMethod("c",       "Matrix", function(x, ...) c.Matrix      (x, ...))
setMethod("c", "sparseVector", function(x, ...) c.sparseVector(x, ...))
}


## METHODS FOR GENERIC: cbind, rbind
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: not yet registered or exported
cbind.Matrix <- function(..., deparse.level = 1)
    .External(R_bind, deparse.level, 1L, substitute(list(...)), ...)
rbind.Matrix <- function(..., deparse.level = 1)
    .External(R_bind, deparse.level, 0L, substitute(list(...)), ...)


## METHODS FOR GENERIC: cbind2, rbind2
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.cbind2 <- function(x, y, ...) cbind.Matrix(x, y, deparse.level = 0L)
.rbind2 <- function(x, y, ...) rbind.Matrix(x, y, deparse.level = 0L)

setMethod("cbind2", signature(x = "Matrix", y = "missing"),
          function(x, y, ...) x)
setMethod("rbind2", signature(x = "Matrix", y = "missing"),
          function(x, y, ...) x)

setMethod("cbind2", signature(x = "Matrix", y = "NULL"), .cbind2)
setMethod("cbind2", signature(x = "NULL", y = "Matrix"), .cbind2)
setMethod("rbind2", signature(x = "Matrix", y = "NULL"), .rbind2)
setMethod("rbind2", signature(x = "NULL", y = "Matrix"), .rbind2)

setMethod("cbind2", signature(x = "Matrix", y = "vector"), .cbind2)
setMethod("cbind2", signature(x = "vector", y = "Matrix"), .cbind2)
setMethod("rbind2", signature(x = "Matrix", y = "vector"), .rbind2)
setMethod("rbind2", signature(x = "vector", y = "Matrix"), .rbind2)

setMethod("cbind2", signature(x = "Matrix", y = "matrix"), .cbind2)
setMethod("cbind2", signature(x = "matrix", y = "Matrix"), .cbind2)
setMethod("rbind2", signature(x = "Matrix", y = "matrix"), .rbind2)
setMethod("rbind2", signature(x = "matrix", y = "Matrix"), .rbind2)

setMethod("cbind2", signature(x = "Matrix", y = "Matrix"), .cbind2)
setMethod("rbind2", signature(x = "Matrix", y = "Matrix"), .rbind2)

rm(.cbind2, .rbind2)
