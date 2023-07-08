## METHODS FOR CLASS: pMatrix
## permutation matrices, i.e., matrices with standard unit vectors
## for all rows _and_ all columns
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: could export without dot
.changeMargin <- function(x) {
    x@margin <- if(x@margin == 1L) 2L else 1L
    x@perm <- invertPerm(x@perm)
    x
}


## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("numeric", "pMatrix",
      function(from) {
          J <- new("pMatrix")
          if((m <- length(from)) == 0L)
              return(J)
          if(m > .Machine$integer.max)
              stop("dimensions cannot exceed 2^31-1")
          from.i <- from
          if(anyNA(r <- range(from)) || any(r != c(1L, m)) ||
             (is.double(from) && any(from != (from.i <- as.integer(from)))) ||
             anyDuplicated.default(from.i))
              stop("'perm' slot must be a permutation of seq_along(perm)")
          nms <- names(from)
          J@Dim <- c(m, m)
          J@Dimnames <- list(nms, nms)
          J@perm <- from.i
          J
      })

setAs("nsparseMatrix", "pMatrix",
      function(from) {
          d <- from@Dim
          if((n <- d[1L]) != d[2L])
              stop("attempt to coerce non-square matrix to pMatrix")
          from <- .sparse2g(from)
          J <- new("pMatrix")
          J@Dim <- d
          J@Dimnames <- from@Dimnames
          from. <- as(from, "RsparseMatrix")
          p <- from.@p
          m <- length(p) - 1L
          if(all(p == 0:m) && !anyDuplicated.default(j <- from.@j)) {
              J@perm <- j + 1L
              return(J)
          }
          from. <- as(from, "CsparseMatrix")
          p <- from.@p
          n <- length(p) - 1L
          if(all(p == 0:n) && !anyDuplicated.default(i <- from.@i)) {
              J@perm <- i + 1L
              J@margin <- 2L
              return(J)
          }
          stop("matrix must have exactly one nonzero element in each row and column")
      })

setAs("Matrix", "pMatrix",
      function(from) as(as(from, "nsparseMatrix"), "pMatrix"))

setAs("matrix", "pMatrix",
      function(from) as(as(from, "nsparseMatrix"), "pMatrix"))


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## NB: pMatrix are orthogonal, hence the transpose and inverse coincide

setMethod("t", signature(x = "pMatrix"),
          function(x) {
              r <- new("pMatrix")
              r@Dim <- x@Dim
              r@Dimnames = x@Dimnames[2:1]
              r@perm <- x@perm
              if(x@margin == 1L)
                  r@margin <- 2L
              r
          })

for(.op in c("%*%", "%&%")) {
setMethod(.op, signature(x = "pMatrix", y = "pMatrix"),
          function(x, y) {
              r <- new("pMatrix")
              r@Dim <- mmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
                  if(y@margin == 1L)
                      y@perm[if(x@margin == 1L) x@perm else invertPerm(x@perm)]
                  else {
                      r@margin <- 2L
                      (if(x@margin == 1L) invertPerm(x@perm) else x@perm)[y@perm]
                  }
              r
          })

setMethod(.op, signature(x = "pMatrix", y = "indMatrix"),
          function(x, y) {
              r <- new("indMatrix")
              r@Dim <- mmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
                  if(y@margin == 1L)
                      y@perm[if(x@margin == 1L) x@perm else invertPerm(x@perm)]
                  else {
                      r@margin <- 2L
                      (if(x@margin == 1L) invertPerm(x@perm) else x@perm)[y@perm]
                  }
              r
          })

setMethod(.op, signature(x = "indMatrix", y = "pMatrix"),
          function(x, y) {
              r <- new("indMatrix")
              r@Dim <- mmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
                  if(x@margin == 1L)
                      (if(y@margin == 1L) y@perm else invertPerm(y@perm))[x@perm]
                  else {
                      r@margin <- 2L
                      x@perm[if(y@margin == 1L) invertPerm(x@perm) else y@perm]
                  }
              r
          })
}
rm(.op)

setMethod("%*%", signature(x = "pMatrix", y = "matrix"),
          function(x, y) {
              mmultDim(x@Dim, dim(y), type = 1L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .m2ge(y[perm, , drop = FALSE], "d")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%*%", signature(x = "matrix", y = "pMatrix"),
          function(x, y) {
              mmultDim(dim(x), y@Dim, type = 1L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .m2ge(x[, perm, drop = FALSE], "d")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%*%", signature(x = "pMatrix", y = "Matrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- as(y[perm, , drop = FALSE], "dMatrix")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%*%", signature(x = "Matrix", y = "pMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- as(x[, perm, drop = FALSE], "dMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "pMatrix", y = "matrix"),
          function(x, y) {
              mmultDim(x@Dim, dim(y), type = 1L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .m2ge(y[perm, , drop = FALSE], "n")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "matrix", y = "pMatrix"),
          function(x, y) {
              mmultDim(dim(x), y@Dim, type = 1L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .m2ge(x[, perm, drop = FALSE], "n")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "pMatrix", y = "Matrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- as(y[perm, , drop = FALSE], "nMatrix")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "Matrix", y = "pMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- as(x[, perm, drop = FALSE], "nMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("crossprod", signature(x = "pMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <- new(if(isTRUE(boolArith)) "ldiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames[c(2L, 2L)]
              r@diag <- "U"
              r
          })

setMethod("crossprod", signature(x = "matrix", y = "pMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(dim(x), y@Dim, type = 2L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .m2ge(t(x)[, perm, drop = FALSE],
                         if(isTRUE(boolArith)) "n" else "d")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("crossprod", signature(x = "Matrix", y = "pMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- as(t(x)[, perm, drop = FALSE],
                      if(isTRUE(boolArith)) "nMatrix" else "dMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "pMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <- new(if(isTRUE(boolArith)) "ldiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames[c(1L, 1L)]
              r@diag <- "U"
              r
          })

setMethod("tcrossprod", signature(x = "pMatrix", y = "matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, dim(y), type = 3L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .m2ge(t(y)[perm, , drop = FALSE],
                         if(isTRUE(boolArith)) "n" else "d")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 3L)
              r
          })

setMethod("tcrossprod", signature(x = "pMatrix", y = "Matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- as(t(y)[perm, , drop = FALSE],
                      if(isTRUE(boolArith)) "nMatrix" else "dMatrix")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 3L)
              r
          })
