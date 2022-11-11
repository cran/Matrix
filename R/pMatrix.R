#### Permutation Matrices -- Coercion and Methods

### NB "pMatrix" extends "indMatrix" and inherits methods -->  indMatrix.R

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("integer", "pMatrix",
      function(from) {
          if ((n <- length(from)) == 0L)
              return(new("pMatrix"))
          if(anyNA(from))
              stop("'perm' slot cannot contain NA")
          if(min(from) < 1L)
              stop("elements of 'perm' slot must be positive integers")
          nms <- names(from)
          new("pMatrix", Dim = c(n, n), Dimnames = list(nms, nms), perm = from)
      })

setAs("numeric", "pMatrix",
      function(from) {
          if ((n <- length(from)) == 0L)
              return(new("pMatrix"))
          if(anyNA(from))
              stop("'perm' slot cannot contain NA")
          r <- range(from)
          if(r[2L] > .Machine$integer.max)
              stop("elements of 'perm' slot cannot exceed 2^31-1")
          if(r[1L] < 1 || any(from != (from.i <- as.integer(from))))
              stop("elements of 'perm' slot must be positive integers")
          nms <- names(from)
         new("pMatrix", Dim = c(n, n), Dimnames = list(nms, nms), perm = from.i)
      })

setAs("nsparseMatrix", "pMatrix",
      function(from) {
          d <- from@Dim
          if((n <- d[1L]) != d[2L])
              stop("attempt to a coerce a non-square matrix to pMatrix")
          from <- .sparse2g(as(from, "RsparseMatrix"))
          p <- from@p
          if(n > 0L && any(p != 0:n))
              stop("matrix must have exactly one nonzero element in each row")
          new("pMatrix", Dim = from@Dim, Dimnames = from@Dimnames,
              perm = from@j + 1L) # validity method checks 'perm' for duplicates
      })

setAs("Matrix", "pMatrix",
      function(from) as(as(from, "nsparseMatrix"), "pMatrix"))

setAs("matrix", "pMatrix",
      function(from) as(as(from, "nsparseMatrix"), "pMatrix"))


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## NB: pMatrix are orthogonal, hence the transpose and inverse coincide

setMethod("t", signature(x = "pMatrix"),
          function(x) {
              x@perm <- invPerm(x@perm)
              x@Dimnames <- x@Dimnames[2:1]
              x
          })

## setMethod("%*%", signature(x = "pMatrix", y = "matrix"), .) # inherited

## setMethod("%*%", signature(x = "pMatrix", y = "Matrix"), .) # inherited

## setMethod("%*%", signature(x = "pMatrix", y = "indMatrix"), .) # inherited

setMethod("%*%", signature(x = "matrix", y = "pMatrix"),
	  function(x, y) {
              mmultDim(dim(x), y@Dim, type = 1L)
              r <- .m2ge(x[, invPerm(y@perm), drop = FALSE], "d")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%*%", signature(x = "Matrix", y = "pMatrix"),
	  function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- as(x[, invPerm(y@perm), drop = FALSE], "dMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%*%", signature(x = "indMatrix", y = "pMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              x@perm <- y@perm[x@perm]
              x@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              x
          })

## setMethod("%&%", signature(x = "pMatrix", y = "matrix"), .) # inherited

## setMethod("%&%", signature(x = "pMatrix", y = "Matrix"), .) # inherited

## setMethod("%&%", signature(x = "pMatrix", y = "indMatrix"), .) # inherited

setMethod("%&%", signature(x = "matrix", y = "pMatrix"),
	  function(x, y) {
              mmultDim(dim(x), y@Dim, type = 1L)
              r <- .m2ge(x[, invPerm(y@perm), drop = FALSE], "n")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "Matrix", y = "pMatrix"),
	  function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- as(x[, invPerm(y@perm), drop = FALSE], "nMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "indMatrix", y = "pMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              x@perm <- y@perm[x@perm]
              x@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              x
          })

setMethod("crossprod", signature(x = "pMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <- new(if(isTRUE(boolArith)) "ldiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames[c(2L, 2L)]
              r@diag <- "U"
              r
          })

setMethod("crossprod", signature(x = "pMatrix", y = "matrix"),
	  function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, dim(y), type = 2L)
              r <- .m2ge(y[invPerm(x@perm), , drop = FALSE],
                         if(isTRUE(boolArith)) "n" else "d")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 2L)
              r
          })

setMethod("crossprod", signature(x = "pMatrix", y = "Matrix"),
	  function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- as(y[invPerm(x@perm), , drop = FALSE],
                      if(isTRUE(boolArith)) "nMatrix" else "dMatrix")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 2L)
              r
          })

setMethod("crossprod", signature(x = "pMatrix", y = "indMatrix"),
	  function(x, y = NULL, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
	      y@perm <- y@perm[invPerm(x@perm)]
              y@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 2L)
	      y
	  })

setMethod("crossprod", signature(x = "matrix", y = "pMatrix"),
	  function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(dim(x), y@Dim, type = 2L)
              r <- .m2ge(t(x)[, invPerm(y@perm), drop = FALSE],
                         if(isTRUE(boolArith)) "n" else "d")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("crossprod", signature(x = "Matrix", y = "pMatrix"),
	  function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- as(t(x)[, invPerm(y@perm), drop = FALSE],
                      if(isTRUE(boolArith)) "nMatrix" else "dMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

## setMethod("crossprod", signature(x = "indMatrix", y = "pMatrix"), .) # inherited

setMethod("tcrossprod", signature(x = "pMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <- new(if(isTRUE(boolArith)) "ldiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames[c(1L, 1L)]
              r@diag <- "U"
              r
          })

## setMethod("tcrossprod", signature(x = "pMatrix", y = "matrix"), .) # inherited

## setMethod("tcrossprod", signature(x = "pMatrix", y = "Matrix"), .) # inherited

## setMethod("tcrossprod", signature(x = "pMatrix", y = "indMatrix"), .) # inherited

setMethod("tcrossprod", signature(x = "matrix", y = "pMatrix"),
	  function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(dim(x), y@Dim, type = 3L)
              r <- .m2ge(x[, y@perm, drop = FALSE],
                         if(isTRUE(boolArith)) "n" else "d")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 3L)
              r
          })

setMethod("tcrossprod", signature(x = "Matrix", y = "pMatrix"),
	  function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- as(t(x)[, y@perm, drop = FALSE],
                      if(isTRUE(boolArith)) "nMatrix" else "dMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 3L)
              r
          })

setMethod("tcrossprod", signature(x = "indMatrix", y = "pMatrix"),
	  function(x, y = NULL, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              x@perm <- invPerm(y@perm)[x@perm]
	      x@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 3L)
	      x
	  })




