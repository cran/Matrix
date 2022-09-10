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

## pMatrix are orthogonal, hence the transpose and inverse coincide
setMethod("t", signature(x = "pMatrix"),
          function(x) {
              x@perm <- invPerm(x@perm)
              x@Dimnames <- x@Dimnames[2:1]
              x
          })

## FIXME: none of these handle '[dD]imnames' consistently with base products

setMethod("%*%", signature(x = "matrix", y = "pMatrix"),
	  function(x, y) { mmultCheck(x,y); x[, invPerm(y@perm)] })
setMethod("%*%", signature(x = "Matrix", y = "pMatrix"),
	  function(x, y) { mmultCheck(x,y); x[, invPerm(y@perm)] })

setMethod("%*%", signature(x = "pMatrix", y = "pMatrix"),
	  function(x, y) {
              stopifnot(identical(x@Dim, y@Dim))
	      x@perm <- x@perm[y@perm]
	      x
	  })

setMethod("crossprod", signature(x = "pMatrix", y = "matrix"),
	  function(x, y) { mmultCheck(x,y, 2L); y[invPerm(x@perm) ,]})
setMethod("crossprod", signature(x = "pMatrix", y = "Matrix"),
	  function(x, y) { mmultCheck(x,y, 2L); y[invPerm(x@perm) ,]})
setMethod("crossprod", signature(x = "pMatrix", y = "pMatrix"),
	  function(x, y) {
	      stopifnot(identical(x@Dim, y@Dim))
	      x@perm <- invPerm(x@perm)[y@perm]
	      x
	  })

setMethod("tcrossprod", signature(x = "pMatrix", y = "pMatrix"),
	  function(x, y) {
	      stopifnot(identical(x@Dim, y@Dim))
	      x@perm <- x@perm[invPerm(y@perm)]
	      x
	  })


setMethod("crossprod", signature(x = "pMatrix", y = "missing"),
          function(x, y=NULL) Diagonal(nrow(x)))
setMethod("tcrossprod", signature(x = "pMatrix", y = "missing"),
          function(x, y=NULL) Diagonal(nrow(x)))

