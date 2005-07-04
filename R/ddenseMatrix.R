### Define Methods that can be inherited for all subclasses

## -- see also ./Matrix.R  e.g., for a show() method

## These methods are the 'fallback' methods for all dense numeric
## matrices in that they simply coerce the ddenseMatrix to a
## dgeMatrix. Methods for special forms override these.

setMethod("norm", signature(x = "ddenseMatrix", type = "missing"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("norm", signature(x = "ddenseMatrix", type = "character"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix"), type))

setMethod("rcond", signature(x = "ddenseMatrix", type = "missing"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("rcond", signature(x = "ddenseMatrix", type = "character"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix"), type))

setMethod("t", signature(x = "ddenseMatrix"),
	  function(x) callGeneric(as(x, "dgeMatrix")))

setMethod("tcrossprod", signature(x = "ddenseMatrix"),
	  function(x) callGeneric(as(x, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "missing"),
          function(x, y = NULL) callGeneric(as(x, "dgeMatrix")))

setMethod("diag", signature(x = "ddenseMatrix"),
          function(x = 1, nrow, ncol = n) callGeneric(as(x, "dgeMatrix")))

setMethod("solve", signature(a = "ddenseMatrix", b = "missing"),
          function(a, b, ...) callGeneric(as(a, "dgeMatrix")))

setMethod("solve", signature(a = "ddenseMatrix", b = "ANY"),
          function(a, b, ...) callGeneric(as(a, "dgeMatrix"), b))

setMethod("lu", signature(x = "ddenseMatrix"),
          function(x, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "missing"),
          function(x, logarithm, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
          callGeneric(as(x, "dgeMatrix"), logarithm))

setMethod("expm", signature(x = "ddenseMatrix"),
          function(x) callGeneric(as(x, "dgeMatrix")))

setMethod("Schur", signature(x = "ddenseMatrix", vectors = "missing"),
          function(x, vectors, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("Schur", signature(x = "ddenseMatrix", vectors = "logical"),
          function(x, vectors, ...) callGeneric(as(x, "dgeMatrix"), vectors))


