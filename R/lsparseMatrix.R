#### Superclass Methods for all sparse logical matrices

setAs("CsparseMatrix", "lsparseMatrix",
      function(from) as(.Call(Csparse_to_nz_pattern, from,
			      is(from, "triangularMatrix")), "lsparseMatrix"))

setAs("lsparseMatrix", "matrix",
      function(from) as(as(from, "ldenseMatrix"), "matrix"))


###------- Work via  as(*, lgC) : ------------

## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a lsparseMatrix.

setMethod("%*%", signature(x = "lsparseMatrix", y = "ldenseMatrix"),
          function(x, y) callGeneric(x, as(y, "lsparseMatrix")))

setMethod("%*%", signature(x = "ldenseMatrix", y = "lsparseMatrix"),
          function(x, y) callGeneric(as(x, "lsparseMatrix"), y))

setMethod("crossprod", signature(x = "lsparseMatrix", y = "ldenseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "lsparseMatrix")))

setMethod("crossprod", signature(x = "ldenseMatrix", y = "lsparseMatrix"),
          function(x, y = NULL) callGeneric(as(x, "lsparseMatrix"), y))

## and coerce lsparse* to lgC*
setMethod("%*%", signature(x = "lsparseMatrix", y = "lsparseMatrix"),
          function(x, y) callGeneric(as(x, "lgCMatrix"), as(y, "lgCMatrix")))

setMethod("crossprod", signature(x = "lsparseMatrix", y = "lsparseMatrix"),
          function(x, y = NULL)
          callGeneric(as(x, "lgCMatrix"), as(y, "lgCMatrix")))

setMethod("!", "lsparseMatrix",
          ## turns FALSE to TRUE --> dense matrix
          function(e1) !as(e1, "denseMatrix"))# was "lgeMatrix"


setMethod("all", signature(x = "lsparseMatrix"),
	  function(x, ..., na.rm = TRUE)
	  !is(x, "triangularMatrix") && all(x@x, ..., na.rm = na.rm))

setMethod("any", signature(x = "lsparseMatrix"),
	  function(x, ..., na.rm = TRUE)
	  ## logical unit-triangular has TRUE diagonal:
	  (is(x, "triangularMatrix") && x@diag == "U") ||
	  any(x@x, ..., na.rm = na.rm))
