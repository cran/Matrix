#### Superclass Methods for all sparse logical matrices

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
          function(e1) !as(e1, "lgeMatrix"))

