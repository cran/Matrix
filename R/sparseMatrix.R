### Define Methods that can be inherited for all subclasses

## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a dgeMatrix.

setMethod("%*%", signature(x = "sparseMatrix", y = "ddenseMatrix"),
          function(x, y) callGeneric(x, as(y, "dgeMatrix")))

setMethod("%*%", signature(x = "ddenseMatrix", y = "sparseMatrix"),
          function(x, y) callGeneric(as(x, "dgeMatrix"), y))

setMethod("crossprod", signature(x = "sparseMatrix", y = "ddenseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "sparseMatrix"),
          function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y))

