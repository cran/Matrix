## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a dgeMatrix.

## setMethod("%*%", signature(x = "dsparseMatrix", y = "ddenseMatrix"),
##           function(x, y) callGeneric(x, as(y, "dgeMatrix")))

## setMethod("%*%", signature(x = "ddenseMatrix", y = "dsparseMatrix"),
##           function(x, y) callGeneric(as(x, "dgeMatrix"), y))

setMethod("crossprod", signature(x = "dsparseMatrix", y = "ddenseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "dsparseMatrix"),
          function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y))

## and coerce dsparse* to dgC*
## setMethod("%*%", signature(x = "dsparseMatrix", y = "dgeMatrix"),
##           function(x, y) callGeneric(as(x, "dgCMatrix"), y))

## setMethod("%*%", signature(x = "dgeMatrix", y = "dsparseMatrix"),
##           function(x, y) callGeneric(x, as(y, "dgCMatrix")))

setMethod("crossprod", signature(x = "dsparseMatrix", y = "dgeMatrix"),
## NB: using   callGeneric(.) here, leads to infinite recursion :
          function(x, y = NULL) .Call(Csparse_dense_crossprod, as(x, "dgCMatrix"), y))

## NB: there's already
##     ("CsparseMatrix", "missing") and ("TsparseMatrix", "missing") methods

setMethod("crossprod", signature(x = "dgeMatrix", y = "dsparseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "dgCMatrix")))

setMethod("image", "dsparseMatrix",
          function(x, ...) image(as(x, "dgTMatrix"), ...))

setMethod("kronecker", signature(X = "dsparseMatrix", Y = "dsparseMatrix"),
          function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
          callGeneric(as(X, "dgTMatrix"), as(Y, "dgTMatrix")))

setMethod("chol", signature(x = "dsparseMatrix", pivot = "ANY"),
           function(x, pivot, ...) {
               px <- as(x, "dsCMatrix")
               if (isTRUE(validObject(px, test=TRUE))) chol(px, pivot)
               else stop("'x' is not positive definite -- chol() undefined.")
           })

setMethod("lu", signature(x = "dsparseMatrix"),
	  function(x, ...) callGeneric(as(x, "dgCMatrix")))


## Group Methods, see ?Arith (e.g.)
## -----

##-> now moved to ./Csparse.R (and 'up' to ./sparseMatrix.R):
##  "Math2" is in ./dMatrix.R


