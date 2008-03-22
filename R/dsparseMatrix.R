### d(ouble)sparseMatrix methods :

## Note that '%*%' are now handled via "sparse*" , "Csparse*" etc

setMethod("crossprod", signature(x = "dsparseMatrix", y = "ddenseMatrix"),
	  function(x, y = NULL) crossprod(x, as(y, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "dsparseMatrix"),
	  function(x, y = NULL) crossprod(as(x, "dgeMatrix"), y))

setMethod("crossprod", signature(x = "dsparseMatrix", y = "dgeMatrix"),
## NB: using   callGeneric(.) here, leads to infinite recursion :
          function(x, y = NULL) .Call(Csparse_dense_crossprod, as(x, "dgCMatrix"), y))

## NB: there's already
##     ("CsparseMatrix", "missing") and ("TsparseMatrix", "missing") methods

setMethod("crossprod", signature(x = "dgeMatrix", y = "dsparseMatrix"),
	  function(x, y = NULL) crossprod(x, as(y, "dgCMatrix")))

setMethod("image", "dsparseMatrix",
	  function(x, ...) image(as(x, "dgTMatrix"), ...))

setMethod("chol", signature(x = "dsparseMatrix"),
	   function(x, pivot=FALSE, ...) {
	       px <- as(x, "dsCMatrix")
	       if (isTRUE(validObject(px, test=TRUE))) chol(px, pivot, ...)
	       else stop("'x' is not positive definite -- chol() undefined.")
	   })

setMethod("lu", signature(x = "dsparseMatrix"),
	  function(x, ...) lu(as(x, "dgCMatrix")))


## Group Methods, see ?Arith (e.g.)
## -----

##-> now moved to ./Csparse.R (and 'up' to ./sparseMatrix.R):
##  "Math2" is in ./dMatrix.R


