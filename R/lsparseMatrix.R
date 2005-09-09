#### Superclass Methods for all sparse logical matrices

###------- Work via  as(*, lgC) : ------------

setMethod("%*%", signature(x = "lsparseMatrix", y = "lsparseMatrix"),
          function(x, y) callGeneric(as(x, "lgCMatrix"), as(y, "lgCMatrix")))

setMethod("crossprod", signature(x = "lsparseMatrix", y = "missing"),
	  function(x, y = NULL)
          .Call("lgCMatrix_crossprod", as(x, "lgCMatrix"), TRUE, NULL),
	  valueClass = "lsCMatrix")

setMethod("tcrossprod", signature(x = "lsparseMatrix"),
	  function(x)
          .Call("lgCMatrix_crossprod", as(x, "lgCMatrix"), FALSE, NULL),
	  valueClass = "lsCMatrix")

