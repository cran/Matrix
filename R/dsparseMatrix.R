## METHODS FOR CLASS: dsparseMatrix (virtual)
## sparse matrices with 'x' slot of type "double"
## ... but _excluding_ ddiMatrix (FIXME?)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... now inherited from [CRT]sparseMatrix
if(FALSE) {
setMethod("symmpart", signature(x = "dsparseMatrix"),
          function(x) forceSymmetric(x + t(x)) / 2)
setMethod("skewpart", signature(x = "dsparseMatrix"),
          function(x) symmetrizeDimnames(x - t(x)) / 2)
} ## MJ

## MJ: no longer needed ... now inherited from Matrix
if(FALSE) {
setMethod("image", "dsparseMatrix",
	  function(x, ...) image(as(x, "dgTMatrix"), ...))
} ## MJ
