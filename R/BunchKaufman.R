## METHODS FOR GENERIC: BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("BunchKaufman", signature(x = "dsyMatrix"),
	  function(x, ...) .Call(dsyMatrix_trf, x))

setMethod("BunchKaufman", signature(x = "dspMatrix"),
	  function(x, ...) .Call(dspMatrix_trf, x))

setMethod("BunchKaufman", signature(x = "matrix"),
	  function(x, uplo = NULL, ...) .Call(matrix_trf, x, uplo))
