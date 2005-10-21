setMethod("crossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL)
	  .Call("Csparse_crossprod", x, trans = FALSE, triplet = FALSE,
		PACKAGE = "Matrix"))

setMethod("t", signature(x = "CsparseMatrix"),
	  function(x)
	  .Call("Csparse_transpose", x, PACKAGE = "Matrix"))

setMethod("tcrossprod", signature(x = "CsparseMatrix"),
	  function(x)
	  .Call("Csparse_crossprod", x, trans = TRUE, triplet = FALSE,
		PACKAGE = "Matrix"))
