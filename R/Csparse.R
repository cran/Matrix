setMethod("crossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL)
	  .Call("Csparse_crossprod", x, trans = FALSE, triplet = FALSE,
		PACKAGE = "Matrix"))

setMethod("t", signature(x = "CsparseMatrix"),
	  function(x)
	  .Call("Csparse_transpose", x, PACKAGE = "Matrix"))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL)
	  .Call("Csparse_crossprod", x, trans = TRUE, triplet = FALSE,
		PACKAGE = "Matrix"))
## FIXME (TODO):
## setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
## 	  function(x, y)
## 	  .Call("Csparse_crossprod_2", x, y, trans = TRUE, triplet = FALSE,
## 		PACKAGE = "Matrix"))


setMethod("%*%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y) .Call("Csparse_Csparse_prod", x, y, PACKAGE = "Matrix"))

setMethod("%*%", signature(x = "CsparseMatrix", y = "denseMatrix"),
          function(x, y) .Call("Csparse_dense_prod", x, y, PACKAGE = "Matrix"))


## FIXME: the is(*,"generalMatrix") test at least makes these work,
##        but they are still ``wrong'', since triangularity is lost  :

setAs("CsparseMatrix", "TsparseMatrix",
      function(from) {
	  if(!is(from, "generalMatrix")) { ## e.g. for triangular | symmetric
	      if     (is(from, "dMatrix")) from <- as(from, "dgCMatrix")
	      else if(is(from, "lMatrix")) from <- as(from, "lgCMatrix")
	      else if(is(from, "zMatrix")) from <- as(from, "zgCMatrix")
	      else stop("undefined method for class ", class(from))
	  }
          .Call("Csparse_to_Tsparse", from, PACKAGE = "Matrix")
      })

setAs("CsparseMatrix", "denseMatrix",
      function(from) {
	  if(!is(from, "generalMatrix")) { ## e.g. for triangular | symmetric
	      if     (is(from, "dMatrix")) from <- as(from, "dgCMatrix")
	      else if(is(from, "lMatrix")) from <- as(from, "lgCMatrix")
	      else if(is(from, "zMatrix")) from <- as(from, "zgCMatrix")
	      else stop("undefined method for class ", class(from))
	  }
          .Call("Csparse_to_dense", from, PACKAGE = "Matrix")
      })
