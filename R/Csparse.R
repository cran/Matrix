setMethod("crossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      a <- .Call(Csparse_crossprod, x, trans = FALSE, triplet = FALSE,
			 PACKAGE = "Matrix")
	      switch(substr(class(a)[1], 1, 1),
		     "d" ={ new("dsCMatrix", i = a@i, p = a@p, x = a@x,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "U",
				factors = list()) },
		     "l" ={ new("lsCMatrix", i = a@i, p = a@p,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "U",
				factors = list()) })
	  })


setMethod("t", signature(x = "CsparseMatrix"),
	  function(x)
	  .Call(Csparse_transpose, x))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      a <- .Call(Csparse_crossprod, x, trans = TRUE, triplet = FALSE,
			 PACKAGE = "Matrix")
	      switch(substr(class(a)[1], 1, 1),
		     "d" ={ new("dsCMatrix", i = a@i, p = a@p, x = a@x,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "L",
				factors = list()) },
		     "l" ={ new("lsCMatrix", i = a@i, p = a@p,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "L",
				factors = list()) })
	  })

## FIXME (TODO):
## setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
## 	  function(x, y)
## 	  .Call(Csparse_crossprod_2, x, y, trans = TRUE, triplet = FALSE,
## 		PACKAGE = "Matrix"))


setMethod("%*%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y) .Call(Csparse_Csparse_prod, x, y))

setMethod("%*%", signature(x = "CsparseMatrix", y = "denseMatrix"),
          function(x, y) .Call(Csparse_dense_prod, x, y))


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
          .Call(Csparse_to_Tsparse, from)
      })

setAs("CsparseMatrix", "denseMatrix",
      function(from) {
	  if(!is(from, "generalMatrix")) { ## e.g. for triangular | symmetric
	      if     (is(from, "dMatrix")) from <- as(from, "dgCMatrix")
	      else if(is(from, "lMatrix")) from <- as(from, "lgCMatrix")
	      else if(is(from, "zMatrix")) from <- as(from, "zgCMatrix")
	      else stop("undefined method for class ", class(from))
	  }
          .Call(Csparse_to_dense, from)
      })

setMethod("tril", "CsparseMatrix",
          function(x, k = 0, ...) {
              k <- as.integer(k[1])
              dd <- dim(x)
              stopifnot(k >= -dd[1], k <= 0)
              .Call(Csparse_band, x, -(dim(x)[1]), k)
          })

setMethod("triu", "CsparseMatrix",
          function(x, k = 0, ...) {
              k <- as.integer(k[1])
              dd <- dim(x)
              stopifnot(k >= 0, k <= dd[1])
              .Call(Csparse_band, x, k, dd[2])
          })

setMethod("band", "CsparseMatrix",
          function(x, k1, k2, ...) {
              k1 <- as.integer(k1[1])
              k2 <- as.integer(k2[1])
              dd <- dim(x)
              stopifnot(k1 <= k2, k1 >= -dd[1], k2 <= dd[1])
              .Call(Csparse_band, x, k1, k2)
          })

setMethod("colSums", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("colMeans", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("rowSums", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("rowMeans", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
