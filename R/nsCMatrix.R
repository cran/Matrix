#### Logical Symmetric Sparse Matrices in Compressed column-oriented format

### contains = "nsparseMatrix"

setAs("nsCMatrix", "matrix",
      function(from) as(as(from, "ngCMatrix"), "matrix"))

setAs("nsCMatrix", "ngCMatrix",
      function(from) .Call(Csparse_symmetric_to_general, from))

## Specific conversions, should they be necessary.  Better to convert as
## as(x, "TsparseMatrix") or as(x, "denseMatrix")
setAs("nsCMatrix", "nsTMatrix",
      function(from) .Call(Csparse_to_Tsparse, from, FALSE))

.nsC2d <- function(from)
    new("dsCMatrix", i = from@i, p = from@p,
	x = rep.int(1, length(from@i)), uplo = from@uplo,
	Dim = from@Dim, Dimnames = from@Dimnames)

.nsC2l <- function(from)
    new("lsCMatrix", i = from@i, p = from@p,
	x = rep.int(TRUE, length(from@i)), uplo = from@uplo,
	Dim = from@Dim, Dimnames = from@Dimnames)

setAs("nsCMatrix", "dsCMatrix", .nsC2d)
setAs("nsCMatrix", "dsparseMatrix", .nsC2d)

setAs("nsCMatrix", "lsCMatrix", .nsC2l)
setAs("nsCMatrix", "lsparseMatrix", .nsC2l)

rm(.nsC2d,.nsC2l) # don't even keep "hidden"

setAs("nsCMatrix", "dgTMatrix",
      function(from) as(as(x, "dsCMatrix"), "dgTMatrix"))

## have rather tril() and triu() methods than
## setAs("nsCMatrix", "ntCMatrix", ....)
setMethod("tril", "nsCMatrix",
	  function(x, k = 0, ...) {
	      if(x@uplo == "L" && k == 0)
		  ## same internal structure + diag
		  new("ntCMatrix", uplo = x@uplo, i = x@i, p = x@p,
		      Dim = x@Dim, Dimnames = x@Dimnames)
	      else tril(as(x, "ngCMatrix"), k = k, ...)
	  })
setMethod("triu", "nsCMatrix",
	  function(x, k = 0, ...) {
	      if(x@uplo == "U" && k == 0)
		  new("ntCMatrix", uplo = x@uplo, i = x@i, p = x@p,
		      Dim = x@Dim, Dimnames = x@Dimnames)
	      else triu(as(x, "ngCMatrix"), k = k, ...)
	  })

## FIXME: generalize to "nsparseMatrix" or (class union)  "symmetric sparse"
setMethod("image", "nsCMatrix",
	  function(x, ...) image(as(as(x, "dsCMatrix"), "dgTMatrix"), ...))


setMethod("chol", signature(x = "nsCMatrix", pivot = "missing"),
	  function(x, pivot, ...) chol(x, pivot = FALSE))
##          .Call(nsCMatrix_chol, x, FALSE))

setMethod("chol", signature(x = "nsCMatrix", pivot = "logical"),
	  function(x, pivot, ...) stop("temporarily disabled"))## FIXME
##          .Call(nsCMatrix_chol, x, pivot))

## Use more general method from CsparseMatrix class
## setMethod("t", signature(x = "nsCMatrix"),
##           function(x)
##           .Call(nsCMatrix_trans, x),
##           valueClass = "nsCMatrix")
