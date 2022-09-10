#### Logical Symmetric Sparse Matrices in Compressed column-oriented format

### contains = "nsparseMatrix"

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
setAs("nsCMatrix", "matrix",
      function(from) as(as(from, "ngCMatrix"), "matrix"))

setAs("nsCMatrix", "ngCMatrix",
      function(from) .Call(Csparse_symmetric_to_general, from))

setAs("nsCMatrix", "nsTMatrix",
      function(from) .Call(Csparse_to_Tsparse, from, FALSE))

setAs("nsCMatrix", "dMatrix", .nC2d)
setAs("nsCMatrix", "dsparseMatrix", .nC2d)
setAs("nsCMatrix", "dsCMatrix", .nC2d)

setAs("nsCMatrix", "lMatrix", .nC2l)
setAs("nsCMatrix", "lsparseMatrix", .nC2l)
setAs("nsCMatrix", "lsCMatrix", .nC2l)
} ## MJ

## MJ: no longer needed ... methods now inherited from CsparseMatrix
if(FALSE) {
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
setMethod("t", signature(x = "nsCMatrix"),
          function(x) .Call(nsCMatrix_trans, x),
          valueClass = "nsCMatrix")
} ## MJ
