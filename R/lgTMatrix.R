#### Logical Sparse Matrices in triplet format

### contains = "lsparseMatrix"
###             ============= ---> superclass methods in ./lsparseMatrix.R

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("matrix", "lgTMatrix",
      function(from) {
	  stopifnot(is.logical(from))
	  dn <- dimnames(from)
	  if(is.null.DN(dn))
	      dn <- list(NULL,NULL)
	  else dimnames(from) <- NULL
	  TorNA <- is.na(from) | from
	  ij <- which(TorNA, arr.ind = TRUE, useNames = FALSE) - 1L
	  if(length(ij) == 0) ij <- matrix(ij, 0, 2)
	  new("lgTMatrix",
	      i = ij[,1],
	      j = ij[,2],
	      x = from[TorNA],
	      Dim = as.integer(dim(from)),
	      Dimnames = dn)
	  })
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
setAs("lgTMatrix", "lgeMatrix",
      function(from) .Call(lgTMatrix_to_lgeMatrix, from))

setAs("lgTMatrix", "matrix",
      function(from) .Call(lgTMatrix_to_matrix, from))

setAs("lgTMatrix", "dgTMatrix",
      function(from)
      ## more efficient than
      ## as(as(as(sM, "lgCMatrix"), "dgCMatrix"), "dgTMatrix")
      new("dgTMatrix", i = from@i, j = from@j,
          x = as.double(from@x),
          ## cannot copy factors, but can we use them?
          Dim = from@Dim, Dimnames= from@Dimnames))

setAs("lgTMatrix", "ltTMatrix",
      function(from) check.gT2tT(from, toClass = "ltTMatrix", do.n=FALSE))
} ## MJ

## MJ: no longer needed ... methods now inherited from Matrix
if(FALSE) {
setAs("lgTMatrix", "triangularMatrix",
      function(from) check.gT2tT(from, toClass = "ltTMatrix", do.n=FALSE))
## Favor coercion to superclasses: here "symmetricMatrix" not "lsTMatrix"
setAs("lgTMatrix", "symmetricMatrix",
      function(from) check.gT2sT(from, toClass = "lsTMatrix", do.n=FALSE))
} ## MJ

## MJ: no longer needed ... method now inherited from TsparseMatrix
if(FALSE) {
setMethod("t", signature(x = "lgTMatrix"),
	  function(x) new("lgTMatrix", i = x@j, j = x@i, x = x@x,
			  Dim = x@Dim[2:1],
			  Dimnames= x@Dimnames[2:1]),
	  valueClass = "lgTMatrix")
} ## MJ
