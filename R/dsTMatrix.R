### Coercion and Methods for Symmetric Triplet Matrices

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
## Now in ./Tsparse.R
## setAs("dsTMatrix", "dsCMatrix",
##       function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("dsTMatrix", "dgTMatrix",
      function(from) .Call(dsTMatrix_as_dgTMatrix, from))

setAs("dsTMatrix", "lsTMatrix",
      function(from) new("lsTMatrix", x = as.logical(from@x),
                         Dim = from@Dim, Dimnames = from@Dimnames,
                         uplo = from@uplo, i = from@i, j = from@j))

## Conversion <--> dense storage is via dsyMatrix :
setAs("dsTMatrix", "dsyMatrix",
      function(from) .Call(dsTMatrix_as_dsyMatrix, from))

setAs("dsTMatrix", "dgeMatrix",
      function(from) as(as(from, "dsyMatrix"), "dgeMatrix"))
setAs("dsTMatrix", "matrix",
      function(from) as(as(from, "dsyMatrix"), "matrix"))
} ## MJ

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
to_dsT <- function(from) as(as(from, "dsyMatrix"), "TsparseMatrix")
setAs("dgeMatrix", "dsTMatrix", to_dsT)
setAs("matrix", "dsTMatrix", to_dsT)
rm(to_dsT)
} ## MJ

## MJ: no longer needed ... method now inherited from TsparseMatrix
if(FALSE) {
setMethod("t", "dsTMatrix",
	  function(x)
	  new("dsTMatrix", Dim = x@Dim, Dimnames = x@Dimnames[2:1],
	      i = x@j, j = x@i, x = x@x,
	      uplo = if (x@uplo == "U") "L" else "U"))
} ## MJ

## setMethod("writeHB", signature(obj = "dsTMatrix"),
##           function(obj, file, ...) callGeneric(as(obj, "CsparseMatrix"), file, ...))
