### Coercion and Methods for Triangular Triplet Matrices

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("matrix", "dtTMatrix",
      function(from) as(as(from, "dtpMatrix"), "dtTMatrix"))
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
setAs("dtTMatrix", "dgTMatrix",
      function(from) tT2gT(from, cl = "dtTMatrix", toClass = "dgTMatrix"))
setAs("dtTMatrix", "generalMatrix",
      function(from) tT2gT(from, cl = "dtTMatrix", toClass = "dgTMatrix"))

setAs("dtTMatrix", "ltTMatrix",
      function(from) new("ltTMatrix", i = from@i, j = from@j,
                         x = as.logical(from@x),
                         uplo = from@uplo, diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))
setAs("dtTMatrix", "ntTMatrix",
      function(from) new("ntTMatrix", i = from@i, j = from@j,
                         uplo = from@uplo, diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

## Conversion to dense storage is first to a dtrMatrix
setAs("dtTMatrix", "dtrMatrix",
      function(from) .Call(dtTMatrix_as_dtrMatrix, from))

setAs("dtTMatrix", "matrix",
      function(from) as(as(from, "dtrMatrix"), "matrix"))

setAs("dtTMatrix", "dgeMatrix",
      function(from) as(as(from, "dtrMatrix"), "dgeMatrix"))
} ## MJ

## MJ: no longer needed ... method now inherited from TsparseMatrix
if(FALSE) {
setMethod("t", "dtTMatrix",
	  function(x)
	  new("dtTMatrix", Dim = x@Dim[2:1], Dimnames = x@Dimnames[2:1],
	      i = x@j, j = x@i, x = x@x, diag = x@diag,
	      uplo = if (x@uplo == "U") "L" else "U"))
} ## MJ
