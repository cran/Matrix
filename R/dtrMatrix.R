## METHODS FOR CLASS: dtrMatrix
## dense (unpacked) triangular matrices with 'x' slot of type "double"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("dgeMatrix", "dtrMatrix", function(from) asTri(from, "dtrMatrix"))
setAs("dtrMatrix", "dtpMatrix",
      dtr2dtp <- function(from) .Call(dtrMatrix_as_dtpMatrix, from))
setAs("matrix", "dtrMatrix",
      function(from) as(..2dge(from), "dtrMatrix"))

.dtr2mat <- function(from, keep.dimnames=TRUE)
    .Call(dtrMatrix_as_matrix, from, keep.dimnames)
## needed for t() method
setAs("dtrMatrix", "matrix",
      function(from) .Call(dtrMatrix_as_matrix, from, TRUE))
setAs("Cholesky", "lMatrix",
      function(from) as(as(from, "dtrMatrix"), "lMatrix"))
setAs("BunchKaufman", "lMatrix",
      function(from) as(as(from, "dtrMatrix"), "lMatrix"))
setAs("dtrMatrix", "sparseMatrix",
      function(from) .dense2C(from, kind = "tri", uplo = from@uplo))
setAs("dtrMatrix", "CsparseMatrix",
      function(from) .dense2C(from, kind = "tri", uplo = from@uplo))
} ## MJ

## MJ: no longer needed ... replacement in ./unpackedMatrix.R
if (FALSE) {
setMethod("t", signature(x = "dtrMatrix"), t_trMatrix)
setMethod("diag", signature(x = "dtrMatrix"),
          .mkSpec.diag(quote(dtrMatrix_getDiag)),
          valueClass = "numeric")
setMethod("diag<-", signature(x = "dtrMatrix"),
	  function(x, value) {
	      .Call(dtrMatrix_setDiag,
		    if(x@diag == "U") .dense.diagU2N(x, "d", isPacked=FALSE) else x,
		    value)
	  })
} ## MJ
