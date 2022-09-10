#### Logical Sparse Triangular Matrices in Compressed column-oriented format

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("matrix", "ltCMatrix",
      function(from) {
	  if(!is.logical(from)) storage.mode(from) <- "logical"
	  .Call(matrix_to_Csparse, from, "ltCMatrix")
      })
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
setAs("ltCMatrix", "matrix",
      function(from) as(as(from, "lgCMatrix"), "matrix"))

setAs("ltCMatrix", "lgCMatrix",
      function(from) copyClass(diagU2N(from), "lgCMatrix",
			       c("i", "p", "x", "Dim", "Dimnames")))

setAs("ltCMatrix", "ltTMatrix",
      function(from) .Call(Csparse_to_Tsparse, from, TRUE))

setAs("ltCMatrix", "dMatrix", # < instead of "dtCMatrix"
      function(from) new("dtCMatrix", i = from@i, p = from@p,
                         x = as.double(from@x), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "ltCMatrix", # to triangular {needed in triu() }
      function(from) as(as(as(from, "lgTMatrix"), "ltTMatrix"), "ltCMatrix"))

## setAs("ltCMatrix", "generalMatrix",
##       function(from) ......)
} ## MJ

## MJ: no longer needed ... method now inherited from CsparseMatrix
if(FALSE) {
setMethod("t", signature(x = "ltCMatrix"),
          function(x) .Call(ltCMatrix_trans, x),
          valueClass = "ltCMatrix")
} ## MJ
