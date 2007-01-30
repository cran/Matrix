#### Logical Sparse Triangular Matrices in Triplet format

### contains = "lsparseMatrix"

setAs("ltTMatrix", "lgTMatrix",
      function(from) new("lgTMatrix", i = from@i, j = from@j, x = from@x,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ltTMatrix", "ltCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, TRUE))
setAs("ltTMatrix", "lgCMatrix",
      function(from) as(.Call(Tsparse_to_Csparse, from, TRUE), "lgCMatrix"))

setAs("ltTMatrix", "dtTMatrix",
      function(from) new("dtTMatrix", i = from@i, j = from@j,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ltTMatrix", "ltrMatrix",
      function(from) .Call(ltTMatrix_as_ltrMatrix, from))

setAs("ltTMatrix", "matrix",
      function(from) as(as(from, "ltrMatrix"), "matrix"))



## untested:
setMethod("image", "ltTMatrix",
          function(x, ...) {
              x <- as(as(x, "dtTMatrix"), "dgTMatrix")
              callGeneric()
          })

## FIXME
## setMethod("t", signature(x = "ltTMatrix"),
##           function(x) .Call(ltTMatrix_trans, x),
##           valueClass = "ltTMatrix")
