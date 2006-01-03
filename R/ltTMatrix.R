#### Logical Sparse Triangular Matrices in Triplet format

### contains = "lsparseMatrix"

setAs("ltTMatrix", "matrix",
      function(from) as(as(from, "dtTMatrix"), "matrix"))

setAs("ltTMatrix", "dtTMatrix",
      function(from) new("dtTMatrix", i = from@i, j = from@j,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

## untested:
setMethod("image", "ltTMatrix",
          function(x, ...) {
              x <- as(as(x, "dtTMatrix"), "dgTMatrix")
              callGeneric()
          })

## FIXME
## setMethod("t", signature(x = "ltTMatrix"),
##           function(x) .Call("ltTMatrix_trans", x, PACKAGE = "Matrix"),
##           valueClass = "ltTMatrix")
