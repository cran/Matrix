#### Logical Sparse Triangular Matrices in Triplet format

### contains = "nsparseMatrix"

setAs("ntTMatrix", "matrix",
      function(from) as(as(from, "ngTMatrix"), "matrix"))

setAs("ntTMatrix", "ngTMatrix",
      function(from) new("ngTMatrix", i = from@i, j = from@j, x = from@x,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ntTMatrix", "dtTMatrix",
      function(from) new("dtTMatrix", i = from@i, j = from@j,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

## untested:
setMethod("image", "ntTMatrix",
          function(x, ...) {
              x <- as(as(x, "dtTMatrix"), "dgTMatrix")
              callGeneric()
          })

## FIXME
## setMethod("t", signature(x = "ntTMatrix"),
##           function(x) .Call(ntTMatrix_trans, x),
##           valueClass = "ntTMatrix")
