#### Logical Sparse Triangular Matrices in Triplet format

### contains = "nsparseMatrix"

setAs("matrix", "ntTMatrix",
      function(from) as(as(from, "dtTMatrix"), "ntTMatrix"))

setAs("ntTMatrix", "ngTMatrix",
      function(from) new("ngTMatrix", i = from@i, j = from@j,
                         Dim = from@Dim, Dimnames = from@Dimnames))
setAs("ntTMatrix", "ntCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, TRUE))
setAs("ntTMatrix", "ngCMatrix",
      function(from) as(.Call(Tsparse_to_Csparse, from, TRUE), "ngCMatrix"))


setAs("ntTMatrix", "dtTMatrix",
      function(from) new("dtTMatrix", i = from@i, j = from@j,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ntTMatrix", "ntrMatrix",
      function(from) .Call(ntTMatrix_as_ntrMatrix, from))

setAs("ntTMatrix", "matrix",
      function(from) as(as(from, "ntrMatrix"), "matrix"))


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
