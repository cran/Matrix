#### Sparse Symmetric non-zero pattern Matrices in Triplet format

### contains = "nsparseMatrix"

setAs("nsTMatrix", "matrix",
      function(from) as(as(from, "ngTMatrix"), "matrix"))

setAs("nsTMatrix", "ngCMatrix", # for diag
      function(from) as(as(from, "nsCMatrix"), "ngCMatrix"))

setAs("nsTMatrix", "ngTMatrix",
      function(from) .Call(nsTMatrix_as_ngTMatrix, from))

setAs("nsTMatrix", "dsTMatrix",
      function(from) new("dsTMatrix", i = from@i, j = from@j,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("nsTMatrix", "nsyMatrix",
      function(from) .Call(nsTMatrix_as_nsyMatrix, from))


## untested:
setMethod("image", "nsTMatrix",
          function(x, ...) {
              x <- as(as(x, "dsTMatrix"), "dgTMatrix")
              callGeneric()
          })
