#### Sparse Symmetric non-zero pattern Matrices in Triplet format

### contains = "nsparseMatrix"

setAs("nsTMatrix", "matrix",
      function(from) as(as(from, "ngTMatrix"), "matrix"))

setAs("nsTMatrix", "ngTMatrix",
      function(from) new("ngTMatrix", i = from@i, j = from@j,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("nsTMatrix", "dsTMatrix",
      function(from) new("dsTMatrix", i = from@i, j = from@j,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

## untested:
setMethod("image", "nsTMatrix",
          function(x, ...) {
              x <- as(as(x, "dsTMatrix"), "dgTMatrix")
              callGeneric()
          })
