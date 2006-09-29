#### Logical Sparse Symmetric Matrices in Triplet format

### contains = "lsparseMatrix"

setAs("lsTMatrix", "matrix",
      function(from) as(as(from, "lgTMatrix"), "matrix"))

setAs("lsTMatrix", "lgTMatrix",
      function(from) new("lgTMatrix", i = from@i, j = from@j, x = from@x,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lsTMatrix", "dsTMatrix",
      function(from) new("dsTMatrix", i = from@i, j = from@j,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

## untested:
setMethod("image", "lsTMatrix",
          function(x, ...) {
              x <- as(as(x, "dsTMatrix"), "dgTMatrix")
              callGeneric()
          })
