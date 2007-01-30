#### Logical Sparse Symmetric Matrices in Triplet format

### contains = "lsparseMatrix"

setAs("lsTMatrix", "matrix",
      function(from) as(as(from, "lgTMatrix"), "matrix"))

setAs("lsTMatrix", "lgCMatrix", # for diag
      function(from) as(as(from, "lsCMatrix"), "lgCMatrix"))

setAs("lsTMatrix", "lgTMatrix",
      function(from) .Call(lsTMatrix_as_lgTMatrix, from))


setAs("lsTMatrix", "dsTMatrix",
      function(from) new("dsTMatrix", i = from@i, j = from@j,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lsTMatrix", "lsyMatrix",
      function(from) .Call(lsTMatrix_as_lsyMatrix, from))

## untested:
setMethod("image", "lsTMatrix",
          function(x, ...) {
              x <- as(as(x, "dsTMatrix"), "dgTMatrix")
              callGeneric()
          })
