#### Logical Symmetric Sparse Matrices in Compressed column-oriented format

### contains = "lsparseMatrix"

setAs("lsCMatrix", "dsCMatrix",
      function(from) new("dsCMatrix", i = from@i, p = from@p,
                         x = rep(1, length(from@i)), uplo = from@uplo,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setMethod("image", "lsCMatrix",
          function(x, ...) {
              x <- as(as(x, "dsCMatrix"), "dgTMatrix")
              callGeneric()
          })

setMethod("chol", signature(x = "lsCMatrix", pivot = "missing"),
          function(x, pivot, LINPACK) .Call("lsCMatrix_chol", x, TRUE))

setMethod("chol", signature(x = "lsCMatrix", pivot = "logical"),
          function(x, pivot, LINPACK) .Call("lsCMatrix_chol", x, pivot))

setMethod("t", signature(x = "lsCMatrix"),
          function(x) .Call("lsCMatrix_trans", x),
          valueClass = "lsCMatrix")
