#### Logical Sparse Triangular Matrices in Compressed column-oriented format

### contains = "lsparseMatrix"

setAs("ltCMatrix", "dtCMatrix",
      function(from) new("dtCMatrix", i = from@i, p = from@p,
                         x = rep(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setMethod("image", "ltCMatrix",
          function(x, ...) {
              x <- as(as(x, "dtCMatrix"), "dgTMatrix")
              callGeneric()
          })

setMethod("t", signature(x = "ltCMatrix"),
          function(x) .Call("ltCMatrix_trans", x),
          valueClass = "ltCMatrix")
