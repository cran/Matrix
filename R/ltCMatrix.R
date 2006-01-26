#### Logical Sparse Triangular Matrices in Compressed column-oriented format

setAs("ltCMatrix", "matrix",
      function(from) as(as(from, "lgCMatrix"), "matrix"))
setAs("matrix", "ltCMatrix",
      function(from) as(as(from, "dtCMatrix"), "ltCMatrix"))

setAs("ltCMatrix", "lgCMatrix",
      function(from) new("lgCMatrix", i = from@i, p = from@p,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ltCMatrix", "dMatrix", # < instead of "dtCMatrix"
      function(from) new("dtCMatrix", i = from@i, p = from@p,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

## setAs("ltCMatrix", "generalMatrix",
##       function(from) ......)

setMethod("image", "ltCMatrix",
          function(x, ...) {
              x <- as(as(x, "dtCMatrix"), "dgTMatrix")
              callGeneric()
          })

setMethod("t", signature(x = "ltCMatrix"),
          function(x) .Call("ltCMatrix_trans", x, PACKAGE = "Matrix"),
          valueClass = "ltCMatrix")
