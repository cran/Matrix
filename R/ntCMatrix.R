#### Logical Sparse Triangular Matrices in Compressed column-oriented format

setAs("ntCMatrix", "matrix",
      function(from) as(as(from, "ngCMatrix"), "matrix"))
setAs("matrix", "ntCMatrix",
      function(from) as(as(from, "dtCMatrix"), "ntCMatrix"))

setAs("ntCMatrix", "ngCMatrix",
      function(from) new("ngCMatrix", i = from@i, p = from@p,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ntCMatrix", "dMatrix", # < instead of "dtCMatrix"
      function(from) new("dtCMatrix", i = from@i, p = from@p,
                         x = rep.int(1, length(from@i)), uplo = from@uplo,
                         diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

## setAs("ntCMatrix", "generalMatrix",
##       function(from) ......)

setMethod("image", "ntCMatrix",
          function(x, ...) {
              x <- as(as(x, "dtCMatrix"), "dgTMatrix")
              callGeneric()
          })

## setMethod("t", signature(x = "ntCMatrix"),
##           function(x) .Call(ntCMatrix_trans, x),
##           valueClass = "ntCMatrix")
