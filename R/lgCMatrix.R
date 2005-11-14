#### Logical Sparse Matrices in Compressed column-oriented format

### contains = "lsparseMatrix"

setMethod("%*%", signature(x = "lgCMatrix", y = "lgCMatrix"),
          function(x, y)
          .Call("lgCMatrix_lgCMatrix_mm", x, y, PACKAGE = "Matrix"),
          valueClass = "lgCMatrix")

setMethod("t", signature(x = "lgCMatrix"),
          function(x) .Call("lgCMatrix_trans", x, PACKAGE = "Matrix"),
          valueClass = "lgCMatrix")

setMethod("crossprod", signature(x = "lgCMatrix", y = "missing"),
	  function(x, y = NULL)
          .Call("lgCMatrix_crossprod", x, TRUE, NULL, PACKAGE = "Matrix"),
	  valueClass = "lsCMatrix")

setMethod("tcrossprod", signature(x = "lgCMatrix"),
	  function(x)
          .Call("lgCMatrix_crossprod", x, FALSE, NULL, PACKAGE = "Matrix"),
	  valueClass = "lsCMatrix")

setAs("lgCMatrix", "dgCMatrix",
      function(from) new("dgCMatrix", i = from@i, p = from@p,
                         x = rep(1, length(from@i)),
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "lgTMatrix",
      function(from) new("lgTMatrix", i = from@i,
                         j = .Call("Matrix_expand_pointers", from@p, PACKAGE = "Matrix"),
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "matrix",
      function(from) .Call("lcsc_to_matrix", from, PACKAGE = "Matrix"))


setMethod("image", "lgCMatrix",
          function(x, ...) {
              x <- as(as(x, "dgCMatrix"), "dgTMatrix")
              callGeneric()
          })

