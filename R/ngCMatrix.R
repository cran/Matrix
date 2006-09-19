#### Logical Sparse Matrices in Compressed column-oriented format

### contains = "nsparseMatrix"

## Can use CsparseMatrix methods for all of these

## setMethod("%*%", signature(x = "ngCMatrix", y = "ngCMatrix"),
##           function(x, y)
##           .Call(ngCMatrix_ngCMatrix_mm, x, y),
##           valueClass = "ngCMatrix")

## setMethod("t", signature(x = "ngCMatrix"),
##           function(x) .Call(ngCMatrix_trans, x),
##           valueClass = "ngCMatrix")


## setMethod("diag", signature(x = "ngCMatrix"),
## 	  function(x, nrow, ncol = n) .Call(ngCMatrix_diag, x))

## setMethod("crossprod", signature(x = "ngCMatrix", y = "missing"),
## 	  function(x, y = NULL)
##           .Call(ngCMatrix_crossprod, x, TRUE, NULL),
## 	  valueClass = "nsCMatrix")

## setMethod("tcrossprod", signature(x = "ngCMatrix", y = "missing"),
## 	  function(x, y = NULL)
##           .Call(ngCMatrix_crossprod, x, FALSE, NULL),
## 	  valueClass = "nsCMatrix")

setAs("ngCMatrix", "dgCMatrix",
      function(from) new("dgCMatrix", i = from@i, p = from@p,
                         x = rep(1, length(from@i)),
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ngCMatrix", "ngTMatrix",
      function(from) new("ngTMatrix", i = from@i,
                         j = .Call(Matrix_expand_pointers, from@p),
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ngCMatrix", "ngeMatrix",
      function(from)
	  new("ngeMatrix", x = c(as(from, "matrix")), # is fast,
	      Dim = from@Dim, Dimnames = from@Dimnames))

setAs("ngCMatrix", "matrix",
      function(from) .Call(ncsc_to_matrix, from))
## TODO (maybe): write  matrix_to_lcsc()  in ../src/ngCMatrix.c
setAs("matrix", "ngCMatrix",
      function(from) as(as(from, "ngTMatrix"), "ngCMatrix"))


setMethod("image", "ngCMatrix",
          function(x, ...) {
              x <- as(as(x, "dgCMatrix"), "dgTMatrix")
              callGeneric()
          })
