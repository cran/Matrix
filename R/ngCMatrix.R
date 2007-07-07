#### Logical Sparse Matrices in Compressed column-oriented format

### contains = "nsparseMatrix"

.ngC2d <- function(from)
    new("dgCMatrix", i = from@i, p = from@p, x = rep.int(1, length(from@i)),
	Dim = from@Dim, Dimnames = from@Dimnames)

.ngC2l <- function(from)
    new("lgCMatrix", i = from@i, p = from@p, x = rep.int(TRUE, length(from@i)),
	Dim = from@Dim, Dimnames = from@Dimnames)

setAs("ngCMatrix", "dMatrix", .ngC2d)# < instead of "dgCMatrix"
setAs("ngCMatrix", "dsparseMatrix", .ngC2d)
setAs("ngCMatrix", "dgCMatrix", .ngC2d)

setAs("ngCMatrix", "lMatrix", .ngC2l)
setAs("ngCMatrix", "lsparseMatrix", .ngC2l)
setAs("ngCMatrix", "lgCMatrix", .ngC2l)

rm(.ngC2d,.ngC2l) # don't even keep "hidden"

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
## not this: .Call(Csparse_to_matrix, from)), since it goes via dense -> double precision

## TODO (maybe): write  matrix_to_lcsc()  in ../src/ngCMatrix.c
setAs("matrix", "ngCMatrix",
      function(from) as(as(from, "ngTMatrix"), "ngCMatrix"))



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

setMethod("image", "ngCMatrix",
          function(x, ...) {
              x <- as(as(x, "dgCMatrix"), "dgTMatrix")
              callGeneric()
          })
