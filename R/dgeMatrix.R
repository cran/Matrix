setAs("matrix", "dgeMatrix",
      function(from) .Call(dup_mMatrix_as_dgeMatrix, from))
if(FALSE)## the above is MUCH faster than
setAs("matrix", "dgeMatrix",
      function(from) {
	  new("dgeMatrix",
	      x = as.double(from),
	      Dim = as.integer(dim(from)),
	      Dimnames = .M.DN(from))
      })

setAs("dgeMatrix", "matrix",
      function(from) {
	  array(from@x, dim = from@Dim, dimnames = from@Dimnames)
      })

## Group Methods, see ?Math (e.g.)

##  "Arith" is in ./Ops.R

setMethod("Math",
	  signature(x = "dgeMatrix"),
	  function(x) {
	      x@x <- callGeneric(x@x)
	      x
	  })

##  "Math2" is in ./dMatrix.R

## "Summary"

## "Compare"  now happens in ./dMatrix.R

## -- end{group generics} -----------------------


##  "[" settings are "up in"  Matrix.R & denseMatrix.R


setMethod("as.vector", signature(x = "dgeMatrix", mode = "missing"),
          function(x) x@x)

setMethod("norm", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) norm(x, type = "O", ...))

setMethod("norm", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
	  .Call(dgeMatrix_norm, x, type),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) rcond(x, type = "O", ...))

setMethod("rcond", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
	  .Call(dgeMatrix_rcond, x, type),
	  valueClass = "numeric")

setMethod("t", signature(x = "dgeMatrix"), t_geMatrix)

## crossprod(x) & tcrossprod(x) :
setMethod("crossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call(dgeMatrix_crossprod, x, FALSE),
	  valueClass = "dpoMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call(dgeMatrix_crossprod, x, TRUE),
	  valueClass = "dpoMatrix")

if(FALSE) { ## this would mask 'base::tcrossprod'
setMethod("tcrossprod", signature(x = "matrix", y = "missing"),
	  function(x, y = NULL)
          .Call(dgeMatrix_crossprod, as(x, "dgeMatrix"), TRUE),
	  valueClass = "dpoMatrix")

setMethod("tcrossprod", signature(x = "numeric", y = "missing"),
	  function(x, y = NULL) callGeneric(as.matrix(as.double(x))))
}

## crossprod (x,y)
setMethod("crossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call(dgeMatrix_dgeMatrix_crossprod, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call(dgeMatrix_matrix_crossprod, x, y, FALSE),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL)
	  .Call(dgeMatrix_matrix_crossprod, x, as.matrix(as.double(y)), FALSE),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as.matrix(as.double(x)), y),
	  valueClass = "dgeMatrix")

## tcrossprod (x,y)
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call(dgeMatrix_dgeMatrix_crossprod, x, y, TRUE),
	  valueClass = "dgeMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call(dgeMatrix_matrix_crossprod, x, y, TRUE),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL)
	  .Call(dgeMatrix_matrix_crossprod, x, rbind(as.double(y)), TRUE),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(rbind(as.double(x)), y),
	  valueClass = "dgeMatrix")

## %*% methods
setMethod("%*%", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, x, y, FALSE),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, y, x, TRUE),
          valueClass = "dgeMatrix")

## DB: Should we retain these methods?  Does the shortcut save enough
## to justify additional signatures?
## dgeMatrix <-> numeric: conceptually dispatch to "matrix" one, but shortcut
setMethod("%*%", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y)
	  .Call(dgeMatrix_matrix_mm, y, rbind(x), TRUE),
	  valueClass = "dgeMatrix")

setMethod("diag", signature(x = "dgeMatrix"),
	  function(x, nrow, ncol = n) .Call(dgeMatrix_getDiag, x))

setMethod("chol", signature(x = "dgeMatrix", pivot = "ANY"), cholMat)

setMethod("solve", signature(a = "dgeMatrix", b = "missing"),
	  function(a, b, ...) .Call(dgeMatrix_solve, a),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "ddenseMatrix"),
	  function(a, b, ...) .Call(dgeMatrix_matrix_solve, a, b),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "matrix"),
	  function(a, b, ...) .Call(dgeMatrix_matrix_solve, a, b),
          valueClass = "dgeMatrix")

## not needed - method for numeric defined for Matrix class
## setMethod("solve", signature(a = "dgeMatrix", b = "numeric"),
## 	  function(a, b, ...)
## 	  .Call(dgeMatrix_matrix_solve, a, as.matrix(as.double(b))))

setMethod("lu", signature(x = "dgeMatrix"),
	  function(x, ...) .Call(dgeMatrix_LU, x),
	  valueClass = "denseLU")

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "missing"),
	  function(x, logarithm, ...)
	  .Call(dgeMatrix_determinant, x, TRUE))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
	  function(x, logarithm, ...)
	  .Call(dgeMatrix_determinant, x, logarithm))


setMethod("expm", signature(x = "dgeMatrix"),
	  function(x) .Call(dgeMatrix_exp, x),
	  valueClass = "dgeMatrix")

##-> ./colSums.R  for colSums,... rowMeans
