 ### Coercion and Methods for Symmetric Matrices

setAs("dsyMatrix", "dgeMatrix",
      function(from) .Call("dsyMatrix_as_dgeMatrix", from) )

setAs("dsyMatrix", "matrix",
      function(from) .Call("dsyMatrix_as_matrix", from) )

setAs("dsyMatrix", "dspMatrix",
      function(from) .Call("dsyMatrix_as_dspMatrix", from) )

setMethod("rcond", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...)
          .Call("dsyMatrix_rcond", x, type),
          valueClass = "numeric")

setMethod("rcond", signature(x = "dsyMatrix", type = "missing"),
          function(x, type, ...)
          .Call("dsyMatrix_rcond", x, "O"),
          valueClass = "numeric")

setMethod("%*%", signature(x = "dsyMatrix", y = "dgeMatrix"),
          function(x, y) .Call("dsyMatrix_dgeMatrix_mm", x, y) )

setMethod("%*%", signature(x = "dgeMatrix", y = "dsyMatrix"),
          function(x, y) .Call("dsyMatrix_dgeMatrix_mm_R", y, x) )

setMethod("solve", signature(a = "dsyMatrix", b = "missing"),
          function(a, b, ...) .Call("dsyMatrix_solve", a),
          valueClass = "dsyMatrix")

setMethod("solve", signature(a = "dsyMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("dsyMatrix_matrix_solve", a, b),
          valueClass = "matrix")

setMethod("solve", signature(a = "dsyMatrix", b = "dgeMatrix"),
          function(a, b, ...)
          .Call("dsyMatrix_dgeMatrix_solve", a, b),
          valueClass = "dgeMatrix")

setMethod("norm", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...) .Call("dsyMatrix_norm", x, type),
          valueClass = "numeric")

setMethod("norm", signature(x = "dsyMatrix", type = "missing"),
          function(x, type, ...) .Call("dsyMatrix_norm", x, "O"),
          valueClass = "numeric")

## Should this create the opposite storage format - i.e. "U" -> "L"
## and vice-versa?
## MM: I think yes, since the other part can be filled arbitrarily (wrongly)
##WAS setMethod("t", signature(x = "dsyMatrix"), function(x) x)
setMethod("t", signature(x = "dsyMatrix"),
	  function(x) {
	      new("dsyMatrix",
                  Dim = x@Dim[2:1], Dimnames = x@Dimnames[2:1],
                  x = as.vector(t(as(x, "matrix"))),
                  uplo = if (x@uplo == "U") "L" else "U",
                  rcond = x@rcond)
          }, valueClass = "dsyMatrix")

setIs("dsyMatrix", "dpoMatrix",
      test = function(obj)
          "try-error" != class(try(.Call("dpoMatrix_chol", obj), TRUE)),
      replace = function(obj, value) {
          obj@uplo <- value@uplo
          obj@rcond <- value@rcond
          obj@factors <- value@factors
          obj@x <- value@x
          obj@Dim <- value@Dim
          obj@Dimnames <- value@Dimnames
          obj}
      )

## Now that we have "chol", we can define  "determinant" methods,
## exactly like in ./dsCMatrix.R
## DB - Probably figure out how to use the BunchKaufman decomposition instead
## {{FIXME: Shouldn't it be possible to have "determinant" work by
## default automatically for "Matrix"es  when there's a "chol" method available?
## -- not have to define showMethod("determinant", ...) for all classes

