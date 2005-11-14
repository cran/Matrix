 ### Coercion and Methods for Symmetric Matrices

setAs("dspMatrix", "dsyMatrix",
      function(from) .Call("dspMatrix_as_dsyMatrix", from, PACKAGE = "Matrix"))

setAs("dspMatrix", "dgeMatrix",
      function(from) as(as(from, "dsyMatrix"), "dgeMatrix"))

setAs("dspMatrix", "matrix",
      function(from) as(as(from, "dsyMatrix"), "matrix"))

setMethod("rcond", signature(x = "dspMatrix", type = "character"),
          function(x, type, ...)
          .Call("dspMatrix_rcond", x, type, PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("rcond", signature(x = "dspMatrix", type = "missing"),
          function(x, type, ...)
          .Call("dspMatrix_rcond", x, "O", PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("%*%", signature(x = "dspMatrix", y = "dgeMatrix"),
          function(x, y)
          .Call("dspMatrix_matrix_mm", x, y, TRUE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dspMatrix", y = "matrix"),
          function(x, y)
          .Call("dspMatrix_matrix_mm", x, y, FALSE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

##setMethod("%*%", signature(x = "dspMatrix", y = "numeric"),
##          function(x, y)
##          .Call("dspMatrix_matrix_mm", x, as.matrix(y), FALSE, PACKAGE = "Matrix"),
##          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dspMatrix", y = "integer"),
          function(x, y) {
              storage.mode(y) <- "double"
              .Call("dspMatrix_matrix_mm", x, as.matrix(y), FALSE, PACKAGE = "Matrix")
          }, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dspMatrix", b = "missing"),
	  function(a, b, ...) .Call("dspMatrix_solve", a, PACKAGE = "Matrix"),
	  valueClass = "dspMatrix")

setMethod("solve", signature(a = "dspMatrix", b = "matrix"),
	  function(a, b, ...)
	  .Call("dspMatrix_matrix_solve", a, b, FALSE, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dspMatrix", b = "dgeMatrix"),
	  function(a, b, ...)
	  .Call("dspMatrix_matrix_solve", a, as(b,"matrix"), TRUE, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

##setMethod("solve", signature(a = "dspMatrix", b = "numeric"),
##	  function(a, b, ...)
##	  .Call("dspMatrix_matrix_solve", a, as.matrix(b), FALSE, PACKAGE = "Matrix"),
##	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dspMatrix", b = "integer"),
	  function(a, b, ...) {
	      storage.mode(b) <- "double"
	      .Call("dspMatrix_matrix_solve", a, as.matrix(b), FALSE, PACKAGE = "Matrix")
	  }, valueClass = "dgeMatrix")

setMethod("norm", signature(x = "dspMatrix", type = "character"),
          function(x, type, ...) .Call("dspMatrix_norm", x, type, PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("norm", signature(x = "dspMatrix", type = "missing"),
          function(x, type, ...) .Call("dspMatrix_norm", x, "O", PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("t", signature(x = "dspMatrix"),
          function(x) as(t(as(x, "dsyMatrix")), "dspMatrix"),
          valueClass = "dspMatrix")

setMethod("unpack", signature(x = "dspMatrix"),
          function(x, ...) as(x, "dsyMatrix"),
          valueClass = "dsyMatrix")

## The following allows  as(*, "dppMatrix").
## However it *requires* that dppMatrix_chol() gives an error
## for non-positive-semi-definite matrices -- which it does since 2005-10-03
## FIXME: This gives an error for singular pos.SEMI-def. matrices:
setIs("dspMatrix", "dppMatrix",
      test = function(obj)
          "try-error" != class(try(.Call("dppMatrix_chol", obj, PACKAGE = "Matrix"), TRUE)),
      replace = function(obj, value) {
          obj@uplo <- value@uplo
          obj@rcond <- value@rcond
          obj@factors <- value@factors
          obj@x <- value@x
          obj@Dim <- value@Dim
          obj@Dimnames <- value@Dimnames
          obj}
      )

