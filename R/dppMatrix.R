#### Positive-definite Symmetric Packed Matrices -- Coercion and Methods

setAs("dppMatrix", "dpoMatrix",
      function(from) as(as(from, "dsyMatrix"), "dpoMatrix"))

setMethod("chol", signature(x = "dppMatrix"),
          function(x, pivot, LINPACK)
          .Call("dppMatrix_chol", x))

setMethod("rcond", signature(x = "dppMatrix", type = "character"),
          function(x, type, ...)
          .Call("dppMatrix_rcond", x, type),
          valueClass = "numeric")

setMethod("rcond", signature(x = "dppMatrix", type = "missing"),
          function(x, type, ...)
          .Call("dppMatrix_rcond", x, "O"),
          valueClass = "numeric")

setMethod("solve", signature(a = "dppMatrix", b = "missing"),
          function(a, b, ...) .Call("dppMatrix_solve", a),
          valueClass = "dppMatrix")

setMethod("solve", signature(a = "dppMatrix", b = "dgeMatrix"),
          function(a, b, ...)
          .Call("dppMatrix_matrix_solve", a, b, TRUE),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dppMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("dppMatrix_matrix_solve", a, b, FALSE),
          valueClass = "dgeMatrix")

##setMethod("solve", signature(a = "dppMatrix", b = "numeric"),
##          function(a, b, ...)
##          .Call("dppMatrix_matrix_solve", a, as.matrix(b), FALSE),
##          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dppMatrix", b = "integer"),
          function(a, b, ...) {
              storage.mode(b) <- "double"
              .Call("dppMatrix_matrix_solve", a, as.matrix(b), FALSE)
          }, valueClass = "dgeMatrix")

setMethod("t", signature(x = "dppMatrix"),
          function(x) as(t(as(x, "dspMatrix")), "dppMatrix"),
          valueClass = "dppMatrix")

