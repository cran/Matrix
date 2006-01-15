#### Positive-definite Symmetric Matrices -- Coercion and Methods

setAs("dpoMatrix", "dppMatrix",
      function(from) as(as(from, "dspMatrix"), "dppMatrix"))

setAs("dpoMatrix", "correlation",
      function(from) {
          sd <- sqrt(diag(from))
          if (!is.null(nms <- from@Dimnames[[1]])) names(sd) <- nms
          ans <- new("correlation", as(t(from/sd)/sd, "dpoMatrix"), sd = sd)
          ans@Dimnames <- from@Dimnames
          ans
      })

setMethod("chol", signature(x = "dpoMatrix"),
          function(x, pivot, LINPACK)
          .Call("dpoMatrix_chol", x, PACKAGE = "Matrix"))

setMethod("rcond", signature(x = "dpoMatrix", type = "character"),
          function(x, type, ...)
          .Call("dpoMatrix_rcond", x, type, PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("rcond", signature(x = "dpoMatrix", type = "missing"),
          function(x, type, ...)
          .Call("dpoMatrix_rcond", x, "O", PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("solve", signature(a = "dpoMatrix", b = "missing"),
          function(a, b, ...)
          .Call("dpoMatrix_solve", a, PACKAGE = "Matrix"),
          valueClass = "dpoMatrix")

setMethod("solve", signature(a = "dpoMatrix", b = "dgeMatrix"),
          function(a, b, ...)
          .Call("dpoMatrix_dgeMatrix_solve", a, b, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dpoMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("dpoMatrix_matrix_solve", a, b, PACKAGE = "Matrix"),
          valueClass = "matrix")

##setMethod("solve", signature(a = "dpoMatrix", b = "numeric"),
##          function(a, b, ...)
##          as.numeric(.Call("dpoMatrix_matrix_solve",
##                           a, as.matrix(b), PACKAGE = "Matrix")),
##          valueClass = "numeric")

