if (FALSE) {
setMethod("eigen", signature(x = "dgeMatrix", only.values = "missing"),
          function(x, symmetric, only.values, EISPACK) # << must match generic
          .Call("dgeMatrix_eigen", x, FALSE, PACKAGE = "Matrix"))

setMethod("eigen", signature(x = "dgeMatrix", only.values = "logical"),
          function(x, symmetric, only.values, EISPACK)
          .Call("dgeMatrix_eigen", x, only.values, PACKAGE = "Matrix"))
} #not yet

setMethod("Schur", signature(x = "dgeMatrix", vectors = "missing"),
          function(x, vectors, ...)
          .Call("dgeMatrix_Schur", x, TRUE, PACKAGE = "Matrix"))

setMethod("Schur", signature(x = "dgeMatrix", vectors = "logical"),
          function(x, vectors, ...)
          .Call("dgeMatrix_Schur", x, vectors, PACKAGE = "Matrix"))

