if (FALSE) {
setMethod("eigen", signature(x = "dgeMatrix", only.values = "missing"),
          function(x, symmetric, only.values, EISPACK) # << must match generic
          .Call("dgeMatrix_eigen", x, FALSE))

setMethod("eigen", signature(x = "dgeMatrix", only.values = "logical"),
          function(x, symmetric, only.values, EISPACK)
          .Call("dgeMatrix_eigen", x, only.values))
} #not yet

setMethod("Schur", signature(x = "dgeMatrix", vectors = "missing"),
          function(x, vectors, ...) .Call("dgeMatrix_Schur", x, TRUE))

setMethod("Schur", signature(x = "dgeMatrix", vectors = "logical"),
          function(x, vectors, ...) .Call("dgeMatrix_Schur", x, vectors))

