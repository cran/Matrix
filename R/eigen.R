setMethod("eigen", signature(x = "dgeMatrix", only.values = "missing"),
          function(x, symmetric, only.values, EISPACK) {# << must match generic
              nCall <- match.call()
              nCall$only.values <- FALSE
              eval(nCall, parent.frame())
          })

setMethod("eigen", signature(x = "dgeMatrix", only.values = "logical"),
          function(x, symmetric, only.values, EISPACK)
          .Call("dgeMatrix_eigen", x, only.values)
          )

setMethod("Schur", signature(x = "dgeMatrix", vectors = "missing"),
          function(x, vectors, ...) Schur(x, TRUE, ...))

setMethod("Schur", signature(x = "dgeMatrix", vectors = "logical"),
          function(x, vectors, ...) .Call("dgeMatrix_Schur", x, vectors)
          )

