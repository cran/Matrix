setMethod("eigen", signature(x = "geMatrix", only.values = "missing"),
          function(x, symmetric, only.values, EISPACK) {
              nCall = match.call()
              nCall$only.values = FALSE
              eval(nCall, parent.frame())
          })

setMethod("eigen", signature(x = "geMatrix", only.values = "logical"),
          function(x, symmetric, only.values, EISPACK)
          .Call("geMatrix_eigen", x, only.values)
          )

setMethod("Schur", signature(x = "geMatrix", vectors = "missing"),
          function(x, vectors, ...) {
              nCall = match.call()
              nCall$vectors = FALSE
              eval(nCall, parent.frame())
          })
         
setMethod("Schur", signature(x = "geMatrix", vectors = "logical"),
          function(x, vectors, ...)
          .Call("geMatrix_Schur", x, vectors)
          )

