setMethod("chol", signature(x = "poMatrix"),
          function(x, pivot, LINPACK)
          .Call("poMatrix_chol", x, PACKAGE = "Matrix"))

setMethod("rcond", signature(x = "poMatrix", type = "character"),
          function(x, type, ...)
          .Call("poMatrix_rcond", x, type, PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("rcond", signature(x = "poMatrix", type = "missing"),
          function(x, type, ...)
          .Call("poMatrix_rcond", x, "O", PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("solve", signature(a = "poMatrix", b = "missing"),
          function(a, b, ...)
          .Call("poMatrix_solve", a, PACKAGE = "Matrix"),
          valueClass = "poMatrix")

setMethod("solve", signature(a = "poMatrix", b = "geMatrix"),
          function(a, b, ...)
          .Call("poMatrix_geMatrix_solve", a, b, PACKAGE = "Matrix"),
          valueClass = "geMatrix")

setMethod("solve", signature(a = "poMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("poMatrix_matrix_solve", a, b, PACKAGE = "Matrix"),
          valueClass = "matrix")

setMethod("solve", signature(a = "poMatrix", b = "numeric"),
          function(a, b, ...)
          as.numeric(.Call("poMatrix_matrix_solve",
                           a, as.matrix(b), PACKAGE = "Matrix")),
          valueClass = "numeric")
