setMethod("solve", signature(a = "sscMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("sscMatrix_matrix_solve", a, b, TRUE, PACKAGE = "Matrix"))

setMethod("chol", signature(x = "sscMatrix", pivot = "missing"),
          function(x, pivot, LINPACK)
          .Call("sscMatrix_chol", x, TRUE, PACKAGE = "Matrix"))

setMethod("chol", signature(x = "sscMatrix", pivot = "logical"),
          function(x, pivot, LINPACK)
          .Call("sscMatrix_chol", x, pivot, PACKAGE = "Matrix"))

setMethod("t", signature(x = "sscMatrix"),
          function(x) .Call("ssc_transpose", x, PACKAGE = "Matrix"),
          valueClass = "sscMatrix")
