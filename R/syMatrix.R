setAs("syMatrix", "geMatrix",
      function(from) {
          .Call("syMatrix_as_geMatrix", from, PACKAGE="Matrix")
      })

setAs("syMatrix", "matrix",
      function(from) {
          .Call("syMatrix_as_matrix", from, PACKAGE="Matrix")
      })

setMethod("%*%", signature(x = "syMatrix", y = "geMatrix"),
          function(x, y)
          .Call("syMatrix_geMatrix_mm", x, y, PACKAGE = "Matrix"))

setMethod("%*%", signature(x = "geMatrix", y = "syMatrix"),
          function(x, y)
          .Call("syMatrix_geMatrix_mm_R", y, x, PACKAGE = "Matrix"))

setMethod("norm", signature(x = "syMatrix", type = "character"),
          function(x, type, ...)
          .Call("syMatrix_norm", x, type, PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("norm", signature(x = "syMatrix", type = "missing"),
          function(x, type, ...)
          .Call("syMatrix_norm", x, "O", PACKAGE = "Matrix"),
          valueClass = "numeric")

setMethod("t", signature(x = "syMatrix"),
          function(x) x)
