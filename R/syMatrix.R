setAs("syMatrix", "geMatrix",
      function(from) {
          .Call("syMatrix_as_geMatrix", from)
      })

setAs("syMatrix", "matrix",
      function(from) {
          .Call("syMatrix_as_matrix", from)
      })

setMethod("%*%", signature(x = "syMatrix", y = "geMatrix"),
          function(x, y)
          .Call("syMatrix_geMatrix_mm", x, y))

setMethod("%*%", signature(x = "geMatrix", y = "syMatrix"),
          function(x, y)
          .Call("syMatrix_geMatrix_mm_R", y, x))

setMethod("norm", signature(x = "syMatrix", type = "character"),
          function(x, type, ...)
          .Call("syMatrix_norm", x, type),
          valueClass = "numeric")

setMethod("norm", signature(x = "syMatrix", type = "missing"),
          function(x, type, ...)
          .Call("syMatrix_norm", x, "O"),
          valueClass = "numeric")

setMethod("t", signature(x = "syMatrix"),
          function(x) x)
