setMethod("t", signature(x = "tscMatrix"),
          function(x) .Call("tsc_transpose", x, PACKAGE = "Matrix"),
          valueClass = "tscMatrix")
