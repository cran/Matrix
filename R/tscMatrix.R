setMethod("t", signature(x = "tscMatrix"),
          function(x) .Call("tsc_transpose", x, PACKAGE = "Matrix"),
          valueClass = "tscMatrix")

setAs("tscMatrix", "tripletMatrix",
      function(from) .Call("tsc_to_triplet", from, PACKAGE = "Matrix"))

setAs("tscMatrix", "geMatrix",
      function(from) as(as(from, "tripletMatrix"), "geMatrix"))
