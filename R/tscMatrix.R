setMethod("t", signature(x = "tscMatrix"),
          function(x) .Call("tsc_transpose", x),
          valueClass = "tscMatrix")

setAs("tscMatrix", "tripletMatrix",
      function(from) .Call("tsc_to_triplet", from))

setAs("tscMatrix", "geMatrix",
      function(from) as(as(from, "tripletMatrix"), "geMatrix"))
