setAs("cscBlocked", "cscMatrix",
      function(from) {
          if (any(dim(from@x)[1:2] != 1))
              stop("code not yet written for nr > 1 or nc > 1")
          new("cscMatrix", p = from@p, i = from@i, x = as.vector(from@x),
                 Dim = c(max(from@i) + 1:1, length(from@p) - 1:1))
      })

setAs("cscBlocked", "tripletMatrix",
      function(from) as(as(from, "cscMatrix"), "tripletMatrix"))
