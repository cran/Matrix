setMethod("solve", signature(a = "sscMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("sscMatrix_matrix_solve", a, b))

setMethod("chol", signature(x = "sscMatrix", pivot = "missing"),
          function(x, pivot, LINPACK)
          .Call("sscMatrix_chol", x, TRUE))

setMethod("chol", signature(x = "sscMatrix", pivot = "logical"),
          function(x, pivot, LINPACK)
          .Call("sscMatrix_chol", x, pivot))

setMethod("t", signature(x = "sscMatrix"),
          function(x) .Call("ssc_transpose", x),
          valueClass = "sscMatrix")

setMethod("determinant", signature(x = "sscMatrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "sscMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
      {
          ldet <- sum(log(chol(x)@D))
          modulus <- if (logarithm) ldet else exp(ldet)
          attr(modulus, "logarithm") <- logarithm
          val <- list(modulus = modulus, sign = as.integer(1))
          class(val) <- "det"
          val
      })

setAs("sscMatrix", "tripletMatrix",
      function(from) .Call("sscMatrix_to_triplet", from))

setAs("sscMatrix", "geMatrix",
      function(from) as(as(from, "tripletMatrix"), "geMatrix"))

setAs("sscMatrix", "matrix",
      function(from) as(as(from, "tripletMatrix"), "matrix"))

