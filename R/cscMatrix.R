setMethod("crossprod", signature(x = "cscMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("csc_crossprod", x))

setMethod("crossprod", signature(x = "cscMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", x, y))

setMethod("crossprod", signature(x = "cscMatrix", y = "numeric"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", x, as.matrix(y)))

setMethod("dim", signature(x = "cscMatrix"),
          function(x) x@Dim, valueClass = "integer")

setMethod("diag", signature(x = "cscMatrix"),
          function(x = 1, nrow, ncol = n)
          .Call("csc_getDiag", x))

setAs("cscMatrix", "tripletMatrix",
      function(from)
      .Call("csc_to_triplet", from)
      )

setAs("cscMatrix", "matrix",
      function(from)
      .Call("csc_to_matrix", from))

setAs("cscMatrix", "geMatrix",
      function(from)
      .Call("csc_to_geMatrix", from))

setAs("matrix", "cscMatrix",
      function(from) {
          storage.mode(from) = "double"
          .Call("matrix_to_csc", from)
      })

setMethod("t", signature(x = "cscMatrix"),
          function(x) .Call("csc_transpose", x),
          valueClass = "cscMatrix")

setMethod("image", "cscMatrix",
          function(x, ...) {
              x = as(x, "tripletMatrix")
              callGeneric()
          })
