setMethod("crossprod", signature(x = "cscMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("csc_crossprod", x, PACKAGE = "Matrix"))

setMethod("crossprod", signature(x = "cscMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", x, y, PACKAGE = "Matrix"))

setMethod("crossprod", signature(x = "cscMatrix", y = "numeric"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", x, as.matrix(y), PACKAGE = "Matrix"))

setMethod("dim", signature(x = "cscMatrix"),
          function(x) x@Dim, valueClass = "integer")

setMethod("diag", signature(x = "cscMatrix"),
          function(x = 1, nrow, ncol = n)
          .Call("csc_getDiag", x, PACKAGE = "Matrix"))

setAs("cscMatrix", "tripletMatrix",
      function(from)
      .Call("csc_to_triplet", from, PACKAGE = "Matrix")
      )

setAs("cscMatrix", "matrix",
      function(from)
      .Call("csc_to_matrix", from, PACKAGE = "Matrix"))

setAs("cscMatrix", "geMatrix",
      function(from)
      .Call("csc_to_geMatrix", from, PACKAGE = "Matrix"))

setAs("matrix", "cscMatrix",
      function(from) {
          storage.mode(from) = "double"
          .Call("matrix_to_csc", from, PACKAGE = "Matrix")
      })

setMethod("t", signature(x = "cscMatrix"),
          function(x) .Call("csc_transpose", x, PACKAGE = "Matrix"),
          valueClass = "cscMatrix")

setMethod("image", "cscMatrix",
          function(x, ...) {
              x = as(x, "tripletMatrix")
              callGeneric()
          })
