#### Sparse Matrices in Compressed column-oriented format

### contains = "dMatrix"

setAs("dgCMatrix", "dgTMatrix",
      function(from) .Call("compressed_to_dgTMatrix", from, TRUE))

setAs("dgCMatrix", "matrix",
      function(from) .Call("csc_to_matrix", from))

setAs("dgCMatrix", "dgeMatrix",
      function(from) .Call("csc_to_dgeMatrix", from))

setAs("dgCMatrix", "dgBCMatrix",
      function(from) new("dgBCMatrix", p = from@p, i = from@i,
                         x = array(from@x, c(1, 1, length(from@x)))))

setAs("matrix", "dgCMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .Call("matrix_to_csc", from)
      })

setMethod("crossprod", signature(x = "dgCMatrix", y = "missing"),
          function(x, y = NULL) .Call("csc_crossprod", x))

setMethod("crossprod", signature(x = "dgCMatrix", y = "matrix"),
          function(x, y = NULL) .Call("csc_matrix_crossprod", x, y))

setMethod("crossprod", signature(x = "dgCMatrix", y = "numeric"),
          function(x, y = NULL) .Call("csc_matrix_crossprod", x, as.matrix(y)))

setMethod("tcrossprod", signature(x = "dgCMatrix"),
          function(x) .Call("csc_tcrossprod", x))

setMethod("diag", signature(x = "dgCMatrix"),
          function(x = 1, nrow, ncol = n) .Call("csc_getDiag", x))

## try to define for "Matrix" -- once and for all -- but that fails -- why?
setMethod("dim", signature(x = "dgCMatrix"),
          function(x) x@Dim, valueClass = "integer")

setMethod("t", signature(x = "dgCMatrix"),
          function(x) .Call("csc_transpose", x),
          valueClass = "dgCMatrix")

setMethod("image", "dgCMatrix",
          function(x, ...) {
              x = as(x, "dgTMatrix")
              callGeneric()
          })
