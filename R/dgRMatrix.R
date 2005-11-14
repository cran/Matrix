#### Sparse Matrices in Compressed row-oriented format

### contains = "dMatrix"

setAs("dgRMatrix", "dgTMatrix",
      function(from) .Call("compressed_to_dgTMatrix", from, FALSE, PACKAGE = "Matrix"))

##setAs("dgRMatrix", "matrix",
##      function(from) .Call("csc_to_matrix", from, PACKAGE = "Matrix"))

##setAs("dgRMatrix", "dgeMatrix",
##      function(from) .Call("csc_to_dgeMatrix", from, PACKAGE = "Matrix"))

##setAs("matrix", "dgRMatrix",
##      function(from) {
##          storage.mode(from) <- "double"
##          .Call("matrix_to_csc", from, PACKAGE = "Matrix")
##      })

##setMethod("diag", signature(x = "dgRMatrix"),
##          function(x = 1, nrow, ncol = n) .Call("csc_getDiag", x, PACKAGE = "Matrix"))

## try to define for "Matrix" -- once and for all -- but that fails -- why?
setMethod("dim", signature(x = "dgRMatrix"),
          function(x) x@Dim, valueClass = "integer")

##setMethod("t", signature(x = "dgRMatrix"),
##          function(x) .Call("csc_transpose", x, PACKAGE = "Matrix"),
##          valueClass = "dgRMatrix")

setMethod("image", "dgRMatrix",
          function(x, ...) {
              x = as(x, "dgTMatrix")
              callGeneric()
          })
