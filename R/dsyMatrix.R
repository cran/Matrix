 ### Coercion and Methods for Symmetric Matrices

setAs("dsyMatrix", "dgeMatrix",
      function(from) .Call("dsyMatrix_as_dgeMatrix", from) )

setAs("dsyMatrix", "matrix",
      function(from) .Call("dsyMatrix_as_matrix", from) )

setMethod("%*%", signature(x = "dsyMatrix", y = "dgeMatrix"),
          function(x, y) .Call("dsyMatrix_dgeMatrix_mm", x, y) )

setMethod("%*%", signature(x = "dgeMatrix", y = "dsyMatrix"),
          function(x, y) .Call("dsyMatrix_dgeMatrix_mm_R", y, x) )

setMethod("norm", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...) .Call("dsyMatrix_norm", x, type),
          valueClass = "numeric")

setMethod("norm", signature(x = "dsyMatrix", type = "missing"),
          function(x, type, ...) .Call("dsyMatrix_norm", x, "O"),
          valueClass = "numeric")

## Should this create the opposite storage format - i.e. "U" -> "L"
## and vice-versa?
setMethod("t", signature(x = "dsyMatrix"), function(x) x)
