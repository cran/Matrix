#### Sparse Matrices in Compressed row-oriented format
####                               --- "R"

### contains = "dMatrix"

setAs("dgRMatrix", "dgTMatrix",
      function(from) .Call(compressed_to_dgTMatrix, from, FALSE))

### FIXME: Activate the following with "little" change in ../src/dgCMatrix.c
### -----  similar to compressed_to_dgTMatrix above       ~~~~~~~~~~~~~~~~~~

##setAs("dgRMatrix", "matrix",
##      function(from) .Call(csc_to_matrix, from))
## easy "cheap" alternative:
setAs("dgRMatrix", "matrix",
      function(from) as(as(from, "dgTMatrix"), "matrix"))

##setAs("dgRMatrix", "dgeMatrix",
##      function(from) .Call(csc_to_dgeMatrix, from))

##setAs("matrix", "dgRMatrix",
##      function(from) {
##          storage.mode(from) <- "double"
##          .Call(matrix_to_csc, from)
##      })


##setMethod("diag", signature(x = "dgRMatrix"),
##          function(x = 1, nrow, ncol = n) .Call(csc_getDiag, x))

setAs("dgRMatrix", "dgCMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgCMatrix"))

## **VERY** cheap substitutes:  work via dgC and t(.)
.to.dgR <- function(from) {
    m <- as(t(from), "dgCMatrix")
    new("dgRMatrix", Dim = dim(from), Dimnames = .M.DN(from),
        p = m@p, j = m@i, x = m@x)
}

setAs("matrix",    "dgRMatrix", .to.dgR)
setAs("dgCMatrix", "dgRMatrix", .to.dgR)
setAs("dgTMatrix", "dgRMatrix", .to.dgR)



##setAs("dgRMatrix", "dgeMatrix",
##      function(from) .Call(csc_to_dgeMatrix, from))

##setAs("matrix", "dgRMatrix",
##      function(from) {
##          storage.mode(from) <- "double"
##          .Call(matrix_to_csc, from)
##      })


## try to define for "Matrix" -- once and for all -- but that fails -- why? __ FIXME __
## setMethod("dim", signature(x = "dgRMatrix"),
##           function(x) x@Dim, valueClass = "integer")

##setMethod("t", signature(x = "dgRMatrix"),
##          function(x) .Call(csc_transpose, x),
##          valueClass = "dgRMatrix")

setMethod("image", "dgRMatrix",
          function(x, ...) {
              x <- as(x, "dgTMatrix")
              callGeneric()
          })


## Want tril(), triu(), band() --- just as "indexing" ---
## return a "close" class:
setMethod("tril", "RsparseMatrix",
	  function(x, k = 0, ...) as_Rsparse(tril(as_Csparse(x), k = k, ...)))
setMethod("triu", "RsparseMatrix",
	  function(x, k = 0, ...) as_Rsparse(triu(as_Csparse(x), k = k, ...)))
setMethod("band", "RsparseMatrix",
	  function(x, k1, k2, ...)
	  as_Rsparse(band(as_Csparse(x), k1 = k1, k2 = k2, ...)))
