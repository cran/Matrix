#### Symmetric Sparse Matrices in compressed column-oriented format

### contains = "dgCMatrix"

setAs("dsCMatrix", "dgTMatrix",
      function(from) .Call("dsCMatrix_to_dgTMatrix", from))

setAs("dsCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))
setAs("dsCMatrix", "matrix",
      function(from) as(as(from, "dgTMatrix"), "matrix"))

if(FALSE) { ## << FIXME
setAs("dsCMatrix", "dsTMatrix",
      function(from) { }) ## < follow "csc_to_dgTMatrix" in ../src/dgTMatrix.c
setAs("dsCMatrix", "dsyMatrix",
      function(from) as(as(from, "dsTMatrix"), "dsyMatrix"))
}# not yet


setMethod("solve", signature(a = "dsCMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("dsCMatrix_matrix_solve", a, b))

setMethod("chol", signature(x = "dsCMatrix", pivot = "missing"),
          function(x, pivot, LINPACK) .Call("dsCMatrix_chol", x, TRUE))

setMethod("chol", signature(x = "dsCMatrix", pivot = "logical"),
          function(x, pivot, LINPACK) .Call("dsCMatrix_chol", x, pivot))

setMethod("t", signature(x = "dsCMatrix"),
          function(x) .Call("ssc_transpose", x),
          valueClass = "dsCMatrix")

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
      {
          ldet <- sum(log(chol(x)@D))
          modulus <- if (logarithm) ldet else exp(ldet)
          attr(modulus, "logarithm") <- logarithm
          val <- list(modulus = modulus, sign = as.integer(1))
          class(val) <- "det"
          val
      })
