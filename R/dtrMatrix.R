#### Triangular Matrices -- Coercion and Methods

setAs("dtrMatrix", "dgeMatrix",
      function(from) .Call("dtrMatrix_as_dgeMatrix", from) )

## needed for t() method
setAs("dtrMatrix", "matrix",
      function(from) .Call("dtrMatrix_as_matrix", from) )

setMethod("%*%", signature(x = "dtrMatrix", y = "dgeMatrix"),
	  function(x, y) .Call("dtrMatrix_dgeMatrix_mm", x, y))

setMethod("%*%", signature(x = "dgeMatrix", y = "dtrMatrix"),
	  function(x, y) .Call("dtrMatrix_dgeMatrix_mm_R", y, x))

setMethod("crossprod", signature(x = "dtrMatrix", y = "missing"),
	  function(x, y = NULL) crossprod(as(x, "dgeMatrix")),
	  valueClass = "dpoMatrix")

setMethod("determinant", signature(x = "dtrMatrix", logarithm = "missing"),
	  function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "dtrMatrix", logarithm = "logical"),
	  function(x, logarithm, ...) {
	      dg <- diag(x)
	      if (logarithm) {
		  modulus <- sum(log(abs(dg)))
		  sgn <- prod(sign(dg))
	      } else {
		  modulus <- prod(dg)
		  sgn <- sign(modulus)
		  modulus <- abs(modulus)
	      }
	      attr(modulus, "logarithm") <- logarithm
	      val <- list(modulus = modulus, sign = sgn)
	      class(val) <- "det"
	      val
	  })

setMethod("norm", signature(x = "dtrMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dtrMatrix_norm", x, type),
	  valueClass = "numeric")

setMethod("norm", signature(x = "dtrMatrix", type = "missing"),
	  function(x, type, ...)
	  .Call("dtrMatrix_norm", x, "O"),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dtrMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dtrMatrix_rcond", x, type),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dtrMatrix", type = "missing"),
	  function(x, type, ...)
	  .Call("dtrMatrix_rcond", x, "O"),
	  valueClass = "numeric")

setMethod("solve", signature(a = "dtrMatrix", b="missing"),
	  function(a, b, ...)
	  .Call("dtrMatrix_solve", a),
	  valueClass = "dtrMatrix")

setMethod("solve", signature(a = "dtrMatrix", b="matrix"),
	  function(a, b, ...)
	  .Call("dtrMatrix_matrix_solve", a, b),
	  valueClass = "matrix")

setMethod("t", signature(x = "dtrMatrix"),
	  function(x) {
	      val <- new("dtrMatrix",
                         Dim = x@Dim, Dimnames = x@Dimnames[2:1],
                         x = as.vector(t(as(x, "matrix"))))
	      if (x@uplo == "U") val@uplo <- "L"
	      if (x@diag == "U") val@diag <- "U"
	      val
	  }, valueClass = "dtrMatrix")

###

## in construction -> not yet exported
Diagonal <- function(n, x = NULL, ncol = n)
{
    ## Purpose: Constructor of diagonal matrices -- *not* diag() extractor!
    stopifnot(length(n) == 1, n == as.integer(n), n >= 0)
    n <- as.integer(n)
    stopifnot(length(ncol <- as.integer(ncol)) == 1, ncol >= 0)
    if(missing(x)) # unit diagonal matrix
        new("dtrMatrix", Dim = c(n,ncol), diag = "U",
            x = rep(0, .......))
    else {
        x <- as.numeric(x)
        stopifnot(length(x) == min(n, ncol))
        new("dtrMatrix", Dim = c(n,ncol), diag = "N",
            x = ..x..and..many..0s)
    }
}
