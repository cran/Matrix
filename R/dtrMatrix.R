#### Triangular Matrices -- Coercion and Methods

setAs("dtrMatrix", "dgeMatrix",
      function(from) .Call("dtrMatrix_as_dgeMatrix", from) )

setAs("dtrMatrix", "dtpMatrix",
      function(from) .Call("dtrMatrix_as_dtpMatrix", from) )

## needed for t() method
setAs("dtrMatrix", "matrix",
      function(from) .Call("dtrMatrix_as_matrix", from) )

setAs("matrix", "dtrMatrix",
      function(from) as(as(from, "dgeMatrix"), "dtrMatrix"))

## Group Methods:
## TODO: carefully check for the cases where the result remains triangular
## instead : inherit them from "dgeMatrix" via definition in ./dMatrix.R

setMethod("%*%", signature(x = "dtrMatrix", y = "dgeMatrix"),
	  function(x, y) .Call("dtrMatrix_matrix_mm", x, y, TRUE, FALSE),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dtrMatrix", y = "matrix"),
	  function(x, y) .Call("dtrMatrix_matrix_mm", x, y, FALSE, FALSE),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "dtrMatrix"),
	  function(x, y) .Call("dtrMatrix_matrix_mm", y, x, TRUE, TRUE),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dtrMatrix"),
	  function(x, y) .Call("dtrMatrix_matrix_mm", y, x, FALSE, TRUE),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dtrMatrix", y = "dtrMatrix"),
	  function(x, y) callGeneric(x, as(y, "dgeMatrix")),
          valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dtrMatrix", y = "missing"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix")),
	  valueClass = "dpoMatrix")

setMethod("determinant", signature(x = "dtrMatrix", logarithm = "missing"),
	  function(x, logarithm, ...) callGeneric(x, TRUE))

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

setMethod("solve", signature(a = "dtrMatrix", b="dgeMatrix"),
	  function(a, b, ...) .Call("dtrMatrix_matrix_solve", a, b, TRUE),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtrMatrix", b="matrix"),
	  function(a, b, ...) .Call("dtrMatrix_matrix_solve", a, b, FALSE),
	  valueClass = "dgeMatrix")

setMethod("t", signature(x = "dtrMatrix"),
	  function(x) {
	      new("dtrMatrix",
                  Dim = x@Dim[2:1], Dimnames = x@Dimnames[2:1],
                  x = as.vector(t(as(x, "matrix"))),
	          uplo = if (x@uplo == "U") "L" else "U",
                  diag = x@diag)
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
