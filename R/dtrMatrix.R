#### Triangular Matrices -- Coercion and Methods

setAs("dtrMatrix", "dgeMatrix",
      function(from) .Call("dtrMatrix_as_dgeMatrix", from, PACKAGE = "Matrix"))

## or rather setIs() {since test can fail }?
setAs("dgeMatrix", "dtrMatrix",
      function(from) {
          ## FIXME: also check for unit-diagonal: 'diag = "U"'
	  if(isTriangular(from, upper = TRUE))
	      new("dtrMatrix", x = from@x, Dim = from@Dim, uplo = "U",
		  Dimnames = from@Dimnames)
	  else if(isTriangular(from, upper = FALSE))
	      new("dtrMatrix", x = from@x, Dim = from@Dim, uplo = "L",
		  Dimnames = from@Dimnames)
	  else stop("not a triangular matrix")
      })


setAs("dtrMatrix", "dtpMatrix",
      function(from) .Call("dtrMatrix_as_dtpMatrix", from, PACKAGE = "Matrix"))

## needed for t() method
setAs("dtrMatrix", "matrix",
      function(from) .Call("dtrMatrix_as_matrix", from, PACKAGE = "Matrix"))

setAs("matrix", "dtrMatrix",
      function(from) as(as(from, "dgeMatrix"), "dtrMatrix"))

## Group Methods:
## TODO: carefully check for the cases where the result remains triangular
## instead : inherit them from "dgeMatrix" via definition in ./dMatrix.R

## Note: Just *because* we have an explicit  dtr -> dge coercion,
##       show( <ddenseMatrix> ) is not okay, and we need our own:
setMethod("show", "dtrMatrix", function(object) prMatrix(object))


setMethod("%*%", signature(x = "dtrMatrix", y = "dgeMatrix"),
	  function(x, y) .Call("dtrMatrix_matrix_mm", x, y, TRUE, FALSE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dtrMatrix", y = "matrix"),
	  function(x, y) .Call("dtrMatrix_matrix_mm", x, y, FALSE, FALSE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "dtrMatrix"),
	  function(x, y) .Call("dtrMatrix_matrix_mm", y, x, TRUE, TRUE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dtrMatrix"),
	  function(x, y) .Call("dtrMatrix_matrix_mm", y, x, FALSE, TRUE, PACKAGE = "Matrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dtrMatrix", y = "dtrMatrix"),
	  function(x, y) callGeneric(x = x, y = as(y, "dgeMatrix")),
          valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dtrMatrix", y = "missing"),
	  function(x, y = NULL) callGeneric(x = as(x, "dgeMatrix")),
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
	  .Call("dtrMatrix_norm", x, type, PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("norm", signature(x = "dtrMatrix", type = "missing"),
	  function(x, type, ...)
	  .Call("dtrMatrix_norm", x, "O", PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dtrMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dtrMatrix_rcond", x, type, PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dtrMatrix", type = "missing"),
	  function(x, type, ...)
	  .Call("dtrMatrix_rcond", x, "O", PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("solve", signature(a = "dtrMatrix", b="missing"),
	  function(a, b, ...)
	  .Call("dtrMatrix_solve", a, PACKAGE = "Matrix"),
	  valueClass = "dtrMatrix")

setMethod("solve", signature(a = "dtrMatrix", b="dgeMatrix"),
	  function(a, b, ...)
          .Call("dtrMatrix_matrix_solve", a, b, TRUE, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtrMatrix", b="matrix"),
	  function(a, b, ...)
          .Call("dtrMatrix_matrix_solve", a, b, FALSE, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

setMethod("t", signature(x = "dtrMatrix"), t_trMatrix)
