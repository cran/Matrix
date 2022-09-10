#### The "mother" of all  Matrix factorizations

## use a "fits all" bail-out method -- eventually this should not happen anymore
setMethod("expand", "MatrixFactorization",
          function(x, ...) .bail.out.1(.Generic, class(x)))

setMethod("show", "MatrixFactorization",
	  function(object) { ## cheap one -- can have better for sub-classes
	      ## cl <- class(object)
	      ## cat(sprintf("'MatrixFactorization' of class \"%s\"\n", cl))
	      cat("'MatrixFactorization' of ")
	      str(object)
	  })
setMethod("show", "BunchKaufman",
	  function(object) {
	      cat("'Bunch-Kaufman' factorization of ")
	      str(object)
	  })
setMethod("show", "pBunchKaufman",
	  function(object) {
	      cat("packed 'Bunch-Kaufman' factorization of ")
	      str(object)
	  })
## However, these result from chol() and should  *just be* a matrix to the non-expert user:
setMethod("show",  "Cholesky", function(object) prMatrix(object))
setMethod("show", "pCholesky", function(object) prMatrix(object))


setMethod("dim", "MatrixFactorization", function(x) x@Dim)
