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

