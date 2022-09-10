## METHODS FOR GENERIC: rcond
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("rcond", signature(x = "ANY", norm = "missing"),
	  function(x, norm, ...) rcond(x, norm = "O", ...))

## FIXME: Need a version of LAPACK's rcond() algorithm, using sparse arithmetic
setMethod("rcond", signature(x = "sparseMatrix", norm = "character"),
	  function(x, norm, useInv = FALSE, ...) {
              ## As a workaround, allow computing 1/(norm(A) * norm(solve(A)))
              if(!isFALSE(useInv)) {
                  Ix <-
                      if(isS4(useInv) && is(useInv, "Matrix"))
                          useInv
                      else solve(x)
                  return(1 / (norm(x, type = norm) * norm(Ix, type = norm)))
              }
              d <- x@Dim
              ## FIXME: qr.R(qr(.)) warns about differing R (permutation!)
              ##        really fix qr.R() *or* go via dense even in those cases
	      rcond(if(d[1L] == d[2L]) {
			warning("rcond(.) via sparse -> dense coercion")
			as(x, "denseMatrix")
		    } else qr.R(qr(if(d[1L] < d[2L]) t(x) else x)),
		    norm = norm, ...)
	  })

setMethod("rcond", signature(x = "denseMatrix", norm = "character"),
	  function(x, norm, ...) rcond(..dense2d(x), norm = norm, ...))

setMethod("rcond", signature(x = "dgeMatrix", norm = "character"),
	  function(x, norm, ...) {
              d <- x@Dim
	      if(d[1L] != d[2L])
		  rcond(qr.R(qr(if(d[1L] < d[2L]) t(x) else x)),
                        norm = norm, ...)
              else .Call(dgeMatrix_rcond, x, norm)
	  })

setMethod("rcond", signature(x = "dtrMatrix", norm = "character"),
	  function(x, norm, ...) .Call(dtrMatrix_rcond, x, norm))

setMethod("rcond", signature(x = "dtpMatrix", norm = "character"),
	  function(x, norm, ...) .Call(dtpMatrix_rcond, x, norm))

setMethod("rcond", signature(x = "dsyMatrix", norm = "character"),
          function(x, norm, ...) .Call(dsyMatrix_rcond, x, norm))

setMethod("rcond", signature(x = "dspMatrix", norm = "character"),
          function(x, norm, ...) .Call(dspMatrix_rcond, x, norm))

setMethod("rcond", signature(x = "dpoMatrix", norm = "character"),
          function(x, norm, ...) .Call(dpoMatrix_rcond, x, norm))

setMethod("rcond", signature(x = "dppMatrix", norm = "character"),
	  function(x, norm, ...) .Call(dppMatrix_rcond, x, norm))
