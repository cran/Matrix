## METHODS FOR GENERIC: expm
## the matrix exponential
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: currently going via ddiMatrix or dgeMatrix in all cases

setMethod("expm", signature(x = "Matrix"),
          function(x) {
              d <- x@Dim
              if(d[1L] != d[2L])
                  stop("matrix is not square")
              expm(as(x, "dMatrix"))
          })

setMethod("expm", signature(x = "dsparseMatrix"),
          function(x) {
              d <- x@Dim
              if(d[1L] != d[2L])
                  stop("matrix is not square")
              expm(.sparse2dense(x))
          })

setMethod("expm", signature(x = "ddiMatrix"),
	  function(x) {
	      if(x@diag == "N") {
		  x@x <- exp(x@x)
	      } else {
		  x@diag <- "N"
		  x@x <- rep.int(exp(1), x@Dim[1L])
	      }
	      x
	  })

setMethod("expm", signature(x = "dgeMatrix"),
	  function(x) .Call(dgeMatrix_exp, x))

setMethod("expm", signature(x = "dtrMatrix"),
	  function(x) {
              r <- .Call(dgeMatrix_exp, .dense2g(x))
              if(x@uplo == "U") triu(r) else tril(r)
          })

setMethod("expm", signature(x = "dtpMatrix"),
	  function(x) {
              r <- .Call(dgeMatrix_exp, .dense2g(x))
              ## Pack without checking:
              .Call(unpackedMatrix_pack, r, TRUE, TRUE, x@uplo == "U")
          })

setMethod("expm", signature(x = "dsyMatrix"),
	  function(x) {
              r <- .Call(dgeMatrix_exp, .dense2g(x))
              forceSymmetric(r)
          })

setMethod("expm", signature(x = "dspMatrix"),
	  function(x) {
              r <- .Call(dgeMatrix_exp, .dense2g(x))
              ## Pack without checking:
              .Call(unpackedMatrix_pack, r, TRUE, FALSE, x@uplo == "U")
          })

## Until R supports it:
setMethod("expm", signature(x = "matrix"),
          function(x) {
              d <- dim(x)
              if(d[1L] != d[2L])
                  stop("matrix is not square")
              storage.mode(x) <- "double"
              expm(if(isDiagonal(x)) forceDiagonal(x) else .m2ge(x))
          })
