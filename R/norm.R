## METHODS FOR GENERIC: norm
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## base::norm(x, type = "2")
norm2 <- function(x) if(anyNA(x)) NA_real_ else svd(x, nu = 0L, nv = 0L)$d[1L]

setMethod("norm", signature(x = "ANY", type = "missing"),
	  function(x, type, ...) norm(x, type = "O", ...))

setMethod("norm", signature(x = "sparseMatrix", type = "character"),
	  function(x, type, ...) {
              if(any(x@Dim == 0L))
                  return(0)
              type <- toupper(substr(type[1L], 1L, 1L))
	      switch(type,
		     "O" = ,
		     "1" = max(colSums(abs(x))), # L_1
		     "I" = max(rowSums(abs(x))), # L_Infinity
		     "F" = sqrt(sum(x^2)),       # Frobenius
		     "M" = max(abs(x)),          # maximum modulus
		     "2" = norm2(x),             # maximum singular value
		     stop("invalid 'type'"))
	  })

setMethod("norm", signature(x = "diagonalMatrix", type = "character"),
	  function(x, type, ...) {
              if((n <- x@Dim[1L]) == 0L)
                  return(0)
	      type <- toupper(substr(type[1L], 1L, 1L))
              if(identical(type, "F"))
                  sqrt(if(x@diag == "U") n else sum(x@x^2))
	      else if(x@diag == "U")
                  1
              else max(abs(x@x))
	  })

setMethod("norm", signature(x = "denseMatrix", type = "character"),
	  function(x, type, ...) norm(..dense2d(x), type = type, ...))

setMethod("norm", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dgeMatrix_norm, x, type))

setMethod("norm", signature(x = "dtrMatrix", type = "character"),
	  function(x, type, ...) {
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dtrMatrix_norm, x, type)
          })

setMethod("norm", signature(x = "dtpMatrix", type = "character"),
	  function(x, type, ...)
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dtpMatrix_norm, x, type))

setMethod("norm", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...)
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dsyMatrix_norm, x, type))

setMethod("norm", signature(x = "dspMatrix", type = "character"),
	  function(x, type, ...)
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dspMatrix_norm, x, type))
