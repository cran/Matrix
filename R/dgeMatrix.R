setAs("matrix", "dgeMatrix",
      function(from) .Call(dup_mMatrix_as_dgeMatrix, from))
if(FALSE)## the above is MUCH faster than
setAs("matrix", "dgeMatrix",
      function(from) {
	  new("dgeMatrix",
	      x = as.double(from),
	      Dim = as.integer(dim(from)),
	      Dimnames = .M.DN(from))
      })

setAs("dgeMatrix", "matrix",
      function(from) {
	  array(from@x, dim = from@Dim, dimnames = from@Dimnames)
      })

## Group Methods, see ?Arith (e.g.)
## ----- only work with NAMESPACE importFrom(methods, ..)

setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
	  signature(e1 = "dgeMatrix", e2 = "dgeMatrix"),
	  function(e1, e2) {
	      ## NB:  triangular, symmetric, etc may need own method
	      d1 <- e1@Dim
	      d2 <- e2@Dim
	      eqD <- d1 == d2
	      if (!eqD[1])
		  stop("Matrices must have same number of rows for arithmetic")
	      same.dim <- eqD[2]
	      if (same.dim) {
		  d <- d1
		  dn <- dimNamesCheck(e1, e2)
	      }
	      else { # nrows differ ----> maybe recycling
		  if(d2[2] %% d1[2] == 0) { # nrow(e2) is a multiple
		      e1@x <- rep.int(e1@x, d2[2] %/% d1[2])
		      d <- d2
		      dn <- e2@Dimnames
		  } else if(d1[2] %% d2[2] == 0) { # nrow(e1) is a multiple
		      e2@x <- rep.int(e2@x, d1[2] %/% d2[2])
		      d <- d1
		      dn <- e1@Dimnames
		  } else
		      stop("number of rows are not compatible for ", .Generic)
	      }

	      ## be smart and preserve, e.g., triangular, or symmetric
	      ## but this sucks: For these,
	      ## 'uplo' and 'diag' also must coincide or be dealt with properly

	      ## ==> triangular, symmetric, etc may need own method
	      ##     also since their @x is `non-typical'

##		 if(same.dim) {
##		     if(extends(class(e1), class(e2))) {
##			 e2@x <- callGeneric(e1@x, e2@x)
##			 e2@Dimnames <- dn
##			 e2
##		     }
##		     else if(extends(class(e2), class(e1))) {
##			 e1@x <- callGeneric(e1@x, e2@x)
##			 e1@Dimnames <- dn
##			 e1
##		     }
##		 }
##		 else
		  new("dgeMatrix", Dim = d, Dimnames = dn,
		      x = callGeneric(e1@x, e2@x))
	  })

setMethod("Arith",
	  signature(e1 = "dgeMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      d <- e1@Dim
	      le <- length(e2)
	      if(le == 1 || le == d[1] || prod(d) == le) { # matching dim
		  e1@x <- callGeneric(e1@x, as.vector(e2))
		  e1
	      } else stop ("length of 2nd arg does not match dimension of first")
	  })

setMethod("Arith",
	  signature(e1 = "numeric", e2 = "dgeMatrix"),
	  function(e1, e2) {
	      d <- e2@Dim
	      le <- length(e1)
	      if(le == 1 || le == d[1] || prod(d) == le) { # matching dim
		  e2@x <- callGeneric(as.vector(e1), e2@x)
		  e2
	      } else stop ("length of 1st arg does not match dimension of 2nd")
	  })

setMethod("Math",
	  signature(x = "dgeMatrix"),
	  function(x) {
	      x@x <- callGeneric(x@x)
	      x
	  })

##  "Math2" is in ./dMatrix.R

## "Summary"

## "Compare"  now happens in ./dMatrix.R

## -- end{group generics} -----------------------


##  "[" settings are "up in"  Matrix.R & denseMatrix.R


setMethod("as.vector", signature(x = "dgeMatrix", mode = "missing"),
          function(x) x@x)

setMethod("norm", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) norm(x, type = "O", ...))

setMethod("norm", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
	  .Call(dgeMatrix_norm, x, type),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) rcond(x, type = "O", ...))

setMethod("rcond", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
	  .Call(dgeMatrix_rcond, x, type),
	  valueClass = "numeric")

setMethod("t", signature(x = "dgeMatrix"), t_geMatrix)

## crossprod(x) & tcrossprod(x) :
setMethod("crossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call(dgeMatrix_crossprod, x, FALSE),
	  valueClass = "dpoMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call(dgeMatrix_crossprod, x, TRUE),
	  valueClass = "dpoMatrix")

if(FALSE) { ## this would mask 'base::tcrossprod'
setMethod("tcrossprod", signature(x = "matrix", y = "missing"),
	  function(x, y = NULL)
          .Call(dgeMatrix_crossprod, as(x, "dgeMatrix"), TRUE),
	  valueClass = "dpoMatrix")

setMethod("tcrossprod", signature(x = "numeric", y = "missing"),
	  function(x, y = NULL) callGeneric(as.matrix(as.double(x))))
}

## crossprod (x,y)
setMethod("crossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call(dgeMatrix_dgeMatrix_crossprod, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call(dgeMatrix_matrix_crossprod, x, y, FALSE),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL)
	  .Call(dgeMatrix_matrix_crossprod, x, as.matrix(as.double(y)), FALSE),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as.matrix(as.double(x)), y),
	  valueClass = "dgeMatrix")

## tcrossprod (x,y)
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call(dgeMatrix_dgeMatrix_crossprod, x, y, TRUE),
	  valueClass = "dgeMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call(dgeMatrix_matrix_crossprod, x, y, TRUE),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL)
	  .Call(dgeMatrix_matrix_crossprod, x, rbind(as.double(y)), TRUE),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y),
	  valueClass = "dgeMatrix")
setMethod("tcrossprod", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(rbind(as.double(x)), y),
	  valueClass = "dgeMatrix")

## %*% methods
setMethod("%*%", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, x, y, FALSE),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, y, x, TRUE),
          valueClass = "dgeMatrix")

## DB: Should we retain these methods?  Does the shortcut save enough
## to justify additional signatures?
## dgeMatrix <-> numeric: conceptually dispatch to "matrix" one, but shortcut
setMethod("%*%", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y) .Call(dgeMatrix_matrix_mm, x, y, FALSE),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y)
	  .Call(dgeMatrix_matrix_mm, y, rbind(x), TRUE),
	  valueClass = "dgeMatrix")

setMethod("diag", signature(x = "dgeMatrix"),
	  function(x = 1, nrow, ncol = n)
	  .Call(dgeMatrix_getDiag, x))

setMethod("chol", signature(x = "dgeMatrix", pivot = "ANY"), cholMat)

setMethod("solve", signature(a = "dgeMatrix", b = "missing"),
	  function(a, b, ...) .Call(dgeMatrix_solve, a),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "ddenseMatrix"),
	  function(a, b, ...) .Call(dgeMatrix_matrix_solve, a, b),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "matrix"),
	  function(a, b, ...) .Call(dgeMatrix_matrix_solve, a, b),
          valueClass = "dgeMatrix")

## not needed - method for numeric defined for Matrix class
## setMethod("solve", signature(a = "dgeMatrix", b = "numeric"),
## 	  function(a, b, ...)
## 	  .Call(dgeMatrix_matrix_solve, a, as.matrix(as.double(b))))

setMethod("lu", signature(x = "dgeMatrix"),
	  function(x, ...) .Call(dgeMatrix_LU, x),
	  valueClass = "denseLU")

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "missing"),
	  function(x, logarithm, ...)
	  .Call(dgeMatrix_determinant, x, TRUE))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
	  function(x, logarithm, ...)
	  .Call(dgeMatrix_determinant, x, logarithm))


setMethod("expm", signature(x = "dgeMatrix"),
	  function(x) .Call(dgeMatrix_exp, x),
	  valueClass = "dgeMatrix")

setMethod("colSums", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call(dgeMatrix_colsums, x, na.rm, TRUE, FALSE),
	  valueClass = "numeric")

setMethod("colMeans", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call(dgeMatrix_colsums, x, na.rm, TRUE, TRUE),
	  valueClass = "numeric")

setMethod("rowSums", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call(dgeMatrix_colsums, x, na.rm, FALSE, FALSE),
	  valueClass = "numeric")

setMethod("rowMeans", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call(dgeMatrix_colsums, x, na.rm, FALSE, TRUE),
	  valueClass = "numeric")
