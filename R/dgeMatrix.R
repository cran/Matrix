setAs("matrix", "dgeMatrix",
      function(from) {
	  new("dgeMatrix",
	      x = as.double(from),
	      Dim = as.integer(dim(from)),
	      Dimnames =
	      if(!is.null(dn <- dimnames(from))) dn else list(NULL,NULL)
	      )
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
	      else { # nrows differ
		  if(d2[2] %% d1[2] == 0) { # nrow(e2) is a multiple
		      e1@x <- rep.int(e1@x, d2[2] %/% d1[2])
		      d <- d2
		      dn <- e2@Dimnames
		  } else if(d1[2] %% d2[2] == 0) { # nrow(e1) is a multiple
		      e2@x <- rep.int(e2@x, d1[2] %/% d2[2])
		      d <- d1
		      dn <- e1@Dimnames
		  } else
		      stop("number of rows are not compatible for arithmetic")
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

if(FALSE) ## unneeded with "Math2" in ./dMatrix.R
setMethod("Math2",
	  signature(x = "dgeMatrix", digits = "numeric"),
	  function(x, digits) {
	      x@x <- callGeneric(x@x, digits = digits)
	      x
	  })

## "Summary"


## TODO :  "Compare" -> returning  logical Matrices


## -- end{group generics} -----------------------

setMethod("as.vector", signature(x = "dgeMatrix", mode = "missing"),
          function(x) x@x)

setMethod("norm", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) norm(x, type = "O", ...))

setMethod("norm", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dgeMatrix_norm", x, type, PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) rcond(x, type = "O", ...))

setMethod("rcond", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dgeMatrix_rcond", x, type, PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("t", signature(x = "dgeMatrix"), t_geMatrix)

setMethod("crossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call("dgeMatrix_crossprod", x, FALSE, PACKAGE = "Matrix"),
	  valueClass = "dpoMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call("dgeMatrix_crossprod", x, TRUE, PACKAGE = "Matrix"),
	  valueClass = "dpoMatrix")
setMethod("tcrossprod", signature(x = "matrix", y = "missing"),
	  function(x, y = NULL) .Call("dgeMatrix_crossprod", as(x, "dgeMatrix"), TRUE, PACKAGE = "Matrix"),
	  valueClass = "dpoMatrix")
setMethod("tcrossprod", signature(x = "numeric", y = "missing"),
	  function(x, y = NULL) callGeneric(as.matrix(as.double(x))))


setMethod("crossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call("dgeMatrix_dgeMatrix_crossprod", x, y, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call("dgeMatrix_matrix_crossprod", x, y, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL)
	  .Call("dgeMatrix_matrix_crossprod", x, as.matrix(as.double(y)), PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as.matrix(as.double(x)), y),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y) .Call("dgeMatrix_matrix_mm", x, y, TRUE, FALSE, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

## dgeMatrix <-> matrix ("matrix" dispatches before "numeric" since R 2.1.0)
setMethod("%*%", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y) {
	      storage.mode(y) <- "double"
	      .Call("dgeMatrix_matrix_mm", x, y, FALSE, FALSE, PACKAGE = "Matrix")
	  }, valueClass = "dgeMatrix")


setMethod("%*%", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y) {
	      storage.mode(x) <- "double"
	      .Call("dgeMatrix_matrix_mm", y, x, FALSE, TRUE, PACKAGE = "Matrix")
	  }, valueClass = "dgeMatrix")

## dgeMatrix <-> numeric: conceptually dispatch to "matrix" one, but shortcut
setMethod("%*%", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y)
	  .Call("dgeMatrix_matrix_mm", x, as.matrix(as.double(y)), FALSE, FALSE, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")
setMethod("%*%", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y)
	  .Call("dgeMatrix_matrix_mm", y, rbind(as.double(x)), FALSE, TRUE, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

setMethod("diag", signature(x = "dgeMatrix"),
	  function(x = 1, nrow, ncol = n)
	  .Call("dgeMatrix_getDiag", x, PACKAGE = "Matrix"))


## DB - I don't think this is a good idea without first checking symmetry
#setMethod("chol", signature(x = "dgeMatrix", pivot = "ANY"), cholMat)

setMethod("solve", signature(a = "dgeMatrix", b = "missing"),
	  function(a, b, ...) .Call("dgeMatrix_solve", a, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "dgeMatrix"),
	  function(a, b, ...) .Call("dgeMatrix_matrix_solve", a, b, TRUE, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "matrix"),
	  function(a, b, ...) {
	      storage.mode(b) <- "double"
	      .Call("dgeMatrix_matrix_solve", a, b, FALSE, PACKAGE = "Matrix")
	  }, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "numeric"),
	  function(a, b, ...)
	  .Call("dgeMatrix_matrix_solve", a, as.matrix(as.double(b)), FALSE, PACKAGE = "Matrix"))

setMethod("lu", signature(x = "dgeMatrix"),
	  function(x, ...) .Call("dgeMatrix_LU", x, PACKAGE = "Matrix"),
          valueClass = "LU")

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "missing"),
	  function(x, logarithm, ...)
	  .Call("dgeMatrix_determinant", x, TRUE, PACKAGE = "Matrix"))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
	  function(x, logarithm, ...)
	  .Call("dgeMatrix_determinant", x, logarithm, PACKAGE = "Matrix"))


setMethod("expm", signature(x = "dgeMatrix"),
	  function(x) .Call("dgeMatrix_exp", x, PACKAGE = "Matrix"),
	  valueClass = "dgeMatrix")

setMethod("colSums", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call("dgeMatrix_colsums", x, na.rm, TRUE, FALSE, PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("colMeans", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call("dgeMatrix_colsums", x, na.rm, TRUE, TRUE, PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("rowSums", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call("dgeMatrix_colsums", x, na.rm, FALSE, FALSE, PACKAGE = "Matrix"),
	  valueClass = "numeric")

setMethod("rowMeans", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call("dgeMatrix_colsums", x, na.rm, FALSE, TRUE, PACKAGE = "Matrix"),
	  valueClass = "numeric")

### The following all serve for	 as.Matrix()
### which is not yet exported (nor tested):

## utilities for Matrix.class() :

## FIXME  base::eigen() has a more sensible test for Hermitian/symmetry !
Hermitian.test <- function(x)
{
    ## Includes Symmetry test for non-complex 'x'
    if ((!inherits(x, "Matrix") && !is.matrix(x)) || (nrow(x) != ncol(x)))
        return(Inf)
    if (is.complex(x))
        max(Mod(x - t(Conj(x))))
    else
        max(abs(x - t(x)))
}

LowerTriangular.test <- function(x)
{
    ## return largest |value| in the lower triangle of x
    if ((!inherits(x, "Matrix") && !is.matrix(x))) return(Inf)
    i <- row(x) < col(x)
    if(!any(i)) return(Inf)
    max(if (is.complex(x)) abs(x[i]) else Mod(x[i]))
}

UpperTriangular.test <- function(x)
{
    if ((!inherits(x, "Matrix") && !is.matrix(x))) return(Inf)
    i <- row(x) > col(x)
    if(!any(i)) return(Inf)
    max(if (is.complex(x)) abs(x[i]) else Mod(x[i]))
}

Orthogonal.test <- function(x, byrow = FALSE, normal = TRUE)
{
    if ((!inherits(x, "Matrix") && !is.matrix(x))) return(Inf)
    if (byrow) { x <- t(x) }
    xx <- crossprod(x)
    if (normal) # check for orthonormal
	max(Mod(xx[row(xx) > col(xx)]), Mod(diag(xx) - 1))
    else
        max(Mod(xx[row(xx) > col(xx)]))
}

Orthonormal.test <- function(x, byrow = FALSE)
{ Orthogonal.test(x, byrow, normal = TRUE) }

is.Hermitian <- function(x, tol = 0) { Hermitian.test(x) <= tol }

is.LowerTriangular <- function(x, tol = 0) { LowerTriangular.test(x) <= tol }

is.UpperTriangular <- function(x, tol = 0) { UpperTriangular.test(x) <= tol }

is.ColOrthonormal <- function(x, tol = sqrt(.Machine$double.eps))
{ Orthonormal.test(x, byrow = FALSE) <= tol }

is.RowOrthonormal <- function(x, tol = sqrt(.Machine$double.eps))
{ Orthonormal.test(x, byrow = TRUE) <= tol }

is.Orthonormal <- function(x, tol = sqrt(.Machine$double.eps), byrow = FALSE)
{
    if (byrow)
	is.RowOrthonormal(x, tol)
    else
	is.ColOrthonormal(x, tol)
}


Matrix.class <- function(x, tol = 0, symmetry = TRUE, unit.diagonal = TRUE,
			 triangularity = c(TRUE, TRUE),
			 orthogonality = c(TRUE, TRUE),
			 normality = c(TRUE, TRUE))
{
    ## basic work horse for as.Matrix()

    val <- "Matrix"
    x <- as.matrix(x)
    if (symmetry) {
	if (is.Hermitian(x, tol)) val <- c("Hermitian", val)
    }
    if (triangularity[1]) {
	if (is.LowerTriangular(x, tol)) {
	    val <- c("LowerTriangular", val)
	    if (unit.diagonal)
		if (max(Mod(diag(x) - 1)) <= tol)
		    val <- c("UnitLowerTriangular", val)
	}
    }
    if (triangularity[2]) {
	if (is.UpperTriangular(x, tol)) {
	    val <- c("UpperTriangular", val)
	    if (unit.diagonal)
		if (max(Mod(diag(x) - 1)) <= tol)
		    val <- c("UnitUpperTriangular", val)
	}
    }
    if (orthogonality[1]) {
	if (is.ColOrthonormal(x, tol))
	    val <- c("ColOrthoNormal", "ColOrthogonal", val)
	else if (Orthogonal.test(x, normal = FALSE) <= tol)
	    val <- c("ColOrthogonal", val)
    }

    if (orthogonality[2]) {
	if (normality[2] && is.RowOrthonormal(x, tol))
	    val <- c("RowOrthoNormal", "RowOrthogonal", val)
	else if (Orthogonal.test(x, byrow = TRUE, normal = FALSE) <= tol)
	    val <- c("RowOrthogonal", val)
    }
    val
}


as.Matrix <- function(x, tol = .Machine$double.eps,
                      integer.max = .Machine$integer.max)
{
    if(is(x, "Matrix")) return(x)
    ## else
    if(!is.matrix(x)) x <- as.matrix(x)
    mc <- Matrix.class(x, tol = tol) ## a character *vector*
    xmode <-
        if(is.logical(x)) "l" else if(is.complex(x)) "z"
        else if(is.numeric(x)) {
            if(is.integer(x) || all(abs(x) < integer.max)) "i" else "d"
        }
        else stop("invalid data type")
    ## .... .... fixme
    as(x, smartFunction(mc))
}
