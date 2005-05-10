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
	      if (same.dim)
		  d <- d1
	      else { # nrows differ
		  if(d2[2] %% d1[2] == 0) { # nrow(e2) is a multiple
		      e1@x <- rep.int(e1@x, d2[2] %/% d1[2])
		      d <- d2
		  } else if(d1[2] %% d2[2] == 0) { # nrow(e1) is a multiple
		      e2@x <- rep.int(e2@x, d1[2] %/% d2[2])
		      d <- d1
		  }
		  else
		      stop("number of rows are not compatible for arithmetic")
	      }
	      dn0 <- list(NULL,NULL)
	      if(identical(dn0, dn <- e1@Dimnames))
		  dn <- e2@Dimnames
	      else if(!identical(dn0, e2@Dimnames) &&
		      !identical(dn,  e2@Dimnames)) {
		  dn <- dn0
		  warning("not using incompatible 'Dimnames' in arithmetical result")
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

## help(Math2)	mentions this uglyness:
setGeneric("round",  group="Math2")
setGeneric("signif", group="Math2")

setMethod("Math2",
	  signature(x = "dgeMatrix", digits = "numeric"),
	  function(x, digits) {
	      x@x <- callGeneric(x@x, digits = digits)
	      x
	  })

## TODO :  "Compare" -> returning  logical Matrices


## -- end{group generics} -----------------------


setMethod("norm", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) norm(x, type = "O", ...))

setMethod("norm", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dgeMatrix_norm", x, type),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) rcond(x, type = "O", ...))

setMethod("rcond", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dgeMatrix_rcond", x, type),
	  valueClass = "numeric")

setMethod("t", signature(x = "dgeMatrix"),
	  function(x) {
	      x@x <- as.vector(t(array(x@x, dim = x@Dim)))# no dimnames here!
	      x@Dim <- x@Dim[2:1]
	      x@Dimnames <- x@Dimnames[2:1]
	      x })

setMethod("crossprod", signature(x = "dgeMatrix", y = "missing"),
	  function(x, y = NULL) .Call("dgeMatrix_crossprod", x, FALSE),
	  valueClass = "dpoMatrix")

setMethod("tcrossprod", signature(x = "dgeMatrix"),
	  function(x) .Call("dgeMatrix_crossprod", x, TRUE),
	  valueClass = "dpoMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y = NULL) .Call("dgeMatrix_dgeMatrix_crossprod", x, y),
	  valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y = NULL) .Call("dgeMatrix_matrix_crossprod", x, y),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y = NULL)
	  .Call("dgeMatrix_matrix_crossprod", x, as.matrix(as.double(y))),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y),
	  valueClass = "dgeMatrix")
setMethod("crossprod", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y = NULL) callGeneric(as.matrix(as.double(x)), y),
	  valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "dgeMatrix", y = "dgeMatrix"),
	  function(x, y) .Call("dgeMatrix_matrix_mm", x, y, TRUE, FALSE),
	  valueClass = "dgeMatrix")

## dgeMatrix <-> matrix ("matrix" dispatches before "numeric" since R 2.1.0)
setMethod("%*%", signature(x = "dgeMatrix", y = "matrix"),
	  function(x, y) {
	      storage.mode(y) <- "double"
	      .Call("dgeMatrix_matrix_mm", x, y, FALSE, FALSE)
	  }, valueClass = "dgeMatrix")


setMethod("%*%", signature(x = "matrix", y = "dgeMatrix"),
	  function(x, y) {
	      storage.mode(x) <- "double"
	      .Call("dgeMatrix_matrix_mm", y, x, FALSE, TRUE)
	  }, valueClass = "dgeMatrix")

## dgeMatrix <-> numeric: conceptually dispatch to "matrix" one, but shortcut
setMethod("%*%", signature(x = "dgeMatrix", y = "numeric"),
	  function(x, y)
	  .Call("dgeMatrix_matrix_mm", x, as.matrix(as.double(y)), FALSE, FALSE),
	  valueClass = "dgeMatrix")
setMethod("%*%", signature(x = "numeric", y = "dgeMatrix"),
	  function(x, y)
	  .Call("dgeMatrix_matrix_mm", y, rbind(as.double(x)), FALSE, TRUE),
	  valueClass = "dgeMatrix")

setMethod("diag", signature(x = "dgeMatrix"),
	  function(x = 1, nrow, ncol = n)
	  .Call("dgeMatrix_getDiag", x))


## DB - I don't think this is a good idea without first checking symmetry
#setMethod("chol", signature(x = "dgeMatrix", pivot = "ANY"), cholMat)

setMethod("solve", signature(a = "dgeMatrix", b = "missing"),
	  function(a, b, ...) .Call("dgeMatrix_solve", a),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "dgeMatrix"),
	  function(a, b, ...) .Call("dgeMatrix_matrix_solve", a, b, TRUE),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "matrix"),
	  function(a, b, ...) {
	      storage.mode(b) <- "double"
	      .Call("dgeMatrix_matrix_solve", a, b, FALSE)
	  }, valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgeMatrix", b = "numeric"),
	  function(a, b, ...)
	  .Call("dgeMatrix_matrix_solve", a, as.matrix(as.double(b)), FALSE))

setMethod("lu", signature(x = "dgeMatrix"),
	  function(x, ...) .Call("dgeMatrix_LU", x), valueClass = "LU")

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "missing"),
	  function(x, logarithm, ...)
	  .Call("dgeMatrix_determinant", x, TRUE))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
	  function(x, logarithm, ...)
	  .Call("dgeMatrix_determinant", x, logarithm))


setMethod("expm", signature(x = "dgeMatrix"),
	  function(x) .Call("dgeMatrix_exp", x),
	  valueClass = "dgeMatrix")

setMethod("colSums", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call("dgeMatrix_colsums", x, na.rm, TRUE, FALSE),
	  valueClass = "numeric")

setMethod("colMeans", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call("dgeMatrix_colsums", x, na.rm, TRUE, TRUE),
	  valueClass = "numeric")

setMethod("rowSums", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call("dgeMatrix_colsums", x, na.rm, FALSE, FALSE),
	  valueClass = "numeric")

setMethod("rowMeans", signature(x = "dgeMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          .Call("dgeMatrix_colsums", x, na.rm, FALSE, TRUE),
	  valueClass = "numeric")

### The following all serve for	 as.Matrix()
### which is not yet exported (nor tested):

## utilities for Matrix.class() :

Hermitian.test <- function(x)
{
    if ((!inherits(x, "Matrix") && !is.matrix(x)) ||
	(nrow(x) != ncol(x))) return(Inf)
    if (is.complex(x)) return(max(Mod(x - t(Conj(x)))))
    max(x - t(x))
}

is.Hermitian <- function(x, tol = 0) { Hermitian.test(x) <= tol }

LowerTriangular.test <- function(x)
{
    if ((!inherits(x, "Matrix") && !is.matrix(x))) return(Inf)
    if (is.complex(x)) return(max(Mod(x[row(x) < col(x)])))
    max(abs(x[row(x) < col(x)]))
}

is.LowerTriangular <- function(x, tol = 0) { LowerTriangular.test(x) <= tol }

UpperTriangular.test <- function(x)
{
    if ((!inherits(x, "Matrix") && !is.matrix(x))) return(Inf)
    if (is.complex(x)) return(max(Mod(x[row(x) > col(x)])))
    max(abs(x[row(x) > col(x)]))
}

is.UpperTriangular <- function(x, tol = 0) { UpperTriangular.test(x) <= tol }

Orthogonal.test <- function(x, byrow = FALSE, normal = TRUE)
{
    if ((!inherits(x, "Matrix") && !is.matrix(x))) return(Inf)
    if (byrow) { x <- t(x) }
    xx <- crossprod(x)
    if (normal) {			# check for orthonormal
	return(max(Mod(xx[row(xx) > col(xx)]), Mod(diag(xx) - 1)))
    }
    max(Mod(xx[row(xx) > col(xx)]))
}

Orthonormal.test <- function(x, byrow = FALSE)
{
    Orthogonal.test(x, byrow, normal = TRUE)
}

is.ColOrthonormal <- function(x, tol = sqrt(.Machine$double.eps))
{
    Orthonormal.test(x, byrow = FALSE) <= tol
}

is.RowOrthonormal <- function(x, tol = sqrt(.Machine$double.eps))
{
    Orthonormal.test(x, byrow = TRUE) <= tol
}

is.Orthonormal <- function(x, tol = sqrt(.Machine$double.eps), byrow = FALSE)
{
    if (byrow) return(is.RowOrthonormal(x, tol))
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


as.Matrix <- function(x, tol = .Machine$double.eps)
{
    if(is(x, "Matrix")) return(x)
    ## else
    as(if(is.matrix(x)) x else as.matrix(x),
       Matrix.class(x, tol = tol))
}

