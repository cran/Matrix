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
              d1 <- e1@Dim
              d2 <- e2@Dim
              eqD <- d1 == d2
              if (!eqD[1])
                  stop("Matrices must have same number of rows for arithmetic")
              if (eqD[2])
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
              ##
              new("dgeMatrix", Dim = d, Dimnames = dn,
                  x = callGeneric(e1@x, e2@x))
          })

setMethod("Arith",
          signature(e1 = "dgeMatrix", e2 = "numeric"),
          ## could use this function twice  for   A o B  and   B o A
          ## if callGeneric(.) became smarter:
          function(e1, e2) {
              d <- e1@Dim
              le <- length(e2)
              if(le == 1 || le == d[1] || prod(d) == le) { # matching dim
                  new("dgeMatrix", Dim = d, Dimnames = e1@Dimnames,
                      x = callGeneric(as.vector(e2), e1@x))
              } else stop ("length of 2nd arg does not match dimension of first")
          })

setMethod("Arith",
          signature(e1 = "numeric", e2 = "dgeMatrix"),
          function(e1, e2) {
              d <- e2@Dim
              le <- length(e1)
              if(le == 1 || le == d[1] || prod(d) == le) { # matching dim
                  new("dgeMatrix", Dim = d, Dimnames = e2@Dimnames,
                      x = callGeneric(as.vector(e1), e2@x))
              } else stop ("length of 1st arg does not match dimension of 2nd")
          })

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

setMethod("crossprod", signature(x = "dgeMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("dgeMatrix_crossprod", x),
          valueClass = "dpoMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "dgeMatrix"),
          function(x, y = NULL)
          .Call("dgeMatrix_dgeMatrix_crossprod", x, y),
          valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("dgeMatrix_matrix_crossprod", x, y),
          valueClass = "dgeMatrix")

setMethod("crossprod", signature(x = "dgeMatrix", y = "numeric"),
          function(x, y = NULL)
          .Call("dgeMatrix_matrix_crossprod", x, as.matrix(y)),
          valueClass = "dgeMatrix")

setMethod("diag", signature(x = "dgeMatrix"),
          function(x = 1, nrow, ncol = n)
          .Call("dgeMatrix_getDiag", x))

setMethod("dim", signature(x = "dgeMatrix"),
          function(x) x@Dim, valueClass = "integer")

setMethod("solve", signature(a = "dgeMatrix", b = "missing"),
          function(a, b, ...) .Call("dgeMatrix_solve", a)
          )

setMethod("solve", signature(a = "dgeMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("dgeMatrix_matrix_solve", a, b)
          )

setMethod("lu", signature(x = "dgeMatrix"),
          function(x, ...) .Call("dgeMatrix_LU", x))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "dgeMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
          .Call("dgeMatrix_determinant", x, logarithm))

setMethod("%*%", signature(x = "dgeMatrix", y = "dgeMatrix"),
          function(x, y)
          .Call("dgeMatrix_dgeMatrix_mm", x, y))

setMethod("expm", signature(x = "dgeMatrix"),
          function(x) .Call("dgeMatrix_exp", x),
          valueClass = "dgeMatrix")

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
    if (normal) {                       # check for orthonormal
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

Matrix.class <-
    function(x, tol = 0, symmetry = TRUE, unit.diagonal = TRUE,
             triangularity = c(TRUE, TRUE),
             orthogonality = c(TRUE, TRUE), normality = c(TRUE, TRUE))
{
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
        if (is.ColOrthonormal(x, tol)) {
            val <- c("ColOrthoNormal", "ColOrthogonal", val)
        } else if (Orthogonal.test(x, normal = FALSE) <= tol) {
            val <- c("ColOrthogonal", val)
        }
    }
    if (orthogonality[2]) {
        if (normality[2] && is.RowOrthonormal(x, tol)) {
            val <- c("RowOrthoNormal", "RowOrthogonal", val)
        } else if (Orthogonal.test(x, byrow = TRUE, normal = FALSE) <= tol) {
            val <- c("RowOrthogonal", val)
        }
    }
    val
}

as.Matrix <- function(x, tol = .Machine$double.eps)
{
    asObject(if (inherits(x, "Matrix")) x else as.matrix(x),
	     Matrix.class(x, tol = tol))
}
