setAs("matrix", "geMatrix",
      function(from) {
          new("geMatrix", x = c(from), Dim = as.integer(attr(from, "dim")))
      })

setAs("geMatrix", "matrix",
      function(from) {
          array(from@x, from@Dim)
      })

setMethod("norm", signature(x = "geMatrix", type = "missing"),
          function(x, type, ...) norm(x, type = "O", ...))

setMethod("norm", signature(x = "geMatrix", type = "character"),
          function(x, type, ...)
          .Call("geMatrix_norm", x, type),
          valueClass = "numeric")

setMethod("rcond", signature(x = "geMatrix", type = "missing"),
          function(x, type, ...) rcond(x, type = "O", ...))

setMethod("rcond", signature(x = "geMatrix", type = "character"),
          function(x, type, ...)
          .Call("geMatrix_rcond", x, type),
          valueClass = "numeric")

setMethod("crossprod", signature(x = "geMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("geMatrix_crossprod", x),
          valueClass = "poMatrix")

setMethod("crossprod", signature(x = "geMatrix", y = "geMatrix"),
          function(x, y = NULL)
          .Call("geMatrix_geMatrix_crossprod", x, y),
          valueClass = "geMatrix")

setMethod("crossprod", signature(x = "geMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("geMatrix_matrix_crossprod", x, y),
          valueClass = "geMatrix")

setMethod("crossprod", signature(x = "geMatrix", y = "numeric"),
          function(x, y = NULL)
          .Call("geMatrix_matrix_crossprod", x, as.matrix(y)),
          valueClass = "geMatrix")

setMethod("diag", signature(x = "geMatrix"),
          function(x = 1, nrow, ncol = n)
          .Call("geMatrix_getDiag", x))

setMethod("dim", signature(x = "geMatrix"),
          function(x) x@Dim, valueClass = "integer")

setMethod("solve", signature(a = "geMatrix", b = "missing"),
          function(a, b, ...) .Call("geMatrix_solve", a)
          )

setMethod("solve", signature(a = "geMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("geMatrix_matrix_solve", a, b)
          )

setMethod("lu", signature(x = "geMatrix"),
          function(x, ...) .Call("geMatrix_LU", x))

setMethod("determinant", signature(x = "geMatrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "geMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
          .Call("geMatrix_determinant", x, logarithm))

setMethod("%*%", signature(x = "geMatrix", y = "geMatrix"),
          function(x, y)
          .Call("geMatrix_geMatrix_mm", x, y))

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

Matrix.class <- function(x, tol = 0, symmetry = TRUE, unit.diagonal = TRUE,
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
        } else {
            if (Orthogonal.test(x, normal = FALSE) <= tol)
                val <- c("ColOrthogonal", val)
        }
    }
    if (orthogonality[2]) {
        if (normality[2] && is.RowOrthonormal(x, tol)) {
            val <- c("RowOrthoNormal", "RowOrthogonal", val)
        } else {
            if (Orthogonal.test(x, byrow = TRUE, normal = FALSE) <= tol)
                val <- c("RowOrthogonal", val)
        }
    }
    val
}

as.Matrix <- function(x, tol = .Machine$double.eps)
{
    if (inherits(x, "Matrix")) return(asObject(x, Matrix.class(x, tol = tol)))
    asObject(as.matrix(x), Matrix.class(x, tol = tol))
}

