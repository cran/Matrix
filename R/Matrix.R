Matrix <-
    function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL)
{
    if (inherits(data, "Matrix")) return(data)
    if (is.matrix(data)) { val <- data }
    else {
        if (missing(nrow))
            nrow <- ceiling(length(data)/ncol)
        else if (missing(ncol))
            ncol <- ceiling(length(data)/nrow)
        val <- .Internal(matrix(data, nrow, ncol, byrow))
        dimnames(val) <- dimnames
    }
    class(val) <- "Matrix"
    val
}

print.Matrix <- function(x, ...)
{
    print(unclass(x), ...)
}

as.matrix.Matrix <- function(x)
{
    unclass(unpack(x))
}

solve.Matrix <- function(a, b)
{   ## short version of a solve method
###    if (missing(b)) return(.Call("R_LapackPP_solve1", a))
    if (missing(b)) {
        m <- nrow(a)
        if (m != ncol(a)) stop("only square matrices can be inverted")
        return(.Call("R_LapackPP_solve", a, diag(m)))
    }
    .Call("R_LapackPP_solve", a, b)
}

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
    max(x[row(x) < col(x)])
}

is.LowerTriangular <- function(x, tol = 0) { LowerTriangular.test(x) <= tol }

UpperTriangular.test <- function(x)
{
    if ((!inherits(x, "Matrix") && !is.matrix(x))) return(Inf)
    if (is.complex(x)) return(max(Mod(x[row(x) > col(x)])))
    max(x[row(x) > col(x)])
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
                val <- c("ColOrthogonal")
        }
    }
    if (orthogonality[2]) {
        if (normality[2] && is.RowOrthonormal(x, tol)) {
            val <- c("ColOrthoNormal", "ColOrthogonal", val)
        } else {
            if (Orthogonal.test(x, byrow = TRUE, normal = FALSE) <= tol)
                val <- c("ColOrthogonal")
        }
    }
    val
}
