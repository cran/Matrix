setMethod("show", signature(object = "Matrix"),
          function(object) print(as(object, "matrix")))

Matrix <-
    function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL)
{
    if (is(data, "Matrix")) return(data)
    if (is.matrix(data)) { val <- data }
    else {
        if (missing(nrow))
            nrow <- ceiling(length(data)/ncol)
        else if (missing(ncol))
            ncol <- ceiling(length(data)/nrow)
        val <- .Internal(matrix(data, nrow, ncol, byrow))
        dimnames(val) <- dimnames
    }
    as(val, "geMatrix")
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
