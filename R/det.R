det.default <- function (x, method = c("qr", "eigenvalues"), ...)
{
    ## old version - in case anyone depends on it
    if (!is.matrix(x) || (n <- ncol(x)) != nrow(x))
        stop("x must be a square matrix")
    method <- match.arg(method)
    if (method == "qr") {
        x <- prod(diag(qr(x)$qr))
        if (n%%2 == 1)
            x
        else -x
    }
    else Re(prod(eigen(x, only.values = TRUE)$values))
}

det.Matrix <- function(x, logarithm = TRUE, ...)
{
    .Call("R_LapackPP_det", x, as.logical(logarithm), PACKAGE="Matrix")
}

det.UnitLowerTriangular <- function(x, logarithm = TRUE, ...)
{
    logarithm <- as.logical(logarithm[1])
    asObject(list(modulus =
                 structure(ifelse(logarithm, 0., 1.), logarithm = logarithm),
                 sign = 1),
            call = match.call(),
            c("det.UnitLowerTriangular", "det"))
}

det.UnitUpperTriangular <- function(x, logarithm = TRUE, ...)
{
    logarithm <- as.logical(logarithm[1])
    asObject(list(modulus =
                 structure(ifelse(logarithm, 0., 1.), logarithm = logarithm),
                 sign = 1),
            call = match.call(),
            c("det.UnitUpperTriangular", "det"))
}

## calculate the determinant of a triangular matrix from its diagonal
diagDet <- function(x, logarithm = TRUE, ...)
{
    logarithm <- as.logical(logarithm)[1]
    asObject(list(modulus =
                 structure(if (logarithm) sum(log(abs(x))) else prod(abs(x)),
                           logarithm = logarithm),
                 sign = prod(sign(x))),
            "det")
}

det.LowerTriangular <- function(x, logarithm = TRUE, ...)
    asObject(diagDet(x, logarithm), c("det.LowerTriangular", "det"))

det.UpperTriangular <- function(x, logarithm = TRUE, ...)
    asObject(diagDet(x, logarithm), c("det.UpperTriangular", "det"))
