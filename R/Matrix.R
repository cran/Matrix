prMatrix <-
    ## private function to be used as show() method possibly more than once
    function(object) {
        d <- dim(object)
        cat(paste(d, collapse= " x "), " Matrix of class ",
            sQuote(class(object)),"\n", sep='')
        m <- as(object, "matrix")
        maxp <- getOption("max.print")
        if(prod(d) <= maxp) print(m)
        else { ## d[1] > maxp / d[2] >= nr :
            nr <- maxp %/% d[2]
            n2 <- ceiling(nr / 2)
            print(head(m, max(1, n2)))
            cat("\n ..........\n\n")
            print(tail(m, max(1, nr - n2)))
        }
        ## DEBUG: cat("str(.):\n") ; str(object)
        invisible()
    }

setMethod("show", signature(object = "Matrix"), prMatrix)

if(FALSE) {## FIXME: we should do this here (for all subclasses),
    ##        -----  but it coerces some to "Matrix" {with no @x slot}
setMethod("dim", signature(x = "Matrix"),
          function(x) x@Dim, valueClass = "integer")
setMethod("dimnames", signature(x = "Matrix"), function(x) x@Dimnames)
}# FIXME

Matrix <-
    function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL)
{
    if (is(data, "Matrix")) return(data)
    if (is.matrix(data)) { val <- data }
    else { ## cut & paste from "base::matrix" :
        if (missing(nrow))
            nrow <- ceiling(length(data)/ncol)
        else if (missing(ncol))
            ncol <- ceiling(length(data)/nrow)
        val <- .Internal(matrix(data, nrow, ncol, byrow))
        dimnames(val) <- dimnames
    }
    as(val, "dgeMatrix")
}


if(FALSE) { ##--- not-yet used -- {almost same code also in ./dgeMatrix.R }

## utility for as.Matrix() {which is currently invalid }
Matrix.class <- function(x, tol = 0, symmetry = TRUE, unit.diagonal = TRUE,
                         triangularity = c(TRUE, TRUE),
                         orthogonality = c(TRUE, TRUE),
                         normality = c(TRUE, TRUE))
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
    asObject(if (inherits(x, "Matrix")) x else as.matrix(x),
	     Matrix.class(x, tol = tol))
}

}## not-yet used
