#### Toplevel ``virtual'' class "Matrix"

## probably not needed eventually:
setAs(from = "ddenseMatrix", to = "matrix",
      function(from) {
          if(length(d <- dim(from)) != 2) stop("dim(.) has not length 2")
          array(from@x, dim = d, dimnames = dimnames(from))
      })

## private function to be used as show() method possibly more than once
prMatrix <- function(object) {
    d <- dim(object)
    cl <- class(object)
    cat(sprintf('%d x %d Matrix of class "%s"\n', d[1], d[2], cl))
##- no longer needed: have no objects of virtual classes:
##     if(cl == "Matrix") { ## have no data slot
##         cat("Dim = ", d)
##         if(any(sapply(object@Dimnames,length) > 0)) {
##             cat("; Dimnames = ")
##             str(object@Dimnames)
##         }
##         cat("\n")
##     } else { # not "Matrix", hence have data 'x' slot
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
##    }
    invisible(object)# as print() S3 methods do
}

setMethod("show", signature(object = "ddenseMatrix"), prMatrix)

setMethod("show", signature(object = "sparseMatrix"),
   function(object) {
       d <- dim(object)
       cl <- class(object)
       cat(sprintf('%d x %d sparse Matrix of class "%s"\n', d[1], d[2], cl))

       maxp <- getOption("max.print")
       if(prod(d) <= maxp) print(as(object, "matrix"))
       else { ## d[1] > maxp / d[2] >= nr :
           cat("\n Not printing large sparse matrix -- maybe increase options(max.print)\n")
           if(FALSE) { ### need storage economic "[,]" method for sparse!!
               nr <- maxp %/% d[2]
               n2 <- ceiling(nr / 2)
               print(head(m, max(1, n2)))
               cat("\n ..........\n\n")
               print(tail(m, max(1, nr - n2)))
           }
       }
        ## DEBUG: cat("str(.):\n") ; str(object)
       invisible(object)
   })

## this may go away {since sparse matrices need something better!} :
setMethod("show", signature(object = "Matrix"), prMatrix)

## should propagate to all subclasses:
setMethod("as.matrix", signature(x = "Matrix"), function(x) as(x, "matrix"))

setMethod("dim", signature(x = "Matrix"),
          function(x) x@Dim, valueClass = "integer")
setMethod("dimnames", signature(x = "Matrix"), function(x) x@Dimnames)
## not exported but used more than once for "dimnames<-" method :
## -- or do only once for all "Matrix" classes ??
dimnamesGets <- function (x, value) {
    d <- dim(x)
    if (!is.list(value) || length(value) != 2 ||
        !(is.null(v1 <- value[[1]]) || length(v1) == d[1]) ||
        !(is.null(v2 <- value[[2]]) || length(v2) == d[2]))
        stop(sprintf("invalid dimnames given for '%s' object", class(x)))
    x@Dimnames <- list(if(!is.null(v1)) as.character(v1),
                       if(!is.null(v2)) as.character(v2))
    x
}
setMethod("dimnames<-", signature(x = "Matrix", value = "list"),
          dimnamesGets)

setMethod("unname", signature("Matrix", force="missing"),
          function(obj) { obj@Dimnames <- list(NULL,NULL); obj})

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

## Methods for operations where one argument is numeric

setMethod("%*%", signature(x = "Matrix", y = "numeric"),
          function(x, y) callGeneric(x, array(y, c(length(y), 1))))

setMethod("%*%", signature(x = "numeric", y = "Matrix"),
          function(x, y) callGeneric(array(x, c(1, length(x))), y))

setMethod("crossprod", signature(x = "Matrix", y = "numeric"),
          function(x, y = NULL) callGeneric(x, array(y, c(length(y), 1))))

setMethod("crossprod", signature(x = "numeric", y = "Matrix"),
          function(x, y = NULL)  callGeneric(array(x, c(1, length(x))), y))

setMethod("solve", signature(a = "Matrix", b = "numeric"),
          function(a, b, ...) callGeneric(a, array(b, c(length(b), 1))))

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
