## Purpose: Constructor of diagonal matrices -- ~= diag() ,
##          but *not* diag() extractor!
Diagonal <- function(n, x = NULL)
{
    ## Allow  Diagonal(4)  and  Diagonal(x=1:5)
    if(missing(n))
        n <- length(x)
    else {
        stopifnot(length(n) == 1, n == as.integer(n), n >= 0)
        n <- as.integer(n)
    }

    if(missing(x)) # unit diagonal matrix
        new("ddiMatrix", Dim = c(n,n), diag = "U")
    else {
        stopifnot(length(x) == n)
        if(is.logical(x))
            cl <- "ldiMatrix"
        else {
            cl <- "ddiMatrix"
            x <- as.numeric(x)
        }
        new(cl, Dim = c(n,n), diag = "N", x = x)
    }
}

setAs("diagonalMatrix", "triangularMatrix",
      function(from) {
          n <- from@Dim[1]
          i <- seq(length = n)
          x <- from@x
          new(if(is.numeric(x)) "dtTMatrix" else "ltTMatrix",
              diag = from@diag, Dim = from@Dim, Dimnames = from@Dimnames,
              x = x, i = i, j = i)
          })

setAs("diagonalMatrix", "matrix",
      function(from) {
          n <- from@Dim[1]
	  diag(x = if(from@diag == "U") { if(is.numeric(from@x)) 1. else TRUE
                                     } else from@x,
               nrow = n, ncol = n)
      })

setAs("diagonalMatrix", "generalMatrix",
      function(from) {
          x <- as(from, "matrix")
          as(x,
             if(is.logical(x)) "lgeMatrix"
## Not yet:
##              else if(is.complex(x)) "zgeMatrix"
##              else if(is.integer(x)) "igeMatrix"
             else "dgeMatrix")
      })

setAs("ddiMatrix", "dgTMatrix",
      function(from) {
	  n <- from@Dim[1]
	  i <- seq(length = n) - 1:1
	  new("dgTMatrix", i = i, j = i,
	      x = if(from@diag == "U") rep(1,n) else from@x,
	      Dim = c(n,n), Dimnames = from@Dimnames) })

setAs("ddiMatrix", "dgCMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgCMatrix"))

setAs("ldiMatrix", "lgTMatrix",
      function(from) {
	  n <- from@Dim[1]
	  i <- (if(from@diag == "U") seq(length = n) else which(from@x)) - 1:1
	  new("lgTMatrix", i = i, j = i,
	      Dim = c(n,n), Dimnames = from@Dimnames) })

setAs("ldiMatrix", "lgCMatrix",
      function(from) as(as(from, "lgTMatrix"), "lgCMatrix"))

setAs("diagonalMatrix", "sparseMatrix",
      function(from)
	  as(from, if(is(from, "dMatrix")) "dgCMatrix" else "lgCMatrix"))

if(FALSE) # now have faster  "ddense" -> "dge"
setAs("ddiMatrix", "dgeMatrix",
      function(from) as(as(from, "matrix"), "dgeMatrix"))

setAs("matrix", "diagonalMatrix",
      function(from) {
	  d <- dim(from)
	  if(d[1] != (n <- d[2])) stop("non-square matrix")
	  if(any(from[row(from) != col(from)] != 0))
	      stop("has non-zero off-diagonal entries")
	  x <- diag(from)
	  if(is.logical(x)) {
	      cl <- "ldiMatrix"
	      uni <- all(x)
	  } else {
	      cl <- "ddiMatrix"
	      uni <- all(x == 1)
	      storage.mode(x) <- "double"
	  }
	  new(cl, Dim = c(n,n), diag = if(uni) "U" else "N",
	      x = if(uni) x[FALSE] else x)
      })

## ``generic'' coercion to  diagonalMatrix : build on  isDiagonal() and diag()
setAs("Matrix", "diagonalMatrix",
      function(from) {
          d <- dim(from)
	  if(d[1] != (n <- d[2])) stop("non-square matrix")
          if(!isDiagonal(from)) stop("matrix is not diagonal")
          ## else:
          x <- diag(from)
          if(is.logical(x)) {
              cl <- "ldiMatrix"
              uni <- all(x)
          } else {
              cl <- "ddiMatrix"
              uni <- all(x == 1)
              storage.mode(x) <- "double"
          }
          new(cl, Dim = c(n,n), diag = if(uni) "U" else "N",
              x = if(uni) x[FALSE] else x)
      })

setMethod("t", signature(x = "diagonalMatrix"),
          function(x) { x@Dimnames <- x@Dimnames[2:1] ; x })

setMethod("isDiagonal", signature(object = "diagonalMatrix"),
          function(object) TRUE)
setMethod("isTriangular", signature(object = "diagonalMatrix"),
          function(object) TRUE)
setMethod("isSymmetric", signature(object = "diagonalMatrix"),
          function(object) TRUE)

setMethod("diag", signature(x = "diagonalMatrix"),
	  function(x = 1, nrow, ncol = n) {
             if(x@diag == "U")
                 rep.int(if(is.logical(x@x)) TRUE else 1, x@Dim[1])
             else x@x
          })

setMethod("!", "ldiMatrix", function(e1) {
    if(e1@diag == "N")
	e1@x <- !e1@x
    else { ## "U"
	e1@diag <- "N"
	e1@x <- rep.int(FALSE, e1@Dim[1])
    }
    x
})

## Basic Matrix Multiplication {many more to add}

## FIXME: extend this for 'ldi', i.e. do "diagonalMatrix"
diagdiagprod <- function(x, y) {
    if(any(dim(x) != dim(y))) stop("non-matching dimensions")
    if(x@diag != "U") {
        if(y@diag != "U") x@x <- x@x * y@x
        return(x)
    } else ## x is unit diagonal
    return(y)
}

setMethod("%*%", signature(x = "ddiMatrix", y = "ddiMatrix"),
	  diagdiagprod, valueClass = "ddiMatrix")

formals(diagdiagprod) <- alist(x=, y=NULL)
setMethod("crossprod", signature(x = "ddiMatrix", y = "ddiMatrix"),
	  diagdiagprod, valueClass = "ddiMatrix")
setMethod("tcrossprod", signature(x = "ddiMatrix", y = "ddiMatrix"),
	  diagdiagprod, valueClass = "ddiMatrix")


diagmatprod <- function(x, y) {
    dx <- dim(x)
    dy <- dim(y)
    if(dx[2] != dy[1]) stop("non-matching dimensions")
    n <- dx[1]
    as(if(x@diag == "U") y else x@x * y, "Matrix")
}

setMethod("%*%", signature(x = "diagonalMatrix", y = "matrix"),
          diagmatprod)
formals(diagmatprod) <- alist(x=, y=NULL)
setMethod("crossprod", signature(x = "diagonalMatrix", y = "matrix"),
          diagmatprod)

diagdgeprod <- function(x, y) {
    dx <- dim(x)
    dy <- dim(y)
    if(dx[2] != dy[1]) stop("non-matching dimensions")
    if(x@diag != "U")
        y@x <- x@x * y@x
    y
}
setMethod("%*%", signature(x = "diagonalMatrix", y = "dgeMatrix"),
	  diagdgeprod, valueClass = "dgeMatrix")
formals(diagdgeprod) <- alist(x=, y=NULL)
setMethod("crossprod", signature(x = "diagonalMatrix", y = "dgeMatrix"),
	  diagdgeprod, valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "diagonalMatrix"),
	  function(x, y) {
              dx <- dim(x)
              dy <- dim(y)
              if(dx[2] != dy[1]) stop("non-matching dimensions")
              as(if(y@diag == "U") x else x * rep.int(y@x, dx[1]), "Matrix")
          })

setMethod("%*%", signature(x = "dgeMatrix", y = "diagonalMatrix"),
	  function(x, y) {
              dx <- dim(x)
              dy <- dim(y)
              if(dx[2] != dy[1]) stop("non-matching dimensions")
              if(y@diag == "N")
                  x@x <- x@x * rep.int(y@x, dx[1])
              x
          })

## crossprod {more of these}

## tcrossprod --- all are not yet there: do the dense ones here:

## FIXME:
## setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "denseMatrix"),
## 	  function(x, y = NULL) {
##           })

## setMethod("tcrossprod", signature(x = "denseMatrix", y = "diagonalMatrix"),
## 	  function(x, y = NULL) {
##           })


### ---------------- diagonal  o   sparse  -----------------------------

## These are cheap implementations via coercion

## FIXME?: In theory, this can be done *FASTER*, in some cases, via tapply1()

setMethod("%*%", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y) as(x, "sparseMatrix") %*% y)

setMethod("%*%", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y) x %*% as(y, "sparseMatrix"))

setMethod("crossprod", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y = NULL) { x <- as(x, "sparseMatrix"); callGeneric() })

setMethod("crossprod", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y = NULL) { y <- as(y, "sparseMatrix"); callGeneric() })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y = NULL) { x <- as(x, "sparseMatrix"); callGeneric() })

setMethod("tcrossprod", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y = NULL) { y <- as(y, "sparseMatrix"); callGeneric() })




## similar to prTriang() in ./Auxiliaries.R :
prDiag <-
    function(x, digits = getOption("digits"), justify = "none", right = TRUE)
{
    cf <- array(".", dim = x@Dim, dimnames = x@Dimnames)
    cf[row(cf) == col(cf)] <-
        sapply(diag(x), format, digits = digits, justify = justify)
    print(cf, quote = FALSE, right = right)
    invisible(x)
}

setMethod("show", signature(object = "diagonalMatrix"),
          function(object) prDiag(object))

