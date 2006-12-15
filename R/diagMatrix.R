#### All methods for "diagonalMatrix" and its subclasses,
####  currently "ddiMatrix", "ldiMatrix"

## Purpose: Constructor of diagonal matrices -- ~= diag() ,
##          but *not* diag() extractor!
Diagonal <- function(n, x = NULL)
{
    ## Allow  Diagonal(4)  and	Diagonal(x=1:5)
    if(missing(n))
	n <- length(x)
    else {
	stopifnot(length(n) == 1, n == as.integer(n), n >= 0)
	n <- as.integer(n)
    }

    if(missing(x)) ## unit diagonal matrix
	new("ddiMatrix", Dim = c(n,n), diag = "U")
    else {
	stopifnot(length(x) == n)
	if(is.logical(x))
	    cl <- "ldiMatrix"
	else if(is.numeric(x)) {
	    cl <- "ddiMatrix"
	    x <- as.numeric(x)
	}
	else if(is.complex(x)) {
	    cl <- "zdiMatrix"  # will not yet work
	} else stop("'x' has invalid data type")
	new(cl, Dim = c(n,n), diag = "N", x = x)
    }
}

### This is modified from a post of Bert Gunter to R-help on  1 Sep 2005.
### Bert's code built on a post by Andy Liaw who most probably was influenced
### by earlier posts, notably one by Scott Chasalow on S-news, 16 Jan 2002
### who posted his bdiag() function written in December 1995.

bdiag <- function(...) {
    if(nargs() == 0) return(new("dgCMatrix"))
    ## else :
    mlist <- if (nargs() == 1) as.list(...) else list(...)
    dims <- sapply(mlist, dim)
    ## make sure we had all matrices:
    if(!(is.matrix(dims) && nrow(dims) == 2))
	stop("some arguments are not matrices")
    csdim <- rbind(rep.int(0:0, 2),
                   apply(sapply(mlist, dim), 1, cumsum))
    ret <- new("dgTMatrix", Dim = as.integer(csdim[nrow(csdim),]))
    add1 <- matrix(1:0, 2,2)
    for(i in seq_along(mlist)) {
	indx <- apply(csdim[i:(i+1),] + add1, 2, function(n) n[1]:n[2])
	if(is.null(dim(indx))) ## non-square matrix
	    ret[indx[[1]],indx[[2]]] <- mlist[[i]]
	else ## square matrix
	    ret[indx[,1],indx[,2]] <- mlist[[i]]
    }
    ## slightly debatable if we really should return Csparse.. :
    as(ret, "CsparseMatrix")
}

diag2T <- function(from) {
    i <- if(from@diag == "U") integer(0) else seq_len(from@Dim[1]) - 1:1
    new(paste(.M.kind(from), "tTMatrix", sep=''),
	diag = from@diag, Dim = from@Dim, Dimnames = from@Dimnames,
	x = from@x, # <- ok for diag = "U" and "N" (!)
	i = i, j = i)
}

setAs("diagonalMatrix", "triangularMatrix", diag2T)
setAs("diagonalMatrix", "sparseMatrix", diag2T)
## is better than this:
## setAs("diagonalMatrix", "sparseMatrix",
##       function(from)
## 	  as(from, if(is(from, "dMatrix")) "dgCMatrix" else "lgCMatrix"))
setAs("diagonalMatrix", "CsparseMatrix",
      function(from) as(diag2T(from), "CsparseMatrix"))

setAs("diagonalMatrix", "matrix",
      function(from) {
          n <- from@Dim[1]
	  diag(x = if(from@diag == "U") { if(is.numeric(from@x)) 1. else TRUE
                                     } else from@x,
               nrow = n, ncol = n)
      })

setAs("diagonalMatrix", "generalMatrix", # prefer sparse:
      function(from) as(from, paste(.M.kind(from), "gCMatrix", sep='')))

.diag.x <- function(m) {
    if(m@diag == "U")
	rep.int(if(is.numeric(m@x)) 1. else TRUE,
		m@Dim[1])
    else m@x
}

.diag.2N <- function(m) {
    if(m@diag == "U") m@diag <- "N"
    m
}

## given the above, the following  4  coercions should be all unneeded;
## we prefer triangular to general:
setAs("ddiMatrix", "dgTMatrix",
      function(from) {
	  .Deprecated("as(, \"sparseMatrix\")")
	  n <- from@Dim[1]
	  i <- seq_len(n) - 1:1
	  new("dgTMatrix", i = i, j = i, x = .diag.x(from),
	      Dim = c(n,n), Dimnames = from@Dimnames) })

setAs("ddiMatrix", "dgCMatrix",
      function(from) as(as(from, "sparseMatrix"), "dgCMatrix"))

setAs("ldiMatrix", "lgTMatrix",
      function(from) {
	  .Deprecated("as(, \"sparseMatrix\")")
	  n <- from@Dim[1]
	  if(from@diag == "U") { # unit-diagonal
	      x <- rep.int(TRUE, n)
	      i <- seq_len(n) - 1:1
	  } else { # "normal"
	      nz <- nz.NA(from@x, na. = TRUE)
	      x <- from@x[nz]
	      i <- which(nz) - 1:1
	  }
	  new("lgTMatrix", i = i, j = i, x = x,
	      Dim = c(n,n), Dimnames = from@Dimnames) })

setAs("ldiMatrix", "lgCMatrix",
      function(from) as(as(from, "lgTMatrix"), "lgCMatrix"))


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
	  } ## TODO: complex
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

## When you assign to a diagonalMatrix, the result should be
## diagonal or sparse
setReplaceMethod("[", signature(x = "diagonalMatrix",
				i = "ANY", j = "ANY", value = "ANY"),
		 function(x, i, j, value) {
		     r <- callGeneric(x = as(x,"sparseMatrix"),
				      i=i, j=j, value=value)
		     if(isDiagonal(r)) as(r, "diagonalMatrix") else r
		 })


setMethod("t", signature(x = "diagonalMatrix"),
          function(x) { x@Dimnames <- x@Dimnames[2:1] ; x })

setMethod("isDiagonal", signature(object = "diagonalMatrix"),
          function(object) TRUE)
setMethod("isTriangular", signature(object = "diagonalMatrix"),
          function(object) TRUE)
setMethod("isSymmetric", signature(object = "diagonalMatrix"),
          function(object) TRUE)

setMethod("chol", signature(x = "ddiMatrix"),# pivot = "ANY"
	  function(x, pivot) {
	      if(any(x@x < 0)) stop("chol() is undefined for diagonal matrix with negative entries")
	      x@x <- sqrt(x@x)
	      x
	  })
## chol(L) is L for logical diagonal:
setMethod("chol", signature(x = "ldiMatrix"), function(x, pivot) x)


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
##       ---------------------
## Note that "ldi" logical are treated as numeric
diagdiagprod <- function(x, y) {
    if(any(dim(x) != dim(y))) stop("non-matching dimensions")
    if(x@diag != "U") {
	if(y@diag != "U") {
	    nx <- x@x * y@x
	    if(is.numeric(nx) && !is.numeric(x@x))
		x <- as(x, "dMatrix")
	    x@x <- as.numeric(nx)
	}
	return(x)
    } else ## x is unit diagonal
    return(y)
}

setMethod("%*%", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
	  diagdiagprod, valueClass = "ddiMatrix")

formals(diagdiagprod) <- alist(x=, y=x)
setMethod("crossprod", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
	  diagdiagprod, valueClass = "ddiMatrix")
setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
	  diagdiagprod, valueClass = "ddiMatrix")
setMethod("crossprod", signature(x = "diagonalMatrix", y = "missing"),
	  diagdiagprod, valueClass = "ddiMatrix")
setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "missing"),
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
setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "matrix"),
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
	      as(if(y@diag == "U") x else x * rep(y@x, each = dx[1]), "Matrix")
	  })

setMethod("%*%", signature(x = "dgeMatrix", y = "diagonalMatrix"),
	  function(x, y) {
	      dx <- dim(x)
	      dy <- dim(y)
	      if(dx[2] != dy[1]) stop("non-matching dimensions")
	      if(y@diag == "N")
		  x@x <- x@x * rep(y@x, each = dx[1])
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


### ---------------- diagonal  o  sparse  -----------------------------


## Use function for several signatures, in order to evade
## ambiguous dispatch for "ddi", since there's also Arith(ddense., ddense.)
diagOdiag <- function(e1,e2) { # result should also be diagonal
    r <- callGeneric(.diag.x(e1), .diag.x(e2)) # error if not "compatible"
    if(is.numeric(r)) {
        if(is.numeric(e2@x)) {
            e2@x <- r; return(.diag.2N(e2)) }
        if(!is.numeric(e1@x))
            ## e.g. e1, e2 are logical;
            e1 <- as(e1, "dMatrix")
    }
    else if(is.logical(r))
        e1 <- as(e1, "lMatrix")
    else stop("intermediate 'r' is of type", typeof(r))
    e1@x <- r
    .diag.2N(e1)
}

setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "diagonalMatrix"),
          diagOdiag)
## These two are just for method disambiguation:
setMethod("Ops", signature(e1 = "ddiMatrix", e2 = "diagonalMatrix"),
          diagOdiag)
setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "ddiMatrix"),
          diagOdiag)

## For almost everything else, diag* shall be treated "as sparse" :
## These are cheap implementations via coercion

## for disambiguation
setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "sparseMatrix"),
	  function(e1,e2) callGeneric(as(e1, "sparseMatrix"), e2))
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "diagonalMatrix"),
	  function(e1,e2) callGeneric(e1, as(e2, "sparseMatrix")))
## in general:
setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "ANY"),
          function(e1,e2) callGeneric(as(e1,"sparseMatrix"), e2))
setMethod("Ops", signature(e1 = "ANY", e2 = "diagonalMatrix"),
          function(e1,e2) callGeneric(e1, as(e2,"sparseMatrix")))



## FIXME?: In theory, this can be done *FASTER*, in some cases, via tapply1()
setMethod("%*%", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y) as(x, "sparseMatrix") %*% y)
## NB: The previous is *not* triggering for  "ddi" o "dgC" (= distance 3)
##     since there's a "ddense" o "Csparse" at dist. 2 => triggers first.
## ==> do this:
setMethod("%*%", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
	  function(x, y) as(x, "CsparseMatrix") %*% y)
## NB: this is *not* needed for Tsparse & Rsparse
## TODO: Write tests in ./tests/ which ensure that many "ops" with diagonal*
##       do indeed work by going throug sparse (and *not* ddense)!

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
	  function(object) {
	      d <- dim(object)
	      cl <- class(object)
	      cat(sprintf('%d x %d diagonal matrix of class "%s"\n',
			  d[1], d[2], cl))
	      prDiag(object)
	  })
