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
	lx <- length(x)
	stopifnot(lx == 1 || lx == n) # but keep 'x' short for now
	if(is.logical(x))
	    cl <- "ldiMatrix"
	else if(is.numeric(x)) {
	    cl <- "ddiMatrix"
	    x <- as.numeric(x)
	}
	else if(is.complex(x)) {
	    cl <- "zdiMatrix"  # will not yet work
	} else stop("'x' has invalid data type")
	new(cl, Dim = c(n,n), diag = "N",
	    x = if(lx == 1) rep.int(x,n) else x)
    }
}

## Pkg 'spdep' had (relatively slow) versions of this as_dsCMatrix_I()
.symDiagonal <- function(n, x = rep.int(1,n), uplo = "U") {
    stopifnot(n == (n. <- as.integer(n)), (n <- n.) >= 0)
    if((lx <- length(x)) == 1) x <- rep.int(x, n)
    else if(lx != n) stop("length(x) must be 1 or n")
    cls <-
        if(is.double(x)) "dsCMatrix"
        else if(is.logical(x)) "lsCMatrix"
        else { ## for now
            storage.mode(x) <- "double"
            "dsCMatrix"
        }
    new(cls, Dim = c(n,n), x = x, uplo = uplo,
        i = if(n) 0:(n - 1L) else integer(0), p = 0:n)
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
    csdim <- rbind(rep.int(0L, 2),
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


.diag2tT <- function(from, uplo = "U", kind = .M.kind(from)) {
    ## to triangular Tsparse
    i <- if(from@diag == "U") integer(0) else seq_len(from@Dim[1]) - 1L
    new(paste(kind, "tTMatrix", sep=''),
	diag = from@diag, Dim = from@Dim, Dimnames = from@Dimnames,
	uplo = uplo,
	x = from@x, # <- ok for diag = "U" and "N" (!)
	i = i, j = i)
}

.diag2sT <- function(from, uplo = "U", kind = .M.kind(from)) {
    ## to symmetric Tsparse
    n <- from@Dim[1]
    i <- seq_len(n) - 1L
    new(paste(kind, "sTMatrix", sep=''),
	Dim = from@Dim, Dimnames = from@Dimnames,
	i = i, j = i, uplo = uplo,
	x = if(from@diag == "N") from@x else ## "U"-diag
	rep.int(switch(kind,
		       "d" = 1.,
		       "l" =,
		       "n" = TRUE,
		       ## otherwise
		       stop("'", kind,"' kind not yet implemented")), n))
}

## diagonal -> triangular,  upper / lower depending on "partner":
diag2tT.u <- function(d, x, kind = .M.kind(d))
    .diag2tT(d, uplo = if(is(x,"triangularMatrix")) x@uplo else "U", kind)


## In order to evade method dispatch ambiguity warnings,
## and because we can save a .M.kind() call, we use this explicit
## "hack"  instead of signature  x = "diagonalMatrix" :
##
## ddi*:
diag2tT <- function(from) .diag2tT(from, "U", "d")
setAs("ddiMatrix", "triangularMatrix", diag2tT)
setAs("ddiMatrix", "sparseMatrix", diag2tT)
## needed too (otherwise <dense> -> Tsparse is taken):
setAs("ddiMatrix", "TsparseMatrix", diag2tT)
setAs("ddiMatrix", "CsparseMatrix",
      function(from) as(.diag2tT(from, "U", "d"), "CsparseMatrix"))
setAs("ddiMatrix", "symmetricMatrix",
      function(from) .diag2sT(from, "U", "d"))
##
## ldi*:
diag2tT <- function(from) .diag2tT(from, "U", "l")
setAs("ldiMatrix", "triangularMatrix", diag2tT)
setAs("ldiMatrix", "sparseMatrix", diag2tT)
## needed too (otherwise <dense> -> Tsparse is taken):
setAs("ldiMatrix", "TsparseMatrix", diag2tT)
setAs("ldiMatrix", "CsparseMatrix",
      function(from) as(.diag2tT(from, "U", "l"), "CsparseMatrix"))
setAs("ldiMatrix", "symmetricMatrix",
      function(from) .diag2sT(from, "U", "l"))


setAs("diagonalMatrix", "nMatrix",
      function(from) {
	  n <- from@Dim[1]
	  i <- if(from@diag == "U") integer(0) else which(isN0(from@x)) - 1L
	  new("ntTMatrix", i = i, j = i, diag = from@diag,
	      Dim = from@Dim, Dimnames = from@Dimnames)
      })


setAs("diagonalMatrix", "matrix",
      function(from) {
          n <- from@Dim[1]
	  diag(x = if(from@diag == "U") { if(is.numeric(from@x)) 1. else TRUE
                                     } else from@x,
               nrow = n, ncol = n)
      })

setMethod("as.vector", signature(x = "diagonalMatrix", mode="missing"),
	  function(x, mode) {
	      n <- x@Dim[1]
	      mod <- mode(x@x)
	      r <- vector(mod, length = n^2)
	      if(n)
		  r[1 + 0:(n - 1) * (n + 1)] <-
		      if(x@diag == "U")
			  switch(mod, "integer"= 1L,
				 "numeric"= 1, "logical"= TRUE)
		      else x@x
	      r
	  })

setAs("diagonalMatrix", "generalMatrix", # prefer sparse:
      function(from) as(as(from, "CsparseMatrix"), "generalMatrix"))

.diag.x <- function(m) {
    if(m@diag == "U")
	rep.int(if(is.numeric(m@x)) 1. else TRUE, m@Dim[1])
    else m@x
}

.diag.2N <- function(m) {
    if(m@diag == "U") m@diag <- "N"
    m
}

if(FALSE) {
## given the above, the following  4  coercions should be all unneeded;
## we prefer triangular to general:
setAs("ddiMatrix", "dgTMatrix",
      function(from) {
	  .Deprecated("as(, \"sparseMatrix\")")
	  n <- from@Dim[1]
	  i <- seq_len(n) - 1L
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
	      i <- seq_len(n) - 1L
	  } else { # "normal"
	      nz <- nz.NA(from@x, na. = TRUE)
	      x <- from@x[nz]
	      i <- which(nz) - 1L
	  }
	  new("lgTMatrix", i = i, j = i, x = x,
	      Dim = c(n,n), Dimnames = from@Dimnames) })

setAs("ldiMatrix", "lgCMatrix",
      function(from) as(as(from, "lgTMatrix"), "lgCMatrix"))
}


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


## In order to evade method dispatch ambiguity warnings,
## we use this hack instead of signature  x = "diagonalMatrix" :
diCls <- names(getClass("diagonalMatrix")@subclasses)
for(cls in diCls) {
    setMethod("diag", signature(x = cls),
	      function(x = 1, nrow, ncol) .diag.x(x))
}


subDiag <- function(x, i, j, ..., drop) {
    x <- as(x, "sparseMatrix")
    x <- if(missing(i))
	x[, j, drop=drop]
    else if(missing(j))
	x[i, , drop=drop]
    else
	x[i,j, drop=drop]
    if(isS4(x) && isDiagonal(x)) as(x, "diagonalMatrix") else x
}

setMethod("[", signature(x = "diagonalMatrix", i = "index",
			 j = "index", drop = "logical"), subDiag)
setMethod("[", signature(x = "diagonalMatrix", i = "index",
			j = "missing", drop = "logical"),
	  function(x, i, j, ..., drop) subDiag(x, i=i, drop=drop))
setMethod("[", signature(x = "diagonalMatrix", i = "missing",
			 j = "index", drop = "logical"),
	  function(x, i, j, ..., drop) subDiag(x, j=j, drop=drop))

## When you assign to a diagonalMatrix, the result should be
## diagonal or sparse ---
## FIXME: this now fails because the "denseMatrix" methods come first in dispatch
## Only(?) current bug:  x[i] <- value  is wrong when  i is *vector*
replDiag <- function(x, i, j, ..., value) {
    x <- as(x, "sparseMatrix")
    if(missing(i))
	x[, j] <- value
    else if(missing(j)) { ##  x[i , ] <- v  *OR*   x[i] <- v
        na <- nargs()
##         message("diagnosing replDiag() -- nargs()= ", na)
	if(na == 4)
            x[i, ] <- value
	else if(na == 3)
            x[i] <- value
        else stop("Internal bug: nargs()=",na,"; please report")
    } else
	x[i,j] <- value
    if(isDiagonal(x)) as(x, "diagonalMatrix") else x
}

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index",
				j = "index", value = "replValue"), replDiag)

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index",
				j = "missing", value = "replValue"),
		 function(x,i,j, ..., value) {
                     ## message("before replDiag() -- nargs()= ", nargs())
                     if(nargs() == 3)
                         replDiag(x, i=i, value=value)
                     else ## nargs() == 4 :
                         replDiag(x, i=i, , value=value)
                 })

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "matrix", # 2-col.matrix
				j = "missing", value = "replValue"),
		 function(x,i,j, ..., value) {
		     if(ncol(i) == 2) {
			 if(all((ii <- i[,1]) == i[,2])) { # replace in diagonal only
			     x@x[ii] <- value
			     x
			 } else { ## no longer diagonal, but remain sparse:
			     x <- as(x, "sparseMatrix")
			     x[i] <- value
			     x
			 }
		     }
		     else { # behave as "base R": use as if vector
			 x <- as(x, "matrix")
			 x[i] <- value
			 Matrix(x)
		     }
		 })

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing",
				j = "index", value = "replValue"),
		 function(x,i,j, ..., value) replDiag(x, j=j, value=value))


setMethod("t", signature(x = "diagonalMatrix"),
          function(x) { x@Dimnames <- x@Dimnames[2:1] ; x })

setMethod("isDiagonal", signature(object = "diagonalMatrix"),
          function(object) TRUE)
setMethod("isTriangular", signature(object = "diagonalMatrix"),
          function(object) TRUE)
setMethod("isSymmetric", signature(object = "diagonalMatrix"),
	  function(object, ...) TRUE)

setMethod("symmpart", signature(x = "diagonalMatrix"), function(x) x)
setMethod("skewpart", signature(x = "diagonalMatrix"), setZero)

setMethod("chol", signature(x = "ddiMatrix"),
	  function(x, pivot, ...) {
	      if(any(x@x < 0))
		  stop("chol() is undefined for diagonal matrix with negative entries")
	      x@x <- sqrt(x@x)
	      x
	  })
## chol(L) is L for logical diagonal:
setMethod("chol", signature(x = "ldiMatrix"), function(x, pivot, ...) x)

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

setMethod("crossprod", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y = NULL) crossprod(as(x, "sparseMatrix"), y))

setMethod("crossprod", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y = NULL) crossprod(x, as(y, "sparseMatrix")))

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y = NULL) tcrossprod(as(x, "sparseMatrix"), y))

setMethod("tcrossprod", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y = NULL) tcrossprod(x, as(y, "sparseMatrix")))


## FIXME?: In theory, this can be done *FASTER*, in some cases, via tapply1()
setMethod("%*%", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y) as(x, "sparseMatrix") %*% y)
## NB: The previous is *not* triggering for  "ddi" o "dgC" (= distance 3)
##     since there's a "ddense" o "Csparse" at dist. 2 => triggers first.
## ==> do this:
setMethod("%*%", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
	  function(x, y) as(x, "CsparseMatrix") %*% y)
setMethod("%*%", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
	  function(x, y) x %*% as(y, "CsparseMatrix"))
## NB: this is *not* needed for Tsparse & Rsparse
## TODO: Write tests in ./tests/ which ensure that many "ops" with diagonal*
##       do indeed work by going through sparse (and *not* ddense)!

setMethod("%*%", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y) x %*% as(y, "sparseMatrix"))


setMethod("solve", signature(a = "diagonalMatrix", b = "missing"),
	  function(a, b, ...) {
	      a@x <- 1/ a@x
	      a@Dimnames <- a@Dimnames[2:1]
	      a
	  })

solveDiag <- function(a, b, ...) {
    if((n <- a@Dim[1]) != nrow(b))
        stop("incompatible matrix dimensions")
    ## trivially invert a 'in place' and multiply:
    a@x <- 1/ a@x
    a@Dimnames <- a@Dimnames[2:1]
    a %*% b
}
setMethod("solve", signature(a = "diagonalMatrix", b = "matrix"),
          solveDiag)
setMethod("solve", signature(a = "diagonalMatrix", b = "Matrix"),
          solveDiag)

## Schur()  ---> ./eigen.R



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

## FIXME:    diagonal  o  triangular  |-->  triangular
## -----     diagonal  o  symmetric   |-->  symmetric
##    {also when other is sparse: do these "here" --
##     before conversion to sparse, since that loses "diagonality"}

## For almost everything else, diag* shall be treated "as sparse" :
## These are cheap implementations via coercion

## For disambiguation --- define this for "sparseMatrix" , then for "ANY";
## and because we can save an .M.kind() call, we use this explicit
## "hack" for all diagonalMatrix *subclasses* instead of just "diagonalMatrix" :
##
## ddi*:
setMethod("Ops", signature(e1 = "ddiMatrix", e2 = "sparseMatrix"),
	  function(e1,e2) callGeneric(diag2tT.u(e1,e2, "d"), e2))
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "ddiMatrix"),
	  function(e1,e2) callGeneric(e1, diag2tT.u(e2,e1, "d")))
## ldi*
setMethod("Ops", signature(e1 = "ldiMatrix", e2 = "sparseMatrix"),
	  function(e1,e2) callGeneric(diag2tT.u(e1,e2, "l"), e2))
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "ldiMatrix"),
	  function(e1,e2) callGeneric(e1, diag2tT.u(e2,e1, "l")))

##  other = "numeric" : stay diagonal if possible
## ddi*:
setMethod("Ops", signature(e1 = "ddiMatrix", e2 = "numeric"),
	  function(e1,e2) {
	      n <- e1@Dim[1]; nsq <- n*n
	      f0 <- callGeneric(0, e2)
	      if(all(is0(f0))) { # remain diagonal
		  if(e1@diag == "U" && (r <- callGeneric(1, e2)) != 1)
		      e1@diag <- "N"
		  else
		      r <- callGeneric(e1@x, e2)
		  e1@x <- if(length(e2) == nsq) r else rep(r, length.out = nsq)
		  return(e1)
	      }
	      callGeneric(diag2tT.u(e1,e2, "d"), e2)
	  })

setMethod("Ops", signature(e1 = "numeric", e2 = "ddiMatrix"),
	  function(e1,e2) {
	      n <- e2@Dim[1]; nsq <- n*n
	      f0 <- callGeneric(e1, 0)
	      if(all(is0(f0))) { # remain diagonal
		  if(e2@diag == "U" && (r <- callGeneric(e1, 1)) != 1)
		      e2@diag <- "N"
		  else
		      r <- callGeneric(e1, e2@x)
		  e2@x <- if(length(e1) == nsq) r else rep(r, length.out = nsq)
		  return(e2)
	      }
	      callGeneric(e1, diag2tT.u(e2,e1, "d"))
	  })
## ldi*:
setMethod("Ops", signature(e1 = "ldiMatrix", e2 = "numeric"),
	  function(e1,e2) {
	      n <- e1@Dim[1]; nsq <- n*n
	      f0 <- callGeneric(FALSE, e2)
	      if(all(is0(f0))) { # remain diagonal
		  if(e1@diag == "U" && (r <- callGeneric(TRUE, e2)) != 1)
		      e1@diag <- "N"
		  else
		      r <- callGeneric(e1@x, e2)
		  e1@x <- if(length(e2) == nsq) r else rep(r, length.out = nsq)
		  return(e1)
	      }
	      callGeneric(diag2tT.u(e1,e2, "l"), e2)
	  })

setMethod("Ops", signature(e1 = "numeric", e2 = "ldiMatrix"),
	  function(e1,e2) {
	      n <- e2@Dim[1]; nsq <- n*n
	      f0 <- callGeneric(e1, FALSE)
	      if(all(is0(f0))) { # remain diagonal
		  if(e2@diag == "U" && (r <- callGeneric(e1, TRUE)) != 1)
		      e2@diag <- "N"
		  else
		      r <- callGeneric(e1, e2@x)
		  e2@x <- if(length(e1) == nsq) r else rep(r, length.out = nsq)
		  return(e2)
	      }
	      callGeneric(e1, diag2tT.u(e2,e1, "l"))
	  })



## Not {"sparseMatrix", "numeric} :  {"denseMatrix", "matrix", ... }
## ddi*:
setMethod("Ops", signature(e1 = "ddiMatrix", e2 = "ANY"),
	  function(e1,e2) callGeneric(diag2tT.u(e1,e2, "d"), e2))
setMethod("Ops", signature(e1 = "ANY", e2 = "ddiMatrix"),
	  function(e1,e2) callGeneric(e1, diag2tT.u(e2,e1, "d")))
## ldi*:
setMethod("Ops", signature(e1 = "ldiMatrix", e2 = "ANY"),
	  function(e1,e2) callGeneric(diag2tT.u(e1,e2, "l"), e2))
setMethod("Ops", signature(e1 = "ANY", e2 = "ldiMatrix"),
	  function(e1,e2) callGeneric(e1, diag2tT.u(e2,e1, "l")))



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
