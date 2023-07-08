#### Containing all  cbind2() and rbind2() methods for all our Matrices

###-- General -----------------------------------------------------------

###-- Dense, incl Diagonal ----------------------------------------------

###-- Sparse ------------------------------------------------------------

setMethod("cbind2", signature(x = "sparseMatrix", y = "matrix"),
	  function(x, y, ...) {
              cbind2(x, .Call(R_dense_as_sparse, y, ".gC", NULL, NULL))
          })
setMethod("cbind2", signature(x = "matrix", y = "sparseMatrix"),
	  function(x, y, ...) {
              cbind2(.Call(R_dense_as_sparse, x, ".gC", NULL, NULL), y)
          })
setMethod("rbind2", signature(x = "sparseMatrix", y = "matrix"),
	  function(x, y, ...) {
              rbind2(x, .Call(R_dense_as_sparse, y, ".gC", NULL, NULL))
          })
setMethod("rbind2", signature(x = "matrix", y = "sparseMatrix"),
	  function(x, y, ...) {
              rbind2(.Call(R_dense_as_sparse, x, ".gC", NULL, NULL), y)
          })

## originally from ./Matrix.R : -------------------------------

## The trivial methods :
setMethod("cbind2", signature(x = "Matrix", y = "NULL"),
          function(x, y, ...) x)
setMethod("cbind2", signature(x = "Matrix", y = "missing"),
          function(x, y, ...) x)
setMethod("cbind2", signature(x = "NULL", y="Matrix"),
          function(x, y, ...) y)
## using "atomicVector" not just "numeric"
setMethod("cbind2", signature(x = "Matrix", y = "atomicVector"),
	  function(x, y, ...) cbind2(x, matrix(y, nrow = nrow(x))))
setMethod("cbind2", signature(x = "atomicVector", y = "Matrix"),
	  function(x, y, ...) cbind2(matrix(x, nrow = nrow(y)), y))
setMethod("cbind2", signature(x = "ANY", y = "Matrix"),
	  function(x, y, ...) .bail.out.2(.Generic, class(x), class(y)))
setMethod("cbind2", signature(x = "Matrix", y = "ANY"),
	  function(x, y, ...) .bail.out.2(.Generic, class(x), class(y)))

setMethod("rbind2", signature(x = "Matrix", y = "NULL"),
          function(x, y, ...) x)
setMethod("rbind2", signature(x = "Matrix", y = "missing"),
          function(x, y, ...) x)
setMethod("rbind2", signature(x = "NULL", y="Matrix"),
          function(x, y, ...) y)
setMethod("rbind2", signature(x = "Matrix", y = "atomicVector"),
	  function(x, y, ...) rbind2(x, matrix(y, ncol = ncol(x))))
setMethod("rbind2", signature(x = "atomicVector", y = "Matrix"),
	  function(x, y, ...) rbind2(matrix(x, ncol = ncol(y)), y))
setMethod("rbind2", signature(x = "ANY", y = "Matrix"),
	  function(x, y, ...) .bail.out.2(.Generic, class(x), class(y)))
setMethod("rbind2", signature(x = "Matrix", y = "ANY"),
	  function(x, y, ...) .bail.out.2(.Generic, class(x), class(y)))

## Makes sure one gets x decent error message for the unimplemented cases:
setMethod("cbind2", signature(x = "Matrix", y = "Matrix"),
	  function(x, y, ...) {
	      rowCheck(x,y)
	      .bail.out.2("cbind2", class(x), class(y))
          })

## Use a working fall back {particularly useful for sparse}:
## FIXME: implement rbind2 via "cholmod" for C* and Tsparse ones
setMethod("rbind2", signature(x = "Matrix", y = "Matrix"),
          function(x, y, ...) {
              colCheck(x,y)
              t(cbind2(t(x), t(y)))
          })

## originally from ./denseMatrix.R : -------------------------------

### cbind2
setMethod("cbind2", signature(x = "denseMatrix", y = "numeric"),
	  function(x, y, ...) {
	      d <- dim(x); nr <- d[1]; nc <- d[2]
	      y <- rep_len(y, nr) # 'silent procrustes'
	      ## beware of (packed) triangular, symmetric, ...
	      x <- .dense2g(x)
	      x@x <- c(x@x, as.double(y))
	      x@Dim[2] <- nc + 1L
	      if(is.character(dn <- x@Dimnames[[2]]))
		  x@Dimnames[[2]] <- c(dn, "")
	      x
	  })
## the same, (x,y) <-> (y,x):
setMethod("cbind2", signature(x = "numeric", y = "denseMatrix"),
	  function(x, y, ...) {
	      d <- dim(y); nr <- d[1]; nc <- d[2]
	      x <- rep_len(x, nr)
	      y <- .dense2g(y)
	      y@x <- c(as.double(x), y@x)
	      y@Dim[2] <- nc + 1L
	      if(is.character(dn <- y@Dimnames[[2]]))
		  y@Dimnames[[2]] <- c("", dn)
	      y
	  })


setMethod("cbind2", signature(x = "denseMatrix", y = "matrix"),
	  function(x, y, ...) cbind2(x, .m2ge(y, ".")))
setMethod("cbind2", signature(x = "matrix", y = "denseMatrix"),
	  function(x, y, ...) cbind2(.m2ge(x, "."), y))

cbind2DN <- function(dnx,dny, ncx,ncy) {
    ## R and S+ are different in which names they take
    ## if they differ -- but there's no warning in any case
    rn <- if(!is.null(dnx[[1]])) dnx[[1]] else dny[[1]]
    cx <- dnx[[2]] ; cy <- dny[[2]]
    cn <- if(is.null(cx) && is.null(cy)) NULL
	  else c(if(!is.null(cx)) cx else character(ncx),
		 if(!is.null(cy)) cy else character(ncy))
    list(rn, cn)
}

setMethod("cbind2", signature(x = "denseMatrix", y = "denseMatrix"),
	  function(x, y, ...) {
	      rowCheck(x,y)
	      ncx <- x@Dim[2]
	      ncy <- y@Dim[2]
	      ## beware of (packed) triangular, symmetric, ...
	      hasDN <- !is.null.DN(dnx <- dimnames(x)) | !is.null.DN(dny <- dimnames(y))
	      x <- .dense2g(x)
	      y <- .dense2g(y)
	      xx <- c(x@x, y@x)
	      ## be careful, e.g., if we have an 'n' and 'd'
	      if(identical((tr <- typeof(xx)), typeof(x@x))) {
		  x@x <- xx
                  x@Dim[2] <- ncx + ncy
		  if(hasDN) x@Dimnames <- cbind2DN(dnx,dny, ncx,ncy)
		  x
	      } else if(identical(tr, typeof(y@x))) {
		  y@x <- xx
		  y@Dim[2] <- ncx + ncy
		  if(hasDN) y@Dimnames <- cbind2DN(dnx,dny, ncx,ncy)
		  y
	      } else stop("resulting x-slot has different type than x's or y's")
	  })

### rbind2 -- analogous to cbind2 --- more to do for @x though:

setMethod("rbind2", signature(x = "denseMatrix", y = "numeric"),
	  function(x, y, ...) {
	      if(is.character(dn <- x@Dimnames[[1]])) dn <- c(dn, "")
	      y <- rbind2(as(x,"matrix"), y)
	      new(paste0(.M.kind(y), "geMatrix"), x = c(y),
                  Dim = x@Dim + 1:0, Dimnames = list(dn, x@Dimnames[[2]]))
	  })
## the same, (x,y) <-> (y,x):
setMethod("rbind2", signature(x = "numeric", y = "denseMatrix"),
	  function(x, y, ...) {
	      if(is.character(dn <- y@Dimnames[[1]])) dn <- c("", dn)
	      x <- rbind2(x, as(y,"matrix"))
	      new(paste0(.M.kind(x), "geMatrix"), x = c(x),
                  Dim = y@Dim + 1:0, Dimnames = list(dn, y@Dimnames[[2]]))
	  })

setMethod("rbind2", signature(x = "denseMatrix", y = "matrix"),
	  function(x, y, ...) rbind2(x, .m2ge(y, ".")))
setMethod("rbind2", signature(x = "matrix", y = "denseMatrix"),
	  function(x, y, ...) rbind2(.m2ge(x, "."), y))

rbind2DN <- function(dnx, dny, nrx,nry) {
    if(!is.null.DN(dnx) || !is.null.DN(dny)) {
	## R and S+ are different in which names they take
	## if they differ -- but there's no warning in any case
	list(if(is.null(rx <- dnx[[1]]) & is.null(ry <- dny[[1]]))
	     NULL else
	     c(if(!is.null(rx)) rx else character(nrx),
	       if(!is.null(ry)) ry else character(nry)),
	     if(!is.null(dnx[[2]])) dnx[[2]] else dny[[2]])
    } else list(NULL, NULL)
}

setMethod("rbind2", signature(x = "denseMatrix", y = "denseMatrix"),
	  function(x, y, ...) {
	      colCheck(x,y)
	      nrx <- x@Dim[1]
	      nry <- y@Dim[1]
	      ## beware of (packed) triangular, symmetric, ...
	      hasDN <- !is.null.DN(dnx <- dimnames(x)) | !is.null.DN(dny <- dimnames(y))
	      x <- .dense2g(x)
	      y <- .dense2g(y)
	      xx <- .Call(R_rbind2_vector, x, y)
	      ## be careful, e.g., if we have an 'n' and 'd'
	      if(identical((tr <- typeof(xx)), typeof(x@x))) {
		  x@x <- xx
		  x@Dim[1] <- nrx + nry
		  if(hasDN) x@Dimnames <- rbind2DN(dnx,dny, nrx,nry)
		  x
	      } else if(identical(tr, typeof(y@x))) {
		  y@x <- xx
		  y@Dim[1] <- nrx + nry
		  if(hasDN) y@Dimnames <- rbind2DN(dnx,dny, nrx,nry)
		  y
	      } else stop("resulting x-slot has different type than x's or y's")
	  })

## originally from ./diagMatrix.R : --------------------------------------

## For diagonalMatrix:  preserve sparseness {not always optimal, but "the law"}

setMethod("cbind2", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y, ...)
              cbind2(.diag2sparse(x, ".gC"), as(y, "CsparseMatrix")))
setMethod("cbind2", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y, ...)
              cbind2(as(x, "CsparseMatrix"), .diag2sparse(y, ".gC")))
setMethod("rbind2", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y, ...)
              rbind2(.diag2sparse(x, ".gC"), as(y, "CsparseMatrix")))
setMethod("rbind2", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y, ...)
              rbind2(as(x, "CsparseMatrix"), .diag2sparse(y, ".gC")))

## in order to evade method dispatch ambiguity, but still remain "general"
## we use this hack instead of signature  x = "diagonalMatrix"
for(cls in names(getClassDef("diagonalMatrix")@subclasses)) {

 setMethod("cbind2", signature(x = cls, y = "matrix"),
	   function(x, y, ...)
               cbind2(.diag2sparse(x, ".gC"), .m2sparse(y, ".gC")))
 setMethod("cbind2", signature(x = "matrix", y = cls),
	   function(x, y, ...)
               cbind2(.m2sparse(x, ".gC"), .diag2sparse(y, ".gC")))
 setMethod("rbind2", signature(x = cls, y = "matrix"),
	   function(x, y, ...)
               rbind2(.diag2sparse(x, ".gC"), .m2sparse(y, ".gC")))
 setMethod("rbind2", signature(x = "matrix", y = cls),
	   function(x, y, ...)
               rbind2(.m2sparse(x, ".gC"), .diag2sparse(y, ".gC")))

 ## These are already defined for "Matrix"
 ## -- repeated here for method dispatch disambiguation	 {"design-FIXME" ?}
 setMethod("cbind2", signature(x = cls, y = "atomicVector"),
	   function(x, y, ...) cbind2(x, matrix(y, nrow = nrow(x))))
 setMethod("cbind2", signature(x = "atomicVector", y = cls),
	   function(x, y, ...) cbind2(matrix(x, nrow = nrow(y)), y))
 setMethod("rbind2", signature(x = cls, y = "atomicVector"),
	   function(x, y, ...) rbind2(x, matrix(y, ncol = ncol(x))))
 setMethod("rbind2", signature(x = "atomicVector", y = cls),
	   function(x, y, ...) rbind2(matrix(x, ncol = ncol(y)), y))
}
rm(cls)

## originally from ./dsparseMatrix.R : --------------------------------

## FIXME: dimnames() handling should happen in C code
## ------> ../src/Csparse.c

## Fast - almost non-checking methods
.cbind2Csp <- function(x,y) .Call(Csparse_horzcat, asCspN(x), asCspN(y))
.rbind2Csp <- function(x,y) .Call(Csparse_vertcat, asCspN(x), asCspN(y))

cbind2sparse <- function(x,y) {
    ## beware of (packed) triangular, symmetric, ...
    if(identical(c(dnx <- dimnames(x),
		   dny <- dimnames(y)),
		 list(NULL,NULL,NULL,NULL)))
	## keep empty dimnames
	.cbind2Csp(x,y)
    else {
	## R and S+ are different in which names they take
	## if they differ -- but there's no warning in any case
	rn <- if(!is.null(dnx[[1]])) dnx[[1]] else dny[[1]]
	cx <- dnx[[2]] ; cy <- dny[[2]]
	cn <- if(is.null(cx) && is.null(cy)) NULL
	else c(if(!is.null(cx)) cx else character(ncol(x)),
	       if(!is.null(cy)) cy else character(ncol(y)))
	ans <- .cbind2Csp(x,y)
	ans@Dimnames <- list(rn, cn)
	ans
    }
}
setMethod("cbind2", signature(x = "sparseMatrix", y = "sparseMatrix"),
	  function(x, y, ...) {
	      rowCheck(x,y)
	      cbind2sparse(x,y)
	  })

rbind2sparse <- function(x,y) {
    ## beware of (packed) triangular, symmetric, ...
    if(identical(c(dnx <- dimnames(x),
		   dny <- dimnames(y)),
		 list(NULL,NULL,NULL,NULL)))
	## keep empty dimnames
	.rbind2Csp(x,y)
    else {
	## R and S+ are different in which names they take
	## if they differ -- but there's no warning in any case
	cn <- if(!is.null(dnx[[2]])) dnx[[2]] else dny[[2]]
	rx <- dnx[[1]] ; ry <- dny[[1]]
	rn <- if(is.null(rx) && is.null(ry)) NULL
	else c(if(!is.null(rx)) rx else character(nrow(x)),
	       if(!is.null(ry)) ry else character(nrow(y)))
	ans <- .rbind2Csp(x,y)
	ans@Dimnames <- list(rn, cn)
	ans
    }
}
setMethod("rbind2", signature(x = "sparseMatrix", y = "sparseMatrix"),
	  function(x, y, ...) {
	      colCheck(x,y)
	      rbind2sparse(x,y)
	  })

if(length(formals(cbind2)) >= 3) { ## newer R -- can use optional 'sparse = NA'

setMethod("cbind2", signature(x = "sparseMatrix", y = "denseMatrix"),
	  function(x, y, sparse = NA, ...) {
	      nr <- rowCheck(x,y)
	      if(is.na(sparse)) # result is sparse if "enough zeros" <==> sparseDefault() in Matrix()
		  sparse <- (nnzero(x,na.counted=TRUE)+nnzero(y,na.counted=TRUE)) * 2 <
		      as.double(nr) * (ncol(x)+ncol(y)) # as.double(): avoid integer overflow in '*'
	      if(sparse) cbind2sparse(x,y) else cbind2(as(x, "denseMatrix"), y)
	  })
setMethod("cbind2", signature(x = "denseMatrix", y = "sparseMatrix"),
	  function(x, y, sparse = NA, ...) {
	      nr <- rowCheck(x,y)
	      if(is.na(sparse)) # result is sparse if "enough zeros" <==> sparseDefault() in Matrix()
		  sparse <- (nnzero(x,na.counted=TRUE)+nnzero(y,na.counted=TRUE)) * 2 <
		      as.double(nr) * (ncol(x)+ncol(y))
	      if(sparse) cbind2sparse(x,y) else cbind2(x, as(y, "denseMatrix"))
	  })
setMethod("rbind2", signature(x = "sparseMatrix", y = "denseMatrix"),
	  function(x, y, sparse = NA, ...) {
	      nc <- colCheck(x,y)
	      if(is.na(sparse)) # result is sparse if "enough zeros" <==> sparseDefault() in Matrix()
		  sparse <- (nnzero(x,na.counted=TRUE)+nnzero(y,na.counted=TRUE)) * 2 <
		      (nrow(x)+nrow(y)) * as.double(nc)
	      if(sparse) rbind2sparse(x,y) else rbind2(as(x, "denseMatrix"), y)
	  })
setMethod("rbind2", signature(x = "denseMatrix", y = "sparseMatrix"),
	  function(x, y, sparse = NA, ...) {
	      nc <- colCheck(x,y)
	      if(is.na(sparse)) # result is sparse if "enough zeros" <==> sparseDefault() in Matrix()
		  sparse <- (nnzero(x,na.counted=TRUE)+nnzero(y,na.counted=TRUE)) * 2 <
		(nrow(x)+nrow(y)) * as.double(nc)
	      if(sparse) rbind2sparse(x,y) else rbind2(x, as(y, "denseMatrix"))
	  })

} else { ## older version of R -- cbind2() has no "..."

setMethod("cbind2", signature(x = "sparseMatrix", y = "denseMatrix"),
	  function(x, y, ...) {
	      nr <- rowCheck(x,y)
	      ## result is sparse if "enough zeros" <==> sparseDefault() in Matrix()
	      sparse <- (nnzero(x,na.counted=TRUE)+nnzero(y,na.counted=TRUE)) * 2 <
		  as.double(nr) * (ncol(x)+ncol(y))
	      if(sparse) cbind2sparse(x,y) else cbind2(as(x, "denseMatrix"), y)
	  })
setMethod("cbind2", signature(x = "denseMatrix", y = "sparseMatrix"),
	  function(x, y, ...) {
	      nr <- rowCheck(x,y)
	      ## result is sparse if "enough zeros" <==> sparseDefault() in Matrix()
	      sparse <- (nnzero(x,na.counted=TRUE)+nnzero(y,na.counted=TRUE)) * 2 <
		  as.double(nr) * (ncol(x)+ncol(y))
	      if(sparse) cbind2sparse(x,y) else cbind2(x, as(y, "denseMatrix"))
	  })
setMethod("rbind2", signature(x = "sparseMatrix", y = "denseMatrix"),
	  function(x, y, ...) {
	      nc <- colCheck(x,y)
	      ## result is sparse if "enough zeros" <==> sparseDefault() in Matrix()
	      sparse <- (nnzero(x,na.counted=TRUE)+nnzero(y,na.counted=TRUE)) * 2 <
		  (nrow(x)+nrow(y)) * as.double(nc)
	      if(sparse) rbind2sparse(x,y) else rbind2(as(x, "denseMatrix"), y)
	  })
setMethod("rbind2", signature(x = "denseMatrix", y = "sparseMatrix"),
	  function(x, y, ...) {
	      nc <- colCheck(x,y)
	      ## result is sparse if "enough zeros" <==> sparseDefault() in Matrix()
	      sparse <- (nnzero(x,na.counted=TRUE)+nnzero(y,na.counted=TRUE)) * 2 <
		  (nrow(x)+nrow(y)) * as.double(nc)
	      if(sparse) rbind2sparse(x,y) else rbind2(x, as(y, "denseMatrix"))
	  })
}# older R -- no "sparse = NA"




if(FALSE) {
    ## FIXME
    ##------------- maybe a bit faster --- but too much to maintain
    ## would have to be done for "rbind2" as well ...
setMethod("cbind2", signature(x = "sparseMatrix", y = "numeric"),
          function(x, y, ...) {
              d <- dim(x); nr <- d[1]; nc <- d[2]; cl <- class(x)
              x <- as(x, "CsparseMatrix")
              if(nr > 0) {
		  y <- rep_len(y, nr) # 'silent procrustes'
                  n0y <- y != 0
                  n.e <- length(x@i)
                  x@i <- c(x@i, (0:(nr-1))[n0y])
                  x@p <- c(x@p, n.e + sum(n0y))
                  x@x <- c(x@x, y[n0y])
              } else { ## nr == 0

              }
              x@Dim[2] <- nc + 1L
              if(is.character(dn <- x@Dimnames[[2]]))
                  x@Dimnames[[2]] <- c(dn, "")
              x
          })
## the same, (x,y) <-> (y,x):
setMethod("cbind2", signature(x = "numeric", y = "sparseMatrix"),
          function(x, y, ...) {
              d <- dim(y); nr <- d[1]; nc <- d[2]; cl <- class(y)
              y <-  as(y, "CsparseMatrix")
              if(nr > 0) {
		  x <- rep_len(x, nr) # 'silent procrustes'
                  n0x <- x != 0
                  y@i <- c((0:(nr-1))[n0x], y@i)
                  y@p <- c(0L, sum(n0x) + y@p)
                  y@x <- c(x[n0x], y@x)
              } else { ## nr == 0

              }
              y@Dim[2] <- nc + 1L
              if(is.character(dn <- y@Dimnames[[2]]))
                  y@Dimnames[[2]] <- c(dn, "")
              y
          })

}## -- no longer

setMethod("rbind2", signature(x = "indMatrix", y = "indMatrix"),
	  function(x, y, ...) {
	      dx <- x@Dim
	      dy <- y@Dim
              if((n <- dx[2L]) != dy[2L])
		  stop(gettextf("number of columns of matrices must match in %s",
				deparse(sys.call(sys.parent()))),
		       call. = FALSE, domain = NA)
              if(x@margin != 1L || y@margin != 1L)
                  return(rbind2(as(x, "RsparseMatrix"), as(y, "RsparseMatrix")))
              mx <- dx[1L]
              my <- dy[1L]
              r <- new("indMatrix")
              r@Dim <- c(mx + my, n)
              r@Dimnames <- rbind2DN(x@Dimnames, y@Dimnames, mx, my)
              r@perm <- c(x@perm, y@perm)
              r
	  })

setMethod("cbind2", signature(x = "indMatrix", y = "indMatrix"),
	  function(x, y, ...) {
	      dx <- x@Dim
	      dy <- y@Dim
              if((m <- dx[1L]) != dy[1L])
		  stop(gettextf("number of rows of matrices must match in %s",
				deparse(sys.call(sys.parent()))),
		       call. = FALSE, domain = NA)
              if(x@margin == 1L || y@margin == 1L)
                  return(cbind2(as(x, "CsparseMatrix"), as(y, "CsparseMatrix")))
              nx <- dx[2L]
              ny <- dy[2L]
              r <- new("indMatrix")
              r@Dim <- c(m, nx + ny)
              r@Dimnames <- cbind2DN(x@Dimnames, y@Dimnames, nx, ny)
              r@perm <- c(x@perm, y@perm)
              r@margin <- 2L
              r
	  })
