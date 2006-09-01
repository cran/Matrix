## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a dgeMatrix.

## setMethod("%*%", signature(x = "dsparseMatrix", y = "ddenseMatrix"),
##           function(x, y) callGeneric(x, as(y, "dgeMatrix")))

## setMethod("%*%", signature(x = "ddenseMatrix", y = "dsparseMatrix"),
##           function(x, y) callGeneric(as(x, "dgeMatrix"), y))

setMethod("crossprod", signature(x = "dsparseMatrix", y = "ddenseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "dsparseMatrix"),
          function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y))

## and coerce dsparse* to dgC*
## setMethod("%*%", signature(x = "dsparseMatrix", y = "dgeMatrix"),
##           function(x, y) callGeneric(as(x, "dgCMatrix"), y))

## setMethod("%*%", signature(x = "dgeMatrix", y = "dsparseMatrix"),
##           function(x, y) callGeneric(x, as(y, "dgCMatrix")))

setMethod("crossprod", signature(x = "dsparseMatrix", y = "dgeMatrix"),
          function(x, y = NULL) callGeneric(as(x, "dgCMatrix"), y))

## NB: there's already
##     ("CsparseMatrix", "missing") and ("TsparseMatrix", "missing") methods

setMethod("crossprod", signature(x = "dgeMatrix", y = "dsparseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "dgCMatrix")))

setMethod("image", "dsparseMatrix",
          function(x, ...) image(as(x, "dgTMatrix"), ...))

setMethod("kronecker", signature(X = "dsparseMatrix", Y = "dsparseMatrix"),
          function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
          callGeneric(as(X, "dgTMatrix"), as(Y, "dgTMatrix")))

setMethod("lu", signature(x = "dsparseMatrix"),
	  function(x, ...) callGeneric(as(x, "dgCMatrix")))


## Group Methods, see ?Arith (e.g.)
## -----

##-> now moved to ./Csparse.R (and 'up' to ./sparseMatrix.R):
##  "Math2" is in ./dMatrix.R


### cbind2
    setMethod("cbind2", signature(x = "dsparseMatrix", y = "numeric"),
	      function(x, y) {
		  d <- dim(x); nr <- d[1]; nc <- d[2]; cl <- class(x)
                  x <- as(x, "dgCMatrix")
                  if(nr > 0) {
                      y <- rep(y, length.out = nr) # 'silent procrustes'
                      n0y <- y != 0
                      n.e <- length(x@i)
                      x@i <- c(x@i, (0:(nr-1))[n0y])
                      x@p <- c(x@p, n.e + sum(n0y))
                      x@x <- c(x@x, y[n0y])
                  } else { ## nr == 0

                  }
		  x@Dim[2] <- nc + 1:1
		  if(is.character(dn <- x@Dimnames[[2]]))
		      x@Dimnames[[2]] <- c(dn, "")
		  x
	      })
    ## the same, (x,y) <-> (y,x):
    setMethod("cbind2", signature(x = "numeric", y = "dsparseMatrix"),
	      function(x, y) {
		  d <- dim(y); nr <- d[1]; nc <- d[2]; cl <- class(y)
                  y <- as(y, "dgCMatrix")
                  if(nr > 0) {
                      x <- rep(x, length.out = nr) # 'silent procrustes'
                      n0x <- x != 0
                      y@i <- c((0:(nr-1))[n0x], y@i)
                      y@p <- c(0:0, sum(n0x) + y@p)
                      y@x <- c(x[n0x], y@x)
                  } else { ## nr == 0

                  }
		  y@Dim[2] <- nc + 1:1
		  if(is.character(dn <- y@Dimnames[[2]]))
		      y@Dimnames[[2]] <- c(dn, "")
		  y
	      })


    setMethod("cbind2", signature(x = "dsparseMatrix", y = "matrix"),
	      function(x, y) callGeneric(x, as(y, "dgCMatrix")))
    setMethod("cbind2", signature(x = "matrix", y = "dsparseMatrix"),
	      function(x, y) callGeneric(as(x, "dgCMatrix"), y))

    setMethod("cbind2", signature(x = "dsparseMatrix", y = "dsparseMatrix"),
	      function(x, y) {
		  nr <- rowCheck(x,y)
		  ## beware of (packed) triangular, symmetric, ...
		  hasDN <- !all(lapply(c(dnx <- dimnames(x),
                                         dny <- dimnames(y)), is.null))
                  ans <- .Call(Csparse_horzcat,
                               as(x, "dgCMatrix"), as(y, "dgCMatrix"))
		  if(hasDN) {
		      ## R and S+ are different in which names they take
		      ## if they differ -- but there's no warning in any case
		      rn <- if(!is.null(dnx[[1]])) dnx[[1]] else dny[[1]]
		      cx <- dnx[[2]] ; cy <- dny[[2]]
		      cn <- if(is.null(cx) && is.null(cy)) NULL
		      else c(if(!is.null(cx)) cx else rep.int("", ncol(x)),
			     if(!is.null(cy)) cy else rep.int("", ncol(y)))
		      ans@Dimnames <- list(rn, cn)
		  }
                  ans
	      })

    setMethod("rbind2", signature(x = "dsparseMatrix", y = "dsparseMatrix"),
	      function(x, y) {
		  nr <- colCheck(x,y)
		  ## beware of (packed) triangular, symmetric, ...
		  hasDN <- !all(lapply(c(dnx <- dimnames(x),
                                         dny <- dimnames(y)), is.null))
                  ans <- .Call(Csparse_vertcat,
                               as(x, "dgCMatrix"), as(y, "dgCMatrix"))
		  if(hasDN) {
		      ## R and S+ are different in which names they take
		      ## if they differ -- but there's no warning in any case
		      cn <- if(!is.null(dnx[[2]])) dnx[[2]] else dny[[2]]
		      rx <- dnx[[1]] ; ry <- dny[[1]]
		      rn <- if(is.null(rx) && is.null(ry)) NULL
		      else c(if(!is.null(rx)) rx else rep.int("", nrow(x)),
			     if(!is.null(ry)) ry else rep.int("", nrow(y)))
		      ans@Dimnames <- list(rn, cn)
		  }
		  ans
	      })

