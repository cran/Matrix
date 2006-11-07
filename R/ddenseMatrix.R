### Define Methods that can be inherited for all subclasses

## This replaces many "d..Matrix" -> "dgeMatrix" ones
## >> but << needs all sub(sub(sub)) classes of "ddenseMatrix" listed
##   -----  in  ../src/Mutils.c

## Should this method return 'from' without duplication when it has
## class dgeMatrix?
setAs("ddenseMatrix", "dgeMatrix",
      function(from) {
          if (class(from) != "dgeMatrix")
              from <- .Call(dup_mMatrix_as_dgeMatrix, from)
          from
      })

## d(ouble) to l(ogical):
setAs("dgeMatrix", "lgeMatrix", d2l_Matrix)
setAs("dtrMatrix", "ltrMatrix", d2l_Matrix)
setAs("dtpMatrix", "ltpMatrix", d2l_Matrix)
setAs("dsyMatrix", "lsyMatrix", d2l_Matrix)
setAs("dspMatrix", "lspMatrix", d2l_Matrix)

setAs("ddenseMatrix", "CsparseMatrix",
      function(from) {
          if (class(from) != "dgeMatrix")
              from <- .Call(dup_mMatrix_as_dgeMatrix, from)
          .Call(dense_to_Csparse, from)
      })

## special case
setAs("dgeMatrix", "dgCMatrix",
      function(from) .Call(dense_to_Csparse, from))

setAs("matrix", "CsparseMatrix",
      function(from) {
	    if(is.numeric(from))
		.Call(dense_to_Csparse, .Call(dup_mMatrix_as_dgeMatrix, from))
	    else if(is.logical(from)) ## FIXME: this works, but maybe wastefully
                as(Matrix(from, sparse=TRUE), "CsparseMatrix")
	    else stop('not-yet-implemented coercion to "CsparseMatrix"')
      })


## special case needed in the Matrix function
setAs("matrix", "dgCMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .Call(dense_to_Csparse, from)
      })

setAs("numeric", "CsparseMatrix",
      function(from)
      .Call(dense_to_Csparse, .Call(dup_mMatrix_as_dgeMatrix, from)))

setMethod("as.numeric", signature(x = "ddenseMatrix"),
	  function(x, ...) as(x, "dgeMatrix")@x)

## -- see also ./Matrix.R  e.g., for a show() method

## These methods are the 'fallback' methods for all dense numeric
## matrices in that they simply coerce the ddenseMatrix to a
## dgeMatrix. Methods for special forms override these.

setMethod("norm", signature(x = "ddenseMatrix", type = "missing"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("norm", signature(x = "ddenseMatrix", type = "character"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix"), type))

setMethod("rcond", signature(x = "ddenseMatrix", type = "missing"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("rcond", signature(x = "ddenseMatrix", type = "character"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix"), type))

## Not really useful; now require *identical* class for result:
## setMethod("t", signature(x = "ddenseMatrix"),
## 	  function(x) callGeneric(as(x, "dgeMatrix")))

setMethod("tcrossprod", signature(x = "ddenseMatrix", y = "missing"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "missing"),
          function(x, y = NULL) callGeneric(as(x, "dgeMatrix")))

setMethod("diag", signature(x = "ddenseMatrix"),
          function(x = 1, nrow, ncol = n) callGeneric(as(x, "dgeMatrix")))

## These methods cause an infinite loop in pre-2.4.0
## setMethod("solve", signature(a = "ddenseMatrix", b = "missing"),
##           function(a, b, ...) callGeneric(as(a, "dgeMatrix")))

## setMethod("solve", signature(a = "ddenseMatrix", b = "ANY"),
##           function(a, b, ...) callGeneric(as(a, "dgeMatrix"), b))

## General method for dense matrix multiplication in case specific methods
## have not been defined.
setMethod("%*%", signature(x = "ddenseMatrix", y = "ddenseMatrix"),
          function(x, y) .Call(dgeMatrix_matrix_mm,
                               .Call(dup_mMatrix_as_dgeMatrix, x), y, FALSE),
          valueClass = "dgeMatrix")

setMethod("lu", signature(x = "ddenseMatrix"),
          function(x, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("chol", signature(x = "ddenseMatrix", pivot = "ANY"), cholMat)

setMethod("determinant", signature(x = "ddenseMatrix", logarithm = "missing"),
          function(x, logarithm, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("determinant", signature(x = "ddenseMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
          callGeneric(as(x, "dgeMatrix"), logarithm))

## now done for "dMatrix":
## setMethod("expm", signature(x = "ddenseMatrix"),
##           function(x) callGeneric(as(x, "dgeMatrix")))

setMethod("Schur", signature(x = "ddenseMatrix", vectors = "missing"),
          function(x, vectors, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("Schur", signature(x = "ddenseMatrix", vectors = "logical"),
          function(x, vectors, ...) callGeneric(as(x, "dgeMatrix"), vectors))


## Cheap version: work via "dgeMatrix" and use the group methods there:
## FIXME(?): try to preserve "symmetric", "triangular", ...
setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
          signature(e1 = "ddenseMatrix", e2 = "ddenseMatrix"),
          function(e1, e2) callGeneric(as(e1, "dgeMatrix"),
                                       as(e2, "dgeMatrix")))
setMethod("Arith",
          signature(e1 = "ddenseMatrix", e2 = "numeric"),
          function(e1, e2) callGeneric(as(e1, "dgeMatrix"), e2))
setMethod("Arith",
          signature(e1 = "numeric", e2 = "ddenseMatrix"),
          function(e1, e2) callGeneric(e1, as(e2, "dgeMatrix")))

setMethod("Math",
          signature(x = "ddenseMatrix"),
          function(x) callGeneric(as(x, "dgeMatrix")))



### for R 2.2.x (and later): -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

### cbind2
setMethod("cbind2", signature(x = "ddenseMatrix", y = "numeric"),
	  function(x, y) {
	      d <- dim(x); nr <- d[1]; nc <- d[2]
	      y <- rep(y, length.out = nr) # 'silent procrustes'
	      ## beware of (packed) triangular, symmetric, ...
	      x <- as(x, "dgeMatrix")
	      x@x <- c(x@x, as.double(y))
	      x@Dim[2] <- nc + 1:1
	      if(is.character(dn <- x@Dimnames[[2]]))
		  x@Dimnames[[2]] <- c(dn, "")
	      x
	  })
## the same, (x,y) <-> (y,x):
setMethod("cbind2", signature(x = "numeric", y = "ddenseMatrix"),
	  function(x, y) {
	      d <- dim(y); nr <- d[1]; nc <- d[2]
	      x <- rep(x, length.out = nr)
	      y <- as(y, "dgeMatrix")
	      y@x <- c(as.double(x), y@x)
	      y@Dim[2] <- nc + 1:1
	      if(is.character(dn <- y@Dimnames[[2]]))
		  y@Dimnames[[2]] <- c("", dn)
	      y
	  })

setMethod("cbind2", signature(x = "ddenseMatrix", y = "matrix"),
	  function(x, y) callGeneric(x, as(y, "dgeMatrix")))
setMethod("cbind2", signature(x = "matrix", y = "ddenseMatrix"),
	  function(x, y) callGeneric(as(x, "dgeMatrix"), y))

setMethod("cbind2", signature(x = "ddenseMatrix", y = "ddenseMatrix"),
	  function(x, y) {
	      nr <- rowCheck(x,y)
	      ncx <- x@Dim[2]
	      ncy <- y@Dim[2]
	      ## beware of (packed) triangular, symmetric, ...
	      hasDN <- !is.null(dnx <- dimnames(x)) |
	      !is.null(dny <- dimnames(y))
	      x <- as(x, "dgeMatrix")
	      y <- as(y, "dgeMatrix")
	      x@x <- c(x@x, y@x)
	      x@Dim[2] <- ncx + ncy
	      if(hasDN) {
		  ## R and S+ are different in which names they take
		  ## if they differ -- but there's no warning in any case
		  rn <- if(!is.null(dnx[[1]])) dnx[[1]] else dny[[1]]
		  cx <- dnx[[2]] ; cy <- dny[[2]]
		  cn <- if(is.null(cx) && is.null(cy)) NULL
		  else c(if(!is.null(cx)) cx else rep.int("", ncx),
			 if(!is.null(cy)) cy else rep.int("", ncy))
		  x@Dimnames <- list(rn, cn)
	      }
	      x
	  })

### rbind2 -- analogous to cbind2 --- more to do for @x though:

setMethod("rbind2", signature(x = "ddenseMatrix", y = "numeric"),
	  function(x, y) {
	      if(is.character(dn <- x@Dimnames[[1]])) dn <- c(dn, "")
	      new("dgeMatrix", Dim = x@Dim + 1:0,
		  Dimnames = list(dn, x@Dimnames[[2]]),
		  x = c(rbind2(as(x,"matrix"), y)))
	  })
## the same, (x,y) <-> (y,x):
setMethod("rbind2", signature(x = "numeric", y = "ddenseMatrix"),
	  function(x, y) {
	      if(is.character(dn <- y@Dimnames[[1]])) dn <- c("", dn)
	      new("dgeMatrix", Dim = y@Dim + 1:0,
		  Dimnames = list(dn, y@Dimnames[[2]]),
		  x = c(rbind2(x, as(y,"matrix"))))
	  })

setMethod("rbind2", signature(x = "ddenseMatrix", y = "matrix"),
	  function(x, y) callGeneric(x, as(y, "dgeMatrix")))
setMethod("rbind2", signature(x = "matrix", y = "ddenseMatrix"),
	  function(x, y) callGeneric(as(x, "dgeMatrix"), y))

setMethod("rbind2", signature(x = "ddenseMatrix", y = "ddenseMatrix"),
	  function(x, y) {
	      nc <- colCheck(x,y)
	      nrx <- x@Dim[1]
	      nry <- y@Dim[1]
	      dn <-
		  if(!is.null(dnx <- dimnames(x)) |
		     !is.null(dny <- dimnames(y))) {
		      ## R and S+ are different in which names they take
		      ## if they differ -- but there's no warning in any case
		      list(if(is.null(rx <- dnx[[1]]) && is.null(ry <- dny[[1]]))
			   NULL else
			   c(if(!is.null(rx)) rx else rep.int("", nrx),
			     if(!is.null(ry)) ry else rep.int("", nry)),
			   if(!is.null(dnx[[2]])) dnx[[2]] else dny[[2]])

		  } else list(NULL, NULL)
	      ## beware of (packed) triangular, symmetric, ...
	      new("dgeMatrix", Dim = c(nrx + nry, nc), Dimnames = dn,
		  x = c(rbind2(as(x,"matrix"), as(y,"matrix"))))
	  })

### FIXME: band() et al should be extended from "ddense" to "dense" !
###        However, needs much work to generalize dup_mMatrix_as_dgeMatrix()

## NB: have extra tril(), triu() methods for symmetric ["dsy" and "dsp"] and
##     for triangular ["dtr" and "dtp"]
setMethod("tril", "ddenseMatrix",
	  function(x, k = 0, ...) {
	      k <- as.integer(k[1])
	      dd <- dim(x); sqr <- dd[1] == dd[2]
	      stopifnot(-dd[1] <= k, k <= dd[1]) # had k <= 0
	      ## returns "lower triangular" if k <= 0 && sqr
	      .Call(ddense_band, x, -dd[1], k)
	  })

setMethod("triu", "ddenseMatrix",
	  function(x, k = 0, ...) {
	      k <- as.integer(k[1])
	      dd <- dim(x); sqr <- dd[1] == dd[2]
	      stopifnot(-dd[1] <= k, k <= dd[1]) # had k >= 0
	      ## returns "upper triangular" if k >= 0
	      .Call(ddense_band, x, k, dd[2])
	  })

setMethod("band", "ddenseMatrix",
	  function(x, k1, k2, ...) {
	      k1 <- as.integer(k1[1])
	      k2 <- as.integer(k2[1])
	      dd <- dim(x); sqr <- dd[1] == dd[2]
	      stopifnot(-dd[1] <= k1, k1 <= k2, k2 <= dd[1])
	      r <- .Call(ddense_band, x, k1, k2)
	      if (k1 < 0  &&  k1 == -k2  && isSymmetric(x)) ## symmetric
		  as(r, paste(.M.kind(x), "syMatrix", sep=''))
	      else
		  r
	  })
