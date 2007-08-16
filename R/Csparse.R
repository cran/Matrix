#### Methods for the virtual class 'CsparseMatrix' of sparse matrices stored in
####  "column compressed" format.
#### -- many more specific things are e.g. in ./dgCMatrix.R

setAs("CsparseMatrix", "TsparseMatrix",
      function(from)
          ## |-> cholmod_C -> cholmod_T -> chm_triplet_to_SEXP
          ## modified to support triangular (../src/Csparse.c)
          .Call(Csparse_to_Tsparse, from, is(from, "triangularMatrix")))


## special cases (when a specific "to" class is specified)
setAs("dgCMatrix", "dgTMatrix",
      function(from) .Call(Csparse_to_Tsparse, from, FALSE))

setAs("dsCMatrix", "dsTMatrix",
      function(from) .Call(Csparse_to_Tsparse, from, FALSE))

setAs("dsCMatrix", "dgCMatrix",
      function(from) .Call(Csparse_symmetric_to_general, from))

for(prefix in c("d", "l", "n"))
    setAs(paste(prefix,"sCMatrix",sep=''), "generalMatrix",
	  function(from) .Call(Csparse_symmetric_to_general, from))

setAs("dtCMatrix", "dtTMatrix",
      function(from) .Call(Csparse_to_Tsparse, from, TRUE))

setAs("CsparseMatrix", "denseMatrix",
      function(from) {
	  ## |-> cholmod_C -> cholmod_dense -> chm_dense_to_dense
	  cld <- getClassDef(class(from))
	  if (extends(cld, "generalMatrix"))
	      .Call(Csparse_to_dense, from)
	  else {
	      ## Csparse_to_dense  loses symmetry and triangularity properties.
	      ## With suitable changes to chm_dense_to_SEXP (../src/chm_common.c)
	      ## we could do this in C code -- or do differently in C {FIXME!}
	      if (extends(cld, "triangularMatrix") && from@diag == "U")
		  from <- .Call(Csparse_diagU2N, from)
	      as(.Call(Csparse_to_dense, from), # -> "[dln]geMatrix"
		 paste(.M.kind(from, cld),
		       .dense.prefixes[.M.shape(from, cld)], "Matrix", sep=''))
	  }
      })

## special cases (when a specific "to" class is specified)
setAs("dgCMatrix", "dgeMatrix",
      function(from) .Call(Csparse_to_dense, from))

## cholmod_sparse_to_dense converts symmetric storage to general
## storage so symmetric classes are ok for conversion to matrix.
## unit triangular needs special handling
setAs("CsparseMatrix", "matrix",
      function(from) {
          ## |-> cholmod_C -> cholmod_dense -> chm_dense_to_matrix
          if (is(from, "triangularMatrix") && from@diag == "U")
              from <- .Call(Csparse_diagU2N, from)
          .Call(Csparse_to_matrix, from)
      })

### Some group methods:

setMethod("Math",
	  signature(x = "CsparseMatrix"),
	  function(x) {
	      f0 <- callGeneric(0.)
	      if(is0(f0)) {
		  ## sparseness, symm., triang.,... preserved
                  cl <- class(x)
                  has.x <- !extends(cl, "nsparseMatrix")
                  ## has.x  <==> *not* nonzero-pattern == "nMatrix"
                  if(has.x) {
                      type <- storage.mode(x@x)
                      r <- callGeneric(x@x)
                  } else { ## nsparseMatrix
                      type <- ""
		      r <- rep.int(as.double(callGeneric(TRUE)),
				   switch(.sp.class(cl),
					  CsparseMatrix = length(x@i),
					  TsparseMatrix = length(x@i),
					  RsparseMatrix = length(x@j)))
                  }
		  if(type == storage.mode(r)) {
		      x@x <- r
		      x
		  } else { ## e.g. abs( <lgC> ) --> integer Csparse
		      ## FIXME: when we have 'i*' classes, use them here:
		      rx <- new(sub("^.", "d", cl))
		      rx@x <- as.double(r)
		      ## result is "same"
		      sNams <- slotNames(cl)
		      for(nm in sNams[sNams != "x"])
			  slot(rx, nm) <- slot(x, nm)
		      rx
		  }
	      } else { ## no sparseness
		  callGeneric(as_dense(x))
	      }
	  })



### workhorse for "[<-" -- both for d* and l*  C-sparse matrices :
replCmat <- function (x, i, j, value)
{
    di <- dim(x)
    dn <- dimnames(x)
    i1 <- if(missing(i)) 0:(di[1] - 1L) else .ind.prep2(i, 1, di, dn)
    i2 <- if(missing(j)) 0:(di[2] - 1L) else .ind.prep2(j, 2, di, dn)
    dind <- c(length(i1), length(i2)) # dimension of replacement region
    lenRepl <- prod(dind)
    lenV <- length(value)
    if(lenV == 0) {
        if(lenRepl != 0)
            stop("nothing to replace with")
        else return(x)
    }
    ## else: lenV := length(value)	 is > 0
    if(lenRepl %% lenV != 0)
	stop("number of items to replace is not a multiple of replacement length")
    if(lenV > lenRepl)
	stop("too many replacement values")

    clx <- class(x)
    clDx <- getClassDef(clx) # extends() , is() etc all use the class definition

    ## keep "symmetry" if changed here:
    x.sym <- extends(clDx, "symmetricMatrix")
    if(x.sym) { ## only half the indices are there..
	x.sym <-
	    (dind[1] == dind[2] && all(i1 == i2) &&
	     (lenRepl == 1 || isSymmetric(array(value, dim=dind))))
	## x.sym : result is *still* symmetric
	x <- .Call(Csparse_symmetric_to_general, x) ## but do *not* redefine clx!

    }
    else if((x.tri <- extends(clDx, "triangularMatrix"))) {
        xU <- x@uplo == "U"
	r.tri <- all(if(xU) i1 <= i2 else i2 <= i1)
	if(r.tri) { ## result is *still* triangular
            if(any(i1 == i2)) # diagonal will be changed
                x <- diagU2N(x) # keeps class (!)
	}
	else { # go to "generalMatrix" and continue
	    x <- as(x, paste(.M.kind(x), "gCMatrix", sep='')) ## & do not redefine clx!
	}
    }

    xj <- .Call(Matrix_expand_pointers, x@p)
    sel <- (!is.na(match(x@i, i1)) &
	    !is.na(match( xj, i2)))
    has.x <- "x" %in% slotNames(clDx)# === slotNames(x),
    ## has.x  <==> *not* nonzero-pattern == "nMatrix"

    if(has.x && sum(sel) == lenRepl) { ## all entries to be replaced are non-zero:
        ## need indices instead of just 'sel', for, e.g.,  A[2:1, 2:1] <- v
	non0 <- cbind(match(x@i[sel], i1),
		      match(xj [sel], i2)) - 1L
	iN0 <- 1L + encodeInd(non0, nr = dind[1])

        if(lenV < lenRepl)
            value <- rep(value, length = lenRepl)
	## Ideally we only replace them where value != 0 and drop the value==0
	## ones; but that would have to (?) go through dgT*
	## v0 <- 0 == value
	## if (lenRepl == 1) and v0 is TRUE, the following is not doing anything
	##-  --> ./dgTMatrix.R	and its	 replTmat()
	## x@x[sel[!v0]] <- value[!v0]
        x@x[sel] <- value[iN0]

	return(if(x.sym) as_CspClass(x, clx) else x)
    }
    ## else go via Tsparse.. {FIXME: a waste! - we already have 'xj' ..}
    ## and inside  Tsparse... the above i1, i2,..., sel  are *all* redone!
    x <- as(x, "TsparseMatrix")
    if(missing(i))
	x[ ,j] <- value
    else if(missing(j))
	x[i, ] <- value
    else
	x[i,j] <- value

    if(any(is0(x@x))) ## drop all values that "happen to be 0"
	drop0(x, clx)
    else as_CspClass(x, clx)
}

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "index", j = "missing",
                                value = "replValue"),
                 function (x, i, value) replCmat(x, i=i, value=value))

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "missing", j = "index",
                                value = "replValue"),
                 function (x, j, value) replCmat(x, j=j, value=value))

setReplaceMethod("[", signature(x = "CsparseMatrix", i = "index", j = "index",
				value = "replValue"),
                 replCmat)

## A[ ij ] <- value,  where ij is (i,j) 2-column matrix
setReplaceMethod("[", signature(x = "CsparseMatrix", i = "matrix", j = "missing",
				value = "replValue"),
		 function(x, i, value)
		 ## goto Tsparse modify and convert back:
		 as(.TM.repl.i.2col(as(x, "TsparseMatrix"), i=i, value=value),
		    "CsparseMatrix"))


setMethod("crossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      if (is(x, "symmetricMatrix")) {
		  warning("crossprod(x) calculated as x %*% x for sparse, symmetric x")
		  return(x %*% x)
	      }
	      .Call(Csparse_crossprod, x, trans = FALSE, triplet = FALSE)
	  })

setMethod("crossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
	  function(x, y = NULL)
	  .Call(Csparse_Csparse_crossprod, x, y, trans = FALSE))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
	  function(x, y = NULL)
	  .Call(Csparse_Csparse_crossprod, x, y, trans = TRUE))

## FIXME: Generalize the class of y.  This specific method is to replace one
##        in dgCMatrix.R
setMethod("crossprod", signature(x = "CsparseMatrix", y = "ddenseMatrix"),
	  function(x, y = NULL) .Call(Csparse_dense_crossprod, x, y))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "matrix"),
	  function(x, y = NULL) .Call(Csparse_dense_crossprod, x, y))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "numeric"),
	  function(x, y = NULL) .Call(Csparse_dense_crossprod, x, y))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
              if (is(x, "symmetricMatrix")) {
                  warning("tcrossprod(x) calculated as x %*% x for sparse, symmetric x")
                  return(x %*% x)
              }
              .Call(Csparse_crossprod, x, trans = TRUE, triplet = FALSE)
	  })

setMethod("t", signature(x = "CsparseMatrix"),
	  function(x) .Call(Csparse_transpose, x, is(x, "triangularMatrix")))

setMethod("%*%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y) .Call(Csparse_Csparse_prod, x, y))

setMethod("%*%", signature(x = "CsparseMatrix", y = "ddenseMatrix"),
          function(x, y) .Call(Csparse_dense_prod, x, y))

setMethod("%*%", signature(x = "CsparseMatrix", y = "matrix"),
          function(x, y) .Call(Csparse_dense_prod, x, y))

## Not needed because of c("Matrix", "numeric") method
##setMethod("%*%", signature(x = "CsparseMatrix", y = "numeric"),
##          function(x, y) .Call(Csparse_dense_prod, x, y))

## FIXME(2): These two are sub-optimal : has  2 x  t(<dense>)  :
setMethod("%*%", signature(x = "ddenseMatrix", y = "CsparseMatrix"),
          function(x, y) t(.Call(Csparse_dense_crossprod, y, t(x))),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "matrix", y = "CsparseMatrix"),
          function(x, y) t(.Call(Csparse_dense_crossprod, y, t(x))),
          valueClass = "dgeMatrix")

## Not needed because of c("numeric", "Matrix") method
##setMethod("%*%", signature(x = "numeric", y = "CsparseMatrix"),
##          function(x, y) t(.Call(Csparse_dense_crossprod, y, x)),
##          valueClass = "dgeMatrix")

## NB: have extra tril(), triu() methods for symmetric ["dsC" and "lsC"]
setMethod("tril", "CsparseMatrix",
	  function(x, k = 0, ...) {
	      k <- as.integer(k[1])
	      dd <- dim(x); sqr <- dd[1] == dd[2]
	      stopifnot(-dd[1] <= k, k <= dd[1]) # had k <= 0
	      r <- .Call(Csparse_band, x, -dd[1], k)
	      ## return "lower triangular" if k <= 0
	      if(sqr && k <= 0)
		  as(r, paste(.M.kind(x), "tCMatrix", sep='')) else r
	  })

setMethod("triu", "CsparseMatrix",
	  function(x, k = 0, ...) {
	      k <- as.integer(k[1])
	      dd <- dim(x); sqr <- dd[1] == dd[2]
	      stopifnot(-dd[1] <= k, k <= dd[1]) # had k >= 0
	      r <- .Call(Csparse_band, x, k, dd[2])
	      ## return "upper triangular" if k >= 0
	      if(sqr && k >= 0)
		  as(r, paste(.M.kind(x), "tCMatrix", sep='')) else r
	  })

setMethod("band", "CsparseMatrix",
	  function(x, k1, k2, ...) {
	      k1 <- as.integer(k1[1])
	      k2 <- as.integer(k2[1])
	      dd <- dim(x); sqr <- dd[1] == dd[2]
	      stopifnot(-dd[1] <= k1, k1 <= k2, k2 <= dd[1])
	      r <- .Call(Csparse_band, x, k1, k2)
	      if(sqr && k1 * k2 >= 0) ## triangular
		  as(r, paste(.M.kind(x), "tCMatrix", sep=''))
	      else if (k1 < 0  &&  k1 == -k2  && isSymmetric(x)) ## symmetric
		  as(r, paste(.M.kind(x), "sCMatrix", sep=''))
	      else
		  r
	  })

setMethod("diag", "CsparseMatrix",
	  function(x, nrow, ncol) {
	      dm <- .Call(Csparse_band, x, 0, 0)
	      dlen <- min(dm@Dim)
	      ind1 <- dm@i + 1L	# 1-based index vector
	      if (is(dm, "nMatrix")) {
		  val <- rep.int(FALSE, dlen)
		  val[ind1] <- TRUE
	      }
	      else if (is(dm, "lMatrix")) {
		  val <- rep.int(FALSE, dlen)
		  val[ind1] <- as.logical(dm@x)
	      }
	      else {
		  val <- rep.int(0, dlen)
		  ## cMatrix not yet active but for future expansion
		  if (is(dm, "cMatrix")) val <- as.complex(val)
		  val[ind1] <- dm@x
	      }
	      val
	  })

setMethod("writeMM", "CsparseMatrix",
	  function(obj, file, ...)
          .Call(Csparse_MatrixMarket, obj, as.character(file)))

