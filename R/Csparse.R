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

setAs("dtCMatrix", "dtTMatrix",
      function(from) .Call(Csparse_to_Tsparse, from, TRUE))

## Current code loses symmetry and triangularity properties.  With suitable
## changes to chm_dense_to_SEXP (../src/chm_common.c) we can avoid this.
setAs("CsparseMatrix", "denseMatrix",
      function(from) {
          ## |-> cholmod_C -> cholmod_dense -> chm_dense_to_dense
          if (is(from, "triangularMatrix") && from@diag == "U")
              from <- .Call(Csparse_diagU2N, from)
          .Call(Csparse_to_dense, from)
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

setMethod("Arith",
	  signature(e1 = "CsparseMatrix", e2 = "CsparseMatrix"),
	  function(e1, e2) callGeneric(as(e1, "dgCMatrix"),
				       as(e2, "dgCMatrix")))

setMethod("Arith",
	  signature(e1 = "CsparseMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      if(length(e2) == 1) { ## e.g.,  Mat ^ a
		  f0 <- callGeneric(0, e2)
		  if(is0(f0)) { # remain sparse, symm., tri.,...
		      e1@x <- callGeneric(e1@x, e2)
		      return(e1)
		  }
	      }
	      ## all other (potentially non-sparse) cases: give up symm, tri,..
	      callGeneric(as(e1, paste(.M.kind(e1), "gCMatrix", sep='')), e2)
	  })

## The same,  e1 <-> e2 :
setMethod("Arith",
	  signature(e1 = "numeric", e2 = "CsparseMatrix"),
	  function(e1, e2) {
	      if(length(e1) == 1) {
		  f0 <- callGeneric(e1, 0)
		  if(is0(f0)) {
		      e2@x <- callGeneric(e1, e2@x)
		      return(e2)
		  }
	      }
	      callGeneric(e1, as(e2, paste(.M.kind(e2), "gCMatrix", sep='')))
	  })



setMethod("Math",
	  signature(x = "CsparseMatrix"),
	  function(x) {
	      f0 <- callGeneric(0.)
	      if(is0(f0)) {
		  ## sparseness, symm., triang.,... preserved
		  x@x <- callGeneric(x@x)
		  x
	      } else { ## no sparseness
		  callGeneric(as_dense(x))
	      }
	  })



### workhorse for "[<-" -- both for d* and l*  C-sparse matrices :
replCmat <- function (x, i, j, value)
{
    di <- dim(x)
    dn <- dimnames(x)
    i1 <- if(missing(i)) 0:(di[1] - 1:1) else .ind.prep2(i, 1, di, dn)
    i2 <- if(missing(j)) 0:(di[2] - 1:1) else .ind.prep2(j, 2, di, dn)
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

    if(is(x, "symmetricMatrix")) ## only half the indices are there..
	x <- .Call(Csparse_symmetric_to_general, x)
    clx <- c(class(x))

    xj <- .Call(Matrix_expand_pointers, x@p)
    sel <- (!is.na(match(x@i, i1)) &
            !is.na(match( xj, i2)))
    has.x <- any("x" == slotNames(x)) # i.e. *not* logical
    if(has.x && sum(sel) == lenRepl) { ## all entries to be replaced are non-zero:
        value <- rep(value, length = lenRepl)
        ## Ideally we only replace them where value != 0 and drop the value==0
        ## ones; but that would have to (?) go through dgT*
        ## v0 <- 0 == value
        ## if (lenRepl == 1) and v0 is TRUE, the following is not doing anything
        ##-  --> ./dgTMatrix.R  and its  replTmat()
        ## x@x[sel[!v0]] <- value[!v0]
        x@x[sel] <- value
        return(x)
    }
    ## else go via Tsparse.. {FIXME: a waste! - we already have 'xj' ..}
    x <- as(x, "TsparseMatrix")
    if(missing(i))
	x[ ,j] <- value
    else if(missing(j))
	x[i, ] <- value
    else
	x[i,j] <- value
    as_CspClass(x, clx)
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


setMethod("crossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
              if (is(x, "symmetricMatrix")) {
                  warning("crossprod(x) calculated as x %*% x for sparse, symmetric x")
                  return(x %*% x)
              }
	      .Call(Csparse_crossprod, x, trans = FALSE, triplet = FALSE)
	  })

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

## FIXME (TODO):
## setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
## 	  function(x, y)
## 	  .Call(Csparse_crossprod_2, x, y, trans = TRUE, triplet = FALSE)


setMethod("%*%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y) .Call(Csparse_Csparse_prod, x, y))

setMethod("%*%", signature(x = "CsparseMatrix", y = "ddenseMatrix"),
          function(x, y) .Call(Csparse_dense_prod, x, y))

setMethod("%*%", signature(x = "CsparseMatrix", y = "matrix"),
          function(x, y) .Call(Csparse_dense_prod, x, y))

## Not needed because of c("Matrix", "numeric") method
##setMethod("%*%", signature(x = "CsparseMatrix", y = "numeric"),
##          function(x, y) .Call(Csparse_dense_prod, x, y))

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
	  function(x, nrow, ncol = n) {
	      dm <- .Call(Csparse_band, x, 0, 0)
	      dlen <- min(dm@Dim)
	      ind1 <- dm@i + 1:1	# 1-based index vector
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


setMethod("colSums", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("colMeans", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("rowSums", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("rowMeans", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
