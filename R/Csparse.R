#### Methods for the virtual class 'CsparseMatrix' of sparse matrices stored in
####  "column compressed" format.
#### -- many more specific things are e.g. in ./dgCMatrix.R

## FIXME: the is(*,"generalMatrix") test at least makes these work,
##        but they are still ``wrong'', since triangularity is lost  :

setAs("CsparseMatrix", "TsparseMatrix",
      function(from) {
	  if(!is(from, "generalMatrix")) { ## e.g. for triangular | symmetric
	      if     (is(from, "dMatrix")) from <- as(from, "dgCMatrix")
	      else if(is(from, "lMatrix")) from <- as(from, "lgCMatrix")
	      else if(is(from, "zMatrix")) from <- as(from, "zgCMatrix")
	      else stop("undefined method for class ", class(from))
	  }
          .Call(Csparse_to_Tsparse, from) # ../src/Csparse.c
          ## |-> cholmod_C -> cholmod_T -> chm_triplet_to_SEXP(* , 1)
      })

setAs("CsparseMatrix", "denseMatrix",
      function(from) {
	  if(!is(from, "generalMatrix")) { ## e.g. for triangular | symmetric
	      if     (is(from, "dMatrix")) from <- as(from, "dgCMatrix")
	      else if(is(from, "lMatrix")) from <- as(from, "lgCMatrix")
	      else if(is(from, "zMatrix")) from <- as(from, "zgCMatrix")
	      else stop("undefined method for class ", class(from))
	  }
          .Call(Csparse_to_dense, from)
      })

### Some group methods:

setMethod("Arith",
	  signature(e1 = "CsparseMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      if(length(e2) == 1) { ## e.g.,  Mat ^ a
		  f0 <- callGeneric(0, e2)
		  if(!is.na(f0) && f0 == 0.) { # remain sparse, symm., tri.,...
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
		  if(!is.na(f0) && f0 == 0.) {
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
	      if(!is.na(f0) && f0 == 0.) {
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

    clx <- c(class(x))
    xj <- .Call(Matrix_expand_pointers, x@p)
    sel <- (!is.na(match(x@i, i1)) &
            !is.na(match( xj, i2)))

    has.x <- any("x" == slotNames(x)) # i.e. *not* logical

    if(sum(sel) == lenRepl) { ## all entries to be replaced are non-zero:
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
    ## else go via Tsparse..
    x <- as(x, "TsparseMatrix")
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
	      a <- .Call(Csparse_crossprod, x, trans = FALSE, triplet = FALSE)
	      switch(substr(class(a)[1], 1, 1),
		     "d" ={ new("dsCMatrix", i = a@i, p = a@p, x = a@x,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "U",
				factors = list()) },
		     "l" ={ new("lsCMatrix", i = a@i, p = a@p,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "U",
				factors = list()) })
	  })


setMethod("t", signature(x = "CsparseMatrix"),
	  function(x)
	  .Call(Csparse_transpose, x))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      a <- .Call(Csparse_crossprod, x, trans = TRUE, triplet = FALSE)
	      switch(substr(class(a)[1], 1, 1),
		     "d" ={ new("dsCMatrix", i = a@i, p = a@p, x = a@x,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "L",
				factors = list()) },
		     "l" ={ new("lsCMatrix", i = a@i, p = a@p,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "L",
				factors = list()) })
	  })

## FIXME (TODO):
## setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
## 	  function(x, y)
## 	  .Call(Csparse_crossprod_2, x, y, trans = TRUE, triplet = FALSE)


setMethod("%*%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y) .Call(Csparse_Csparse_prod, x, y))

setMethod("%*%", signature(x = "CsparseMatrix", y = "denseMatrix"),
          function(x, y) .Call(Csparse_dense_prod, x, y))

## NB: have extra tril(), triu() methods for symmetric ["dsC" and "lsC"]
setMethod("tril", "CsparseMatrix",
	  function(x, k = 0, ...) {
	      k <- as.integer(k[1])
	      dd <- dim(x)
	      stopifnot(-dd[1] <= k, k <= dd[1]) # had k <= 0
	      r <- .Call(Csparse_band, x, -dd[1], k)
	      ## return "lower triangular" if k <= 0
	      if(k <= 0) as(r, paste(.M.kind(x), "tCMatrix", sep='')) else r
	  })

setMethod("triu", "CsparseMatrix",
	  function(x, k = 0, ...) {
	      k <- as.integer(k[1])
	      dd <- dim(x)
	      stopifnot(-dd[1] <= k, k <= dd[1]) # had k >= 0
	      r <- .Call(Csparse_band, x, k, dd[2])
	      ## return "upper triangular" if k >= 0
	      if(k >= 0) as(r, paste(.M.kind(x), "tCMatrix", sep='')) else r
	  })

setMethod("band", "CsparseMatrix",
	  function(x, k1, k2, ...) {
	      k1 <- as.integer(k1[1])
	      k2 <- as.integer(k2[1])
	      dd <- dim(x)
	      stopifnot(-dd[1] <= k1, k1 <= k2, k2 <= dd[1])
	      r <- .Call(Csparse_band, x, k1, k2)
	      if(k1 * k2 >= 0) ## triangular
		  as(r, paste(.M.kind(x), "tCMatrix", sep=''))
	      else if (k1 < 0  &&  k1 == -k2  && isSymmetric(x)) ## symmetric
		  as(r, paste(.M.kind(x), "sCMatrix", sep=''))
	      else
		  r
	  })

setMethod("colSums", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("colMeans", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("rowSums", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
setMethod("rowMeans", signature(x = "CsparseMatrix"), .as.dgC.Fun,
	  valueClass = "numeric")
