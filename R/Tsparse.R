#### "TsparseMatrix" : Virtual class of sparse matrices in triplet-format

setAs("TsparseMatrix", "CsparseMatrix",
      ## FIXME: this loses 'triangular' (and 'symmetric' ??)
      function(from)
      .Call(Tsparse_to_Csparse, from) ## ../src/Tsparse.c
      ## |-> cholmod_T -> cholmod_C -> chm_sparse_to_SEXP(* , 1)
      )

### "[" :
### -----

## Want to allow 'numeric', 'logical' and 'character' indices

## Test for numeric/logical/character
## method-*internally* ; this is not strictly OO, but allows to use
## the following utility and hence much more compact code.

.ind.prep <- function(xi, i, margin, di, dn)
{
    ## Purpose: do the ``common things'' for "*gTMatrix" indexing
    ##		for 1 dimension, 'margin' ,
    ##          and return match(.,.) + li = length of corresponding dimension
    ##
    ## i is "index";  xi = "x@i";  margin in {1,2};
    ## di = dim(x)      { used when i is "logical" }
    ## dn = dimnames(x) { used when i is character }

    dn <- dn[[margin]]
    has.dn <- is.character(dn)
    if(is(i, "numeric")) {
	storage.mode(i) <- "integer"
	if(any(ineg <- i < 0:0)) {
	    if(any(i > 0:0))
		stop("you cannot mix negative and positive indices")
	    i0 <- (0:(di[margin]-1:1))[i]
	} else {
	    if(length(i) && max(i) > di[margin])
		stop("indexing out of range 0:",di[margin])
	    if(any(z <- i == 0)) i <- i[!z]
	    i0 <- i - 1:1 # transform to 0-indexing
	}
	if(has.dn) dn <- dn[i]
    }
    else if (is(i, "logical")) {
	i0 <- (0:(di[margin]-1:1))[i]
	if(has.dn) dn <- dn[i]
    } else { ## character
	if(!has.dn)
	    stop(gettextf("no 'dimnames[[%d]]': cannot use character indexing"),
		 margin, domain = NA)
	i0 <- match(i, dn)
	if(any(is.na(i0))) stop("invalid character indexing")
	dn <- dn[i0]
	i0 <- i0 - 1:1
    }
    list(m = match(xi, i0, nomatch=0), li = length(i0), dn = dn)
}


.ind.prep2 <- function(i, margin, di, dn)
{
    ## Purpose: do the ``common things'' for "*gTMatrix" sub-assignment
    ##		for 1 dimension, 'margin' ,
    ##          and return match(.,.) + li = length of corresponding dimension
    ##
    ## i is "index"; margin in {1,2};
    ## di = dim(x)      { used when i is "logical" }

    dn <- dn[[margin]]
    has.dn <- is.character(dn)
    if(is(i, "numeric")) {
        storage.mode(i) <- "integer"
	if(any(ineg <- i < 0:0)) {
	    if(any(i > 0:0))
		stop("you cannot mix negative and positive indices")
	    i0 <- (1:di[margin])[i]
	} else	{
	    if(length(i) && max(i) > di[margin])
		stop("indexing out of range 0:",di[margin])
	    if(any(z <- i == 0)) i <- i[!z]
	    i0 <- i
	}
    }
    else if (is(i, "logical")) {
        i0 <- (1:di[margin])[i]
    } else { ## character
        if(!has.dn)
            stop(gettextf("no 'dimnames[[%d]]': cannot use character indexing"),
                 margin, domain = NA)
        i0 <- match(i, dn)
        if(any(is.na(i0))) stop("invalid character indexing")
    }
    i0 - 1:1  # transform to 0-indexing
}


## Otherwise have to write methods for all possible combinations of
##  (i , j) \in
##  (numeric, logical, character, missing) x (numeric, log., char., miss.)


## Select rows
setMethod("[", signature(x = "TsparseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select rows
	      ip <- .ind.prep(x@i, i, 1, dim(x), dimnames(x))
	      x@Dim[1] <- ip$li
	      if(!is.null(ip$dn)) x@Dimnames[[1]] <- ip$dn
	      sel <- ip$m > 0
	      x@i <- ip$m[sel] - 1:1
	      x@j <- x@j[sel]
	      if (!is(x, "lsparseMatrix")) x@x <- x@x[sel]
	      if (drop && any(x@Dim == 1:1)) drop(as(x,"matrix")) else x
	  })


## Select columns
setMethod("[", signature(x = "TsparseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select columns
	      ip <- .ind.prep(x@j, j, 2, dim(x), dimnames(x))
	      x@Dim[2] <- ip$li
	      if(!is.null(ip$dn)) x@Dimnames[[2]] <- ip$dn
	      sel <- ip$m > 0
	      x@i <- x@i[sel]
	      x@j <- ip$m[sel] - 1:1
	      if (!is(x, "lsparseMatrix")) x@x <- x@x[sel]
	      if (drop && any(x@Dim == 1:1)) drop(as(x,"matrix")) else x
	  })


## [.data.frame has : drop = if (missing(i)) TRUE else length(cols) == 1)

setMethod("[", signature(x = "TsparseMatrix",
			 i = "index", j = "index", drop = "logical"),
	  function (x, i, j, ..., drop)
      {
	  ## (i,j, drop) all specified
          di <- dim(x)
          dn <- dimnames(x)
          ip1 <- .ind.prep(x@i, i, 1, di, dn)
          ip2 <- .ind.prep(x@j, j, 2, di, dn)
          x@Dim <- nd <- c(ip1$li, ip2$li)
          x@Dimnames <- list(ip1$dn, ip2$dn)
          sel <- ip1$m > 0:0  &  ip2$m > 0:0
          x@i <- ip1$m[sel] - 1:1
          x@j <- ip2$m[sel] - 1:1
          if (!is(x, "lsparseMatrix")) x@x <- x@x[sel]
	  if (drop && any(nd == 1)) drop(as(x,"matrix")) else x
      })


## workhorse for "[<-" :
replTmat <- function (x, i, j, value)
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

    ## Note: *T*matrix maybe non-unique: an entry can be split
    ##	  into a *sum* of several ones :
    x <- uniq(x) # -> ./Auxiliaries.R

    sel <- ((m1 <- match(x@i, i1, nomatch=0)) > 0:0 &
            (m2 <- match(x@j, i2, nomatch=0)) > 0:0)

    has.x <- any("x" == slotNames(x)) # i.e. *not* logical

    ## the simplest case: for all Tsparse, even for i or j missing
    if(all(value == 0)) { ## just drop the non-zero entries
	if(any(sel)) { ## non-zero there
	    x@i <- x@i[!sel]
	    x@j <- x@j[!sel]
            if(has.x)
		x@x <- x@x[!sel]
	}
	return(x)
    }

    ## else --  some( value != 0 ) --
    if(lenV > lenRepl)
        stop("too many replacement values")

    ## another simple, typical case:
    if(lenRepl == 1) {
        if(any(sel)) { ## non-zero there
            if(has.x)
                x@x[sel] <- value
        } else { ## new non-zero
            x@i <- c(x@i, i1)
            x@j <- c(x@j, i2)
            if(has.x)
                x@x <- c(x@x, value)
        }
        return(x)
    }

    v0 <- 0 == (value <- rep(value, length = lenRepl))
    ## value[1:lenRepl]:  which are structural 0 now, which not?

    if(any(sel)) {
	## the 0-based indices of non-zero -- WRT to submatrix
	non0 <- cbind(match(x@i[sel], i1),
		      match(x@j[sel], i2)) - 1:1
	iN0 <- 1:1 + encodeInd(non0, nr = dind[1])

	## 1) replace those that are already non-zero (when value != 0)
	vN0 <- !v0[iN0]
        if(has.x)
            x@x[sel[vN0]] <- value[iN0[vN0]]

	iI0 <- (1:lenRepl)[-iN0]	# == complementInd(non0, dind)
    } else iI0 <- 1:lenRepl

    if(length(iI0)) {
	## 2) add those that were structural 0 (where value != 0)
	vN0 <- !v0[iI0]
	ij0 <- decodeInd(iI0[vN0] - 1:1, nr = dind[1])
	x@i <- c(x@i, i1[ij0[,1] + 1:1])
	x@j <- c(x@j, i2[ij0[,2] + 1:1])
        if(has.x)
            x@x <- c(x@x, value[iI0[vN0]])
    }
    x
}

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "index", j = "missing",
                                value = "replValue"),
                 function (x, i, value) replTmat(x, i=i, value=value))

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "missing", j = "index",
                                value = "replValue"),
                 function (x, j, value) replTmat(x, j=j, value=value))

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "index", j = "index",
				value = "replValue"),
                 replTmat)




setMethod("crossprod", signature(x = "TsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      a <- .Call(Csparse_crossprod, x, trans = FALSE, triplet = TRUE)
	      switch(substr(class(a)[1], 1, 1),
		     "d" ={ new("dsCMatrix", i = a@i, p = a@p, x = a@x,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "U",
				factors = list()) },
		     "l" ={ new("lsCMatrix", i = a@i, p = a@p,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "U",
				factors = list())})
	  })

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      a <- .Call(Csparse_crossprod, x, trans = TRUE, triplet = TRUE)
	      switch(substr(class(a)[1], 1, 1),
		     "d" ={ new("dsCMatrix", i = a@i, p = a@p, x = a@x,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "L",
				factors = list()) },
		     "l" ={ new("lsCMatrix", i = a@i, p = a@p,
				Dim = a@Dim, Dimnames = a@Dimnames, uplo = "L",
				factors = list()) })
	  })

setMethod("colSums", signature(x = "TsparseMatrix"), .as.dgT.Fun,
	  valueClass = "numeric")
setMethod("colMeans", signature(x = "TsparseMatrix"), .as.dgT.Fun,
	  valueClass = "numeric")
setMethod("rowSums", signature(x = "TsparseMatrix"), .as.dgT.Fun,
	  valueClass = "numeric")
setMethod("rowMeans", signature(x = "TsparseMatrix"), .as.dgT.Fun,
	  valueClass = "numeric")
