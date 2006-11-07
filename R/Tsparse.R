#### "TsparseMatrix" : Virtual class of sparse matrices in triplet-format

setAs("TsparseMatrix", "CsparseMatrix",
      ## |-> cholmod_T -> cholmod_C -> chm_sparse_to_SEXP
      ## adjusted for triangular matrices not represented in cholmod
      function(from) .Call(Tsparse_to_Csparse, from, ## ../src/Tsparse.c
                           is(from, "triangularMatrix"))
      )

## special cases

setAs("dgTMatrix", "dgCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("dsTMatrix", "dsCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("dtTMatrix", "dtCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, TRUE))

setAs("ngTMatrix", "ngCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))
setAs("lgTMatrix", "lgCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

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
{    ## Purpose: do the ``common things'' for "*gTMatrix" sub-assignment
    ##		for 1 dimension, 'margin' ,
    ##          and return match(.,.) + li = length of corresponding dimension
    ##
    ## i is "index"; margin in {1,2};
    ## di = dim(x)      { used when i is "logical" }

    ## difference to .ind.prep(): use 1-indices; no match(xi,..), no dn at end
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
              if(is(x, "symmetricMatrix"))
		  x <- as(x, paste(.M.kind(x), "geMatrix", sep=''))
	      ip <- .ind.prep(x@i, i, 1, dim(x), dimnames(x))
	      x@Dim[1] <- ip$li
	      if(!is.null(ip$dn)) x@Dimnames[[1]] <- ip$dn
	      sel <- ip$m > 0
	      x@i <- ip$m[sel] - 1:1
	      x@j <- x@j[sel]
	      if (!is(x, "nsparseMatrix")) x@x <- x@x[sel]
	      if (drop && any(x@Dim == 1:1)) drop(as(x,"matrix")) else x
	  })


## Select columns
setMethod("[", signature(x = "TsparseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select columns
              if(is(x, "symmetricMatrix"))
		  x <- as(x, paste(.M.kind(x), "geMatrix", sep=''))
	      ip <- .ind.prep(x@j, j, 2, dim(x), dimnames(x))
	      x@Dim[2] <- ip$li
	      if(!is.null(ip$dn)) x@Dimnames[[2]] <- ip$dn
	      sel <- ip$m > 0
	      x@i <- x@i[sel]
	      x@j <- ip$m[sel] - 1:1
	      if (!is(x, "nsparseMatrix")) x@x <- x@x[sel]
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
	  if(is(x, "symmetricMatrix")) {
	      isSym <- length(i) == length(j) && all(i == j)
	      ## result is *still* symmetric --> keep symmetry!
	      if(!isSym)
		  ## result no longer symmetric -> to "generalMatrix"
		  x <- as(x, paste(.M.kind(x), "gTMatrix", sep=''))
	  } else isSym <- FALSE
	  if(isSym) {
	      offD <- x@i != x@j
	      ip1 <- .ind.prep(c(x@i,x@j[offD]), i, 1, di, dn)
	      ip2 <- .ind.prep(c(x@j,x@i[offD]), j, 2, di, dn)
	  } else {
	      ip1 <- .ind.prep(x@i, i, 1, di, dn)
	      ip2 <- .ind.prep(x@j, j, 2, di, dn)
	  }
	  x@Dim <- nd <- c(ip1$li, ip2$li)
	  x@Dimnames <- list(ip1$dn, ip2$dn)

	  sel <- ip1$m > 0:0  &	 ip2$m > 0:0
	  if(isSym) { # only those corresponding to upper/lower triangle
	      sel <- sel &
	      (if(x@uplo == "U") ip1$m <= ip2$m else ip2$m <= ip1$m)
	  }
	  x@i <- ip1$m[sel] - 1:1
	  x@j <- ip2$m[sel] - 1:1
	  if (!is(x, "nsparseMatrix"))
	      x@x <- c(x@x, if(isSym) x@x[offD])[sel]
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

    get.ind.sel <- function(ii,ij)
	(match(x@i, ii, nomatch = 0) > 0:0 &
	 match(x@j, ij, nomatch = 0) > 0:0)

    if((sym.x <- is(x, "symmetricMatrix"))) {
	r.sym <- dind[1] == dind[2] && i1 == i2 &&
	(lenRepl == 1 || isSymmetric(value <- array(value, dim=dind)))
	if(r.sym) { ## result is *still* symmetric --> keep symmetry!
	    ## now consider only those indices above / below diagonal:
	    xU <- x@uplo == "U"
	    useI <- if(xU) i1 <= i2 else i2 <= i1
	    i1 <- i1[useI]
	    i2 <- i2[useI]
	    ## select also the corresponding triangle
	    if(lenRepl > 1)
		value <- value[(if(xU)upper.tri else lower.tri)(value, diag=TRUE)]
	}
	else { # go to "generalMatrix" and continue
	    x <- as(x, paste(.M.kind(x), "gTMatrix", sep=''))
	}
    }

    sel <- get.ind.sel(i1,i2)
    has.x <- any("x" == slotNames(x)) # i.e. *not* nonzero-pattern

    ## the simplest case: for all Tsparse, even for i or j missing
    if(all0(value)) { ## just drop the non-zero entries
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

    if(sym.x && r.sym)
       lenRepl <- length(value) # shorter (since only "triangle")
    else
       value <- rep(value, length = lenRepl)

    v0 <- is0(value)
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
              if (is(x, "symmetricMatrix")) {
                  x <- as(x, "CsparseMatrix")
                  warning("crossprod(x) calculated as x %*% x for sparse, symmetric x")
                  return(x %*% x)
              }
	      .Call(Csparse_crossprod, x, trans = FALSE, triplet = TRUE)
	  })

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
	      .Call(Csparse_crossprod, x, trans = TRUE, triplet = TRUE)
	  })

## Must define methods for y = "missing" first so they have precedence
## (this will change in R-2.4.0).

setMethod("crossprod", signature(x = "TsparseMatrix", y = "ANY"),
	  function(x, y = NULL) callGeneric(as(x, "CsparseMatrix"), y))

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "ANY"),
	  function(x, y = NULL) callGeneric(as(x, "CsparseMatrix"), y))

setMethod("%*%", signature(x = "TsparseMatrix", y = "ANY"),
          function(x, y) callGeneric(as(x, "CsparseMatrix"), y))

setMethod("%*%", signature(x = "ANY", y = "TsparseMatrix"),
          function(x, y) callGeneric(x, as(y, "CsparseMatrix")))

## Not yet.  Don't have methods for y = "CsparseMatrix" and general x
#setMethod("%*%", signature(x = "ANY", y = "TsparseMatrix"),
#          function(x, y) callGeneric(x, as(y, "CsparseMatrix")))

setMethod("colSums", signature(x = "TsparseMatrix"), .as.dgT.Fun,
	  valueClass = "numeric")
setMethod("colMeans", signature(x = "TsparseMatrix"), .as.dgT.Fun,
	  valueClass = "numeric")
setMethod("rowSums", signature(x = "TsparseMatrix"), .as.dgT.Fun,
	  valueClass = "numeric")
setMethod("rowMeans", signature(x = "TsparseMatrix"), .as.dgT.Fun,
	  valueClass = "numeric")

## Want tril(), triu(), band() --- just as "indexing" ---
## return a "close" class:
setMethod("tril", "TsparseMatrix",
	  function(x, k = 0, ...) as_Tsparse(tril(as_Csparse(x), k = k, ...)))
setMethod("triu", "TsparseMatrix",
	  function(x, k = 0, ...) as_Tsparse(triu(as_Csparse(x), k = k, ...)))
setMethod("band", "TsparseMatrix",
	  function(x, k1, k2, ...)
	  as_Tsparse(band(as_Csparse(x), k1 = k1, k2 = k2, ...)))

setMethod("t", signature(x = "TsparseMatrix"),
	  function(x) {
	      r <- new(class(x))
	      r@i <- x@j
	      r@j <- x@i
	      if(any("x" == slotNames(x)))
		  r@x <- x@x
	      r@Dim <- rev(x@Dim)
	      r@Dimnames <- rev(x@Dimnames)
	      r
      })

