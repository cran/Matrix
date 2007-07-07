#### All Methods in relation with the sparseVector (sub)classes


# atomicVector : classUnion (logical,integer,double,....)
setAs("atomicVector", "sparseVector",
      function(from) {
	  n <- length(from)
	  r <- new(paste(.V.kind(from), "sparseVector", sep=''), length = n)
	  ii <- from != 0
	  r@x <- from[ii]
	  r@i <- seq_len(n)[ii]
	  r
      })


for(T in c("d","i","l","z")) {
    setAs("xsparseVector", paste(T, "sparseVector", sep=''),
          function(from) {
              from@x <- as(from@x, .type.kind[T])
              ## and now "the hack":
              class(from) <- paste(T, "sparseVector", sep='')
              from
          })
}

setAs("sparseVector", "nsparseVector",
      function(from) {
          if(any(is.na(from@x)))
              stop("cannot coerce 'NA's to \"nsparseVector\"")
          new("nsparseVector", i = from@i, length = from@length)
      })

sp2vec <- function(x, mode = .type.kind[substr(cl, 1,1)]) {
    cl <- class(x)
    r <- vector(mode, x@length)
    r[x@i] <-
	if(cl != "nsparseVector") { # cheap test for 'has x slot'
	    if(is(x@x, mode)) x@x else as(x@x, mode)
	} else TRUE
    r
}

setAs("sparseVector", "vector", function(from) sp2vec(from))

setMethod("as.vector", signature(x = "sparseVector", mode = "missing"),
	  sp2vec)
setMethod("as.vector", signature(x = "sparseVector", mode = "character"),
	  sp2vec)

setMethod("as.numeric", "sparseVector", function(x) sp2vec(x, mode = "double"))

## the "catch all remaining" method:
setAs("ANY", "sparseVector",
      function(from) as(as.vector(from), "sparseVector"))

setAs("diagonalMatrix", "sparseVector",
      function(from) {
	  kind <- .M.kind(from) ## currently only "l" and "d" --> have 'x'
	  n <- nrow(from)
	  new(paste(kind, "sparseVector", sep=''),
	      length = n, # 1-based indexing
	      i = as.integer(seq(1L, by = n+1, length.out = n)),
	      x = if(from@diag != "U") from@x else
		  switch(kind, "d" = 1, "l" = TRUE, "i" = 1L, "z" = 1+0i))
	 })

setAs("sparseMatrix", "sparseVector",
      function(from) as(as(from, "TsparseMatrix"), "sparseVector"))

setAs("TsparseMatrix", "sparseVector",
      function(from) {
	  d <- dim(from)
	  n <- d[1] * d[2] # length of vector
	  kind <- .M.kind(from)
	  if(is_duplicatedT(from, nr = d[1]))
	      from <- uniqTsparse(from)
	  r <- new(paste(kind, "sparseVector", sep=''), length = n)
	  r@i <- 1L + from@i + d[1] * from@j
	  if(kind != "n") ## have 'x' slot
	      r@x <- from@x
	  r
      })



## TODO -- also want  (sparseVector, dim) |---> sparseMatrix
##  because of (nrow,ncol) specification can not (?)  use as(.).
##  Hence use  Matrix(.) ?  or my  spMatrix(.) ?

## For now, define this utility function:
spV2M <- function (x, nrow, ncol, byrow = FALSE)
{
    ## Purpose:	 sparseVector --> sparseMatrix	constructor
    ## ----------------------------------------------------------------------
    ## Arguments: x: "sparseVector" object
    ##		nrow, ncol, byrow: as for matrix() or Matrix()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 May 2007

    cx <- class(x)
    stopifnot(extends(cx, "sparseVector"))
    if(!missing(ncol)) { ncol <- as.integer(ncol)
			 if(ncol <= 0) stop("'ncol' must be >= 1") }
    if(!missing(nrow)) { nrow <- as.integer(nrow)
			 if(nrow <= 0) stop("'nrow' must be >= 1") }
    n <- length(x)
    if(missing(nrow)) {
	if(missing(ncol)) { ## both missing: --> (n x 1)
	    ncol <- 1L
	    nrow <- n
	} else {
	    if(n %% ncol != 0) warning("'ncol' is not a factor of length(x)")
	    nrow <- as.integer(ceiling(n / ncol))
	}
    } else {
	if(missing(ncol)) {
	    if(n %% nrow != 0) warning("'nrow' is not a factor of length(x)")
	    ncol <- as.integer(ceiling(n / nrow))
	} else { ## both nrow and ncol specified
	    if(ncol * nrow <  n) stop("nrow * ncol < length(x)")
	    if(ncol * nrow != n) warning("nrow * ncol != length(x)")
	}
    }
    ## now nrow * ncol >= n
    ##	   ~~~~~~~~~~~~~~~~
    cld <- getClassDef(cx)
    kind <- .M.kindC(cld)		# "d", "n", "l", "z", ...
    has.x <- kind != "n"
    r <- new(paste(kind,"gTMatrix", sep=''), Dim = c(nrow, ncol))
    ## now "compute"  the (i,j,x) slots given x@(i,x)
    i0 <- x@i - 1L
    if(byrow) {
	r@j <- i0 %% ncol
	r@i <- i0 %/% ncol
    } else {				# default{byrow = FALSE}
	r@i <- i0 %% nrow
	r@j <- i0 %/% nrow
    }
    if(has.x) r@x <- x@x
    r
}

setMethod("length", "sparseVector", function(x) x@length)

setMethod("show", signature(object = "sparseVector"),
   function(object) {
       n <- object@length
       cl <- class(object)
       cat(sprintf('sparse vector (nnz/length = %d/%d) of class "%s"\n',
		   length(object@i), n, cl))
       maxp <- max(1, getOption("max.print"))
       if(n <= maxp) {
	   prSpVector(object, maxp = maxp)
       } else { # n > maxp : will cut length of what we'll display :
	   ## cannot easily show head(.) & tail(.) because of "[1] .." printing of tail
	   prSpVector(object[seq_len(maxp)], maxp = maxp)
	   cat(" ............................",
	       "\n ........suppressing ", n - maxp,
	       " entries in show(); maybe adjust 'options(max.print= *)'",
	       "\n ............................\n\n", sep='')
       }
       invisible(object)
   })

prSpVector <- function(x, digits = getOption("digits"),
		    maxp = getOption("max.print"), zero.print = ".")
{
    cld <- getClassDef(cl <- class(x))
    stopifnot(extends(cld, "sparseVector"), maxp >= 1)
    if(is.logical(zero.print))
	zero.print <- if(zero.print) "0" else " "
##     kind <- .M.kindC(cld)
##     has.x <- kind != "n"
    n <- x@length
    if(n > maxp) {# n > maxp =: nn : will cut length of what we'll display :
	x <- x[seq_len(maxp)] # need "[" to work ...
	n <- as.integer(maxp)
    }
    xi <- x@i
    logi <- extends(cld, "lsparseVector") || extends(cld, "nsparseVector")
    cx <- if(logi) rep.int("N", n) else character(n)
    cx[ -xi ] <- zero.print
    cx[	 xi ] <- {
	if(logi) "|" else
	## numeric (or --not yet-- complex): 'has.x' in any cases
	format(x@x, digits = digits)
    }
    ## right = TRUE : cheap attempt to get better "." alignment
    print(cx, quote = FALSE, right = TRUE, max = maxp)
    invisible(x) # TODO? in case of n > maxp, "should" return original x
}

## This is a simplified intI() {-> ./Tsparse.R } -- for sparseVector indexing:
intIv <- function(i, n)
{
    ## Purpose: translate numeric | logical index     into  1-based integer
    ## --------------------------------------------------------------------
    ## Arguments: i: index vector (numeric | logical)
    ##		  n: array extent { ==	length(.) }
    if(missing(i))
	seq_len(n)
    else if(is(i, "numeric")) {
	storage.mode(i) <- "integer"
	if(any(i < 0L)) {
	    if(any(i > 0L))
		stop("you cannot mix negative and positive indices")
	    seq_len(n)[i]
	} else {
	    if(length(i) && max(i) > n)
		stop("indexing out of range 0:",n)
	    if(any(z <- i == 0))
		i <- i[!z]
	    i
	}
    }
    else if (is(i, "logical")) {
	seq_len(n)[i]
    } else stop("index must be numeric or logical for 'sparseVector' indexing")
}


setMethod("[", signature(x = "sparseVector", i = "index"),
	  function (x, i, j, drop) {
	      cld <- getClassDef(class(x))
	      has.x <- !extends(cld, "nsparseVector")
	      n <- x@length
	      ii <- intIv(i, n)
	      anyDup <- any(iDup <- duplicated(ii))
	      m <- match(x@i, ii, nomatch = 0)
	      sel <- m > 0L
	      x@length <- length(ii)
	      x@i <- m[sel]
	      if(anyDup) {
		  i.i <- match(ii[iDup], ii)
		  jm <- lapply(i.i, function(.) which(. == m))
		  sel <- c(which(sel), unlist(jm))
		  x@i <- c(x@i, rep.int(which(iDup), sapply(jm, length)))
	      }
	      if (has.x)
		  x@x <- x@x[sel]
	      x
	  })

## This is much analogous to replTmat in ./Tsparse.R:
replSPvec <- function (x, i, value)
{
    n <- x@length
    ii <- intIv(i, n)
    lenRepl <- length(ii)
    lenV <- length(value)
    if(lenV == 0) {
	if(lenRepl != 0)
	    stop("nothing to replace with")
	else return(x)
    }
    ## else: lenV := length(value) > 0
    if(lenRepl %% lenV != 0)
	stop("number of items to replace is not a multiple of replacement length")
    anyDup <- any(duplicated(ii))
    if(anyDup) { ## multiple *replacement* indices: last one wins
	## TODO: in R 2.6.0 use	 duplicate(*, fromLast=TRUE)
	ir <- lenRepl:1
	keep <- match(ii, ii[ir]) == ir
	ii <- ii[keep]
	lenV <- length(value <- rep(value, length = lenRepl)[keep])
	lenRepl <- length(ii)
    }

    cld <- getClassDef(class(x))
    has.x <- !extends(cld, "nsparseVector")
    m <- match(x@i, ii, nomatch = 0)
    sel <- m > 0L

    ## the simplest case
    if(all0(value)) { ## just drop the non-zero entries
	if(any(sel)) { ## non-zero there
	    x@i <- x@i[!sel]
	    if(has.x)
		x@x <- x@x[!sel]
	}
	return(x)

    }
    ## else --	some( value != 0 ) --
    if(lenV > lenRepl)
	stop("too many replacement values")
    else if(lenV < lenRepl)
	value <- rep(value, length = lenRepl)
    ## now:  length(value) == lenRepl

    v0 <- is0(value)
    ## value[1:lenRepl]:  which are structural 0 now, which not?

    if(any(sel)) {
	## indices of non-zero entries -- WRT to subvector
	iN0 <- m[sel] ## == match(x@i[sel], ii)

	## 1a) replace those that are already non-zero with new val.
	vN0 <- !v0[iN0]
	if(any(vN0) && has.x)
	    x@x[sel][vN0] <- value[iN0[vN0]]

	## 1b) replace non-zeros with 0 --> drop entries
	if(any(!vN0)) {
	    i <- which(sel)[!vN0]
	    if(has.x)
		x@x <- x@x[-i]
	    x@i <- x@i[-i]
	}
	iI0 <- if(length(iN0) < lenRepl)
	    seq_len(lenRepl)[-iN0]
    } else iI0 <- seq_len(lenRepl)

    if(length(iI0) && any(vN0 <- !v0[iI0])) {
	## 2) add those that were structural 0 (where value != 0)
	ij0 <- iI0[vN0]
	x@i <- c(x@i, ii[ij0])
	if(has.x)
	    x@x <- c(x@x, value[ij0])
    }
    x

}

setReplaceMethod("[", signature(x = "sparseVector", i = "index", j = "missing",
				value = "replValue"),
		 replSPvec)



## a "method" for c(<(sparse)Vector>, <(sparse)Vector>):
c2v <- function(x, y) {
    ## these as(., "sp..V..") check input implicitly:
    cx <- class(x <- as(x, "sparseVector"))
    cy <- class(y <- as(y, "sparseVector"))
    if(cx != cy) { ## find "common" class; result does have 'x' slot
        cxy <- c(cx,cy)
        commType <- {
            if(all(cxy %in% c("nsparseVector", "lsparseVector")))
                "lsparseVector"
            else { # ==> "numeric" ("integer") or "complex"
                xslot1 <- function(u, cl.u)
                    if(cl.u != "nsparseVector") u@x[1] else TRUE
                switch(typeof(xslot1(x, cx) + xslot1(y, cy)),
                       ## "integer", "double", or "complex"
                       "integer" = "isparseVector",
                       "double" = "dsparseVector",
                       "complex" = "zsparseVector")
            }
        }
        if(cx != commType) x <- as(x, commType)
        if(cy != commType) y <- as(y, commType)
        cx <- commType
    }
    ## now *have* common type -- transform 'x' into result:
    nx <- x@length
    x@length <- nx + y@length
    x@i <- c(x@i, nx + y@i)
    if(cx != "nsparseVector")
        x@x <- c(x@x, y@x)
    x
}


### Group Methods (!)

## o "Ops" , "Arith", "Compare"  :  ---> in ./Ops.R

