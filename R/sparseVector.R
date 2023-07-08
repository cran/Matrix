## METHODS FOR CLASS: sparseVector (virtual)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("nsparseVector", "lsparseVector",
      function(from)
          new("lsparseVector", length = from@length, i = from@i,
              x = rep.int(TRUE, length(from@i))))
setAs("nsparseVector", "isparseVector",
      function(from)
          new("isparseVector", length = from@length, i = from@i,
              x = rep.int(1L, length(from@i))))
setAs("nsparseVector", "dsparseVector",
      function(from)
          new("dsparseVector", length = from@length, i = from@i,
              x = rep.int(1, length(from@i))))
setAs("nsparseVector", "zsparseVector",
      function(from)
          new("zsparseVector", length = from@length, i = from@i,
              x = rep.int(1+0i, length(from@i))))

setAs("sparseVector", "nsparseVector",
      function(from)
	  new("nsparseVector", length = from@length, i = from@i))
setAs("sparseVector", "lsparseVector",
      function(from)
          new("lsparseVector", length = from@length, i = from@i,
              x = as.logical(from@x)))
setAs("sparseVector", "isparseVector",
      function(from)
          new("isparseVector", length = from@length, i = from@i,
              x = as.integer(from@x)))
setAs("sparseVector", "dsparseVector",
      function(from)
          new("dsparseVector", length = from@length, i = from@i,
              x = as.double(from@x)))
setAs("sparseVector", "zsparseVector",
      function(from)
          new("zsparseVector", length = from@length, i = from@i,
              x = as.complex(from@x)))

spV2M <- function(x, nrow, ncol, byrow = FALSE,
                  check = TRUE, symmetric = FALSE) {
    if(check && !is(x, "sparseVector"))
	stop("'x' must inherit from \"sparseVector\"")
    if(!missing(ncol)) { ncol <- as.integer(ncol)
			 if(ncol < 0) stop("'ncol' must be >= 0") }
    if(!missing(nrow)) { nrow <- as.integer(nrow)
			 if(nrow < 0) stop("'nrow' must be >= 0") }
    n <- length(x)
    if(symmetric) {
	if(missing(nrow)) stop("Must specify 'nrow' when 'symmetric' is true")
	if(!missing(ncol) && nrow != ncol)
	    stop("'nrow' and 'ncol' must be the same when 'symmetric' is true")
	## otherwise  ncol will not used at all when (symmetric)
	if(check && as.double(nrow)^2 != n)
	    stop("'x' must have length nrow^2 when 'symmetric' is true")
	## x <- x[indTri(nrow, upper=TRUE, diag=TRUE)]
    } else if(missing(nrow)) {
	nrow <- as.integer(
	    if(missing(ncol)) { ## both missing: --> (n x 1)
		ncol <- 1L
		n
	    } else {
		if(n %% ncol != 0) warning("'ncol' is not a factor of length(x)")
		as.integer(ceiling(n / ncol))
	    })
    } else if(missing(ncol)) {
        ncol <- if(symmetric) nrow else {
            if(n %% nrow != 0) warning("'nrow' is not a factor of length(x)")
            as.integer(ceiling(n / nrow)) }
    } else {                          ## both nrow and ncol specified
        n.n <- as.double(ncol) * nrow # no integer overflow
        if(n.n <  n) stop("nrow * ncol < length(x)", domain = NA)
        if(n.n != n) warning("nrow * ncol != length(x)", domain = NA)
    }
    ## now nrow * ncol >= n  (or 'symmetric')
    ##	   ~~~~~~~~~~~~~~~~
    kind <- .M.kind(x) # "d", "n", "l", "i", "z", ...
    has.x <- kind != "n"
    clStem <- if(symmetric) "sTMatrix" else "gTMatrix"
    ## "careful_new()" :
    cNam <- paste0(kind, clStem)
    chngCl <- is.null(newCl <- getClassDef(cNam))
    if(chngCl) { ## e.g. "igTMatrix" is not yet implemented
	if(kind == "z")
	    stop(gettextf("Class %s is not yet implemented", dQuote(cNam)),
		 domain = NA)
	## coerce to "double":
	newCl <- getClassDef(paste0("d", clStem))
    }
    r <- new(newCl, Dim = c(nrow, ncol))
    ## now "compute"  the (i,j,x) slots given x@(i,x)
    i0 <- x@i - 1L
    if(byrow) { ## need as.integer(.) since <sparseVector> @ i can be double
	j <- as.integer(i0 %% ncol)
	i <- as.integer(i0 %/% ncol)
    } else { ## default{byrow = FALSE}
	i <- as.integer(i0 %% nrow)
	j <- as.integer(i0 %/% nrow)
    }
    if(has.x)
	x <- if(chngCl) as.numeric(x@x) else x@x
    if(symmetric) {  ## using  uplo = "U"
	i0 <- i <= j ## i.e., indTri(nrow, upper=TRUE, diag=TRUE)
	i <- i[i0]
	j <- j[i0]
	if(has.x) x <- x[i0]
    }
    r@j <- j
    r@i <- i
    if(has.x) r@x <- x
    r
}

.sparseV2Mat <- function(from)
    spV2M(from, nrow = from@length, ncol = 1L, check = FALSE)

setAs("sparseVector",        "Matrix", .sparseV2Mat)
setAs("sparseVector",  "sparseMatrix", .sparseV2Mat)
setAs("sparseVector", "TsparseMatrix", .sparseV2Mat)
setAs("sparseVector", "CsparseMatrix", function(from) .T2C(.sparseV2Mat(from)))
setAs("sparseVector", "RsparseMatrix", function(from) .T2R(.sparseV2Mat(from)))

sp2vec <- function(x, mode = .type.kind[.M.kind(x)]) {
    ## sparseVector  ->  vector
    has.x <- .hasSlot(x, "x")## has "x" slot
    m.any <- (mode == "any")
    if(m.any)
	mode <- if(has.x) mode(x@x) else "logical"
    else if(has.x) # is.<mode>() is much faster than inherits() | is():
        xxOk <- switch(mode,
		       "double" = is.double(x@x),
		       "logical" = is.logical(x@x),
		       "integer" = is.integer(x@x),
		       "complex" = is.complex(x@x),
		       ## otherwise (does not happen with default 'mode'):
		       inherits(x@x, mode))
    r <- vector(mode, x@length)
    r[x@i] <-
	if(has.x) {
	    if(m.any || xxOk) x@x else as(x@x, mode)
	} else TRUE
    r
}

## Need 'base' functions calling as.*() to dispatch to our S4 methods:
as.vector.sparseVector <- sp2vec
as.matrix.sparseVector <- function(x, ...) as.matrix.default(sp2vec(x))
 as.array.sparseVector <- function(x, ...)  as.array.default(sp2vec(x))

setAs("sparseVector", "vector",  function(from) sp2vec(from))
setAs("sparseVector", "logical", function(from) sp2vec(from, mode = "logical"))
setAs("sparseVector", "integer", function(from) sp2vec(from, mode = "integer"))
setAs("sparseVector", "numeric", function(from) sp2vec(from, mode = "double"))

setMethod("as.vector",  "sparseVector", sp2vec)
setMethod("as.logical", "sparseVector", function(x) sp2vec(x, mode = "logical"))
setMethod("as.numeric", "sparseVector", function(x) sp2vec(x, mode = "double"))


## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("ANY", "sparseVector",
      function(from) as(as.vector(from), "sparseVector"))

setAs("ANY", "nsparseVector",
      function(from) as(as(from, "sparseVector"), "nsparseVector"))

setAs("atomicVector", "sparseVector",
      function(from) {
	  r <- new(paste0(.V.kind(from), "sparseVector"))
	  r@length <- length(from)
	  r@i <- ii <- which(isN0(from))
          r@x <- from[ii]
	  r
      })

setAs("atomicVector", "dsparseVector",
      function(from) {
	  r <- new("dsparseVector")
	  r@length <- length(from)
          r@i <- ii <- which(isN0(from))
	  r@x <- as.double(from)[ii]
	  r
      })

setAs("CsparseMatrix", "sparseVector",
      function(from) .Call(CR2spV, from))

setAs("RsparseMatrix", "sparseVector",
      function(from) .Call(CR2spV, from))

setAs("TsparseMatrix", "sparseVector",
      function(from) .Call(CR2spV, .T2C(from)))

setAs("diagonalMatrix", "sparseVector",
      function(from) {
          n <- (d <- from@Dim)[1L]
          nn <- prod(d)
          kind <- .M.kind(from)
          to <- new(paste0(kind, "sparseVector"))
          to@length <-
              if(nn <= .Machine$integer.max)
                  as.integer(nn)
              else nn
          to@i <- indDiag(n)
          to@x <-
              if(from@diag == "N")
                  from@x
              else rep.int(switch(kind,
                                  "l" = TRUE,
                                  "i" = 1L,
                                  "d" = 1,
                                  "z" = 1+0i),
                           n)
          to
      })

setAs("indMatrix", "sparseVector",
      function(from) {
          d <- from@Dim
          m <- d[1L]
          n <- d[2L]
          mn <- prod(d)
          perm <- from@perm
          to <- new("nsparseVector")
          if(mn <= .Machine$integer.max) {
              to@length <- as.integer(mn)
              to@i <-
                  if(from@margin == 1L)
                      seq.int(to = 0L, by = 1L, length.out = m) + perm * m
                  else seq.int(from = 0L, by = m, length.out = n) + perm
          } else {
              to@length <- mn
              to@i <-
                  if(from@margin == 1L)
                      seq.int(to = 0, by = 1, length.out = m) + perm * as.double(m)
                  else seq.int(from = 0, by = as.double(m), length.out = n) + as.double(perm)
          }
          to
      })


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##' Construct new sparse vector , *dropping* zeros

##' @param class  character, the sparseVector class
##' @param x      numeric/logical/...:  the 'x' slot -- if missing ==> "nsparseVector"
##' @param i      integer: index of non-zero entries
##' @param length integer: the 'length' slot

##' @return a sparseVector, with 0-dropped 'x' (and 'i')
newSpV <- function(class, x, i, length, drop0 = TRUE, checkSort = TRUE) {
    if(has.x <- !missing(x)) {
	if(length(x) == 1 && (li <- length(i)) != 1) ## recycle x :
	    x <- rep.int(x, li)
	if(drop0 && isTRUE(any(x0 <- x == 0))) {
	    keep <- is.na(x) | !x0
	    x <- x[keep]
	    i <- i[keep]
	}
    }
    if(checkSort && is.unsorted(i)) {
	ii <- sort.list(i)
	if(has.x) x <- x[ii]
	i <- i[ii]
    }
    if(has.x)
	new(class, x = x, i = i, length = length)
    else
	new(class,        i = i, length = length)
}
## a "version" of 'prev' with changed contents:
newSpVec <- function(class, x, prev)
    newSpV(class, x=x, i=prev@i, length=prev@length)

## Exported:
sparseVector <- function(x, i, length) {
    newSpV(class = paste0(if(missing(x)) "n" else .V.kind(x), "sparseVector"),
           x=x, i=i, length=length)
}

setMethod("dim<-", signature(x = "sparseVector"),
	  function(x, value) {
	      if(!is.numeric(value) || length(value) != 2L)
		  stop("dimensions must be numeric of length 2")
              if(anyNA(value))
		  stop("dimensions cannot contain NA")
              if(any(value < 0))
                  stop("dimensions cannot contain negative values")
              if(!is.integer(value)) {
                  if(any(value > .Machine$integer.max))
                      stop("dimensions cannot exceed 2^31-1")
                  value <- as.integer(value)
              }
	      if((p <- prod(value)) != (len <- length(x)))
		  stop(gettextf("assigned dimensions [product %.0f] do not match object length [%.0f]",
                                p, len, domain = NA))
	      spV2M(x, nrow = value[1L], ncol = value[2L])
	  })

setMethod("length", "sparseVector", function(x) x@length)

setMethod("mean", signature(x = "sparseVector"),
	  function(x, trim = 0, na.rm = FALSE, ...) {
              if(is.numeric(trim) && length(trim) == 1L && !is.na(trim) &&
                 trim == 0) {
                  ## Be fast in this special case :
                  if(isTRUE(na.rm))
                      ## FIXME: don't allocate !is.na(x)
                      x <- x[!is.na(x)]
                  sum(x) / length(x)
	      } else {
                  ## FIXME: don't allocate as.numeric(x); need 'sort' method
                  warning("suboptimally using as.numeric(x) to compute trimmed mean of sparseVector 'x'")
                  mean.default(as.numeric(x), trim = trim, na.rm = na.rm, ...)
              }
          })

setMethod("t", "sparseVector",
          function(x) spV2M(x, nrow = 1L, ncol = x@length, check = FALSE))

setMethod("show", signature(object = "sparseVector"),
   function(object) {
       n <- object@length
       cl <- class(object)
       cat(sprintf('sparse vector (nnz/length = %d/%.0f) of class "%s"\n',
		   length(object@i), as.double(n), cl))
       maxp <- max(1, getOption("max.print"))
       if(n <= maxp) {
	   prSpVector(object, maxp = maxp)
       } else { # n > maxp : will cut length of what we'll display :
	   ## cannot easily show head(.) & tail(.) because of "[1] .." printing of tail
	   prSpVector(head(object, maxp), maxp = maxp)
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
    cld <- getClassDef(class(x))
    stopifnot(extends(cld, "sparseVector"), maxp >= 1)
    if(is.logical(zero.print))
	zero.print <- if(zero.print) "0" else " "
##     kind <- .M.kindC(cld)
##     has.x <- kind != "n"
    n <- x@length
    if(n > 0) {
	if(n > maxp) { # n > maxp =: nn : will cut length of what we'll display :
	    x <- head(x, maxp)
	    n <- maxp
	}
        xi <- x@i
        is.n <- extends(cld, "nsparseVector")
        logi <- is.n || extends(cld, "lsparseVector")
        cx <- if(logi) rep.int("N", n) else character(n)
        cx[if(length(xi)) -xi else TRUE] <- zero.print
        cx[ xi] <- {
	    if(is.n) "|" else if(logi) c(":","|")[x@x + 1L] else
	    ## numeric (or --not yet-- complex): 'has.x' in any cases
	    format(x@x, digits = digits)
        }
        ## right = TRUE : cheap attempt to get better "." alignment
        print(cx, quote = FALSE, right = TRUE, max = maxp)
    }
    invisible(x) # TODO? in case of n > maxp, "should" return original x
}

## This is a simplified intI() {-> ./Tsparse.R } -- for sparseVector indexing:
intIv <- function(i, n, cl.i = getClass(class(i)))
{
### Note: undesirable to use this for negative indices;
### ----  using seq_len(n) below means we are  NON-sparse ...
### Fixed, for "x[i] with negative i" at least.

    ## Purpose: translate numeric | logical index     into  1-based integer
    ## --------------------------------------------------------------------
    ## Arguments: i: index vector (numeric | logical) *OR* sparseVector
    ##		  n: array extent { ==	length(.) }
    if(missing(i))
	seq_len(n)
    else if(extends(cl.i, "numeric")) {
        ## not ok, when max(i) > .Machine$integer.max !  storage.mode(i) <- "integer"
        int2i(i,n) ##-> ./Tsparse.R
    }
    else if (extends(cl.i, "logical")) {
	seq_len(n)[i]
    } else if(extends(cl.i, "nsparseVector")) {
	i@i # the indices are already there !
    } else if(extends(cl.i, "lsparseVector")) {
	i@i[i@x] # "drop0", i.e. FALSE; NAs ok
    } else if (extends(cl.i, "sparseVector")) { ## 'i'sparse, 'd'sparse	 (etc)
	as.integer(i@x[i@i])
    }
    else
        stop("index must be numeric, logical or sparseVector for indexing sparseVectors")
} ## intIv()


setMethod("head", signature(x = "sparseVector"),
	  function(x, n = 6, ...) {
	      stopifnot(length(n) == 1)
	      if(n >= (nx <- x@length)) return(x)
	      if(is.integer(x@i)) n <- as.integer(n) else stopifnot(n == round(n))
              if(n < 0) n <- max(0L, n + nx)
	      x@length <- n
	      if(length(x@i)) {
		  ## now be careful *NOT* to use seq_len(n), as this be efficient for huge n
		  ## As we *know* that '@i' is sorted increasingly: [x@i <= n] <==> [1:kk]
		  x@i <- x@i[ii <- seq_len(which.max(x@i > n) - 1L)]
		  if(.hasSlot(x, "x")) ## has.x: has "x" slot
		      x@x <- x@x[ii]
	      }
	      x
	  })
setMethod("tail", signature(x = "sparseVector"),
	  function(x, n = 6, ...) {
	      stopifnot(length(n) == 1)
	      if(n >= (nx <- x@length)) return(x)
	      if(is.integer(x@i)) n <- as.integer(n) else stopifnot(n == round(n))
	      if(n < 0) n <- max(0L, n + nx)
	      x@length <- n
	      if((N <- length(x@i))) {
		  ## now be careful *NOT* to use seq_len(n) ... (see above)
		  n <- nx-n # and keep indices > n
		  ii <- if(any(G <- x@i > n)) which.max(G):N else FALSE
		  x@i <- x@i[ii] - n
		  if(.hasSlot(x, "x")) ## has.x: has "x" slot
		      x@x <- x@x[ii]
              }
	      x
	  })


setMethod("[", signature(x = "sparseVector", i = "index"),
	  function (x, i, j, ..., drop) {
	      has.x <- .hasSlot(x, "x")## has "x" slot
	      n <- x@length
	      if(extends(cl.i <- getClass(class(i)), "numeric") && any(i < 0)) {
		  if(any(i > 0))
		      stop("you cannot mix negative and positive indices")
		  if(any(z <- i == 0)) i <- i[!z]
		  ## all (i < 0), negative indices:
		  ## want to remain sparse --> *not* using intIv()
		  ##
		  ## TODO: more efficient solution would use C ..
		  i <- unique(sort(-i)) # so we need to drop the 'i's
		  nom <- is.na(m <- match(x@i, i))
		  ## eliminate those non-0 which do match:
		  x@i <- x@i[nom]
		  if(has.x) x@x <- x@x[nom]
		  ## now all x@i "appear in 'i' but must be adjusted for the removals:
		  x@i <- x@i - findInterval(x@i, i)
		  x@length <- n - length(i)
              } else { ## i >= 0  or  non-numeric 'i'
                  ii <- intIv(i, n, cl.i=cl.i)
                  m <- match(x@i, ii, nomatch = 0)
                  sel <- m > 0L
                  x@length <- length(ii)
                  x@i <- m[sel]
		  if(any(iDup <- duplicated(ii))) {
                      i.i <- match(ii[iDup], ii)
                      jm <- lapply(i.i, function(.) which(. == m))
                      if (has.x) sel <- c(which(sel), unlist(jm))
                      x@i <- c(x@i, rep.int(which(iDup), lengths(jm)))
                  }
                  if(doSort <- is.unsorted(x@i)) {
                      io <- order(x@i, method="radix")
                      x@i <- x@i[io]
                  }
                  if (has.x)
                      x@x <- if(doSort) x@x[sel][io] else x@x[sel]
              }
	      x
	  })

setMethod("[", signature(x = "sparseVector", i = "lsparseVector"),
	  function (x, i, j, ..., drop) x[sort.int(i@i[i@x])])
setMethod("[", signature(x = "sparseVector", i = "nsparseVector"),
	  function (x, i, j, ..., drop) x[sort.int(i@i)])

##--- Something else:  Allow    v[ <sparseVector> ] -- exactly similarly:
if(FALSE) { ## R_FIXME: Not working, as internal "[" only dispatches on 1st argument
setMethod("[", signature(x = "atomicVector", i = "lsparseVector"),
	  function (x, i, j, ..., drop) x[sort.int(i@i[i@x])])
setMethod("[", signature(x = "atomicVector", i = "nsparseVector"),
	  function (x, i, j, ..., drop) x[sort.int(i@i)])
}

##' Implement   x[i] <- value

##' @param x  a "sparseVector"
##' @param i  an "index" (integer, logical, ..)
##' @param value

##' @return  a "sparseVector" of the same length as 'x'
## This is much analogous to replTmat in ./Tsparse.R:
replSPvec <- function (x, i, value)
{
    n <- x@length
    ii <- intIv(i, n)
    lenRepl <- length(ii)
    if(!lenRepl) return(x)
    ## else:  lenRepl = length(ii) > 0
    lenV <- length(value)
    if(lenV == 0)
        stop("nothing to replace with")
    ## else: lenV := length(value) > 0
    if(lenRepl %% lenV != 0)
	stop("number of items to replace is not a multiple of replacement length")
    if(anyDuplicated(ii)) { ## multiple *replacement* indices: last one wins
	## TODO: in R 2.6.0 use	 duplicate(*, fromLast=TRUE)
	ir <- lenRepl:1
	keep <- match(ii, ii[ir]) == ir
	ii <- ii[keep]
	lenV <- length(value <- rep(value, length.out = lenRepl)[keep])
	lenRepl <- length(ii)
    }

    has.x <- .hasSlot(x, "x")## has "x" slot
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
	value <- rep(value, length.out = lenRepl)
    ## now:  length(value) == lenRepl > 0

    v0 <- is0(value)
    ## value[1:lenRepl]:  which are structural 0 now, which not?
    v.sp <- inherits(value, "sparseVector")

    if(any(sel)) {
	## indices of non-zero entries -- WRT to subvector
	iN0 <- m[sel] ## == match(x@i[sel], ii)

	## 1a) replace those that are already non-zero with new val.
	vN0 <- !v0[iN0]
	if(any(vN0) && has.x) {
	    vs <- value[iN0[vN0]]
	    x@x[sel][vN0] <- if(v.sp) sp2vec(vs, mode=typeof(x@x)) else vs
	}
	## 1b) replace non-zeros with 0 --> drop entries
	if(any(!vN0)) {
	    i <- which(sel)[!vN0]
	    if(has.x)
		x@x <- x@x[-i]
	    x@i <- x@i[-i]
	}
	iI0 <- if(length(iN0) < lenRepl) seq_len(lenRepl)[-iN0] # else NULL
    } else iI0 <- seq_len(lenRepl)

    if(length(iI0) && any(vN0 <- !v0[iI0])) {
	## 2) add those that were structural 0 (where value != 0)
	ij0 <- iI0[vN0]
	ii <- c(x@i, ii[ij0]) # new x@i, must be sorted:
	iInc <- sort.list(ii)
	x@i <- ii[iInc]
	if(has.x) # new @x, sorted along '@i':
	    x@x <- c(x@x, if(v.sp)
			      sp2vec(value[ij0], mode=typeof(x@x))
			  else value[ij0]
		     )[iInc]
    }
    x
}

setReplaceMethod("[", signature(x = "sparseVector", i = "index", j = "missing",
				value = "replValueSp"),
		 replSPvec)

setReplaceMethod("[", signature(x = "sparseVector",
                                i = "sparseVector", j = "missing",
				value = "replValueSp"),
                 ## BTW, the important case: 'i' a *logical* sparseVector
		 replSPvec)
rm(replSPvec)

## Something else:  Also allow	  x[ <sparseVector> ] <- v  e.g. for atomic x :

if(FALSE) { ## R_FIXME: Not working, as internal "[<-" only dispatches on 1st argument
## Now "the work is done" inside  intIv() :
setReplaceMethod("[", signature(x = "atomicVector",
				i = "sparseVector", j = "missing",
				value = "replValue"),
		 function (x, i, value)
		 callGeneric(x, i = intIv(i, x@length), value=value))
}

## MJ: unused
if(FALSE) {
## a "method" for c(<(sparse)Vector>, <(sparse)Vector>):
## FIXME: This is not exported, nor used (nor documented)
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

## sort.default() does
##		x[order(x, na.last = na.last, decreasing = decreasing)]
## but that uses a *dense* integer order vector
## ==> need direct sort() method for "sparseVector" for mean(*,trim), median(),..
sortSparseV <- function(x, decreasing = FALSE, na.last = NA) {
    if(length(ina <- which(is.na(x)))) {
        if(is.na(na.last)) x <- x[-ina]
    }
    ## TODO
    .NotYetImplemented()
}

##' Uniquify sparceVectors, i.e., bring them in "regularized" from,
##' --- similar in spirit (and action!) as  uniqTsparse(.) for "TsparseMatrix"
##' __FIXME__ better name ??  , then export and document!  __TODO__
uniqSpVec <- function(x) {
    ii <- sort.list(x@i, method = "radix")
    x@i <- x@i[ii]
    x@x <- x@x[ii]
    x
}
} ## MJ

all.equal.sparseV <- function(target, current, ...)
{
    if(!is(target, "sparseVector") || !is(current, "sparseVector")) {
	return(paste0("target is ", data.class(target), ", current is ",
		      data.class(current)))
    }
    lt <- length(target)
    lc <- length(current)
    if(lt != lc) {
	return(paste0("sparseVector", ": lengths (", lt, ", ", lc, ") differ"))
    }

    t.has.x <- .hasSlot(target,  "x")## has "x" slot
    c.has.x <- .hasSlot(current, "x")## has "x" slot
    nz.t <- length(i.t <- target @i)
    nz.c <- length(i.c <- current@i)
    t.x <- if(t.has.x)	target@x else rep.int(TRUE, nz.t)
    c.x <- if(c.has.x) current@x else rep.int(TRUE, nz.c)
    if(nz.t != nz.c || any(i.t != i.c)) { ## "work" if indices are not the same
	i1.c <- setdiff(i.t, i.c)# those in i.t, not yet in i.c
	i1.t <- setdiff(i.c, i.t)
	if((n1t <- length(i1.t))) {
	    target@i <- i.t <- c(i.t, i1.t)
	    t.x <- c(t.x, rep.int(if(t.has.x) 0 else 0L, n1t))
	}
	if((n1c <- length(i1.c))) {
	    current@i <- i.c <- c(i.c, i1.c)
	    c.x <- c(c.x, rep.int(if(c.has.x) 0 else 0L, n1c))
	}
    }
    if(is.unsorted(i.t)) {  ## method="quick" {"radix" not ok for large range}
	ii <- sort.list(i.t, method = "quick", na.last=NA)
	target@i <- i.t <- i.t[ii]
	t.x <- t.x[ii]
    }
    if(is.unsorted(i.c)) {
	ii <- sort.list(i.c, method = "quick", na.last=NA)
	current@i <- i.c <- i.c[ii]
	c.x <- c.x[ii]
    }

    ## Now, we have extended both target and current
    ## *and* have sorted the respective i-slot, the i-slots should match!
    stopifnot(all(i.c == i.t))

    if(is.logical(t.x))
        all.equal.raw(t.x, c.x, ...)
    else
        all.equal.numeric(t.x, c.x, ...)
} ## all.equal.sparseV


## For these, we remain sparse:
setMethod("all.equal", c(target = "sparseVector", current = "sparseVector"),
	  all.equal.sparseV)
setMethod("all.equal", c(target = "sparseVector", current = "sparseMatrix"),
	  function(target, current, ...)
	  all.equal.sparseV(target, as(current, "sparseVector"), ...))
setMethod("all.equal", c(target = "sparseMatrix", current = "sparseVector"),
	  function(target, current, ...)
	  all.equal.sparseV(as(target, "sparseVector"), current, ...))
## For the others, where one is "dense", "go to" dense rather now than later:
setMethod("all.equal", c(target = "ANY", current = "sparseVector"),
	  function(target, current, ...)
	  all.equal(target, as.vector(current), ...))
setMethod("all.equal", c(target = "sparseVector", current = "ANY"),
	  function(target, current, ...)
	  all.equal(as.vector(target), current, ...))


## S3 method for 'c' [but only for dispatch on 1st arg, hence also exported as fn]
c.sparseVector <- function(...) {
    svl <- lapply(list(...), as, Class = "sparseVector")
    ## cls <- unique(unlist(lapply(svl, is)))
    ns <- vapply(svl, slot, 1, "length")
    if((N <- sum(ns)) < .Machine$integer.max) { # some purism ..
	ns <- as.integer(ns)
	N <- as.integer(N)
    }
    narg <- length(ns)
    iss <- lapply(svl, slot, "i")
    ## new 'i' slot:
    ii <- unlist(iss) + rep(cumsum(c(0L, ns[-narg])), lengths(iss))
    ## result must have 'x' slot if we have any
    has.x <- any(have.x <- vapply(svl, .hasSlot, logical(1L), name = "x"))
    if(has.x) {
	cls <- if     (any(vapply(svl, is, NA, "zsparseVector"))) "zsparseVector"
	       else if(any(vapply(svl, is, NA, "dsparseVector"))) "dsparseVector"
	       else if(any(vapply(svl, is, NA, "isparseVector"))) "isparseVector"
	       else "lsparseVector"
	if(!(all.x <- all(have.x)))
	    one <- if     (identical(cls, "lsparseVector")) TRUE
		   else if(identical(cls, "isparseVector")) 1L else 1.
	xx <- unlist(if(all.x) lapply(svl, slot, "x")
		     else lapply(seq_len(narg), function(i) {
			 if(have.x[[i]]) svl[[i]]@x
			 else rep_len(one, length(iss[[i]]))
		     }))
	new(cls, x = xx,     i = ii, length = N)
    } else ## no "x" slot
	new("nsparseVector", i = ii, length = N)
}

### rep(x, ...) -- rep() is primitive with internal default method with these args:
### -----------
### till R 2.3.1, it had  rep.default()  which we use as 'model' here.

repSpV <- function(x, times) {
    ## == rep.int(<sparseVector>, times)"
    times <- as.integer(times)# truncating as rep.default()
    n <- x@length
    has.x <- .hasSlot(x, "x")## has "x" slot
    ## just assign new correct slots:
    if(times <= 1) { ## be quick for {0, 1} times
        if(times < 0) stop("'times >= 0' is required")
        if(times == 0) {
            x@length <- 0L
            x@i <- integer(0)
            if(has.x) x@x <- rep.int(x@x, 0)
        }
        return(x)
    }
    n. <- as.double(n)
    if(n. * times >= .Machine$integer.max)
        n <- n. # so won't have overflow in subsequent multiplys
    x@length <- n * times
    x@i <- rep.int(x@i, times) + n * rep(0:(times-1L), each=length(x@i))
    ## := outer(x@i, 0:(times-1) * n, "+")   but a bit faster
    if(has.x) x@x <- rep.int(x@x, times)
    x
}

setMethod("rep", "sparseVector",
	  function(x, times, length.out, each, ...) {
	      if (length(x) == 0)
		  return(if(missing(length.out)) x else head(x, length.out))
	      if (!missing(each)) {
		  tm <- rep.int(each, length(x))
		  x <- rep(x, tm) # "recursively"
		  if(missing(length.out) && missing(times))
		      return(x)
	      } ## else :
	      if (!missing(length.out)) # takes precedence over times
		  times <- ceiling(length.out/length(x))
	      r <- repSpV(x, times)
	      if (!missing(length.out) && length(r) != length.out) {
		  if(length.out > 0) head(r, length.out) else r[integer(0)]
	      }
	      else r
	  })


### Group Methods (!)
## "Ops" : ["Arith", "Compare", "Logic"]:  ---> in ./Ops.R
##						     -----
## "Summary"  ---> ./Summary.R
##		     ---------
## "Math", "Math2": ./Math.R
##		     -------


##' indices of vector x[] to construct  Toeplitz matrix
##' FIXME: write in C, port to  R('stats' package), and use in stats::toeplitz()
ind4toeplitz <- function(n) {
    A <- matrix(raw(), n, n)
    abs(as.vector(col(A) - row(A))) + 1L
}

.toeplitz.spV <-  function(x, symmetric=TRUE, repr = c("C","T","R"), giveCsparse) {
    ## semantically "identical" to stats::toeplitz
    n <- length(x)
    r <- spV2M(x[ind4toeplitz(n)], n,n, symmetric = symmetric, check = FALSE)
    ##   ^^^^^ returning TsparseMatrix
    if(!missing(giveCsparse)) {
	if(missing(repr)) {
	    repr <- if(giveCsparse) "C" else "T"
	    warning(gettextf("'giveCsparse' has been deprecated; setting 'repr = \"%s\"' for you",
                             repr),
                    domain = NA)
	} else ## !missing(repr)
            Matrix.msg("'giveCsparse' has been deprecated; will use 'repr' instead")
    }
    switch(match.arg(repr), "C" = .T2C(r), "T" = r, "R" = .T2R(r))
}
setMethod("toeplitz", "sparseVector", .toeplitz.spV)
rm(.toeplitz.spV)
