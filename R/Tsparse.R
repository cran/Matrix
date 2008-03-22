#### "TsparseMatrix" : Virtual class of sparse matrices in triplet-format

## more efficient than going via Csparse:
setAs("matrix", "TsparseMatrix",
      function(from)
      if(is.numeric(from)) mat2dgT(from)
      else if(is.logical(from)) as(Matrix(from, sparse=TRUE), "TsparseMatrix")
      else stop("not-yet-implemented coercion to \"TsparseMatrix\""))

## in ../src/Tsparse.c :  |-> cholmod_T -> cholmod_C -> chm_sparse_to_SEXP
## adjusted for triangular matrices not represented in cholmod
.T.2.C <- function(from) .Call(Tsparse_to_Csparse, from, ##
			       is(from, "triangularMatrix"))

setAs("TsparseMatrix", "CsparseMatrix", .T.2.C)

.T.2.n <- function(from) {
    if(any(is0(from@x))) ## 0 or FALSE -- the following should have drop0Tsp(.)
	from <- as(drop0(from), "TsparseMatrix")
    if(is(from, "triangularMatrix")) # i.e. ?tTMatrix
	new("ntTMatrix", i = from@i, j = from@j,
	    uplo = from@uplo, diag = from@diag,
	    Dim = from@Dim, Dimnames = from@Dimnames)
    else if(is(from, "symmetricMatrix")) # i.e. ?sTMatrix
	new("nsTMatrix", i = from@i, j = from@j, uplo = from@uplo,
	    Dim = from@Dim, Dimnames = from@Dimnames)
    else
	new("ngTMatrix", i = from@i, j = from@j,
	    Dim = from@Dim, Dimnames = from@Dimnames)
}

setAs("TsparseMatrix", "nsparseMatrix", .T.2.n)
setAs("TsparseMatrix", "nMatrix", .T.2.n)

.T.2.l <- function(from) {
    cld <- getClassDef(class(from))
    xx <- if(extends(cld, "nMatrix"))
	rep.int(TRUE, length(from@i)) else as.logical(from@x)
    if(extends(cld, "triangularMatrix")) # i.e. ?tTMatrix
	new("ltTMatrix", i = from@i, j = from@j, x = xx,
	    uplo = from@uplo, diag = from@diag,
	    Dim = from@Dim, Dimnames = from@Dimnames)
    else if(extends(cld, "symmetricMatrix")) # i.e. ?sTMatrix
	new("lsTMatrix", i = from@i, j = from@j, x = xx, uplo = from@uplo,
	    Dim = from@Dim, Dimnames = from@Dimnames)
    else
	new("lgTMatrix", i = from@i, j = from@j, x = xx,
	    Dim = from@Dim, Dimnames = from@Dimnames)
}

setAs("TsparseMatrix", "lsparseMatrix", .T.2.l)
setAs("TsparseMatrix", "lMatrix", .T.2.l)



## Special cases   ("d", "l", "n")  %o%  ("g", "s", "t") :
## used e.g. in triu()

setAs("dgTMatrix", "dgCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("dsTMatrix", "dsCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("dtTMatrix", "dtCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, TRUE))


setAs("lgTMatrix", "lgCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("lsTMatrix", "lsCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("ltTMatrix", "ltCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, TRUE))


setAs("ngTMatrix", "ngCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("nsTMatrix", "nsCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, FALSE))

setAs("ntTMatrix", "ntCMatrix",
      function(from) .Call(Tsparse_to_Csparse, from, TRUE))

### "[" :
### -----

## Test for numeric/logical/character
## method-*internally* ; this is not strictly OO, but allows to use
## the following utility and hence much more compact code.

## Otherwise have to write methods for all possible combinations of
##  (i , j) \in
##  (numeric, logical, character, missing) x (numeric, log., char., miss.)


intI <- function(i, n, dn, give.dn = TRUE)
{
    ## Purpose: translate numeric | logical | character index
    ##		into 0-based integer
    ## ----------------------------------------------------------------------
    ## Arguments: i: index vector (numeric | logical | character)
    ##		  n: array extent                   { ==  dim(.) [margin] }
    ##		 dn: character col/rownames or NULL { == dimnames(.)[[margin]] }
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 23 Apr 2007

    has.dn <- is.character(dn)
    DN <- has.dn && give.dn
    if(is(i, "numeric")) {
	storage.mode(i) <- "integer"
	if(any(i < 0L)) {
	    if(any(i > 0L))
		stop("you cannot mix negative and positive indices")
	    i0 <- (0:(n - 1L))[i]
	} else {
	    if(length(i) && max(i) > n)
		stop("indexing out of range 0:",n)
	    if(any(z <- i == 0)) i <- i[!z]
	    i0 <- i - 1L		# transform to 0-indexing
	}
	if(DN) dn <- dn[i]
    }
    else if (is(i, "logical")) {
	i0 <- (0:(n - 1L))[i]
	if(DN) dn <- dn[i]
    } else { ## character
	if(!has.dn)
	    stop("no 'dimnames[[.]]': cannot use character indexing")
	i0 <- match(i, dn)
	if(any(is.na(i0))) stop("invalid character indexing")
	if(DN) dn <- dn[i0]
	i0 <- i0 - 1L
    }
    if(!give.dn) i0 else list(i0 = i0, dn = dn)
}

.ind.prep <- function(xi, intIlist, iDup = duplicated(i0), anyDup = any(iDup))
{
    ## Purpose: do the ``common things'' for "*gTMatrix" indexing for 1 dim.
    ##          and return match(.,.) + li = length of corresponding dimension
    ##
    ## xi = "x@i" ; intIlist = intI(i, dim(x)[margin], ....)

    i0 <- intIlist$i0
    stopifnot(is.numeric(i0))# cheap fast check (i0 may have length 0 !)

    m <- match(xi, i0, nomatch=0)
    if(anyDup) { # assuming   anyDup <- any(iDup <- duplicated(i0))
        ## i0i: where in (non-duplicated) i0 are the duplicated ones
	i0i <- match(i0[iDup], i0)
        i.x <- which(iDup) - 1L
	jm <- lapply(i0i, function(.) which(. == m))
    }

    c(list(m = m, li = length(i0),
	   i0 = i0, anyDup = anyDup, dn = intIlist$dn),
      ## actually,  iDup  is rarely needed in calling code
      if(anyDup) list(iDup = iDup, i0i = i0i, i.x = i.x,
                      jm = unlist(jm), i.xtra = rep.int(i.x, sapply(jm, length))))
}


.ind.prep2 <- function(i, margin, di, dn)
{    ## Purpose: do the ``common things'' for "*gTMatrix" sub-assignment
    ##		for 1 dimension, 'margin' ,
    ##          and return match(.,.) + li = length of corresponding dimension
    ##
    ## i is "index"; margin in {1,2};
    ## di = dim(x)      { used when i is not character }

    ## difference to .ind.prep(): use 1-indices; no match(xi,..), no dn at end

    intI(i, n = di[margin], dn = dn[[margin]], give.dn = FALSE)
}


## Select rows
setMethod("[", signature(x = "TsparseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select rows
	      clx <- getClassDef(class(x))
	      has.x <- !extends(clx, "nsparseMatrix")
	      x.sym <- extends(clx, "symmetricMatrix")
              x.tri <- extends(clx, "triangularMatrix")
	      gDo <- (x.sym || (x.tri && x@diag == "U"))
	      if(gDo)
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))

	      ip <- .ind.prep(x@i, intI(i, n = dim(x)[1], dimnames(x)[[1]]))
	      Di1 <- ip$li
	      drop.it <- drop && (Di1 == 1L || x@Dim[2] == 1L)
	      if(x.tri && !drop.it && !gDo) # triangular, result not
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))

	      x@Dim[1] <- Di1
	      if(!is.null(ip$dn)) x@Dimnames[[1]] <- ip$dn
	      sel <- ip$m > 0L
	      x@i <- ip$m[sel] - 1L
              if(ip$anyDup) { ## duplicated rows selected: extend sel
                  sel <- c(which(sel), ip$jm)
                  x@i <- c(x@i, ip$i.xtra)
              }
	      x@j <- x@j[sel]
	      if (has.x) x@x <- x@x[sel]
	      if (drop.it) drop(as(x,"matrix")) else x
	  })

## Select columns
setMethod("[", signature(x = "TsparseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select columns
	      clx <- getClassDef(class(x))
	      has.x <- !extends(clx, "nsparseMatrix")
	      x.sym <- extends(clx, "symmetricMatrix")
              x.tri <- extends(clx, "triangularMatrix")
	      gDo <- (x.sym || (x.tri && x@diag == "U"))
	      if(gDo)
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))

	      ip <- .ind.prep(x@j, intI(j, n = dim(x)[2], dimnames(x)[[2]]))
	      Di2 <- ip$li
	      drop.it <- drop && (x@Dim[1] == 1L || Di2 == 1L)
	      if(x.tri && !drop.it && !gDo) # triangular, result not
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))

	      x@Dim[2] <- Di2
	      if(!is.null(ip$dn)) x@Dimnames[[2]] <- ip$dn
	      sel <- ip$m > 0L
	      x@j <- ip$m[sel] - 1L
	      if(ip$anyDup) { ## duplicated columns selected: extend sel
		  sel <- c(which(sel), ip$jm)
		  x@j <- c(x@j, ip$i.xtra)
	      }
	      x@i <- x@i[sel]
	      if (has.x) x@x <- x@x[sel]
	      if (drop.it) drop(as(x,"matrix")) else x
	  })


## [.data.frame has : drop = if (missing(i)) TRUE else length(cols) == 1)


setMethod("[", signature(x = "TsparseMatrix",
			 i = "index", j = "index", drop = "logical"),
	  function (x, i, j, ..., drop)
      {
	  ## (i,j, drop) all specified
	  di <- dim(x)
	  dn <- dimnames(x)
	  clx <- getClassDef(class(x))
	  has.x <- !extends(clx, "nsparseMatrix")
	  isSym <- extends(clx, "symmetricMatrix")

	  if(isSym) {
	      isSym <- length(i) == length(j) && mode(i) == mode(j) && all(i == j)
	      ## result will *still* be symmetric --> keep symmetry!
	      if(!isSym)
		  ## result no longer symmetric -> to "generalMatrix"
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))
	  } else if(extends(clx, "triangularMatrix") && x@diag == "U") {
	      x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))
	  }
	  if(isSym) { ## has only stored "half" of the indices,
	      ## OTOH,	i === j, so only need one intI() call
	      ip1 <- intI(i, n=di[1], dn[[1]]) # -> (i0, dn_1)
	      anyDup <- any(iDup <- duplicated(ip1$i0))
	      ip1 <- .ind.prep(x@i, ip1, iDup=iDup, anyDup=anyDup)
	      ip2 <- {
		  if(anyDup) list(m = match(x@j, ip1$i0, nomatch=0),
				  li = ip1$li)
		  else	     .ind.prep(x@j, ip1, iDup=iDup, anyDup=anyDup)
	      }
	      if(!is.null(dn[[2]]))	# fix result colnames
		  ip2$dn <- dn[[2]][ip1$i0 + 1L]

	  } else {
	      ip1 <- .ind.prep(x@i, intI(i, n = di[1], dn= dn[[1]]))
	      ij <- intI(j, n = di[2], dn= dn[[2]])
	      ip2 <- .ind.prep(x@j, ij)
	  }
	  nd <- c(ip1$li, ip2$li)
	  x@Dim <- nd
	  x@Dimnames <- list(ip1$dn, ip2$dn)

	  if(isSym) {
	      sel <- ip1$m  &  ip2$m
	      ii <- ip1$m[sel] - 1L
	      jj <- ip2$m[sel] - 1L
	      if(anyDup) { ## careful algorithm --- TODO: in C
		  sel <- which(sel)
		  ## keep non-duplicated and "increment" for duplicated ones
		  ij <- pmin(ii, jj)
		  jj <- pmax(ii, jj) ; ii <- ij
		  ## length(ii) == length(jj) == length(sel)  _and_  ii <= jj
		  ix <- ip1$i.x
		  for(k in seq_along(ix)) {
		      ## "recursively" add 1 row+column corresp. iDup[k]
		      i0 <- ip1$i0i[k] -1L # the (0-ind)column we want to repl
		      i.x <- ix[k]	# < i0
		      e1 <- ii == i0
		      e2 <- jj == i0
		      j1 <- jj[e1]	# >= ii[e1] == i0
		      j2 <- ii[e2]	# <= jj[e2] == i0 < i.x
		      ## now the "diagonal special":
		      if(any(e1 & e2)) {
			  ## (e1 & e2)[m] = TRUE  <==>	ii[m] == jj[m] == i0
			  isD <- e1[e2] # logical of same length as j2
			  stopifnot(sum(isD) == 1)
			  j2[isD] <- i.x # instead of i0
		      }
		      l1x <- j1 < i.x	# & j12 <- j1[l1x]
		      j11 <- j1[!l1x]	# those >= i.x
		      s1 <- sel[e1]
		      ii  <- c(ii,  j2, j1[l1x],  rep.int(i.x, length(j11)))
		      jj  <- c(jj,  rep.int(ix[k], length(j2)+sum(l1x)), j11)
		      sel <- c(sel, sel[e2], s1[l1x], s1[!l1x])
		      stopifnot(ii <= jj, length(sel) == length(ii))
		  }
		  if(x@uplo == "U") { ## i <= j : upper triangle
		      x@i <- ii
		      x@j <- jj
		  } else { ## i >= j : lower left triangle
		      x@i <- jj
		      x@j <- ii
		  }
	      }
	      else { ## not any Dup

		  if(x@uplo == "U") { ## i <= j : upper triangle
		      x@i <- pmin(ii, jj)
		      x@j <- pmax(ii, jj)
		  } else { ## i >= j : lower left triangle
		      x@i <- pmax(ii, jj)
		      x@j <- pmin(ii, jj)
		  }
	      }

	  }
	  else if(!ip1$anyDup && !ip2$anyDup) {
	      ## "normal case":	 no duplicated indices (and not symmetric)

	      sel <- ip1$m  &  ip2$m
	      x@i <- ip1$m[sel] - 1L
	      x@j <- ip2$m[sel] - 1L
	  }
	  else { ## not Sym   &&  (ip1$anyDup || ip2$anyDup) :
	      ## duplicated rows or columns -- currently the cheap solution,
	      ## Basically  implement  X[i,j] as  X[i,] [,j] :
	      ## FIXME: we are recomputing ip2 here

	      ## - i - ------------------------------
	      sel <- ip1$m > 0L
	      x@i <- ip1$m[sel] - 1L
	      if(ip1$anyDup) { ## duplicated rows selected: extend sel
		  sel <- c(which(sel), ip1$jm)
		  x@i <- c(x@i, ip1$i.xtra)
	      }
	      x@j <- x@j[sel]
	      if (has.x) x@x <- x@x[sel]

	      ## - j - ------------------------------
	      ## ip2 <- .ind.prep(x@j, intI(j, n = di[2], dn = dn[[2]]))
	      ## FIXME can we do better: current x@j is original x@j[sel]
	      ip2 <- .ind.prep(x@j, ij)
	      sel <- ip2$m > 0L
	      x@j <- ip2$m[sel] - 1L
	      if(ip2$anyDup) { ## duplicated columns selected: extend sel
		  sel <- c(which(sel), ip2$jm)
		  x@j <- c(x@j, ip2$i.xtra)
	      }
	      x@i <- x@i[sel]
	  }

	  if (has.x)
	      x@x <- x@x[sel]

	  if (drop && any(nd == 1)) drop(as(x, "matrix")) else x
      })

## FIXME: Learn from .TM... below or rather  .M.sub.i.2col(.) in ./Matrix.R
## ------ the following should be much more efficient than the   ./Matrix.R code :
if(FALSE)
## A[ ij ]  where ij is (i,j) 2-column matrix :
setMethod("[", signature(x = "TsparseMatrix",
			 i = "matrix", j = "missing"),# drop="ANY"
	  function (x, i, j, ..., drop)
      {
	  di <- dim(x)
	  dn <- dimnames(x)
          ## TODO check  i (= 2-column matrix of indices) ---
          ##      as in  .M.sub.i.2col() in ./Matrix.R
          j <- i[,2]
          i <- i[,1]
	  if(is(x, "symmetricMatrix")) {
	      isSym <- all(i == j)
	      if(!isSym)
		  x <- as(x, paste(.M.kind(x), "gTMatrix", sep=''))
	  } else isSym <- FALSE

	  if(isSym) {
	      offD <- x@i != x@j
	      ip1 <- .ind.prep(c(x@i,x@j[offD]), intI(i, n= di[1], dn=dn[[1]]))
	      ip2 <- .ind.prep(c(x@j,x@i[offD]), intI(j, n= di[2], dn=dn[[2]]))
	  } else {
	      ip1 <- .ind.prep(x@i, intI(i, n = di[1], dn = dn[[1]]))
	      ip2 <- .ind.prep(x@j, intI(j, n = di[2], dn = dn[[2]]))
	  }

          stop("FIXME: NOT YET FINISHED IMPLEMENTATION")

          ## The M[i_vec, j_vec] had -- we need "its diagonal" :
          sel <- ip1$m  &  ip2$m
	  if(isSym) { # only those corresponding to upper/lower triangle
	      sel <- sel &
	      (if(x@uplo == "U") ip1$m <= ip2$m else ip2$m <= ip1$m)
	  }
	  x@i <- ip1$m[sel] - 1L
	  x@j <- ip2$m[sel] - 1L
	  if (!is(x, "nsparseMatrix"))
	      x@x <- c(x@x, if(isSym) x@x[offD])[sel]
	  if (drop && any(nd == 1)) drop(as(x,"matrix")) else x

      })


###========= Sub-Assignment aka *Replace*Methods =========================

### FIXME: make this `very fast'  for the very very common case of
### -----   M[i,j] <- v  with   i,j = length-1-numeric;  v= length-1 number
###                            *and* M[i,j] == 0 previously

## ---------     ----- FIXME(2): keep in sync with replCmat() in ./Csparse.R
## workhorse for "[<-" :
replTmat <- function (x, i, j, ..., value)
{
## NOTE:  need '...', i.e., exact signature such that setMethod()
##	  does not use .local() such that nargs() will work correctly:
    di <- dim(x)
    dn <- dimnames(x)
    iMi <- missing(i)
    jMi <- missing(j)
    na <- nargs()
    if(na == 3) { ## "vector (or 2-col) indexing"  M[i] <- v
	if(getOption("verbose"))
	   message(sprintf(paste("diagnosing replTmat(x,i,j,v): nargs()= %d;",
				 "missing(i,j)= (%d,%d)."), na, iMi,jMi))
	if(iMi) stop("internal bug: missing 'i' in replTmat(): please report")
	if(is.character(i))
	    stop("[ <character> ] indexing not allowed: forgot a \",\" ?")
	if(is.matrix(i))
	    stop("internal bug: matrix 'i' in replTmat(): please report")
	## Now: have  M[i] <- v	 with vector logical or "integer" i :
	## Tmatrix maybe non-unique, have an entry split into a sum of several ones:

	if(!is(x,"generalMatrix")) {
	    cl <- class(x)
	    x <- as(x, paste(.M.kind(x), "gTMatrix", sep=''))
            if(getOption("verbose")) ## hmm, only sub-optimal in *some* cases
	    message("'sub-optimal sparse 'x[i] <- v' assignment: Coercing class ",
		    cl," to ",class(x))
	}
	nr <- di[1]
	x.i <- encodeInd2(x@i, x@j, nr)
	if(any(duplicated(x.i))) { ## == if(is_duplicatedT(x, nr = di[1]))
	    x <- uniqTsparse(x)
	    x.i <- encodeInd2(x@i, x@j, nr)
	}

	if(is.logical(i)) { # full-size logical indexing
	    n <- prod(di)
	    if(n) {
		if(length(i) < n) i <- rep(i, length.out = n)
		i <- (0:(n-1))[i] # -> 0-based index vector as well {maybe LARGE!}
	    } else i <- integer(0)
	} else i <- as.integer(i) - 1L ## 0-based indices

        clx <- class(x)
        clDx <- getClassDef(clx) # extends(), is() etc all use the class definition
        has.x <- "x" %in% slotNames(clDx) # === slotNames(x)

	## now have 0-based indices   x.i (entries) and	 i (new entries)

	## the simplest case:
	if(all0(value)) { ## just drop the non-zero entries
	    sel <- is.na(match(x.i, i))
	    if(any(!sel)) { ## non-zero there
		x@i <- x@i[sel]
		x@j <- x@j[sel]
		if(has.x)
		    x@x <- x@x[sel]
	    }
	    return(x)
	}

	m <- length(i)
	if(length(value) != m) { ## use recycling rules
	    if(m %% length(value) != 0)
		warning("number of items to replace is not a multiple of replacement length")
	    value <- rep(value, length.out = m)
	}

	## matching existing non-zeros and new entries
	isE <- !is.na(mi <- match(i, x.i)) ## use  which(isE) , mi[isE]
	## 1) Change the matching non-zero entries
	if(has.x)
	    x@x[mi[isE]] <- value[isE]
	## 2) add the new non-zero entries
	i <- i[!isE]
	x@i <- c(x@i, i %%  nr)
	x@j <- c(x@j, i %/% nr)
	if(has.x)
	    x@x <- c(x@x, value[!isE])
	return(x)
    }
    ## nargs() == 4 :

    i1 <- if(iMi) 0:(di[1] - 1L) else .ind.prep2(i, 1, di, dn)
    i2 <- if(jMi) 0:(di[2] - 1L) else .ind.prep2(j, 2, di, dn)
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
    if(!iMi && any(duplicated(i1))) {
        ## a bit faster than  keep <- !rev(duplicated(rev(i1))) :
        ir <- dind[1]:1 ; keep <- match(i1, i1[ir]) == ir
        i1 <- i1[keep]
        lenV <- length(value <- rep(value, length = lenRepl)[keep])
        dind[1] <- length(i1)
        lenRepl <- dind[1] * dind[2]
    }

    if(!jMi && any(duplicated(i2))) {
        ## a bit faster than  keep <- !rev(duplicated(rev(i2))) :
        ir <- dind[2]:1 ; keep <- match(i2, i2[ir]) == ir
        i2 <- i2[keep]
        lenV <- length(value <- rep(value, length = lenRepl)[keep])
        dind[2] <- length(i2)
        lenRepl <- dind[1] * dind[2]
    }
    clx <- class(x)
    clDx <- getClassDef(clx) # extends() , is() etc all use the class definition
    stopifnot(extends(clDx, "TsparseMatrix"))
    ## Tmatrix maybe non-unique, have an entry split into a sum of several ones:
    if(is_duplicatedT(x, nr = di[1]))
	x <- uniqTsparse(x)

    toGeneral <- FALSE
    if((sym.x <- extends(clDx, "symmetricMatrix"))) {
	r.sym <- (dind[1] == dind[2]) && all(i1 == i2) &&
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
	else toGeneral <- TRUE
    }
    else if((tri.x <- extends(clDx, "triangularMatrix"))) {
        xU <- x@uplo == "U"
	r.tri <- ((any(dind == 1) || dind[1] == dind[2]) &&
		  all(if(xU) i1 <= i2 else i2 <= i1))
	if(r.tri) { ## result is *still* triangular
            if(any(i1 == i2)) # diagonal will be changed
                x <- diagU2N(x) # keeps class (!)
	}
	else toGeneral <- TRUE
    }
    if(toGeneral) { # go to "generalMatrix" and continue
        x <- as(x, paste(.M.kind(x), "gTMatrix", sep=''))
        clDx <- getClassDef(clx <- class(x))
    }

    get.ind.sel <- function(ii,ij)
	(match(x@i, ii, nomatch = 0) & match(x@j, ij, nomatch = 0))

    ## sel[k] := TRUE iff k-th non-zero entry (typically x@x[k]) is to be replaced
    sel <- get.ind.sel(i1,i2)
    has.x <- "x" %in% slotNames(clDx) # === slotNames(x)

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
    ## now have  lenV <= lenRepl

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

    if(sym.x && r.sym) # value already adjusted, see above
       lenRepl <- length(value) # shorter (since only "triangle")
    else if(lenV < lenRepl)
       value <- rep(value, length = lenRepl)

    ## now:  length(value) == lenRepl

    v0 <- is0(value)
    ## value[1:lenRepl]:  which are structural 0 now, which not?

    if(any(sel)) {
	## the 0-based indices of non-zero entries -- WRT to submatrix
	non0 <- cbind(match(x@i[sel], i1),
		      match(x@j[sel], i2)) - 1L
	iN0 <- 1L + encodeInd(non0, nr = dind[1])

	## 1a) replace those that are already non-zero with non-0 values
	vN0 <- !v0[iN0]
	if(any(vN0) && has.x)
	    x@x[sel][vN0] <- value[iN0[vN0]]

	## 1b) replace non-zeros with 0 --> drop entries
	if(any(!vN0)) {
	    ii <- which(sel)[!vN0]
	    if(has.x)
		x@x <- x@x[-ii]
	    x@i <- x@i[-ii]
	    x@j <- x@j[-ii]
	}
	iI0 <- if(length(iN0) < lenRepl)
	    seq_len(lenRepl)[-iN0] # == complementInd(non0, dind)
    } else iI0 <- seq_len(lenRepl)

    if(length(iI0) && any(vN0 <- !v0[iI0])) {
	## 2) add those that were structural 0 (where value != 0)
	ij0 <- decodeInd(iI0[vN0] - 1L, nr = dind[1])
	x@i <- c(x@i, i1[ij0[,1] + 1L])
	x@j <- c(x@j, i2[ij0[,2] + 1L])
        if(has.x)
            x@x <- c(x@x, value[iI0[vN0]])
    }
    x
} ## end{replTmat}

## A[ ij ] <- value,  where ij is (i,j) 2-column matrix :
## ----------------   ./Matrix.R has a general cheap method
## This one should become as fast as possible -- is also used from Csparse.R --
.TM.repl.i.2col <- function (x, i, j, ..., value)
{
    nA <- nargs()
    if(nA != 3) stop("nargs() = ", nA, " should never happen; please report.")

    ## else: nA == 3  i.e.,  M [ cbind(ii,jj) ] <- value or M [ Lmat ] <- value
    if(is.logical(i)) {
	if(getOption("verbose"))
	    message(".TM.repl.i.2col(): drop 'matrix' case ...")
	## c(i) : drop "matrix" to logical vector
	x[c(i)] <- value
	return(x)
    } else if(!is.numeric(i) || ncol(i) != 2)
	stop("such indexing must be by logical or 2-column numeric matrix")
    if(!is.integer(i)) storage.mode(i) <- "integer"
    if(any(i < 0))
	stop("negative values are not allowed in a matrix subscript")
    if(any(is.na(i)))
	stop("NAs are not allowed in subscripted assignments")
    if(any(i0 <- (i == 0))) # remove them
	i <- i[ - which(i0, arr.ind = TRUE)[,"row"], ]
    if(length(attributes(i)) > 1) # more than just 'dim'; simplify: will use identical
	attributes(i) <- list(dim = dim(i))
    ## now have integer i >= 1
    m <- nrow(i)
    if(m == 0)
	return(x)
    if(length(value) == 0)
	stop("nothing to replace with")
    ## mod.x <- .type.kind[.M.kind(x)]
    if(length(value) != m) { ## use recycling rules
	if(m %% length(value) != 0)
	    warning("number of items to replace is not a multiple of replacement length")
	value <- rep(value, length = m)
    }
    clx <- class(x)
    clDx <- getClassDef(clx) # extends() , is() etc all use the class definition
    stopifnot(extends(clDx, "TsparseMatrix"))

    di <- dim(x)
    nr <- di[1]
    nc <- di[2]
    i1 <- i[,1]
    i2 <- i[,2]
    if(any(i1 > nr)) stop("row indices must be <= nrow(.) which is ", nr)
    if(any(i2 > nc)) stop("column indices must be <= ncol(.) which is ", nc)

    ## Tmatrix maybe non-unique, have an entry split into a sum of several ones:
    if(is_duplicatedT(x, nr = nr))
	x <- uniqTsparse(x)

    toGeneral <- FALSE
    if((sym.x <- extends(clDx, "symmetricMatrix"))) {
	## Tests to see if the assignments are symmetric as well
	r.sym <- all(i1 == i2)
	if(!r.sym) { # do have *some* Lower or Upper entries
	    iL <- i1 > i2
	    iU <- i1 < i2
	    r.sym <- sum(iL) == sum(iU) # same number
	    if(r.sym) {
		iLord <- order(i1[iL], i2[iL])
		iUord <- order(i2[iU], i1[iU]) # row <-> col. !
		r.sym <- {
		    identical(i[iL,    ][iLord,],
			      i[iU, 2:1][iUord,]) &&
		    all(value[iL][iLord] ==
			value[iU][iUord])
		}
	    }
	}
	if(r.sym) { ## result is *still* symmetric --> keep symmetry!
	    if(getOption("verbose"))
		message("keeping Tsparse matrix *symmetric* in sub-assignment")
	    ## now consider only those indices above / below diagonal:
	    xU <- x@uplo == "U"
	    useI <- if(xU) i1 <= i2 else i2 <= i1
	    i1 <- i1[useI]
	    i2 <- i2[useI]
	    value <- value[useI]
	}
	else toGeneral <- TRUE
    }
    else if((tri.x <- extends(clDx, "triangularMatrix"))) {
	xU <- x@uplo == "U"
	r.tri <- all(if(xU) i1 <= i2 else i2 <= i1)
	if(r.tri) { ## result is *still* triangular
	    if(any(i1 == i2)) # diagonal will be changed
		x <- diagU2N(x) # keeps class (!)
	}
	else toGeneral <- TRUE
    }
    if(toGeneral) { # go to "generalMatrix" and continue
	x <- as(x, paste(.M.kind(x), "gTMatrix", sep=''))
	clDx <- getClassDef(clx <- class(x))
    }

    i <- i - 1L # 0-indexing
    ii.v <- encodeInd (i, nr)
    if(any(d <- duplicated(rev(ii.v)))) { # reverse: "last" duplicated FALSE
	warning("duplicate ij-entries in 'Matrix[ ij ] <- value'; using last")
	nd <- !rev(d)
	## i  <- i    [nd, , drop=FALSE]
	ii.v  <- ii.v [nd]
	value <- value[nd]
    }
    ii.x <- encodeInd2(x@i, x@j, nr)
    m1 <- match(ii.v, ii.x)
    i.repl <- !is.na(m1) # those that need to be *replaced*

    if(isN <- extends(clDx, "nMatrix")) { ## no 'x' slot
	isN <- all(value %in% c(FALSE, TRUE)) # will result remain  "nMatrix" ?
	if(!isN)
	    x <- as(x, paste(if(extends(clDx, "lMatrix")) "l" else "d",
			     .sparse.prefixes[.M.shape(x)], "TMatrix", sep=''))
    }
    has.x <- !isN ## isN  <===> "remains pattern matrix" <===> has no 'x' slot

    if(any(i.repl)) { ## some to replace at matching (@i, @j)
	if(has.x)
	    x@x[m1[i.repl]] <- value[i.repl]
	else { # nMatrix ; eliminate entries that are set to FALSE; keep others
	    if(any(isF <- !value[i.repl]))  {
		ii <- m1[i.repl][isF]
		x@i <- x@i[ -ii]
		x@j <- x@j[ -ii]
	    }
	}
    }
    if(!all(i.repl)) { ## some new entries
	i.j <- decodeInd(ii.v[!i.repl], nr)
	x@i <- c(x@i, i.j[,1])
	x@j <- c(x@j, i.j[,2])
	if(has.x)
	    x@x <- c(x@x, value[!i.repl])
    }

    x
} ## end{.TM.repl.i.2col}

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "index", j = "missing",
				value = "replValue"),
		 replTmat)

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "missing", j = "index",
				value = "replValue"),
		 replTmat)

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "index", j = "index",
				value = "replValue"),
		 replTmat)

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "matrix", j = "missing",
				value = "replValue"),
		 .TM.repl.i.2col)


setMethod("crossprod", signature(x = "TsparseMatrix", y = "missing"),
	  function(x, y = NULL) {
              if (is(x, "symmetricMatrix")) {
                  x <- .T.2.C(x)
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
	  function(x, y = NULL) crossprod(.T.2.C(x), y))

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "ANY"),
	  function(x, y = NULL) tcrossprod(.T.2.C(x), y))

setMethod("%*%", signature(x = "TsparseMatrix", y = "ANY"),
	  function(x, y) .T.2.C(x) %*% y)

setMethod("%*%", signature(x = "ANY", y = "TsparseMatrix"),
	  function(x, y) x %*% .T.2.C(y))

## Not yet.  Don't have methods for y = "CsparseMatrix" and general x
#setMethod("%*%", signature(x = "ANY", y = "TsparseMatrix"),
#          function(x, y) callGeneric(x, as(y, "CsparseMatrix")))

setMethod("solve", signature(a = "TsparseMatrix", b = "ANY"),
	  function(a, b) solve(as(a, "CsparseMatrix"), b))
setMethod("solve", signature(a = "TsparseMatrix", b = "missing"),
	  function(a, b) solve(as(a, "CsparseMatrix")))


## Want tril(), triu(), band() --- just as "indexing" ---
## return a "close" class:
setMethod("tril", "TsparseMatrix",
	  function(x, k = 0, ...)
	  as(tril(.T.2.C(x), k = k, ...), "TsparseMatrix"))
setMethod("triu", "TsparseMatrix",
	  function(x, k = 0, ...)
	  as(triu(.T.2.C(x), k = k, ...), "TsparseMatrix"))
setMethod("band", "TsparseMatrix",
	  function(x, k1, k2, ...)
	  as(band(.T.2.C(x), k1 = k1, k2 = k2, ...), "TsparseMatrix"))


## For the "general" T ones (triangular & symmetric have special methods):
setMethod("t", signature(x = "TsparseMatrix"),
	  function(x) {
	      r <- new(class(x))
	      r@i <- x@j
	      r@j <- x@i
	      if(any("x" == slotNames(x)))
		  r@x <- x@x
	      r@Dim <- x@Dim[2:1]
	      r@Dimnames <- x@Dimnames[2:1]
	      r
      })


setMethod("writeMM", "TsparseMatrix",
	  function(obj, file, ...)
          .Call(Csparse_MatrixMarket, as(obj, "CsparseMatrix"),
                as.character(file)))
