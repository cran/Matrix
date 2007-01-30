#### "TsparseMatrix" : Virtual class of sparse matrices in triplet-format

setAs("TsparseMatrix", "CsparseMatrix",
      ## |-> cholmod_T -> cholmod_C -> chm_sparse_to_SEXP
      ## adjusted for triangular matrices not represented in cholmod
      function(from) .Call(Tsparse_to_Csparse, from, ## ../src/Tsparse.c
                           is(from, "triangularMatrix"))
      )

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
    ## di = dim(x)      { used when i is not character }
    ## dn = dimnames(x) { used when i is character }

    dn <- dn[[margin]]
    has.dn <- is.character(dn)
    if(is(i, "numeric")) {
	storage.mode(i) <- "integer"
        n <- di[margin]
	if(any(i < 0:0)) {
	    if(any(i > 0:0))
		stop("you cannot mix negative and positive indices")
	    i0 <- (0:(n - 1:1))[i]
	} else {
	    if(length(i) && max(i) > n)
		stop("indexing out of range 0:",n)
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
    ## di = dim(x)      { used when i is not character }

    ## difference to .ind.prep(): use 1-indices; no match(xi,..), no dn at end
    dn <- dn[[margin]]
    has.dn <- is.character(dn)
    if(is(i, "numeric")) {
        storage.mode(i) <- "integer"
        n <- di[margin]
	if(any(i < 0:0)) {
	    if(any(i > 0:0))
		stop("you cannot mix negative and positive indices")
	    i0 <- seq_len(n)[i]
	} else	{
	    if(length(i) && max(i) > n)
		stop("indexing out of range 0:",n)
	    if(any(z <- i == 0)) i <- i[!z]
	    i0 <- i
	}
    }
    else if (is(i, "logical")) {
        i0 <- seq_len(di[margin])[i]
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
	      clx <- getClassDef(class(x))
	      has.x <- !extends(clx, "nsparseMatrix")
	      x.sym <- extends(clx, "symmetricMatrix")
	      if(x.sym)
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))
	      ip <- .ind.prep(x@i, i, 1, dim(x), dimnames(x))
	      Di1 <- ip$li
	      drop.it <- drop && (Di1 == 1:1 || x@Dim[2] == 1:1)
	      if(!drop.it && !x.sym && extends(clx, "triangularMatrix"))
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))
	      x@Dim[1] <- Di1
	      if(!is.null(ip$dn)) x@Dimnames[[1]] <- ip$dn
	      sel <- ip$m > 0
	      x@i <- ip$m[sel] - 1:1
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
	      if(x.sym)
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))
	      ip <- .ind.prep(x@j, j, 2, dim(x), dimnames(x))
	      Di2 <- ip$li
	      drop.it <- drop && (x@Dim[1] == 1:1 || Di2 == 1:1)
	      if(!drop.it && !x.sym && extends(clx, "triangularMatrix"))
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))
	      x@Dim[2] <- Di2
	      if(!is.null(ip$dn)) x@Dimnames[[2]] <- ip$dn
	      sel <- ip$m > 0
	      x@i <- x@i[sel]
	      x@j <- ip$m[sel] - 1:1
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
	      isSym <- length(i) == length(j) && all(i == j)
	      ## result is *still* symmetric --> keep symmetry!
	      if(!isSym)
		  ## result no longer symmetric -> to "generalMatrix"
		  x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))
	  }
	  if(isSym) {
	      offD <- x@i != x@j
	      ip1 <- .ind.prep(c(x@i,x@j[offD]), i, 1, di, dn)
	      ip2 <- .ind.prep(c(x@j,x@i[offD]), j, 2, di, dn)
	  } else {
	      ip1 <- .ind.prep(x@i, i, 1, di, dn)
	      ip2 <- .ind.prep(x@j, j, 2, di, dn)
	  }
	  nd <- c(ip1$li, ip2$li)
	  drop.it <- drop && any(nd == 1)
	  if(!drop.it && !isSym && extends(clx, "triangularMatrix"))
	      x <- as(x, paste(.M.kind(x, clx), "gTMatrix", sep=''))
	  x@Dim <- nd
	  x@Dimnames <- list(ip1$dn, ip2$dn)
	  sel <- ip1$m > 0:0  &	 ip2$m > 0:0
	  if(isSym) { # only those corresponding to upper/lower triangle
	      sel <- sel &
	      (if(x@uplo == "U") ip1$m <= ip2$m else ip2$m <= ip1$m)
	  }
	  x@i <- ip1$m[sel] - 1:1
	  x@j <- ip2$m[sel] - 1:1
	  if (has.x)
	      x@x <- c(x@x, if(isSym) x@x[offD])[sel]
	  if (drop.it) drop(as(x,"matrix")) else x
      })


## FIXME: Learn from .TM... below or rather  .M.sub.i.2col(.) in ./Matrix.R
## ------ the following should be much more efficient than the   ./Matrix.R code :
if(FALSE)
## A[ ij ]  where ij is (i,j) 2-column matrix :
setMethod("[", signature(x = "TsparseMatrix",
			 i = "matrix", j = "missing"),# drop="ANY"
	  function (x, i, j, drop)
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
	      ip1 <- .ind.prep(c(x@i,x@j[offD]), i, 1, di, dn)
	      ip2 <- .ind.prep(c(x@j,x@i[offD]), j, 2, di, dn)
	  } else {
	      ip1 <- .ind.prep(x@i, i, 1, di, dn)
	      ip2 <- .ind.prep(x@j, j, 2, di, dn)
	  }

          stop("FIXME: NOT YET FINISHED IMPLEMENTATION")

          ## The M[i_vec, j_vec] had -- we need "its diagonal" :
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


###========= Sub-Assignment aka *Replace*Methods =========================

### FIXME: make this `very fast'  for the very very common case of
### -----   M[i,j] <- v  with   i,j = length-1-numeric;  v= length-1 number
###                            *and* M[i,j] == 0 previously

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

    clx <- class(x)
    clDx <- getClassDef(clx) # extends() , is() etc all use the class definition
    stopifnot(extends(clDx, "TsparseMatrix"))
    ## Tmatrix maybe non-unique, have an entry split into a sum of several ones:
    if(is_duplicatedT(x, nr = di[1]))
	x <- uniqTsparse(x)

    get.ind.sel <- function(ii,ij)
	(match(x@i, ii, nomatch = 0) > 0:0 &
	 match(x@j, ij, nomatch = 0) > 0:0)

    toGeneral <- FALSE
    if((sym.x <- extends(clDx, "symmetricMatrix"))) {
	r.sym <- dind[1] == dind[2] && all(i1 == i2) &&
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
    else if(lenV < lenRepl)
       value <- rep(value, length = lenRepl)

    v0 <- is0(value)
    ## value[1:lenRepl]:  which are structural 0 now, which not?

    if(any(sel)) {
	## the 0-based indices of non-zero entries -- WRT to submatrix
	non0 <- cbind(match(x@i[sel], i1),
		      match(x@j[sel], i2)) - 1:1
	iN0 <- 1:1 + encodeInd(non0, nr = dind[1])

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
	ij0 <- decodeInd(iI0[vN0] - 1:1, nr = dind[1])
	x@i <- c(x@i, i1[ij0[,1] + 1:1])
	x@j <- c(x@j, i2[ij0[,2] + 1:1])
        if(has.x)
            x@x <- c(x@x, value[iI0[vN0]])
    }
    x
}


## A[ ij ] <- value,  where ij is (i,j) 2-column matrix :
## ----------------   ./Matrix.R has a general cheap method
## This one should become as fast as possible:
.TM.repl.i.2col <- function (x, i, value)
{
    nA <- nargs()
    if(nA != 3) stop("nargs() = ", nA, " should never happen; please report.")

    ## else: nA == 3  i.e.,  M [ cbind(ii,jj) ] <- value
    if(is.logical(i)) {
	message(".TM.repl.i.2col(): drop 'matrix' case ...")
	i <- c(i) # drop "matrix"
	return( callNextMethod() )
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
	    ## message("keeping Tsparse matrix *symmetric* in sub-assignment")
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

    i <- i - 1:1 # 0-indexing
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

setReplaceMethod("[", signature(x = "TsparseMatrix", i = "matrix", j = "missing",
				value = "replValue"),
		 .TM.repl.i.2col)


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

## Not needed, have identical "sparseMatrix"
## setMethod("colSums", signature(x = "TsparseMatrix"), .as.dgT.Fun,
## 	  valueClass = "numeric")
## setMethod("colMeans", signature(x = "TsparseMatrix"), .as.dgT.Fun,
## 	  valueClass = "numeric")
##
## Here, "sparseMatrix" uses .as.dgC.Fun:
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

