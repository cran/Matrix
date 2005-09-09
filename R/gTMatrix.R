#### "gTMatrix" : Virtual class of general sparse matrices in triplet-format

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
        i0 <- i - 1:1 # transform to 0-indexing
        if(has.dn) dn <- dn[i]
    }
    else if (is(i, "logical")) {
        i0 <- (0:(di[margin]-1:1))[i]
        if(has.dn) dn <- dn[i]
    } else { ## character
        if(!has.dn)
            stop(gettextf(
                          "no 'dimnames[[%d]]': cannot use character indexing"),
                 margin, domain = NA)
        i0 <- match(i, dn, nomatch=0)
        dn <- dn[i0]
        i0 <- i0 - 1:1
    }
    list(m = match(xi, i0, nomatch=0), li = length(i0), dn = dn)
}


## Otherwise have to write methods for all possible combinations of
##  (i , j) \in
##  (numeric, logical, character, missing) x (numeric, log., char., miss.)


## Select rows
setMethod("[", signature(x = "gTMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select rows
              ip <- .ind.prep(x@i, i, 1, dim(x), dimnames(x))
	      x@Dim[1] <- ip$li
              x@Dimnames[1] <- ip$dn
	      sel <- ip$m > 0
	      x@i <- ip$m[sel] - 1:1
	      x@j <- x@j[sel]
	      x@x <- x@x[sel]
	      if (drop && any(x@Dim == 1:1)) drop(as(x,"matrix")) else x
	  })


## Select columns
setMethod("[", signature(x = "gTMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select columns
              ip <- .ind.prep(x@j, j, 2, dim(x), dimnames(x))
	      x@Dim[2] <- ip$li
              x@Dimnames[2] <- ip$dn
	      sel <- ip$m > 0
	      x@i <- x@i[sel]
	      x@j <- ip$m[sel] - 1:1
	      x@x <- x@x[sel]
	      if (drop && any(x@Dim == 1:1)) drop(as(x,"matrix")) else x
	  })


## [.data.frame has : drop = if (missing(i)) TRUE else length(cols) == 1)

setMethod("[", signature(x = "gTMatrix",
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
          x@x <- x@x[sel]
	  if (drop && any(nd == 1)) drop(as(x,"matrix")) else x
      })
