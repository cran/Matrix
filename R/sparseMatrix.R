### Define Methods that can be inherited for all subclasses

### Idea: Coercion between *VIRTUAL* classes -- as() chooses "closest" classes
### ----  should also work e.g. for  dense-triangular --> sparse-triangular !

##-> see als ./dMatrix.R, ./ddenseMatrix.R  and  ./lMatrix.R

setAs("ANY", "sparseMatrix", function(from) as(from, "CsparseMatrix"))


## "graph" coercions -- this needs the graph package which is currently
##  -----               *not* required on purpose
## Note: 'undirected' graph <==> 'symmetric' matrix

## Add some utils that may no longer be needed in future versions of the 'graph' package
graph.has.weights <- function(g) "weight" %in% names(edgeDataDefaults(g))

graph.wgtMatrix <- function(g)
{
    ## Purpose: work around "graph" package's  as(g, "matrix") bug
    ## ----------------------------------------------------------------------
    ## Arguments: g: an object inheriting from (S4) class "graph"
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, based on Seth Falcon's code;  Date: 12 May 2006

    ## MM: another buglet for the case of  "no edges":
    if(numEdges(g) == 0) {
      p <- length(nd <- nodes(g))
      return( matrix(0, p,p, dimnames = list(nd, nd)) )
    }
    ## Usual case, when there are edges:
    has.w <- "weight" %in% names(edgeDataDefaults(g))
    if(has.w) {
        w <- unlist(edgeData(g, attr = "weight"))
        has.w <- any(w != 1)
    } ## now 'has.w' is TRUE  iff  there are weights != 1
    m <- as(g, "matrix")
    ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
    if(has.w) { ## fix it if needed
        tm <- t(m)
        tm[tm != 0] <- w
        t(tm)
    }
    else m
}


setAs("graphAM", "sparseMatrix",
      function(from) {
	  symm <- edgemode(from) == "undirected" && isSymmetric(from@adjMat)
	  ## This is only ok if there are no weights...
	  if(graph.has.weights(from)) {
	      as(graph.wgtMatrix(from),
		 if(symm) "dsTMatrix" else "dgTMatrix")
	  }
	  else { ## no weights: 0/1 matrix -> logical
	      as(as(from, "matrix"),
		 if(symm) "nsTMatrix" else "ngTMatrix")
	  }
      })

setAs("graph", "CsparseMatrix",
      function(from) as(as(from, "graphNEL"), "CsparseMatrix"))

setAs("graphNEL", "CsparseMatrix",
      function(from) as(as(from, "TsparseMatrix"), "CsparseMatrix"))

setAs("graphNEL", "TsparseMatrix",
      function(from) {
          nd <- nodes(from)
          dm <- rep.int(length(nd), 2)
	  symm <- edgemode(from) == "undirected"

 	  if(graph.has.weights(from)) {
	      eWts <- edgeWeights(from)
	      lens <- unlist(lapply(eWts, length))
	      i <- rep.int(0:(dm[1]-1), lens) # column indices (0-based)
	      To <- unlist(lapply(eWts, names))
	      j <- as.integer(match(To,nd) - 1:1) # row indices (0-based)
	      ## symm <- symm && <weights must also be symmetric>: improbable
	      ## if(symm) new("dsTMatrix", .....) else
	      new("dgTMatrix", i = i, j = j, x = unlist(eWts),
		  Dim = dm, Dimnames = list(nd, nd))
	  }
 	  else { ## no weights: 0/1 matrix -> logical
              edges <- lapply(from@edgeL[nd], "[[", "edges")
              lens <- unlist(lapply(edges, length))
              ## nnz <- sum(unlist(lens))  # number of non-zeros
              i <- rep.int(0:(dm[1]-1), lens) # column indices (0-based)
              j <- as.integer(unlist(edges) - 1) # row indices (0-based)
              if(symm) {            # symmetric: ensure upper triangle
                  tmp <- i
                  flip <- i > j
                  i[flip] <- j[flip]
                  j[flip] <- tmp[flip]
                  new("nsTMatrix", i = i, j = j, Dim = dm,
                      Dimnames = list(nd, nd), uplo = "U")
              } else {
                  new("ngTMatrix", i = i, j = j, Dim = dm,
                      Dimnames = list(nd, nd))
              }
          }
      })

setAs("sparseMatrix", "graph", function(from) as(from, "graphNEL"))
setAs("sparseMatrix", "graphNEL",
      function(from) as(as(from, "TsparseMatrix"), "graphNEL"))

Tsp2grNEL <- function(from) {
    d <- dim(from)
    if(d[1] != d[2])
	stop("only square matrices can be used as incidence matrices for graphs")
    n <- d[1]
    if(n == 0) return(new("graphNEL"))
    if(is.null(rn <- dimnames(from)[[1]]))
	rn <- as.character(1:n)
    from <- uniq(from) ## Need to 'uniquify' the triplets!

    if(isSymmetric(from)) { # either "symmetricMatrix" or otherwise
	##-> undirected graph: every edge only once!
	if(!is(from, "symmetricMatrix")) {
	    ## a general matrix which happens to be symmetric
	    ## ==> remove the double indices
	    from <- tril(from)
	}
        eMode <- "undirected"
    } else {
        eMode <- "directed"
    }
    ## every edge is there only once, either upper or lower triangle
    ft1 <- cbind(rn[from@i + 1:1], rn[from@j + 1:1])
    ## not yet: graph::ftM2graphNEL(.........)
    ftM2graphNEL(ft1, W = from@x, V= rn, edgemode= eMode)

}
setAs("TsparseMatrix", "graphNEL", Tsp2grNEL)


### Subsetting -- basic things (drop = "missing") are done in ./Matrix.R

### FIXME : we defer to the "*gT" -- conveniently, but not efficient for gC !

## [dl]sparse -> [dl]gT   -- treat both in one via superclass
##                        -- more useful when have "z" (complex) and even more

setMethod("[", signature(x = "sparseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, drop) {
              cl <- class(x)
              viaCl <- paste(.M.kind(x,cl), "gTMatrix", sep='')
              x <- callGeneric(x = as(x, viaCl), i=i, drop=drop)
              ## try_as(x, c(cl, sub("T","C", viaCl)))
              if(is(x, "Matrix") && extends(cl, "CsparseMatrix"))
                  as(x, sub("T","C", viaCl)) else x
          })

setMethod("[", signature(x = "sparseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, drop) {
              cl <- class(x)
              viaCl <- paste(.M.kind(x,cl), "gTMatrix", sep='')
              x <- callGeneric(x = as(x, viaCl), j=j, drop=drop)
              ## try_as(x, c(cl, sub("T","C", viaCl)))
              if(is(x, "Matrix") && extends(cl, "CsparseMatrix"))
                  as(x, sub("T","C", viaCl)) else x
          })

setMethod("[", signature(x = "sparseMatrix",
			 i = "index", j = "index", drop = "logical"),
	  function (x, i, j, drop) {
              cl <- class(x)
              viaCl <- paste(.M.kind(x,cl), "gTMatrix", sep='')
              x <- callGeneric(x = as(x, viaCl), i=i, j=j, drop=drop)
              ## try_as(x, c(cl, sub("T","C", viaCl)))
              if(is(x, "Matrix") && extends(cl, "CsparseMatrix"))
                  as(x, sub("T","C", viaCl)) else x
          })


## setReplaceMethod("[", signature(x = "sparseMatrix", i = "index", j = "missing",
##                                 value = "numeric"),
##                  function (x, i, value) {

##                      stop("NOT YET")

##                      as(r, class(x))
##                  })

## setReplaceMethod("[", signature(x = "sparseMatrix", i = "missing", j = "index",
##                                 value = "numeric"),
##                  function (x, j, value) {

##                      stop("NOT YET")

##                      as(r, class(x))
##                  })

## setReplaceMethod("[", signature(x = "sparseMatrix", i = "index", j = "index",
##                                 value = "numeric"),

##                      stop("NOT YET")

##                      as(r, class(x))
##                  })



setMethod("-", signature(e1 = "sparseMatrix", e2 = "missing"),
          function(e1) { e1@x <- -e1@x ; e1 })
## with the following exceptions:
setMethod("-", signature(e1 = "nsparseMatrix", e2 = "missing"),
          function(e1) callGeneric(as(e1, "dgCMatrix")))
setMethod("-", signature(e1 = "pMatrix", e2 = "missing"),
          function(e1) callGeneric(as(e1, "ngTMatrix")))

## Group method  "Arith"

## have CsparseMatrix methods (-> ./Csparse.R )
## which may preserve "symmetric", "triangular" -- simply defer to those:

setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
	  signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
				       as(e2, "CsparseMatrix")))
setMethod("Arith",
	  signature(e1 = "sparseMatrix", e2 = "numeric"),
	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"), e2))
setMethod("Arith",
	  signature(e1 = "numeric", e2 = "sparseMatrix"),
	  function(e1, e2) callGeneric(e1, as(e2, "CsparseMatrix")))

setMethod("Math",
	  signature(x = "sparseMatrix"),
	  function(x) callGeneric(as(x, "CsparseMatrix")))



### --- show() method ---

## FIXME(?) -- ``merge this'' (at least ``synchronize'') with
## - - -   prMatrix() from ./Auxiliaries.R
prSpMatrix <- function(object, digits = getOption("digits"),
                       maxp = getOption("max.print"), zero.print = ".")
{
    stopifnot(is(object, "sparseMatrix"))
    d <- dim(object)
    if(prod(d) > maxp) { # "Large" => will be "cut"
        ## only coerce to dense that part which won't be cut :
        nr <- maxp %/% d[2]
	m <- as(object[1:max(1, nr), ,drop=FALSE], "Matrix")
    } else {
        m <- as(object, "matrix")
    }
    logi <- is(object,"lsparseMatrix") || is(object,"nsparseMatrix")
    if(logi)
	x <- array("N", # or as.character(NA),
		   dim(m), dimnames=dimnames(m))
    else {
	x <- apply(m, 2, format)
	if(is.null(dim(x))) {# e.g. in	1 x 1 case
	    dim(x) <- dim(m)
	    dimnames(x) <- dimnames(m)
	}
    }
    x <- emptyColnames(x)
    if(is.logical(zero.print))
	zero.print <- if(zero.print) "0" else " "
    if(logi) {
	x[!m] <- zero.print
	x[m] <- "|"
    } else { # non logical
	## show only "structural" zeros as 'zero.print', not all of them..
	## -> cannot use 'm'
	iN0 <- 1:1 + encodeInd(non0ind(object), nr = nrow(x))
	if(length(iN0)) x[-iN0] <- zero.print else x[] <- zero.print
    }
    print(x, quote = FALSE, max = maxp)
    invisible(object)
}

setMethod("show", signature(object = "sparseMatrix"),
   function(object) {
       d <- dim(object)
       cl <- class(object)
       cat(sprintf('%d x %d sparse Matrix of class "%s"\n', d[1], d[2], cl))
       maxp <- getOption("max.print")
       if(prod(d) <= maxp)
	   prSpMatrix(object, maxp = maxp)
       else { ## d[1] > maxp / d[2] >= nr : -- this needs [,] working:
	   nr <- maxp %/% d[2]
	   n2 <- ceiling(nr / 2)
	   nR <- d[1] # nrow
	   prSpMatrix(object[seq(length = min(nR, max(1, n2))), drop = FALSE])
	   cat("\n ..........\n\n")
	   prSpMatrix(object[seq(to = nR, length = min(max(1, nr-n2), nR)),
                             drop = FALSE])
           invisible(object)
       }
   })


## not exported:
setMethod("isSymmetric", signature(object = "sparseMatrix"),
	  function(object, tol = 100*.Machine$double.eps) {
	      ## pretest: is it square?
	      d <- dim(object)
	      if(d[1] != d[2]) return(FALSE)
	      ## else slower test
	      if (is(object, "dMatrix"))
		  ## use gC; "T" (triplet) is *not* unique!
		  isTRUE(all.equal(as(object, "dgCMatrix"),
				   as(t(object), "dgCMatrix"), tol = tol))
	      else if (is(object, "lMatrix"))
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(object, "lgCMatrix"),
			    as(t(object), "lgCMatrix"))
	      else if (is(object, "nMatrix"))
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(object, "ngCMatrix"),
			    as(t(object), "ngCMatrix"))
	      else stop("not yet implemented")
	  })


setMethod("isTriangular", signature(object = "sparseMatrix"),
	  function(object, upper = NA)
              isTriC(as(object, "CsparseMatrix"), upper))

setMethod("isDiagonal", signature(object = "sparseMatrix"),
	  function(object) {
	      gT <- as(object, "TsparseMatrix")
	      all(gT@i == gT@j)
	  })


setMethod("diag", signature(x = "sparseMatrix"),
	  function(x, nrow, ncol = n) diag(as(x, "CsparseMatrix")))

## .as.dgT.Fun
setMethod("colSums",  signature(x = "sparseMatrix"), .as.dgT.Fun)
setMethod("colMeans", signature(x = "sparseMatrix"), .as.dgT.Fun)
## .as.dgC.Fun
setMethod("rowSums", signature(x = "sparseMatrix"), .as.dgC.Fun)
setMethod("rowMeans", signature(x = "sparseMatrix"),.as.dgC.Fun)
