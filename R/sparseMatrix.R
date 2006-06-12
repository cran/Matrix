### Define Methods that can be inherited for all subclasses

### Idea: Coercion between *VIRTUAL* classes -- as() chooses "closest" classes
### ----  should also work e.g. for  dense-triangular --> sparse-triangular !

##-> see  ./dMatrix.R  and  ./lMatrix.R


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
		 if(symm) "lsTMatrix" else "lgTMatrix")
	  }
      })

## FIXME: in the case of NEL or other sparse graphs, we really should *NOT* go
##        via a *dense* adjacency matrix as we do here :
setAs("graph", "sparseMatrix",
      function(from) as(as(from, "graphAM"), "sparseMatrix"))
## but rather
if(FALSE) { #------------------------- NOT YET -----------------
setAs("graph", "sparseMatrix",
      function(from) as(as(from, "graphNEL"), "sparseMatrix"))

setAs("graphNEL", "sparseMatrix",
      function(from) {
	  nd <- nodes(from)
	  symm <- edgemode(from) == "undirected"
	  if(graph.has.weights(from)) {
	      ## symm <- symm && <weights must also be symmetric>: improbable
	      ## if(symm) new("dsTMatrix", .....) else
	      new("dgTMatrix", .....)
	  }
	  else { ## no weights: 0/1 matrix -> logical
	      if(symm) new("lsTMatrix", .....)
	      else     new("lgTMatrix", .....)
	  }
      })
}# not yet

setAs("sparseMatrix", "graph", function(from) as(from, "graphNEL"))
setAs("sparseMatrix", "graphNEL",
      function(from) as(as(from, "TsparseMatrix"), "graphNEL"))
setAs("TsparseMatrix", "graphNEL",
      function(from) {
          d <- dim(from)
          if(d[1] != d[2])
              stop("only square matrices can be used as incidence matrices for grphs")
          n <- d[1]
          if(n == 0) return(new("graphNEL"))
          if(is.null(rn <- dimnames(from)[[1]]))
              rn <- as.character(1:n)
          from <- uniq(from) ## Need to 'uniquify' the triplets!
          if(isSymmetric(from)) { # because it's "dsTMatrix" or otherwise
              upper <- from@i <= from@j
              ft1 <- cbind(from@i + 1:1, from@j + 1:1)
              graph::ftM2graphNEL(rbind(ft1, ft1[, 2:1]),
                                  W = from@x, V=rn, edgemode="undirected")

          } else { ## not symmetric

              graph::ftM2graphNEL(cbind(from@i + 1:1, from@j + 1:1),
                                  W = from@x, V=rn, edgemode="directed")
          }
          ## stop("'dgTMatrix -> 'graphNEL' method is not yet implemented")
      })




### Subsetting -- basic things (drop = "missing") are done in ./Matrix.R

### FIXME : we defer to the "*gT" -- conveniently, but not efficient for gC !

## [dl]sparse -> [dl]gT   -- treat both in one via superclass
##                        -- more useful when have "z" (complex) and even more

setMethod("[", signature(x = "sparseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, drop) {
              cl <- class(x)
              viaCl <- if(is(x,"dMatrix")) "dgTMatrix" else "lgTMatrix"
              x <- callGeneric(x = as(x, viaCl), i=i, drop=drop)
              ## try_as(x, c(cl, sub("T","C", viaCl)))
              if(is(x, "Matrix") && extends(cl, "CsparseMatrix"))
                  as(x, sub("T","C", viaCl)) else x
          })

setMethod("[", signature(x = "sparseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, drop) {
              cl <- class(x)
              viaCl <- if(is(x,"dMatrix")) "dgTMatrix" else "lgTMatrix"
              x <- callGeneric(x = as(x, viaCl), j=j, drop=drop)
              ## try_as(x, c(cl, sub("T","C", viaCl)))
              if(is(x, "Matrix") && extends(cl, "CsparseMatrix"))
                  as(x, sub("T","C", viaCl)) else x
          })

setMethod("[", signature(x = "sparseMatrix",
			 i = "index", j = "index", drop = "logical"),
	  function (x, i, j, drop) {
              cl <- class(x)
              viaCl <- if(is(x,"dMatrix")) "dgTMatrix" else "lgTMatrix"
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
setMethod("-", signature(e1 = "lsparseMatrix", e2 = "missing"),
          function(e1) callGeneric(as(e1, "dgCMatrix")))
setMethod("-", signature(e1 = "pMatrix", e2 = "missing"),
          function(e1) callGeneric(as(e1, "lgTMatrix")))

### --- show() method ---

prSpMatrix <- function(object, zero.print = ".")
{
    stopifnot(is(object, "sparseMatrix"))
    m <- as(object, "matrix")
    x <- apply(m, 2, format)
    if(is.null(dim(x))) {# e.g. in  1 x 1 case
        dim(x) <- dim(m)
        dimnames(x) <- dimnames(m)
    }
    x <- emptyColnames(x)
    if(is.logical(zero.print))
	zero.print <- if(zero.print) "0" else " "
    ## FIXME: show only "structural" zeros as 'zero.print', not all of them..
    x[m == 0.] <- zero.print
    if(is(object,"lsparseMatrix"))
        x[m] <- "|"
    print(noquote(x))
    invisible(object)
}

setMethod("show", signature(object = "sparseMatrix"),
   function(object) {
       d <- dim(object)
       cl <- class(object)
       cat(sprintf('%d x %d sparse Matrix of class "%s"\n', d[1], d[2], cl))
       maxp <- getOption("max.print")
       if(prod(d) <= maxp)
	   prSpMatrix(object)
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


## .as.dgT.Fun is in ./Tsparse.R
setMethod("colSums",  signature(x = "sparseMatrix"), .as.dgT.Fun)
setMethod("colMeans", signature(x = "sparseMatrix"), .as.dgT.Fun)
## .as.dgC.Fun is in ./Csparse.R
setMethod("rowSums", signature(x = "sparseMatrix"), .as.dgC.Fun)
setMethod("rowMeans", signature(x = "sparseMatrix"),.as.dgC.Fun)

