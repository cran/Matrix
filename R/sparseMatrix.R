### Define Methods that can be inherited for all subclasses

## An idea: Coercion between *VIRTUAL* classes
## -- making sure that result is *actual*!

## setAs("denseMatrix", "sparseMatrix",
##       function(from) {
##       })

## setAs("dMatrix", "lMatrix",
##       function(from) {
##       })


## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a dgeMatrix.

setMethod("%*%", signature(x = "sparseMatrix", y = "ddenseMatrix"),
          function(x, y) callGeneric(x, as(y, "dgeMatrix")))

setMethod("%*%", signature(x = "ddenseMatrix", y = "sparseMatrix"),
          function(x, y) callGeneric(as(x, "dgeMatrix"), y))

setMethod("crossprod", signature(x = "sparseMatrix", y = "ddenseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "sparseMatrix"),
          function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y))

## "graph" coercions -- this needs the graph package which is currently
##  -----               *not* required on purpose
## Note: 'undirected' graph <==> 'symmetric' matrix

setAs("graphNEL", "sparseMatrix",
      function(from) {
          .Call("graphNEL_as_dgTMatrix",
                from,
                symmetric = (from@edgemode == "undirected"))
      })
setAs("graph", "sparseMatrix",
      function(from) as(as(from,"graphNEL"), "sparseMatrix"))

##! if(FALSE) {##--- not yet

setAs("sparseMatrix", "graph", function(from) as(from, "graphNEL"))
setAs("sparseMatrix", "graphNEL",
      function(from) as(as(from, "dgTMatrix"), "graphNEL"))
setAs("dgTMatrix", "graphNEL",
      function(from) {
          d <- dim(from)
          if(d[1] != d[2])
              stop("only square matrices can be used as incidence matrices for grphs")
          n <- d[1]
          if(n == 0) return(new("graphNEL"))
          if(is.null(rn <- dimnames(from)[[1]]))
              rn <- as.character(1:n)
          if(isSymmetric(from)) { # because it's "dsTMatrix" or otherwise
              ## Need to 'uniquify' the triplets!
              upper <- from@i <= from@j
              graph::ftM2graphNEL(cbind(from@i + 1:1, from@j + 1:1),
                                  W = from@x, V=rn, edgemode="undirected")

          } else { ## not symmetric

              graph::ftM2graphNEL(cbind(from@i + 1:1, from@j + 1:1),
                                  W = from@x, V=rn, edgemode="directed")
          }
          stop("'dgTMatrix -> 'graphNEL' method is not yet implemented")
          ## new("graphNEL", nodes = paste(1:n) , edgeL = ...)
      })

##! }#--not_yet



### Subsetting -- basic things (drop = "missing") are done in ./Matrix.R

## 1)  dsparse -> dgT
setMethod("[", signature(x = "dsparseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "dgTMatrix"), i=i, drop=drop))

setMethod("[", signature(x = "dsparseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "dgTMatrix"), j=j, drop=drop))

setMethod("[", signature(x = "dsparseMatrix",
			 i = "index", j = "index", drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "dgTMatrix"), i=i, j=j, drop=drop))

## 2)  lsparse -> lgT
setMethod("[", signature(x = "lsparseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "lgTMatrix"), i=i, drop=drop))

setMethod("[", signature(x = "lsparseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "lgTMatrix"), j=j, drop=drop))

setMethod("[", signature(x = "lsparseMatrix",
			 i = "index", j = "index", drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "lgTMatrix"), i=i, j=j, drop=drop))




### --- show() method ---

emptyColnames <- function(x)
{
    ## Useful for compact printing of (parts) of sparse matrices
    ## possibly  dimnames(x) "==" NULL :
    dimnames(x) <- list(dimnames(x)[[1]], rep("", dim(x)[2]))
    x
}

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
	  function(object, ...) {
	      ## pretest: is it square?
	      d <- dim(object)
	      if(d[1] != d[2]) return(FALSE)
	      ## else slower test
	      if (is(object("dMatrix")))
		  ## use gC; "T" (triplet) is *not* unique!
		  isTRUE(all.equal(as(object, "dgCMatrix"),
				   as(t(object), "dgCMatrix"), ...))
	      else if (is(object("lMatrix")))
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(object, "lgCMatrix"),
			    as(t(object), "lgCMatrix"))
	      else stop("not yet implemented")
	  })
