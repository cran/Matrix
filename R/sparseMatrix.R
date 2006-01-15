### Define Methods that can be inherited for all subclasses

### Idea: Coercion between *VIRTUAL* classes -- as() chooses "closest" classes
### ----  should also work e.g. for  dense-triangular --> sparse-triangular !
## setAs("denseMatrix", "sparseMatrix",
##        function(from) {
##            as(as(from, "dgeMatrix")
##        })

## setAs("dMatrix", "lMatrix",
##       function(from) {
##       })


## "graph" coercions -- this needs the graph package which is currently
##  -----               *not* required on purpose
## Note: 'undirected' graph <==> 'symmetric' matrix

setAs("graphNEL", "sparseMatrix",
      function(from) {
          .Call("graphNEL_as_dgTMatrix",
                from,
                symmetric = (from@edgemode == "undirected"),
                PACKAGE = "Matrix")
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
	  function(object, upper) {
	      ## pretest: is it square?
	      d <- dim(object)
	      if(d[1] != d[2]) return(FALSE)
	      ## else slower test
              object <- as(object, "TsparseMatrix")
              i <- object@i
              j <- object@j
              if(upper)
                  all(i < j)## FIXME or "0" that are not structural..
              else
                  all(i > j)## FIXME or "0" that are not structural..
          })

setMethod("isDiagonal", signature(object = "sparseMatrix"),
	  function(object) {
	      gT <- as(object, "TsparseMatrix")
	      all(gT@i == gT@j)
	  })

