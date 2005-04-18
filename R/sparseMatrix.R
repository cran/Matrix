### Define Methods that can be inherited for all subclasses

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



setMethod("[", signature(x = "dsparseMatrix",
			 i = "missing", j = "missing", drop = "ANY"), ## "x[]"
	  function (x, i, j, drop) x)

setMethod("[", signature(x = "dsparseMatrix", i = "numeric", j = "missing",
			 drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "dgTMatrix"), i=i, drop=drop))

setMethod("[", signature(x = "dsparseMatrix", i = "missing", j = "numeric",
			 drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "dgTMatrix"), j=j, drop=drop))

setMethod("[", signature(x = "dsparseMatrix", i = "numeric", j = "missing",
			 drop = "missing"),
	  function(x,i,j, drop) callGeneric(x, i=i, drop= TRUE))

setMethod("[", signature(x = "dsparseMatrix", i = "missing", j = "numeric",
			 drop = "missing"),
	  function(x,i,j, drop) callGeneric(x, j=j, drop= TRUE))

setMethod("[", signature(x = "dsparseMatrix",
			 i = "numeric", j = "numeric", drop = "logical"),
	  function (x, i, j, drop)
          callGeneric(x = as(x, "dgTMatrix"), i=i, j=j, drop=drop))

setMethod("[", signature(x = "dsparseMatrix",
			 i = "numeric", j = "numeric", drop = "missing"),
	  function(x,i,j, drop) callGeneric(x, i=i, j=j, drop= TRUE))

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
