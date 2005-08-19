#### Toplevel ``virtual'' class "Matrix"

## probably not needed eventually:
setAs(from = "ddenseMatrix", to = "matrix",
      function(from) {
	  if(length(d <- dim(from)) != 2) stop("dim(.) has not length 2")
	  array(from@x, dim = d, dimnames = dimnames(from))
      })

## private function to be used as show() method possibly more than once
prMatrix <- function(object) {
    d <- dim(object)
    cl <- class(object)
    cat(sprintf('%d x %d Matrix of class "%s"\n', d[1], d[2], cl))
    m <- as(object, "matrix")
    maxp <- getOption("max.print")
    if(prod(d) <= maxp) print(m)
    else { ## d[1] > maxp / d[2] >= nr :
	nr <- maxp %/% d[2]
	n2 <- ceiling(nr / 2)
	print(head(m, max(1, n2)))
	cat("\n ..........\n\n")
	print(tail(m, max(1, nr - n2)))
    }
    ## DEBUG: cat("str(.):\n") ; str(object)
    invisible(object)# as print() S3 methods do
}

setMethod("show", signature(object = "ddenseMatrix"), prMatrix)

##- ## FIXME: The following is only for the "dMatrix" objects that are not
##- ##	      "dense" nor "sparse" -- i.e. "packed" ones :
##- ## But these could be printed better -- "." for structural zeros.
##- setMethod("show", signature(object = "dMatrix"), prMatrix)
##- ## and improve this as well:
##- setMethod("show", signature(object = "pMatrix"), prMatrix)
##- ## this should now be superfluous [keep for safety for the moment]:
setMethod("show", signature(object = "Matrix"), prMatrix)

## should propagate to all subclasses:
setMethod("as.matrix", signature(x = "Matrix"), function(x) as(x, "matrix"))

setMethod("dim", signature(x = "Matrix"),
	  function(x) x@Dim, valueClass = "integer")
setMethod("dimnames", signature(x = "Matrix"), function(x) x@Dimnames)
## not exported but used more than once for "dimnames<-" method :
## -- or do only once for all "Matrix" classes ??
dimnamesGets <- function (x, value) {
    d <- dim(x)
    if (!is.list(value) || length(value) != 2 ||
	!(is.null(v1 <- value[[1]]) || length(v1) == d[1]) ||
	!(is.null(v2 <- value[[2]]) || length(v2) == d[2]))
	stop(sprintf("invalid dimnames given for '%s' object", class(x)))
    x@Dimnames <- list(if(!is.null(v1)) as.character(v1),
		       if(!is.null(v2)) as.character(v2))
    x
}
setMethod("dimnames<-", signature(x = "Matrix", value = "list"),
	  dimnamesGets)

setMethod("unname", signature("Matrix", force="missing"),
	  function(obj) { obj@Dimnames <- list(NULL,NULL); obj})

Matrix <-
    function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL)
{
    if (is(data, "Matrix")) return(data)
    if (is.matrix(data)) { val <- data }
    else { ## cut & paste from "base::matrix" :
	if (missing(nrow))
	    nrow <- ceiling(length(data)/ncol)
	else if (missing(ncol))
	    ncol <- ceiling(length(data)/nrow)
	val <- .Internal(matrix(data, nrow, ncol, byrow))
	dimnames(val) <- dimnames
    }
    as(val, "dgeMatrix")
}

## Methods for operations where one argument is numeric

## Using as.matrix() and rbind()
## in order to get dimnames from names {at least potentially}:

setMethod("%*%", signature(x = "Matrix", y = "numeric"),
	  function(x, y) callGeneric(x, as.matrix(y)))

setMethod("%*%", signature(x = "numeric", y = "Matrix"),
	  function(x, y) callGeneric(rbind(x), y))

setMethod("crossprod", signature(x = "Matrix", y = "numeric"),
	  function(x, y = NULL) callGeneric(x, as.matrix(y)))

setMethod("crossprod", signature(x = "numeric", y = "Matrix"),
	  function(x, y = NULL)	 callGeneric(rbind(x), y))

setMethod("solve", signature(a = "Matrix", b = "numeric"),
	  function(a, b, ...) callGeneric(a, as.matrix(b)))

## Subsetting : The "missing" cases can be dealt with here, "at the top":

## "x[]":
setMethod("[", signature(x = "Matrix",
			 i = "missing", j = "missing", drop = "ANY"),
	  function (x, i, j, drop) x)

## missing 'drop' --> 'drop = TRUE'
##                     -----------
## select rows
setMethod("[", signature(x = "Matrix", i = "numeric", j = "missing",
			 drop = "missing"),
	  function(x,i,j, drop) callGeneric(x, i=i, drop= TRUE))
## select columns
setMethod("[", signature(x = "Matrix", i = "missing", j = "numeric",
			 drop = "missing"),
	  function(x,i,j, drop) callGeneric(x, j=j, drop= TRUE))
setMethod("[", signature(x = "Matrix", i = "numeric", j = "numeric",
                         drop = "missing"),
	  function(x,i,j, drop) callGeneric(x, i=i, j=j, drop= TRUE))


## "FIXME:"
## How can we get at   A[ ij ]	where ij is (i,j) 2-column matrix?
##  and                A[ LL ]	where LL is a logical *vector*




if(FALSE) ## The following can't work as long as cbind is function(..., *)
setMethod("cbind", signature(a = "Matrix", b = "Matrix"),
          function(a, b, ...) {
              da <- Dim(a)
              db <- Dim(b)
              if(da[1] != db[1])
                  stop("Matrices must have same number of rows for cbind()ing")
          })



