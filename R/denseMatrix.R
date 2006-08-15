### Simple fallback methods for all dense matrices
### These are "cheap" to program, but potentially far from efficient;
### Methods for specific subclasses will overwrite these:

## dense to sparse:
setAs("denseMatrix", "dsparseMatrix",
## MM thought that  as() will take the ``closest'' match; but that fails!
##      function(from) as(as(from, "dgeMatrix"), "dsparseMatrix"))
      function(from) as(as(from, "dgeMatrix"), "dgCMatrix"))

setAs("denseMatrix", "CsparseMatrix",
      function(from) {
          cl <- class(from)
	  notGen <- !is(from, "generalMatrix")
	  if (notGen) { ## e.g. for triangular | symmetric
              ## FIXME: this is a *waste* in the case of packed matrices!
	      if     (extends(cl, "dMatrix")) from <- as(from, "dgeMatrix")
	      else if(extends(cl, "lMatrix")) from <- as(from, "lgeMatrix")
	      else if(extends(cl, "zMatrix")) from <- as(from, "zgeMatrix")
	      else stop("undefined method for class ", cl)
	  }
	  .Call(dense_to_Csparse, from)
      })

setAs("denseMatrix", "TsparseMatrix",
      function(from) as(as(from, "CsparseMatrix"), "TsparseMatrix"))


setMethod("show", signature(object = "denseMatrix"),
          function(object) prMatrix(object))
##- ## FIXME: The following is only for the "dMatrix" objects that are not
##- ##	      "dense" nor "sparse" -- i.e. "packed" ones :
##- ## But these could be printed better -- "." for structural zeros.
##- setMethod("show", signature(object = "dMatrix"), prMatrix)
##- ## and improve this as well:
##- setMethod("show", signature(object = "pMatrix"), prMatrix)
##- ## this should now be superfluous [keep for safety for the moment]:

## Using "index" for indices should allow
## integer (numeric), logical, or character (names!) indices :

## use geClass() when 'i' or 'j' are missing:
## since  symmetric, triangular, .. will not be preserved anyway:
setMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, drop) {
	      r <- as(x, "matrix")[i, , drop=drop]
	      if(is.null(dim(r))) r else as(r, geClass(x))
	  })

setMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, j, drop) {
	      r <- as(x, "matrix")[, j, drop=drop]
	      if(is.null(dim(r))) r else as(r, geClass(x))
	  })

setMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
			 drop = "logical"),
	  function (x, i, j, drop) {
	      r <- callGeneric(x = as(x, "matrix"), i=i, j=j, drop=drop)
	      if(is.null(dim(r))) r else as_geClass(r, class(x))
	  })

## Now the "[<-" ones --- see also those in ./Matrix.R
## It's recommended to use setReplaceMethod() rather than setMethod("[<-",.)
## even though the former is currently just a wrapper for the latter

## FIXME: 1) These are far from efficient
## -----  2) value = "numeric" is only ok for "ddense*"
setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
				value = "numeric"),
		 function (x, i, value) {
		     r <- as(x, "matrix")
		     r[i, ] <- value
		     as(r, geClass(x))
		 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
				value = "numeric"),
		 function (x, j, value) {
		     r <- as(x, "matrix")
		     r[, j] <- value
		     as(r, geClass(x))
		 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
				value = "numeric"),
		 function (x, i, j, value) {
		     r <- as(x, "matrix")
		     r[i, j] <- value
		     as_geClass(r, class(x)) ## was as(r, class(x))
		 })


setMethod("isSymmetric", signature(object = "denseMatrix"),
	  function(object, tol = 100*.Machine$double.eps) {
	      ## pretest: is it square?
	      d <- dim(object)
	      if(d[1] != d[2]) return(FALSE)
	      ## else slower test
	      if (is(object,"dMatrix"))
		  isTRUE(all.equal(as(object, "dgeMatrix"),
				   as(t(object), "dgeMatrix"), tol = tol))
	      else if (is(object, "lMatrix"))# not possible currently
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(object, "lgeMatrix"),
			    as(t(object), "lgeMatrix"))
	      else if (is(object, "zMatrix"))
                  stop("'zMatrix' not yet implemented")
	      else if (is(object, "iMatrix"))
                  stop("'iMatrix' not yet implemented")
	  })

setMethod("isTriangular", signature(object = "triangularMatrix"),
	  function(object, ...) TRUE)

setMethod("isTriangular", signature(object = "denseMatrix"), isTriMat)

setMethod("isDiagonal", signature(object = "denseMatrix"), .is.diagonal)

.as.dge.Fun <- function(x, na.rm = FALSE, dims = 1) {
    x <- as(x, "dgeMatrix")
    callGeneric()
}
setMethod("colSums",  signature(x = "denseMatrix"), .as.dge.Fun)
setMethod("colMeans", signature(x = "denseMatrix"), .as.dge.Fun)
setMethod("rowSums",  signature(x = "denseMatrix"), .as.dge.Fun)
setMethod("rowMeans", signature(x = "denseMatrix"), .as.dge.Fun)
