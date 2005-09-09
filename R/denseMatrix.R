### Simple fallback methods for all dense matrices:
### These are "cheap" to program, but potentially far from efficient;
### Methods for specific subclasses will overwrite these:

## Using "index" for indices should allow
## integer (numeric), logical, or character (names!) indices :

setMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
			 drop = "logical"),
          function (x, i, drop) {
              r <- as(x, "matrix")[i, , drop=drop]
              if(is.null(dim(r))) r else as(r, class(x))
          })

setMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
			 drop = "logical"),
          function (x, j, drop) {
              r <- as(x, "matrix")[, j, drop=drop]
              if(is.null(dim(r))) r else as(r, class(x))
          })

setMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
			 drop = "logical"),
          function (x, i, j, drop) {
              r <- callGeneric(x = as(x, "matrix"), i=i, j=j, drop=drop)
              if(is.null(dim(r))) r else as(r, class(x))
          })

## Now the "[<-" ones :

## see also those in ./Matrix.R
## It's recommended to use setReplaceMethod() rather than setMethod("[<-",.)
## even though the former is currently just a wrapper for the latter

setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
                                value = "numeric"),
                 function (x, i, value) {
                     r <- as(x, "matrix")
                     r[i, ] <- value
                     as(r, class(x))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
                                value = "numeric"),
                 function (x, j, value) {
                     r <- as(x, "matrix")
                     r[, j] <- value
                     as(r, class(x))
                 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
                                value = "numeric"),
                 function (x, i, j, value) {
                     r <- as(x, "matrix")
                     r[i, j] <- value
                     as(r, class(x))
                 })



## not exported:
setMethod("isSymmetric", signature(object = "denseMatrix"),
	  function(object, ...) {
	      ## pretest: is it square?
	      d <- dim(object)
	      if(d[1] != d[2]) return(FALSE)
	      ## else slower test
	      if (is(object("dMatrix")))
		  isTRUE(all.equal(as(object, "dgeMatrix"),
				   as(t(object), "dgeMatrix"), ...))
	      else if (is(object("lMatrix")))# not possible currently
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(object, "lgeMatrix"),
			    as(t(object), "lgeMatrix"))
	      else stop("not yet implemented")
	  })
