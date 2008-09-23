#### Superclass Methods for all sparse logical matrices


C2l <- function(from) {
    if(extends(cld <- getClassDef(cl <- class(from)), "lsparseMatrix"))
	return(from)
    ## else
    is.n <- extends(cld, "nsparseMatrix")
    r <- as(.Call(Csparse_to_nz_pattern, if(is.n) from else .Call(Csparse_drop, from, 0),
		  extends(cld, "triangularMatrix")),
	    "lsparseMatrix")
    if(!is.n && any(ina <- is.na(from@x)))
	## NAs must remain NA
	is.na(r@x) <- ina    # strong assumption: r@x "matches" from@x
    r
}

setAs("CsparseMatrix", "lMatrix", C2l)
setAs("CsparseMatrix", "lsparseMatrix", C2l)

setAs("lsparseMatrix", "matrix",
      function(from) as(as(from, "ldenseMatrix"), "matrix"))

setAs("lsparseMatrix", "dsparseMatrix", function(from) as(from, "dMatrix"))


###------- Work via  as(*, lgC) : ------------

## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a lsparseMatrix.

setMethod("%*%", signature(x = "lsparseMatrix", y = "ldenseMatrix"),
	  function(x, y) x %*% as(y, "lsparseMatrix"))

setMethod("%*%", signature(x = "ldenseMatrix", y = "lsparseMatrix"),
	  function(x, y) as(x, "lsparseMatrix") %*% y)

setMethod("crossprod", signature(x = "lsparseMatrix", y = "ldenseMatrix"),
	  function(x, y = NULL) crossprod(x, as(y, "lsparseMatrix")))

setMethod("crossprod", signature(x = "ldenseMatrix", y = "lsparseMatrix"),
	  function(x, y = NULL) crossprod(as(x, "lsparseMatrix"), y))

## and coerce lsparse* to lgC*
setMethod("%*%", signature(x = "lsparseMatrix", y = "lsparseMatrix"),
	  function(x, y) as(x, "lgCMatrix") %*% as(y, "lgCMatrix"))

setMethod("crossprod", signature(x = "lsparseMatrix", y = "lsparseMatrix"),
	  function(x, y = NULL)
	  crossprod(as(x, "lgCMatrix"), as(y, "lgCMatrix")))


setMethod("all", signature(x = "lsparseMatrix"),
	  function(x, ..., na.rm = FALSE) {
	      d <- x@Dim
	      l.x <- length(x@x)
	      if(l.x == prod(d)) ## fully non-zero
		  all(x@x, ..., na.rm = na.rm)
	      else if(is(x, "symmetricMatrix") && l.x == choose(d[1]+1, 2)) {
		  if(.Generic %in% summGener1)
		      all(x@x, ..., na.rm = na.rm)
		  else all(as(x, "generalMatrix")@x, ..., na.rm = na.rm)
	      }
	      else FALSE ## has at least one structural 0
	  })

## setMethod("any", ) ---> ./lMatrix.R

setMethod("image", "lsparseMatrix", function(x, ...) image(as(x,"dMatrix")))
