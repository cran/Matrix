#### Toplevel ``virtual'' class "Matrix"


### Virtual coercions -- via smart "helpers" (-> ./Auxiliaries.R)

setAs("Matrix", "sparseMatrix", function(from) as(from, "CsparseMatrix"))
setAs("Matrix", "CsparseMatrix", function(from) as_Csparse(from))
setAs("Matrix", "denseMatrix",  function(from) as_dense(from))

## Maybe TODO:
## setAs("Matrix", "nMatrix", function(from) ....)

## Most of these work; this is a last resort:
setAs(from = "Matrix", to = "matrix", # do *not* call base::as.matrix() here:
      function(from) .bail.out.2("coerce", class(from), class(to)))
setAs(from = "matrix", to = "Matrix", function(from) Matrix(from))

## ## probably not needed eventually:
## setAs(from = "ddenseMatrix", to = "matrix",
##       function(from) {
## 	  if(length(d <- dim(from)) != 2) stop("dim(.) has not length 2")
## 	  array(from@x, dim = d, dimnames = dimnames(from))
##       })

## should propagate to all subclasses:
setMethod("as.matrix", signature(x = "Matrix"), function(x) as(x, "matrix"))
## for 'Matrix' objects, as.array() should be equivalent:
setMethod("as.array",  signature(x = "Matrix"), function(x) as(x, "matrix"))

## head and tail apply to all Matrix objects for which subscripting is allowed:
setMethod("head", signature(x = "Matrix"), utils::head.matrix)
setMethod("tail", signature(x = "Matrix"), utils::tail.matrix)

setMethod("drop", signature(x = "Matrix"),
	  function(x) if(all(dim(x) != 1)) x else drop(as(x, "matrix")))

## slow "fall back" method {subclasses should have faster ones}:
setMethod("as.vector", signature(x = "Matrix", mode = "missing"),
	  function(x) as.vector(as(x, "matrix")))

## mainly need these for "dMatrix" or "lMatrix" respectively, but why not general:
setMethod("as.numeric", signature(x = "Matrix"),
	  function(x, ...) as.numeric(as.vector(x)))
setMethod("as.logical", signature(x = "Matrix"),
	  function(x, ...) as.logical(as.vector(x)))


## "base" has an isSymmetric() S3-generic since R 2.3.0
setMethod("isSymmetric", signature(object = "symmetricMatrix"),
          function(object,tol) TRUE)
setMethod("isSymmetric", signature(object = "triangularMatrix"),
          ## TRUE iff diagonal:
          function(object,tol) isDiagonal(object))

setMethod("isTriangular", signature(object = "triangularMatrix"),
          function(object, ...) TRUE)

setMethod("isTriangular", signature(object = "matrix"), isTriMat)

setMethod("isDiagonal", signature(object = "matrix"), .is.diagonal)



setMethod("dim", signature(x = "Matrix"),
	  function(x) x@Dim, valueClass = "integer")

setMethod("length", "Matrix", function(x) prod(dim(x)))

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

setMethod("all", signature(x = "Matrix"),
          function(x, ..., na.rm) { x <- as(x, "lMatrix"); callGeneric()})
setMethod("any", signature(x = "Matrix"),
          function(x, ..., na.rm) { x <- as(x, "lMatrix"); callGeneric()})

## NOTE:  "&" and "|"  are now in group "Logic" c "Ops" --> ./Ops.R
##        "!" is in ./not.R


Matrix <-
    function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL,
	      sparse = NULL, forceCheck = FALSE)
{
    sparseDefault <- function(m) prod(dim(m)) > 2*sum(isN0(as(m, "matrix")))

    i.M <- is(data, "Matrix")

    if(is.null(sparse1 <- sparse) && (i.M || is(data, "matrix")))
	sparse <- sparseDefault(data)

    doDN <- TRUE
    if (i.M) {
        if(!missing(nrow) || !missing(ncol)|| !missing(byrow))
            warning("'nrow', 'ncol', etc, are disregarded when 'data' is \"Matrix\" already")
	sM <- is(data,"sparseMatrix")
	if(!forceCheck && ((sparse && sM) || (!sparse && !sM)))
	    return(data)
	## else : convert  dense <-> sparse -> at end
    }
    else if (!is.matrix(data)) { ## cut & paste from "base::matrix" :
	if (missing(nrow))
	    nrow <- ceiling(length(data)/ncol)
	else if (missing(ncol))
	    ncol <- ceiling(length(data)/nrow)
	if(length(data) == 1 && is0(data) && !identical(sparse, FALSE)) {
            ## Matrix(0, ...) : always sparse unless "sparse = FALSE":
	    if(is.null(sparse)) sparse1 <- sparse <- TRUE
            i.M <- sM <- TRUE
	    ## will be sparse: do NOT construct full matrix!
	    data <- new(if(is.numeric(data)) "dgTMatrix" else
			if(is.logical(data)) "lgTMatrix" else
			stop("invalid 'data'"),
			Dim = as.integer(c(nrow,ncol)),
			Dimnames = if(is.null(dimnames)) list(NULL,NULL)
			else dimnames)
	} else { ## normal case
	    data <- .Internal(matrix(data, nrow, ncol, byrow))
	    if(is.null(sparse))
		sparse <- sparseDefault(data)
	    dimnames(data) <- dimnames
	}
        doDN <- FALSE
    } else if(!missing(nrow) || !missing(ncol)|| !missing(byrow))
	warning("'nrow', 'ncol', etc, are disregarded for matrix 'data'")

    ## 'data' is now a "matrix" or "Matrix"
    if (doDN && !is.null(dimnames))
	dimnames(data) <- dimnames

    ## check for symmetric / triangular / diagonal :
    isSym <- isSymmetric(data)
    if((isTri <- !isSym))
	isTri <- isTriangular(data)
    isDiag <- isSym # cannot be diagonal if it isn't symmetric
    if(isDiag)
	isDiag <- isDiagonal(data)

    ## Find proper matrix class 'cl'
    cl <-
	if(isDiag && !isTRUE(sparse1))
	    "diagonalMatrix" # -> will automatically check for type
	else {
	    ## consider it's type
	    ctype <-
		if(is(data,"Matrix")) class(data)
		else {
		    if("complex" == (ctype <- typeof(data)))
			"z" else ctype
		}
	    ctype <- substr(ctype, 1,1) # "d", "l", "i" or "z"
	    if(ctype == "z")
		stop("complex matrices not yet implemented in Matrix package")
	    if(ctype == "i") {
		warning("integer matrices not yet implemented in 'Matrix'; ",
			"using 'double' ones'")
		ctype <- "d"
	    }
	    paste(ctype,
		  if(sparse) {
		      if(isSym) "sCMatrix" else
		      if(isTri) "tCMatrix" else "gCMatrix"
		  } else { ## dense
		      if(isSym) "syMatrix" else
		      if(isTri) "trMatrix" else "geMatrix"
		  }, sep="")
	}

    ## Can we coerce and be done?
    if(!canCoerce(data,cl)) { ## try to coerce ``via'' virtual classes
	if(sparse && !sM)
	    data <- as(data, "sparseMatrix")
	else if(!sparse && !is(data, "denseMatrix"))
	    data <- as(data, "denseMatrix")
	if(isTri && !is(data, "triangularMatrix"))
	    data <- as(data, "triangularMatrix")
	else if(isSym && !is(data, "symmetricMatrix"))
	    data <- as(data, "symmetricMatrix")
    }
    ## now coerce in any case .. maybe producing sensible error message:
    as(data, cl)
}

## Methods for operations where one argument is numeric

## Using as.matrix() and rbind()
## in order to get dimnames from names {at least potentially}:

setMethod("%*%", signature(x = "Matrix", y = "numeric"),
	  function(x, y) callGeneric(x, as.matrix(y)))

setMethod("%*%", signature(x = "numeric", y = "Matrix"),
	  function(x, y) callGeneric(matrix(x, nrow = 1, byrow=TRUE), y))

setMethod("crossprod", signature(x = "Matrix", y = "numeric"),
	  function(x, y = NULL) callGeneric(x, as.matrix(y)))
setMethod("crossprod", signature(x = "numeric", y = "Matrix"),
	  function(x, y = NULL)	 callGeneric(as.matrix(x), y))

## The as.matrix() promotion seems illogical to MM,
## but is according to help(tcrossprod, package = "base") :
setMethod("tcrossprod", signature(x = "Matrix", y = "numeric"),
	  function(x, y = NULL) callGeneric(x, as.matrix(y)))
setMethod("tcrossprod", signature(x = "numeric", y = "Matrix"),
	  function(x, y = NULL)	 callGeneric(as.matrix(x), y))

## maybe not optimal
setMethod("solve", signature(a = "Matrix", b = "missing"),
	  function(a, b, ...) solve(a, Diagonal(nrow(a))))

setMethod("solve", signature(a = "Matrix", b = "numeric"),
	  function(a, b, ...) callGeneric(a, as.matrix(b)))
## when no sub-class method is found, bail out
setMethod("solve", signature(a = "Matrix", b = "matrix"),
	  function(a, b, ...) .bail.out.2("solve", class(a), "matrix"))
setMethod("solve", signature(a = "Matrix", b = "Matrix"),
	  function(a, b, ...) .bail.out.2("solve", class(a), class(b)))

## bail-out methods in order to get better error messages
setMethod("%*%", signature(x = "Matrix", y = "Matrix"),
	  function (x, y)
          stop(gettextf('not-yet-implemented method for <%s> %%*%% <%s>',
                        class(x), class(y))))

setMethod("crossprod", signature(x = "Matrix", y = "ANY"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))
setMethod("crossprod", signature(x = "ANY", y = "Matrix"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))
setMethod("tcrossprod", signature(x = "Matrix", y = "ANY"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))
setMethod("tcrossprod", signature(x = "ANY", y = "Matrix"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))

## cheap fallbacks
setMethod("crossprod", signature(x = "Matrix", y = "Matrix"),
	  function(x, y = NULL) t(x) %*% y)
setMethod("tcrossprod", signature(x = "Matrix", y = "Matrix"),
	  function(x, y = NULL) x %*% t(y))

## There are special sparse methods; this is a "fall back":
setMethod("kronecker", signature(X = "Matrix", Y = "ANY",
				 FUN = "ANY", make.dimnames = "ANY"),
	  function(X, Y, FUN, make.dimnames, ...) {
	      if(is(X, "sparseMatrix"))
		  warning("using slow kronecker() method")
	      X <- as(X, "matrix") ; Matrix(callGeneric()) })

setMethod("kronecker", signature(X = "ANY", Y = "Matrix",
				 FUN = "ANY", make.dimnames = "ANY"),
	  function(X, Y, FUN, make.dimnames, ...) {
	      if(is(Y, "sparseMatrix"))
		  warning("using slow kronecker() method")
	      Y <- as(Y, "matrix") ; Matrix(callGeneric()) })


## FIXME: All of these should never be called
setMethod("chol", signature(x = "Matrix"),
	  function(x, pivot = FALSE) .bail.out.1(.Generic, class(x)))
setMethod("determinant", signature(x = "Matrix"),
	  function(x, logarithm = TRUE) .bail.out.1(.Generic, class(x)))

setMethod("diag", signature(x = "Matrix"),
	  function(x, nrow, ncol) .bail.out.1(.Generic, class(x)))
setMethod("t", signature(x = "Matrix"),
	  function(x) .bail.out.1(.Generic, class(x)))

setMethod("dim<-", signature(x = "Matrix", value = "ANY"),
	  function(x, value) {
	      if(!is.numeric(value) || length(value) != 2)
		  stop("dim(.) value must be numeric of length 2")
	      if(prod(dim(x)) != prod(value <- as.integer(value)))
		  stop("dimensions don't match the number of cells")
	      clx <- class(x)
	      if(substring(clx,2) == "geMatrix") {
		  x@Dim <- value
		  if(length(x@factors) > 0)
		      x@factors <- list()
		  x
	      } else if(extends(clx, "denseMatrix")) {
		  x <- as_geSimpl2(x, clx)
		  dim(x) <- value
	      } else { ## FIXME: this is very inefficient for large sparse x
		  Matrix(as.vector(x), value[1], value[2])
	      }
	  })

## MM: More or less "Cut & paste" from
## --- diff.default() from  R/src/library/base/R/diff.R :
setMethod("diff", signature(x = "Matrix"),
	  function(x, lag = 1, differences = 1, ...) {
	      if (length(lag) > 1 || length(differences) > 1 ||
		  lag < 1 || differences < 1)
		  stop("'lag' and 'differences' must be integers >= 1")
	      xlen <- nrow(x)
	      if (lag * differences >= xlen)
		  return(x[,FALSE][0])	# empty of proper mode

	      i1 <- -1:-lag
	      for (i in 1:differences)
		  x <- x[i1, , drop = FALSE] -
		      x[-nrow(x):-(nrow(x)-lag+1), , drop = FALSE]
	      x
	  })

setMethod("image", "Matrix",
	  function(x, ...) { # coercing to sparse is not inefficient,
	      ##	       since we need 'i' and 'j' for levelplot()
	      x <- as(as(x, "sparseMatrix"), "dMatrix")
	      callGeneric()
	  })


## Group Methods

##-> see ./Ops.R

### --------------------------------------------------------------------------
###
### Subsetting "["  and
### SubAssign  "[<-" : The "missing" cases can be dealt with here, "at the top":

## Using "index" for indices should allow
## integer (numeric), logical, or character (names!) indices :

## "x[]":
setMethod("[", signature(x = "Matrix",
			 i = "missing", j = "missing", drop = "ANY"),
	  function (x, i, j, drop) x)

## missing 'drop' --> 'drop = TRUE'
##                     -----------
## select rows
setMethod("[", signature(x = "Matrix", i = "index", j = "missing",
			 drop = "missing"),
	  function(x,i,j, drop) {
	      if(nargs() == 1) { ## e.g. M[0] , M[TRUE],  M[1:2]
		  if(any(i)) as.vector(x)[i] else as.vector(x[1,1])[FALSE]
	      } else callGeneric(x, i=i, drop= TRUE)})

## select columns
setMethod("[", signature(x = "Matrix", i = "missing", j = "index",
			 drop = "missing"),
	  function(x,i,j, drop) callGeneric(x, j=j, drop= TRUE))
setMethod("[", signature(x = "Matrix", i = "index", j = "index",
                         drop = "missing"),
	  function(x,i,j, drop) callGeneric(x, i=i, j=j, drop= TRUE))

## bail out if any of (i,j,drop) is "non-sense"
setMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY", drop = "ANY"),
	  function(x,i,j, drop)
          stop("invalid or not-yet-implemented 'Matrix' subsetting"))

## logical indexing, such as M[ M >= 7 ] *BUT* also M[ M[,1] >= 3,],
## The following is *both* for    M [ <logical>   ]
##                 and also for   M [ <logical> , ]
.M.sub.i.logical <- function (x, i, j, drop)
{
    nA <- nargs()
    if(nA == 2) { ##  M [ M >= 7 ]
        ## FIXME: when both 'x' and 'i' are sparse, this can be very inefficient
	as(x, geClass(x))@x[as.vector(i)]
	## -> error when lengths don't match
    } else if(nA == 3) { ##  M [ M[,1, drop=FALSE] >= 7, ]
	stop("not-yet-implemented 'Matrix' subsetting") ## FIXME

    } else stop("nargs() = ", nA,
		".  Extraneous illegal arguments inside '[ .. ]' ?")
}
setMethod("[", signature(x = "Matrix", i = "lMatrix", j = "missing",
			 drop = "ANY"),
	  .M.sub.i.logical)
setMethod("[", signature(x = "Matrix", i = "logical", j = "missing",
			 drop = "ANY"),
	  .M.sub.i.logical)


## A[ ij ]  where ij is (i,j) 2-column matrix :
.M.sub.i.2col <- function (x, i, j, drop)
{
    nA <- nargs()
    if(nA == 2) { ##  M [ cbind(ii,jj) ]
	if(!is.integer(nc <- ncol(i)))
	    stop("'i' has no integer column number",
		 " should never happen; please report")
	if(is.logical(i))
	    return(.M.sub.i.logical(x,i,j,drop))
	else if(!is.numeric(i) || nc != 2)
	    stop("such indexing must be by logical or 2-column numeric matrix")
	m <- nrow(i)
        if(m == 0) return(vector(mode = .type.kind[.M.kind(x)]))
        ## else
	i1 <- i[,1]
	i2 <- i[,2]
	## potentially inefficient -- FIXME --
	unlist(lapply(seq_len(m), function(j) x[i1[j], i2[j]]))

    } else stop("nargs() = ", nA,
		".  Extraneous illegal arguments inside '[ .. ]' ?")
}
setMethod("[", signature(x = "Matrix", i = "matrix", j = "missing"),# drop="ANY"
	  .M.sub.i.2col)


### "[<-" : -----------------

## x[] <- value :
setReplaceMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                                value = "ANY"),## double/logical/...
	  function (x, value) {
	      ## Fails for 'nMatrix' ... FIXME : make sure have method there
	      x@x <- rep(value, length = length(x@x))
	      validObject(x)# check if type and lengths above match
	      x
          })

## A[ ij ] <- value,  where ij is (i,j) 2-column matrix :
## ----------------   The cheap general method --- FIXME: provide special ones
.M.repl.i.2col <- function (x, i, j, value)
{
    nA <- nargs()
    if(nA == 3) { ##  M [ cbind(ii,jj) ] <- value
	if(!is.integer(nc <- ncol(i)))
	    stop("'i' has no integer column number",
		 " should never happen; please report")
	else if(!is.numeric(i) || nc != 2)
	    stop("such indexing must be by logical or 2-column numeric matrix")
	if(is.logical(i)) {
	    message(".M.repl.i.2col(): drop 'matrix' case ...")
	    i <- c(i) # drop "matrix"
	    return( callNextMethod() )
        }
	if(!is.integer(i)) storage.mode(i) <- "integer"
	if(any(i < 0))
	    stop("negative values are not allowed in a matrix subscript")
	if(any(is.na(i)))
	    stop("NAs are not allowed in subscripted assignments")
	if(any(i0 <- (i == 0))) # remove them
            i <- i[ - which(i0, arr.ind = TRUE)[,"row"], ]
        ## now have integer i >= 1
	m <- nrow(i)
	## mod.x <- .type.kind[.M.kind(x)]
	if(length(value) > 0 && m %% length(value) != 0)
	    warning("number of items to replace is not a multiple of replacement length")
	## recycle:
	value <- rep(value, length = m)
	i1 <- i[,1]
	i2 <- i[,2]
	## inefficient -- FIXME -- (also loses "symmetry" unnecessarily)
	for(k in seq_len(m))
	    x[i1[k], i2[k]] <- value[k]

	x
    } else stop("nargs() = ", nA,
		".  Extraneous illegal arguments inside '[ .. ]' ?")
}

setReplaceMethod("[", signature(x = "Matrix", i = "matrix", j = "missing",
				value = "replValue"),
	  .M.repl.i.2col)


setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
				value = "Matrix"),
		 function (x, i, j, value) {
### *TEMPORARY* diagnostic output:
##                  cat("<Matrix1>[i,j] <- <Matrix1>:\n<Matrix1> = x :")
##                  str(x)
##                  cat("<Matrix2> = value :")
##                  str(value)
##                  cat("i :"); if(!missing(i)) str(i) else cat("<missing>\n")
##                  cat("j :"); if(!missing(j)) str(j) else cat("<missing>\n")

                     callGeneric(x=x, i=i, j=j, value = as.vector(value))
                 })
setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
				value = "Matrix"),
		 function (x, i, j, value)
		 callGeneric(x=x, i=i, j=j, value = as.vector(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
				value = "matrix"),
		 function (x, i, j, value)
		 callGeneric(x=x, i=i, j=j, value = c(value)))

## (ANY,ANY,ANY) is used when no `real method' is implemented :
setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
                                value = "ANY"),
	  function (x, i, j, value) {
              if(!is.atomic(value))
		  stop(sprintf("RHS 'value' (class %s) matches 'ANY', but must match matrix class %s",
			       class(value),class(x)))
              else stop("not-yet-implemented 'Matrix[<-' method")
          })
