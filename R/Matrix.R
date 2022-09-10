#### Toplevel ``virtual'' class "Matrix"

## MJ: no longer needed ... have methods for all relevant subclasses
if(FALSE) {
### Virtual coercions -- via smart "helpers" (-> ./Auxiliaries.R)
setAs("Matrix", "CsparseMatrix", function(from) as_Csparse(from))
setAs("Matrix", "sparseMatrix", function(from) as(from, "CsparseMatrix"))
setAs("Matrix", "denseMatrix",  function(from) as_dense(from))
} ## MJ

## Anything: we build on  as.matrix(.) :
## ---       authors can always provide their own specific  setAs(*, "Matrix")
setAs("ANY", "Matrix", function(from) Matrix(as.matrix(from)))

## Most of these work; this is a last resort:
setAs("Matrix", "matrix", # do *not* call base::as.matrix() here:
      function(from) .bail.out.2("coerce", class(from), class(to)))
setAs("matrix", "Matrix", function(from) Matrix(from))

## not better; just for those hating typing
.asmatrix <- function(x) as(x, "matrix")

## Such that also base functions dispatch properly on our classes:
if(.Matrix.avoiding.as.matrix) {
    as.matrix.Matrix <- function(x, ...) {
        if(nonTRUEoption("Matrix.quiet.as.matrix") &&
           nonTRUEoption("Matrix.quiet"))
            warning("as.matrix(<Matrix>) is deprecated (to become a no-op in the future).\nUse  as(x, \"matrix\")  or  .asmatrix(x)  instead.")
        as(x, "matrix")
    }
    as.array.Matrix <- function(x, ...) {
        warning("as.array(<Matrix>) is deprecated. Use  as(x, \"matrix\")  or .asmatrix(x) instead.")
        as(x, "matrix")
    }
} else { ## regularly -- documented since 2005 that this works
    as.array.Matrix <- as.matrix.Matrix <- function(x, ...) as(x, "matrix")
}

## should propagate to all subclasses:
setMethod("as.matrix", signature(x = "Matrix"),
          function(x, ...) as(x, "matrix"))
## for 'Matrix' objects, as.array() should be equivalent:
setMethod("as.array",  signature(x = "Matrix"),
          function(x, ...) as(x, "matrix"))

## head and tail apply to all Matrix objects for which subscripting is allowed:
setMethod("head", signature(x = "Matrix"), utils::head.matrix)
setMethod("tail", signature(x = "Matrix"), utils::tail.matrix)

setMethod("drop", signature(x = "Matrix"),
	  function(x) if(all(x@Dim != 1L)) x else drop(as(x, "matrix")))

## slow "fall back" method {subclasses should have faster ones}:
setMethod("as.vector", "Matrix",
	  function(x, mode) as.vector(as(x, "matrix"), mode))
## so base functions calling as.vector() work too:
## S3 dispatch works for base::as.vector(), but S4 dispatch does not
as.vector.Matrix <- function(x, mode) as.vector(as(x, "matrix"), mode)

if(FALSE) { ## still does not work for c(1, Matrix(2))
## For the same reason (and just in case) also do both S3 and S4 here:
c.Matrix <- function(...) unlist(lapply(list(...), as.vector))
## NB: Must use   signature  '(x, ..., recursive = FALSE)' :
setMethod("c", "Matrix", function(x, ..., recursive) c.Matrix(x, ...))
## The above is not sufficient for  c(NA, 3:2, <Matrix>, <matrix>)
setMethod("c", "numMatrixLike", function(x, ..., recursive) c.Matrix(x, ...))
}# not yet

setAs("Matrix", "vector",  function(from) as.vector (as(from, "matrix")))
setAs("Matrix", "numeric", function(from) as.numeric(as(from, "matrix")))
setAs("Matrix", "logical", function(from) as.logical(as(from, "matrix")))
setAs("Matrix", "integer", function(from) as.integer(as(from, "matrix")))
setAs("Matrix", "complex", function(from) as.complex(as(from, "matrix")))

## mainly need these for "dMatrix" or "lMatrix" respectively, but why not general:
setMethod("as.numeric", signature(x = "Matrix"),
	  function(x, ...) as.numeric(as.vector(x)))
setMethod("as.logical", signature(x = "Matrix"),
	  function(x, ...) as.logical(as.vector(x)))

setMethod("mean", signature(x = "sparseMatrix"),
	  function(x, ...) mean(as(x,"sparseVector"), ...))
setMethod("mean", signature(x = "sparseVector"),
	  function(x, trim = 0, na.rm = FALSE, ...)
      {
	  if (na.rm) # remove NAs such that new length() is ok
	      x <- x[!is.na(x)] # remains sparse!
	  if(is0(trim)) sum(x) / length(x)
	  else {
	      ## fast trimmed mean for sparseVector:
	      ## ---> we'd need fast & sparse  sort(<sparseV>).
	      ##      Normally this means to define a xtfrm() method;
	      ##      however, that plus  x[order(x, ..)]  will NOT be sparse
	      ## TODO: sortSparseVector(.)
	      warning("trimmed mean of 'sparseVector' -- suboptimally using as.numeric(.)")
	      mean(as.numeric(x), trim=trim)
	  }
      })
## for the non-"sparseMatrix" ones:
setMethod("mean", signature(x = "Matrix"),
	  function(x, trim = 0, na.rm = FALSE, ...)
      {
	  if (na.rm)
	      x <- x[!is.na(x)]
	  if(is0(trim)) sum(x) / length(x)
	  else mean(as.numeric(x), trim=trim)
      })


## for non-"sparseMatrix" :
setMethod("cov2cor", signature(V = "Matrix"),
	  function(V) { ## was as(cov2cor(as(V, "matrix")), "dpoMatrix"))
	      r <- V
	      p <- (d <- dim(V))[1]
	      if(p != d[2]) stop("'V' is not a square matrix")
	      Is <- sqrt(1/diag(V)) # diag( 1/sigma_i )
	      if(any(!is.finite(Is)))
		  warning("diag(.) had 0 or NA entries; non-finite result is doubtful")
              Is <- Diagonal(x = Is)
              r <- Is %*% V %*% Is
	      r[cbind(1:p,1:p)] <- 1 # exact in diagonal
	      as(forceSymmetric(r), "dpoMatrix")
          })

## MJ: no longer needed ... replacement in ./unpackedMatrix.R
if(FALSE) {
## "base" has an isSymmetric() S3-generic since R 2.3.0
setMethod("isTriangular", signature(object = "matrix"), isTriMat)
setMethod("isDiagonal", signature(object = "matrix"), .is.diagonal)

setMethod("symmpart", signature(x = "matrix"),
          function(x) symmetrizeDimnames(x + t(x)) / 2)
setMethod("skewpart", signature(x = "matrix"),
          function(x) symmetrizeDimnames(x - t(x)) / 2)
} ## MJ

setMethod("symmpart", signature(x = "Matrix"),
	  function(x) symmpart(as(x, "dMatrix")))

setMethod("skewpart", signature(x = "Matrix"),
	  function(x) skewpart(as(x, "dMatrix")))

setMethod("dim", signature(x = "Matrix"), function(x) x@Dim)

setMethod("length", "Matrix", function(x) prod(x@Dim))

setMethod("dimnames", signature(x = "Matrix"), function(x) x@Dimnames)

## not exported but used more than once for "dimnames<-" method :
## -- or do only once for all "Matrix" classes ??
dimnamesGets <- function(x, value) {
    if (is.character(s <- validDN(value, dim(x)))) stop(s)
    x@Dimnames <- fixupDN(value)
    x
}
dimnamesGetsNULL <- function(x) {
    message("dimnames(.) <- NULL translated to\ndimnames(.) <- list(NULL,NULL)")
    x@Dimnames <- list(NULL, NULL)
    x
}

setMethod("dimnames<-", signature(x = "Matrix", value = "list"),
          dimnamesGets)

setMethod("dimnames<-", signature(x = "Matrix", value = "NULL"),
	  function(x, value) dimnamesGetsNULL(x))

setMethod("dimnames<-", signature(x = "compMatrix", value = "list"),
          function(x, value) {
              if(length(x@factors)) x@factors <- list()
              dimnamesGets(x, value)
          })

setMethod("dimnames<-", signature(x = "compMatrix", value = "NULL"),
          function(x, value) {
              if(length(x@factors)) x@factors <- list()
              dimnamesGetsNULL(x)
          })

setMethod("unname", signature("Matrix", force = "missing"),
	  function(obj) { obj@Dimnames <- list(NULL, NULL); obj})

if(FALSE) {
## This version was used in Matrix <= 1.4-1 but has several "bugs":
## * Matrix(`dimnames<-`(diag(1), list("A", "B")), doDiag = TRUE)
##   is not a diagonalMatrix
## * Matrix(<diagonalMatrix>, doDiag = FALSE) is a diagonalMatrix
##   rather than a symmetricMatrix
## * Matrix(<dgeMatrix>, sparse =  TRUE, forceCheck = FALSE),
##   Matrix(<dgCMatrix>, sparse = FALSE, forceCheck = FALSE)
##   _do_ check for structure, possibly resulting in symmetricMatrix
##   or triangularMatrix
## * Matrix(sparseVector(1, 1L, 10L)), etc. fail because spV2M()
##   does not "see" the missingness of 'nrow', 'ncol'
## * Matrix(table(1)), Matrix(table(1, 1, 1)), etc. give bad errors
## * Matrix(structure(matrix(1, 2, 2), class = "zzz")) throws a bad
##   error rather than unclassing as it essentially does for "table"
## * Matrix(x, ...), where 'x' is a vector of type other than logical,
##   integer, and double, evaluates .External(Mmatrix, ...) before
##   throwing an error, hence allocating unnecessarily
## * Matrix(0, n, n) constructs an intermediate dsCMatrix instead
##   of returning a ddiMatrix right away
## * Matrix(0, nrow, ), Matrix(0, , ncol) are wrong for 'nrow', 'ncol'
##   in the interval (0, 1) and throw bad errors for 'nrow', 'ncol'
##   equal to 0
## * Calls is(data, .) much more than necessary ...
Matrix <- function (data = NA, nrow = 1, ncol = 1, byrow = FALSE,
                    dimnames = NULL, sparse = NULL,
                    doDiag = TRUE, forceCheck = FALSE)
{
    i.M <- is(data, "Matrix")
    sM <- FALSE
    if(i.M) {
	if(is(data, "diagonalMatrix")) return(data) # in all cases
	sV <- FALSE
    } else if(inherits(data, "table")) # special treatment
	class(data) <- "matrix" # "matrix" first for S4 dispatch
    else if(is(data, "sparseVector")) {
	data <- spV2M(data, nrow, ncol, byrow=byrow)
	i.M <- sparse <- forceCheck <- sM <- sV <- TRUE
    }
    if(is.null(sparse) && (i.M || is(data, "matrix")))
	sparse <- sparseDefault(data)
    doDN <- TRUE # by default
    if (i.M) {
	if (!sV) {
	    if(!missing(nrow) || !missing(ncol)|| !missing(byrow))
		warning("'nrow', 'ncol', etc, are disregarded when 'data' is \"Matrix\" already")
	    sM <- is(data,"sparseMatrix")
	    if(!forceCheck && ((sparse && sM) || (!sparse && !sM)))
		return(data)
	    ## else : convert  dense <-> sparse -> at end
	}
    }
    else if(!is.matrix(data)) { ## cut & paste from "base::matrix" :
	## avoid copying to strip attributes in simple cases
	if (is.object(data) || !is.atomic(data)) data <- as.vector(data)
	if(length(data) == 1 && is0(data) && !identical(sparse, FALSE)) {
	    ## Matrix(0, ...) : always sparse unless "sparse = FALSE":
	    if(is.null(sparse)) sparse <- TRUE
	    i.M <- sM <- TRUE
	    if (missing(nrow)) nrow <- ceiling(1/ncol) else
	    if (missing(ncol)) ncol <- ceiling(1/nrow)
            isSym <- nrow == ncol
	    ## will be sparse: do NOT construct full matrix!
	    data <- new(paste0(if(is.numeric(data)) "d" else
                               if(is.logical(data)) "l" else
                               stop("invalid 'data'"),
                               if(isSym) "s" else "g", "CMatrix"),
			p = rep.int(0L, ncol+1L),
			Dim = as.integer(c(nrow,ncol)),
			Dimnames = if(is.null.DN(dimnames)) list(NULL,NULL)
			else dimnames)
	} else { ## normal case
	    data <- .External(Mmatrix,
			      data, nrow, ncol, byrow, dimnames,
			      missing(nrow), missing(ncol))
	    if(is.null(sparse))
		sparse <- sparseDefault(data)
	}
        doDN <- FALSE # .. set above
    } else if(!missing(nrow) || !missing(ncol)|| !missing(byrow)) ## i.m == is.matrix(.)
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
	isDiag <- doDiag && isDiagonal(data)

    ## try to coerce ``via'' virtual classes
    if(isDiag) { ## diagonal is preferred to sparse !
	data <- as(data, "diagonalMatrix")
	isSym <- FALSE
    } else if(sparse && !sM)
	data <- as(data, "sparseMatrix")
    else if(!sparse) {
	if(i.M) { ## data is 'Matrix'
	    if(!is(data, "denseMatrix"))
		data <- as(data, "denseMatrix")
	} else { ## data is "matrix" (and result "dense" -> go via "general"
	    ctype <- typeof(data)
	    if (ctype == "complex")
		stop("complex matrices not yet implemented in Matrix package")
	    if (ctype == "integer") ## integer Matrices not yet implemented
		storage.mode(data) <- "double"
	    data <- new(paste0(.M.kind(data), "geMatrix"),
			Dim = dim(data),
			Dimnames = .M.DN(data),
			x = c(data))
	}
    }

    if(isTri && !is(data, "triangularMatrix")) {
	if(attr(isTri,"kind") == "L") tril(data) else triu(data)
    } else if(isSym && !is(data, "symmetricMatrix"))
	forceSymmetric(data)
    else
	data
}
} else {
## This version fixes the issues documented above
Matrix <- function(data = NA, nrow = 1, ncol = 1, byrow = FALSE,
                   dimnames = NULL, sparse = NULL,
                   doDiag = TRUE, forceCheck = FALSE)
{
    i.M <- i.sM <- i.dM <- i.sV <- i.m <- FALSE
    mnrow <- missing(nrow)
    mncol <- missing(ncol)
    if(isS4(data)) {
        cld <- getClassDef(class(data))
        i.M <- extends(cld, "Matrix")
        if(i.M) {
            i.sM <- extends(cld, "sparseMatrix")
            i.dM <- i.sM && extends(cld, "diagonalMatrix")
        } else if(extends(cld, "sparseVector")) {
            ## need to transmit missingness to 'spV2M'
            call. <- quote(spV2M(x = data, nrow =, ncol =, byrow = byrow))
            if(!mnrow)
                call.[[3L]] <- quote(nrow)
            if(!mncol)
                call.[[4L]] <- quote(ncol)
            data <- eval(call.)
            i.M <- i.sM <- i.sV <- forceCheck <- TRUE
        }
    } else {
        i.m <- is.matrix(data)
    }
    if(!i.M) {
        ## validate non-Matrix 'data', throwing type errors _early_
        if(is.object(data)) {
            if(i.m)
                class(data) <- NULL # retaining 'dim'
            else
                data <- as.vector(data)
        }
        mode. <- mode(data)
        kind <- switch(mode., numeric = "d", logical = "l",
                       stop("invalid 'data'"))
    }
    if(i.M || i.m) {
        ## 'data' is a Matrix or a numeric or logical matrix
        ## without a 'class' attribute
        if(!i.sV && !(mnrow && mncol && missing(byrow)))
            warning("'nrow', 'ncol', 'byrow' disregarded for [mM]atrix 'data'")
        if(!is.null(dimnames))
            dimnames(data) <- dimnames
        if(is.null(sparse))
            sparse <- sparseDefault(data)
        if(i.M) {
            ## return early in these cases:
            if(i.dM)
                ## !doDiag has been documented to result in a coercion to
                ## symmetricMatrix; we must use diag2*() below because the
                ## "usual" as(<diagonalMatrix>, "(Csparse|unpacked)Matrix")
                ## inherits from triangularMatrix, _not_ symmetricMatrix
                return(if(doDiag)
                           data
                       else if(sparse)
                           .diag2sparse(data, ".sC", "U", FALSE)
                       else .diag2dense(data, ".sy", "U"))
            if(!forceCheck)
                return(if(i.sM == sparse)
                           data
                       else if(sparse)
                           as(data, "CsparseMatrix")
                       else as(data, "unpackedMatrix"))
        }
    } else {
        ## 'data' is a numeric or logical vector or non-matrix array
        ## without a 'class' attribute
        if(length(data) == 1L && !is.na(data) && data == 0 &&
           (is.null(sparse) || sparse)) {
            ## Matrix(0, ...): sparseMatrix unless sparse=FALSE
            ## MJ: we should _try_ to behave as R's do_matrix()
            ##     in the edge cases ... integer overflow is "OK"
            ##     since anyNA(Dim) is caught by validity methods
            if(mnrow == mncol) {
                nrow <- as.integer(nrow)
                ncol <- as.integer(ncol)
            } else if(mnrow) {
                ncol <- as.integer(ncol)
                if(ncol == 0L)
                    stop("data is too long")
                nrow <- as.integer(ceiling(1 / ncol))
            } else {
                nrow <- as.integer(nrow)
                if(nrow == 0L)
                    stop("data is too long")
                ncol <- as.integer(ceiling(1 / nrow))
            }
            square <- nrow == ncol
            if(is.null(dimnames))
                dimnames <- list(NULL, NULL)
            if(square && doDiag)
                return(new(paste0(kind, "diMatrix"),
                           Dim = c(nrow, ncol),
                           Dimnames = dimnames,
                           x = vector(mode., nrow)))
            data <- new(paste0(kind, if(square) "s" else "g", "CMatrix"),
                        Dim = c(nrow, ncol),
			Dimnames = dimnames,
                        p = integer(ncol + 1))
            i.M <- i.sM <- sparse <- TRUE
	} else {
            ## usual case: vector|array->matrix
	    data <- .External(Mmatrix,
			      data, nrow, ncol, byrow, dimnames, mnrow, mncol)
            if(is.null(sparse))
                sparse <- sparseDefault(data)
            i.m <- TRUE
	}
    }

    ## 'data' is a Matrix (but _not_ a diagonalMatrix) or a
    ## numeric or logical matrix without a 'class' attribute
    if(doDiag && isDiagonal(data))
        ## as(<[mM]atrix>, "diagonalMatrix") uses check = TRUE (a waste)
        return(.M2diag(data, check = FALSE))
    if(i.m || i.sM != sparse) {
        data <- as(data, if(sparse) "CsparseMatrix" else "unpackedMatrix")
        if(i.m)
            ## as(<matrix>, "CsparseMatrix"), as(<matrix>, "unpackedMatrix")
            ## already check for symmetric, triangular structure
            return(data)
    }
    if(!is(data, "generalMatrix"))
        data
    else if(isSymmetric(data))
        forceSymmetric(data)
    else if(!(it <- isTriangular(data)))
        data
    else if(attr(it, "kind") == "U")
        triu(data)
    else tril(data)
}
}

## There are special sparse methods in  ./kronecker.R  ; this is a "fall back":
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

## The ``Right Thing'' to do :
## base::det() calls [base::]determinant();
## our det() should call our determinant() :
det <- base::det
environment(det) <- environment()## == asNamespace("Matrix")

## MJ: no longer needed as methods are available for all subclasses
if(FALSE) {
setMethod("diag", signature(x = "Matrix"),
	  function(x, nrow, ncol, names=TRUE) .bail.out.1("diag", class(x)))
setMethod("t", signature(x = "Matrix"),
	  function(x) .bail.out.1(.Generic, class(x)))
setMethod("chol", signature(x = "Matrix"),
	  function(x, pivot, ...) .bail.out.1("chol", class(x)))
} ## MJ

## TODO: activate later
if(FALSE)
setMethod("diag<-", signature(x = "Matrix"),
	  function(x, value) .bail.out.1("diag", class(x)))

## NB: "sparseMatrix" works via "sparseVector"
setMethod("rep", "Matrix", function(x, ...) rep(as(x, "matrix"), ...))

## We want to use all.equal.numeric() *and* make sure that uses
## not just base::as.vector but the generic with our methods:
all.equal_num <- base::all.equal.numeric ## from <R>/src/library/base/R/all.equal.R
environment(all.equal_num) <- environment()## == as.environment("Matrix")

all.equal_Mat <- function(target, current, check.attributes = TRUE,
                          factorsCheck = FALSE, ...)
{
    msg <- attr.all_Mat(target, current, check.attributes=check.attributes,
                        factorsCheck=factorsCheck, ...)
    if(is.list(msg)) msg[[1]]
    else .a.e.comb(msg,
		   all.equal_num(as.vector(target), as.vector(current),
				 check.attributes=check.attributes, ...))
}
## The all.equal() methods for dense matrices (and fallback):
setMethod("all.equal", c(target = "Matrix", current = "Matrix"),
	  all.equal_Mat)
setMethod("all.equal", c(target = "Matrix", current = "ANY"),
	  all.equal_Mat)
setMethod("all.equal", c(target = "ANY", current = "Matrix"),
	  all.equal_Mat)
## -> ./sparseMatrix.R, ./sparseVector.R  have specific methods



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

## Group Methods

## NOTE:  "&" and "|"  are now in group "Logic" c "Ops" --> ./Ops.R
##        "!" is in ./not.R

## Further, see ./Ops.R
##                ~~~~~


### --------------------------------------------------------------------------
###
### Subsetting "["  and
### SubAssign  "[<-" : The "missing" cases can be dealt with here, "at the top":

## Using "index" for indices should allow
## integer (numeric), logical, or character (names!) indices :

## "x[]" *or* x[,] *or* x[, , drop=<TRUE/FALSE>]
setMethod("[", signature(x = "Matrix",
			 i = "missing", j = "missing", drop = "missing"),
	  function(x,i,j, ..., drop) {
	      Matrix.msg("M[m,m,m] : nargs()=",nargs(), .M.level = 2)
	      if(nargs() == 2) { ##  M[]
		  x
	      } else { # nargs() = 3
		  callGeneric(x, , , drop=TRUE)
	      }
	  })
setMethod("[", signature(x = "Matrix",
			 i = "missing", j = "missing", drop = "logical"),
	  function (x, i, j, ..., drop) {
	      Matrix.msg(sprintf("M[m,m, %s] : nargs()=%d", format(drop), nargs()), .M.level = 2)
	      if(isTRUE(drop) && any(x@Dim == 1L)) drop(as(x, "matrix")) else x
	  })
## otherwise:  drop is neither missing nor logical
setMethod("[", signature(x = "Matrix",
			 i = "missing", j = "missing", drop = "ANY"),
	  function (x, i, j, ..., drop) stop("invalid 'drop' in Matrix subsetting"))


## missing 'drop' --> 'drop = TRUE'
##                     -----------
## select rows __ or __ vector indexing:
setMethod("[", signature(x = "Matrix", i = "index", j = "missing", drop = "missing"),
	  function(x,i,j, ..., drop) {
	      Matrix.msg("M[i,m,m] : nargs()=",nargs(), .M.level = 2)
	      if(nargs() == 2) { ## e.g. M[0] , M[TRUE], M[1:2], M[-7]
                  .M.vectorSub(x,i)
	      } else {
		  callGeneric(x, i=i, , drop=TRUE)
		  ##		      ^^
	      }
	  })
## select columns
setMethod("[", signature(x = "Matrix", i = "missing", j = "index", drop = "missing"),
	  function(x,i,j, ..., drop) {
	      Matrix.msg("M[m,i,m] : nargs()=",nargs(), .M.level = 2)
	      callGeneric(x, , j=j, drop= TRUE)
	  })
## select both rows *and* columns
setMethod("[", signature(x = "Matrix", i = "index", j = "index", drop = "missing"),
	  function(x,i,j, ..., drop) {
	      Matrix.msg("M[i,i,m] : nargs()=",nargs(), .M.level = 2)
	      callGeneric(x, i=i, j=j, drop= TRUE)
	  })

## bail out if any of (i,j,drop) is "non-sense"
setMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY", drop = "ANY"),
	  function(x,i,j, ..., drop)
	  stop("invalid or not-yet-implemented 'Matrix' subsetting"))

## logical indexing, such as M[ M >= 7 ] *BUT* also M[ M[,1] >= 3,],
## The following is *both* for    M [ <logical>   ]
##                 and also for   M [ <logical> , ]
.M.sub.i.logical <- function (x, i, j, ..., drop)
{
    nA <- nargs() # counts 'M[i]' as 2 arguments,  'M[i,]' as 3
    Matrix.msg("M[logi,m,m] : nargs()=", nA, .M.level = 2)
    if(nA == 2) { ##  M [ M >= 7 ]
	## FIXME: when both 'x' and 'i' are sparse, this can be very inefficient
	if(is(x, "sparseMatrix"))
	    message("<sparse>[ <logic> ]: .M.sub.i.logical() maybe inefficient")
	as(as(x, "generalMatrix"), "denseMatrix")@x[as.vector(i)]
	## -> error when lengths don't match
    }
    else if(nA == 3) { ## M[ <logic>, ]  e.g.,  M [ M[,1, drop=FALSE] >= 7, ]  or M[TRUE,]
	if(length(i) && x@Dim[1L] && !anyNA(i) && all(i)) ## select everything
	    x
	else ## not selecting all -> result is *NOT* diagonal/triangular/symmetric/..
	    ## keep j missing, but  drop = "logical"
	    callGeneric(as(x,"generalMatrix"), i = i, , drop = TRUE)

    } else stop(gettextf(
		"nargs() = %d.  Extraneous illegal arguments inside '[ .. ]' (i.logical)?",
			 nA), domain=NA)
}

## instead of using 'drop = "ANY"' {against ambiguity notices}:
for(ii in c("nMatrix", "lMatrix", "logical"))
    setMethod("[", signature(x = "Matrix", i = ii, j = "missing", drop = "missing"),
	      .M.sub.i.logical)
rm(ii)

##' x[ ij ]  where ij is (i,j) 2-column matrix
##' @note only called from  .M.sub.i.2col(x, i) below
subset.ij <- function(x, ij) {
    m <- nrow(ij)
    if(m > 3) {
        cld <- getClassDef(class(x))
	sym.x <- extends(cld, "symmetricMatrix")
	if(sym.x) {
	    W <- if(x@uplo == "U") # stored only [i,j] with i <= j
		ij[,1] > ij[,2] else ij[,1] < ij[,2]
	    if(any(W))
		ij[W,] <- ij[W, 2:1]
        }
        if(extends(cld, "sparseMatrix")) {
	    ## do something smarter:
	    di <- dim(x)
	    if(!extends(cld, "CsparseMatrix")) {
		x <- as(x, "CsparseMatrix") # simpler; our standard
		cld <- getClassDef(class(x))
	    }
	    tri.x <- extends(cld, "triangularMatrix")
	    if(tri.x) {
		## need these for the 'x' slot in any case
		if (x@diag == "U") x <- .Call(Csparse_diagU2N, x)
		## slightly more efficient than non0.i() or non0ind():
		ij.x <- .Call(compressed_non_0_ij, x, isC=TRUE)
	    } else { ## symmetric / general : for symmetric, only "existing" part
		ij.x <- non0.i(x, cld)
	    }

	    m1 <- .Call(m_encodeInd, ij.x, di, orig1=FALSE, checkBounds=FALSE)
            m2 <- .Call(m_encodeInd, ij,   di, orig1= TRUE, checkBounds= TRUE)
	    mi <- match(m2, m1, nomatch=0)
	    mmi <- mi != 0L ## == (m2 %in% m1)
	    ## Result: all FALSE or 0  apart from where we match non-zero entries
	    ans <- vector(mode = .type.kind[.M.kind(x)], length = m)
	    ## those that are *not* zero:
	    ans[mmi] <- if(extends(cld, "nsparseMatrix")) TRUE else x@x[mi[mmi]]
	    if(any(ina <- is.na(m2))) # has one or two NA in that (i,j) row
		is.na(ans) <- ina
	    ans

        } else { ## non-sparse : dense
            ##---- NEVER happens:  'denseMatrix' has its own setMethod(.) !
            message("m[<ij-matrix>]: inefficiently indexing single elements - should not happen, please report!")
            i1 <- ij[,1]
            i2 <- ij[,2]
            ## very inefficient for large m
            unlist(lapply(seq_len(m), function(j) x[i1[j], i2[j]]))
        }
    } else { # 1 <= m <= 3
        i1 <- ij[,1]
        i2 <- ij[,2]
        unlist(lapply(seq_len(m), function(j) x[i1[j], i2[j]]))
    }
}

## A[ ij ]  where ij is (i,j) 2-column matrix -- but also when that is logical mat!
.M.sub.i.2col <- function (x, i, j, ..., drop)
{
    nA <- nargs()
    if(nA != 2)
        stop(domain=NA, gettextf(
            "nargs() = %d.  Extraneous illegal arguments inside '[ .. ]' (i.2col)?", nA))
    ## else: (nA == 2):	 M [ cbind(ii,jj) ] or M [ <logical matrix> ]
    if(!is.integer(nc <- ncol(i)))
        stop(".M.sub.i.2col(): 'i' has no integer column number;\n should never happen; please report")
    if(is.logical(i))
        return(.M.sub.i.logical(x, i=i)) # call with 2 args!
    else if(!is.numeric(i) || nc != 2)
        stop("such indexing must be by logical or 2-column numeric matrix")
    if(!nrow(i)) return(vector(mode = .type.kind[.M.kind(x)]))
    ## else
    subset.ij(x, i)

}
setMethod("[", signature(x = "Matrix", i = "matrix", j = "missing"),# drop="ANY"
	  .M.sub.i.2col)
## just against ambiguity notices :
setMethod("[", signature(x = "Matrix", i = "matrix", j = "missing", drop="missing"),
	  .M.sub.i.2col)


### "[<-" : -----------------

## A[ ij ] <- value,  where ij is (i,j) 2-column matrix :
## ----------------
## The cheap general method, now only used for "pMatrix","indMatrix"
## sparse all use  .TM.repl.i.mat()
## NOTE:  need '...' below such that setMethod() does
##	  not use .local() such that nargs() will work correctly:
.M.repl.i.2col <- function (x, i, j, ..., value)
{
    nA <- nargs()
    if(nA == 3) { ##  M [ cbind(ii,jj) ] <- value  or M [ Lmat ] <- value
	if(!is.integer(nc <- ncol(i)))
	    stop(".M.repl.i.2col(): 'i' has no integer column number;\n should never happen; please report")
	else if(!is.numeric(i) || nc != 2)
	    stop("such indexing must be by logical or 2-column numeric matrix")
	if(is.logical(i)) {
	    message(".M.repl.i.2col(): drop 'matrix' case ...")
	    ## c(i) : drop "matrix" to logical vector
	    return( callGeneric(x, i=c(i), value=value) )
	}
	if(!is.integer(i)) storage.mode(i) <- "integer"
	if(any(i < 0))
	    stop("negative values are not allowed in a matrix subscript")
	if(anyNA(i))
	    stop("NAs are not allowed in subscripted assignments")
	if(any(i0 <- (i == 0))) # remove them
            i <- i[ - which(i0, arr.ind = TRUE)[,"row"], ]
        ## now have integer i >= 1
	m <- nrow(i)
	## mod.x <- .type.kind[.M.kind(x)]
	if(length(value) > 0 && m %% length(value) != 0)
	    warning("number of items to replace is not a multiple of replacement length")
	## recycle:
	value <- rep_len(value, m)
	i1 <- i[,1]
	i2 <- i[,2]
	if(m > 2)
	    message("m[ <ij-matrix> ] <- v: inefficiently treating single elements")
	## inefficient -- FIXME -- (also loses "symmetry" unnecessarily)
	for(k in seq_len(m))
	    x[i1[k], i2[k]] <- value[k]

	x
    } else stop(gettextf(
		"nargs() = %d.  Extraneous illegal arguments inside '[ .. ]' ?",
			 nA), domain=NA)
}

setReplaceMethod("[", signature(x = "Matrix", i = "matrix", j = "missing",
				value = "replValue"),
		 .M.repl.i.2col)

## Three catch-all methods ... would be very inefficient for sparse*
## --> extra methods in ./sparseMatrix.R
setReplaceMethod("[", signature(x = "Matrix", i = "missing", j = "ANY",
				value = "Matrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, , j=j, value = as.vector(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "missing",
				value = "Matrix"),
		 function (x, i, j, ..., value)
		     if(nargs() == 3)
			 callGeneric(x=x, i=i, value = as.vector(value))
		     else
			 callGeneric(x=x, i=i, , value = as.vector(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
				value = "Matrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, i=i, j=j, value = as.vector(value)))


setReplaceMethod("[", signature(x = "Matrix", i = "missing", j = "ANY",
				value = "matrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, , j=j, value = c(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "missing",
				value = "matrix"),
		 function (x, i, j, ..., value)
		     if(nargs() == 3)
			 callGeneric(x=x, i=i, value = c(value))
		     else
			 callGeneric(x=x, i=i, , value = c(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
				value = "matrix"),
		 function (x, i, j, value)
		 callGeneric(x=x, i=i, j=j, value = c(value)))

##  M [ <lMatrix> ] <- value; used notably for x = "CsparseMatrix"  -------------------
.repl.i.lDMat <- function (x, i, j, ..., value)
{
    ## nA <- nargs()
    ## if(nA != 3) stop(gettextf("nargs() = %d should never happen; please report.", nA), domain=NA)
    ## else: nA == 3  i.e.,  M [ Lmat ] <- value
    ## x[i] <- value ; return(x)
    `[<-`(x, i=which(as.vector(i)), value=value)
}
setReplaceMethod("[", signature(x = "Matrix", i = "ldenseMatrix", j = "missing",
				value = "replValue"), .repl.i.lDMat)
setReplaceMethod("[", signature(x = "Matrix", i = "ndenseMatrix", j = "missing",
				value = "replValue"), .repl.i.lDMat)
.repl.i.lSMat <- function (x, i, j, ..., value)
{
    ## nA <- nargs()
    ## if(nA != 3) stop(gettextf("nargs() = %d should never happen; please report.", nA), domain=NA)
    ## else: nA == 3  i.e.,  M [ Lmat ] <- value
    ## x[i] <- value ; return(x)
    `[<-`(x, i=which(as(i, "sparseVector")), value=value)
}
setReplaceMethod("[", signature(x = "Matrix", i = "lsparseMatrix", j = "missing",
				value = "replValue"), .repl.i.lSMat)
setReplaceMethod("[", signature(x = "Matrix", i = "nsparseMatrix", j = "missing",
				value = "replValue"), .repl.i.lSMat)

## (ANY,ANY,ANY) is used when no `real method' is implemented :
setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
                                value = "ANY"),
	  function (x, i, j, value) {
              if(!is.atomic(value))
		  stop(gettextf(
		"RHS 'value' (class %s) matches 'ANY', but must match matrix class %s",
			       class(value), class(x)), domain=NA)
              else stop("not-yet-implemented 'Matrix[<-' method")
          })
