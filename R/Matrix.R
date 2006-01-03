#### Toplevel ``virtual'' class "Matrix"

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

## slow "fall back" method {subclasses should have faster ones}:
setMethod("as.vector", signature(x = "Matrix", mode = "missing"),
          function(x) as.vector(as(x, "matrix")))


## Note that isSymmetric is *not* exported ---
### but also note that "base" eigen may get an isSymmetric() that *would* be exported!
setMethod("isSymmetric", signature(object = "symmetricMatrix"),
          function(object,tol) TRUE)
setMethod("isSymmetric", signature(object = "triangularMatrix"),
          ## FIXME: 'TRUE' if *diagonal*, i.e. return(isDiagonal(object))
          function(object,tol) FALSE)

setMethod("isDiagonal", signature(object = "sparseMatrix"),
          function(object) {
              gT <- as(object, "TsparseMatrix")
              all(gT@i == gT@j)
          })

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
    function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL,
              sparse = NULL)
{
    sparseDefault <- function(m)
        prod(dim(m)) > 2*sum(as(m, "matrix") != 0)

    i.M <- is(data, "Matrix")
    if(is.null(sparse) && (i.M || is(data, "matrix")))
        sparse <- sparseDefault(data)

    if (i.M) {
        sM <- is(data,"sparseMatrix")
        if((sparse && sM) || (!sparse && !sM))
            return(data)
        ## else : convert  dense <-> sparse -> at end
    }
    else if (!is.matrix(data)) { ## cut & paste from "base::matrix" :
	if (missing(nrow))
	    nrow <- ceiling(length(data)/ncol)
	else if (missing(ncol))
	    ncol <- ceiling(length(data)/nrow)
	data <- .Internal(matrix(data, nrow, ncol, byrow))
        if(is.null(sparse))
            sparse <- sparseDefault(data)
	dimnames(data) <- dimnames
    }

    ## 'data' is now a "matrix" or "Matrix"
    ## FIXME: consider it's type (logical,....)
    ## ctype <- substr(class(data), 1,1) # "d", "l", ...
    ## FIXME(2): check for symmetric / triangular / ...
### TODO: Compare with as.Matrix() and its tests in ./dgeMatrix.R
    if(sparse)
        as(data, "dgCMatrix")
    else
        as(data, "dgeMatrix")
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

## bail-out methods in order to get better error messages
setMethod("%*%", signature(x = "Matrix", y = "Matrix"),
	  function (x, y)
          stop(gettextf('not-yet-implemented method for <%s> %%*%% <%s>',
                        class(x), class(y))))

setMethod("crossprod", signature(x = "Matrix", y = "ANY"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))
setMethod("crossprod", signature(x = "ANY", y = "Matrix"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))

## There are special sparse methods; this is a "fall back":
setMethod("kronecker", signature(X = "Matrix", Y = "ANY",
                                 FUN = "ANY", make.dimnames = "ANY"),
          function(X, Y, FUN, make.dimnames, ...) {
              X <- as(X, "matrix") ; Matrix(callGeneric()) })
setMethod("kronecker", signature(X = "ANY", Y = "Matrix",
                                 FUN = "ANY", make.dimnames = "ANY"),
          function(X, Y, FUN, make.dimnames, ...) {
              Y <- as(Y, "matrix") ; Matrix(callGeneric()) })


setMethod("t", signature(x = "Matrix"),
	  function(x) .bail.out.1(.Generic, class(x)))

## Group Methods
setMethod("+", signature(e1 = "Matrix", e2 = "missing"), function(e1) e1)
## "fallback":
setMethod("-", signature(e1 = "Matrix", e2 = "missing"),
          function(e1) {
              warning("inefficient method used for \"- e1\"")
              0-e1
          })

## bail-outs:
setMethod("Compare", signature(e1 = "Matrix", e2 = "Matrix"),
          function(e1, e2) {
              d <- dimCheck(e1,e2)
              .bail.out.2(.Generic, class(e1), class(e2))
          })
setMethod("Compare", signature(e1 = "Matrix", e2 = "ANY"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))
setMethod("Compare", signature(e1 = "ANY", e2 = "Matrix"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))



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
	  function(x,i,j, drop) callGeneric(x, i=i, drop= TRUE))
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

## "FIXME:"
## How can we get at   A[ ij ]	where ij is (i,j) 2-column matrix?
##  and                A[ LL ]	where LL is a logical *vector*
## -> [.data.frame uses nargs() - can we do this in the *generic* ?


### "[<-" : -----------------

## x[] <- value :
setReplaceMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                                value = "index"),##  double/logical/...
	  function (x, value) { x@x <- value ; validObject(x); x })

## Otherwise (value is not "index"): bail out
setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
                                value = "ANY"),
	  function (x, i, j, value)
                 if(!is(value,"index"))
                 stop("RHS 'value' must be of class \"index\"")
                 else stop("not-yet-implemented 'Matrix[<-' method"))



## NOTE: the following only works for R 2.2.x (and later) ---
## ----  *and* 'Matrix' must have been *installed* by R >= 2.2.x

if(paste(R.version$major, R.version$minor, sep=".") >= "2.2") {

    ## The trivial methods :
    setMethod("cbind2", signature(x = "Matrix", y = "NULL"),
	      function(x, y) x)
    setMethod("cbind2", signature(x = "Matrix", y = "missing"),
	      function(x, y) x)
    setMethod("cbind2", signature(x = "NULL", y="Matrix"),
	      function(x, y) x)

    setMethod("rbind2", signature(x = "Matrix", y = "NULL"),
	      function(x, y) x)
    setMethod("rbind2", signature(x = "Matrix", y = "missing"),
	      function(x, y) x)
    setMethod("rbind2", signature(x = "NULL", y="Matrix"),
	      function(x, y) x)

    ## Makes sure one gets x decent error message for the unimplemented cases:
    setMethod("cbind2", signature(x = "Matrix", y = "Matrix"),
              function(x, y) {
                  rowCheck(x,y)
                  stop(gettextf("cbind2() method for (%s,%s) not-yet defined",
                                class(x), class(y)))
              })

    ## Use a working fall back {particularly useful for sparse}:
    ## FIXME: implement rbind2 via "cholmod" for C* and Tsparse ones
    setMethod("rbind2", signature(x = "Matrix", y = "Matrix"),
              function(x, y) {
                  colCheck(x,y)
                  t(cbind2(t(x), t(y)))
              })

}## R-2.2.x and newer
