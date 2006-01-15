#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)

## For %*% (M = Matrix; v = vector (double or integer {complex maybe?}):
.M.v <- function(x, y) callGeneric(x, as.matrix(y))
.v.M <- function(x, y) callGeneric(rbind(x), y)

.has.DN <- ## has non-trivial Dimnames slot?
    function(x) !identical(list(NULL,NULL), x@Dimnames)

.bail.out.1 <- function(fun, cl) {
    stop(gettextf('not-yet-implemented method for %s(<%s>)', fun, cl),
         call. = FALSE)
}
.bail.out.2 <- function(fun, cl1, cl2) {
    stop(gettextf('not-yet-implemented method for %s(<%s>, <%s>)',
                  fun, cl1, cl2), call. = FALSE)
}

## chol() via "dpoMatrix"
cholMat <- function(x, pivot, LINPACK) {
    px <- as(x, "dpoMatrix")
    if (isTRUE(validObject(px, test=TRUE))) chol(px)
    else stop("'x' is not positive definite -- chol() undefined.")
}

dimCheck <- function(a, b) {
    da <- dim(a)
    db <- dim(b)
    if(any(da != db))
	stop(gettextf("Matrices must have same dimensions in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE)
    da
}

dimNamesCheck <- function(a, b) {
    ## assume dimCheck() has happened before
    nullDN <- list(NULL,NULL)
    h.a <- !identical(nullDN, dna <- dimnames(a))
    h.b <- !identical(nullDN, dnb <- dimnames(b))
    if(h.a || h.b) {
	if (!h.b) dna
	else if(!h.a) dnb
	else { ## both have non-trivial dimnames
	    r <- dna # "default" result
	    for(j in 1:2) {
		dn <- dnb[[j]]
		if(is.null(r[[j]]))
		    r[[j]] <- dn
		else if (!is.null(dn) && any(r[[j]] != dn))
		    warning(gettextf("dimnames [%d] mismatch in %s", j,
				     deparse(sys.call(sys.parent()))),
			    call. = FALSE)
	    }
	    r
	}
    }
    else
	nullDN
}

rowCheck <- function(a, b) {
    da <- dim(a)
    db <- dim(b)
    if(da[1] != db[1])
	stop(gettextf("Matrices must have same number of rows in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE)
    ## return the common nrow()
    da[1]
}

colCheck <- function(a, b) {
    da <- dim(a)
    db <- dim(b)
    if(da[2] != db[2])
	stop(gettextf("Matrices must have same number of columns in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE)
    ## return the common ncol()
    da[2]
}

emptyColnames <- function(x)
{
    ## Useful for compact printing of (parts) of sparse matrices
    ## possibly  dimnames(x) "==" NULL :
    dimnames(x) <- list(dimnames(x)[[1]], rep("", dim(x)[2]))
    x
}

prTriang <- function(x, digits = getOption("digits"),
                     justify = "none", right = TRUE)
{
    ## modeled along stats:::print.dist
    diag <- TRUE
    upper <- x@uplo == "U"

    m <- as(x, "matrix")
    cf <- format(m, digits = digits, justify = justify)
    if(upper)
        cf[row(cf) > col(cf)] <- "."
    else
        cf[row(cf) < col(cf)] <- "."
    print(cf, quote = FALSE, right = right)
    invisible(x)
}

prMatrix <- function(x, digits = getOption("digits")) {
    d <- dim(x)
    cl <- class(x)
    cat(sprintf('%d x %d Matrix of class "%s"\n', d[1], d[2], cl))
    maxp <- getOption("max.print")
    if(prod(d) <= maxp) {
        if(is(x, "triangularMatrix"))
            prTriang(x, digits = digits)
        else
            print(as(x, "matrix"), digits = digits)
    }
    else { ## d[1] > maxp / d[2] >= nr :
        m <- as(x, "matrix")
	nr <- maxp %/% d[2]
	n2 <- ceiling(nr / 2)
	print(head(m, max(1, n2)))
	cat("\n ..........\n\n")
	print(tail(m, max(1, nr - n2)))
    }
    ## DEBUG: cat("str(.):\n") ; str(x)
    invisible(x)# as print() S3 methods do
}

## For sparseness handling
non0ind <- function(x) {
    if(is.numeric(x))
        return(if((n <- length(x))) (0:(n-1))[x != 0] else integer(0))

    ## else return a (i,j) matrix of non-zero-indices

    stopifnot(is(x, "sparseMatrix"))
    if(is(x, "TsparseMatrix"))
	return(unique(cbind(x@i,x@j)))

    isCol <- function(M) any("i" == slotNames(M))
    .Call("compressed_non_0_ij", x, isCol(x), PACKAGE = "Matrix")
}

### There is a test on this in ../tests/dgTMatrix.R !
uniq <- function(x) {
    if(is(x, "TsparseMatrix")) {
        ## Purpose: produce a *unique* triplet representation:
        ##		by having (i,j) sorted and unique
        ## -----------------------------------------------------------
        ## The following is *not* efficient {but easy to program}:
        if(is(x, "dgTMatrix")) as(as(x, "dgCMatrix"), "dgTMatrix")
        else if(is(x, "lgTMatrix")) as(as(x, "lgCMatrix"), "lgTMatrix")
        else stop("not implemented for class", class(x))

    } else x      # not 'gT' ; i.e. "uniquely" represented in any case
}

if(FALSE) ## try an "efficient" version
uniq_gT <- function(x)
{
    ## Purpose: produce a *unique* triplet representation:
    ##		by having (i,j) sorted and unique
    ## ----------------------------------------------------------------------
    ## Arguments: a "gT" Matrix
    stopifnot(is(x, "gTMatrix"))
    if((n <- length(x@i)) == 0) return(x)
    ii <- order(x@i, x@j)
    if(any(ii != 1:n)) {
        x@i <- x@i[ii]
        x@j <- x@j[ii]
        x@x <- x@x[ii]
    }
    ij <- x@i + nrow(x) * x@j
    if(any(dup <- duplicated(ij))) {

    }
    ### We should use a .Call() based utility for this!

}

t_geMatrix <- function(x) {
    x@x <- as.vector(t(array(x@x, dim = x@Dim))) # no dimnames here
    x@Dim <- x@Dim[2:1]
    x@Dimnames <- x@Dimnames[2:1]
    ## FIXME: how to set factors?
    x
}

## t( [dl]trMatrix ) and  t( [dl]syMatrix ) :
t_trMatrix <- function(x) {
    x@x <- as.vector(t(as(x, "matrix")))
    x@Dim <- x@Dim[2:1]
    x@Dimnames <- x@Dimnames[2:1]
    x@uplo <- if (x@uplo == "U") "L" else "U"
    # and keep x@diag
    x
}

fixupDense <- function(m, from) {
    if(is(m, "triangularMatrix")) {
        m@uplo <- from@uplo
        m@diag <- from@diag
    } else if(is(m, "symmetricMatrix")) {
        m@uplo <- from@uplo
    }
    m
}

## -> ./ldenseMatrix.R :
l2d_Matrix <- function(from) {
    stopifnot(is(from, "lMatrix"))
    fixupDense(new(sub("^l", "d", class(from)),
                   x = as.double(from@x),
                   Dim = from@Dim, Dimnames = from@Dimnames,
                   factors = list()), ## FIXME: treat 'factors' smartly
               from)
}

if(FALSE)# unused
l2d_meth <- function(x) {
    cl <- class(x)
    as(callGeneric(as(x, sub("^l", "d", cl))), cl)
}

## -> ./ddenseMatrix.R :
d2l_Matrix <- function(from) {
    stopifnot(is(from, "dMatrix"))
    fixupDense(new(sub("^d", "l", class(from)),
                   Dim = from@Dim, Dimnames = from@Dimnames,
                   factors = list()), ## FIXME: treat 'factors' smartly
               from)
}


try_as <- function(x, classes, tryAnyway = FALSE) {
    if(!tryAnyway && !is(x, "Matrix"))
        return(x)
    ## else
    ok <- canCoerce(x, classes[1])
    while(!ok && length(classes <- classes[-1])) {
        ok <- canCoerce(x, classes[1])
    }
    if(ok) as(x, classes[1]) else x
}

if(paste(R.version$major, R.version$minor, sep=".") < "2.3")
    ## This will be in R 2.3.0
canCoerce <- function(object, Class) {
  ## Purpose:  test if 'object' is coercable to 'Class', i.e.,
  ##	       as(object, Class) will {typically} work
  ## ----------------------------------------------------------------------
  ## Author: John Chambers, Date:  6 Oct 2005
   is(object, Class) ||
   !is.null(selectMethod("coerce", c(class(object), Class),
			 optional = TRUE,
			 useInherited = c(from = TRUE, to = FALSE)))
}

.is.triangular <- function(object, upper = TRUE) {
    ## pretest: is it square?
    d <- dim(object)
    if(d[1] != d[2]) return(FALSE)
    ## else slower test
    if(!is.matrix(object))
        object <- as(object,"matrix")
    ## == 0 even works for logical & complex:
    if(upper)
        all(object[lower.tri(object)] == 0)
    else
        all(object[upper.tri(object)] == 0)
}

.is.diagonal <- function(object) {
    d <- dim(object)
    if(d[1] != (n <- d[2])) FALSE
    else all(object[rep(c(FALSE, rep.int(TRUE,n)), length = n^2)] == 0)
}
