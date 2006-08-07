#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)

## For %*% (M = Matrix; v = vector (double or integer {complex maybe?}):
.M.v <- function(x, y) callGeneric(x, as.matrix(y))
.v.M <- function(x, y) callGeneric(rbind(x), y)

.M.DN <- function(x) if(!is.null(dn <- dimnames(x))) dn else list(NULL,NULL)

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

## Note: !isPacked(.)  i.e. `full' still contains
## ----  "*sy" and "*tr" which have "undefined" lower or upper part
isPacked <- function(x)
{
    ## Is 'x' a packed (dense) matrix ?
    is(x,"Matrix") && !is.null(x@x) && length(x@x) < prod(dim(x))
}

emptyColnames <- function(x)
{
    ## Useful for compact printing of (parts) of sparse matrices
    ## possibly	 dimnames(x) "==" NULL :
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
    ## else
    stopifnot(is(x, "sparseMatrix"))
    ## return a 2-column (i,j) matrix of
    ## 0-based indices of non-zero entries  :
    non0.i <- function(M) {
	if(is(M, "TsparseMatrix"))
	    return(unique(cbind(M@i,M@j)))
	if(is(M, "pMatrix"))
	    return(cbind(seq(length=nrow(M)), M@perm) - 1:1)
	## else:
	isC <- any("i" == slotNames(M)) # is Csparse (not Rsparse)
	.Call(compressed_non_0_ij, M, isC)
    }

    if(is(x, "symmetricMatrix")) { # also get "other" triangle
	ij <- non0.i(x)
	notdiag <- ij[,1] != ij[,2]# but not the diagonals again
	rbind(ij, ij[notdiag, 2:1])
    }
    else
	non0.i(x)
}

## nr= nrow: since  i in {0,1,.., nrow-1}  these are 1:1 "decimal" encodings:
## Further, these map to and from the usual "Fortran-indexing" (but 0-based)
encodeInd <- function(ij, nr) ij[,1] + ij[,2] * nr
decodeInd <- function(code, nr) cbind(code %% nr, code %/% nr)

complementInd <- function(ij, dim)
{
    ## Purpose: Compute the complement of the 2-column 0-based ij-matrix
    ##		but as 1-based indices
    n <- prod(dim)
    if(n == 0) return(integer(0))
    ii <- 1:n
    ii[-(1 + encodeInd(ij, nr = dim[1]))]
}

unionInd <- function(ij1, ij2) unique(rbind(ij1, ij2))

intersectInd <- function(ij1, ij2, nrow) {
    ## from 2-column (i,j) matrices where i in {0,.., nrow-1},
    ## return only the *common* entries
    decodeInd(intersect(encodeInd(ij1, nrow),
			encodeInd(ij2, nrow)), nrow)
}

WhichintersectInd <- function(ij1, ij2, nrow) {
    ## from 2-column (i,j) matrices where i \in {0,.., nrow-1},
    ## find *where*  common entries are in ij1 & ij2
    m1 <- match(encodeInd(ij1, nrow), encodeInd(ij2, nrow))
    ni <- !is.na(m1)
    list(which(ni), m1[ni])
}


### There is a test on this in ../tests/dgTMatrix.R !

uniqTsparse <- function(x, class.x = c(class(x))) {
    ## Purpose: produce a *unique* triplet representation:
    ##		by having (i,j) sorted and unique
    ## -----------------------------------------------------------
    ## The following is not quite efficient {but easy to program,
    ## and as() are based on C code  (all of them?)
    ##
    ## FIXME: Do it fast for the case where 'x' is already 'uniq'

    switch(class.x,
	   "dgTMatrix" = as(as(x, "dgCMatrix"), "dgTMatrix"),
	   "dsTMatrix" = as(as(x, "dsCMatrix"), "dsTMatrix"),
	   "dtTMatrix" = as(as(x, "dtCMatrix"), "dtTMatrix"),
	   ## do we need this for "logical" ones, there's no sum() there!
	   "lgTMatrix" = as(as(x, "lgCMatrix"), "lgTMatrix"),
	   "lsTMatrix" = as(as(x, "lsCMatrix"), "lsTMatrix"),
	   "ltTMatrix" = as(as(x, "ltCMatrix"), "ltTMatrix"),
	   ## otherwise:
	   stop("not yet implemented for class ", clx))
}

## Note: maybe, using
## ----    xj <- .Call(Matrix_expand_pointers, x@p)
## would be slightly more efficient than as( <dgC> , "dgTMatrix")
## but really efficient would be to use only one .Call(.) for uniq(.) !

uniq <- function(x) {
    if(is(x, "TsparseMatrix")) uniqTsparse(x) else x
    ## else:  not 'Tsparse', i.e. "uniquely" represented in any case
}

if(FALSE) ## try an "efficient" version
uniq_gT <- function(x)
{
    ## Purpose: produce a *unique* triplet representation:
    ##		by having (i,j) sorted and unique
    ## ------------------------------------------------------------------
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
		   Dim = from@Dim, Dimnames = from@Dimnames),
	       from)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}

if(FALSE)# unused
l2d_meth <- function(x) {
    cl <- class(x)
    as(callGeneric(as(x, sub("^l", "d", cl))), cl)
}

## return "d" or "l" or "z"
.M.kind <- function(x, clx = class(x)) {
    if(is.matrix(x)) { ## 'old style matrix'
	if     (is.numeric(x)) "d"
	else if(is.logical(x)) "l"
	else if(is.complex(x)) "z"
	else stop("not yet implemented for matrix w/ typeof ", typeof(x))
    }
    else if(extends(clx, "dMatrix")) "d"
    else if(extends(clx, "lMatrix")) "l"
    else if(extends(clx, "zMatrix")) "z"
    else stop(" not yet be implemented for ", clx)
}

.M.shape <- function(x, clx = class(x)) {
    if(is.matrix(x)) { ## 'old style matrix'
	if     (isDiagonal  (x)) "d"
	else if(isTriangular(x)) "t"
	else if(isSymmetric (x)) "s"
	else "g" # general
    }
    else if(extends(clx, "diagonalMatrix"))  "d"
    else if(extends(clx, "triangularMatrix"))"t"
    else if(extends(clx, "symmetricMatrix")) "s"
    else "g"
}


class2 <- function(cl, kind = "l", do.sub = TRUE) {
    ## Find "corresponding" class; since pos.def. matrices have no pendant:
    if	   (cl == "dpoMatrix") paste(kind, "syMatrix", sep='')
    else if(cl == "dppMatrix") paste(kind, "spMatrix", sep='')
    else if(do.sub) sub("^d", kind, cl)
    else cl
}

geClass <- function(x) {
    if     (is(x, "dMatrix")) "dgeMatrix"
    else if(is(x, "lMatrix")) "lgeMatrix"
    else if(is(x, "zMatrix")) "zgeMatrix"
    else stop("general Matrix class not yet implemented for ",
	      class(x))
}

.dense.prefixes <- c("d" = "di",
                     "t" = "tr",
                     "s" = "sy",
                     "g" = "ge")

.sparse.prefixes <- c("d" = "t", ## map diagonal to triangular
                      "t" = "t",
                      "s" = "s",
                      "g" = "g")

## Used, e.g. after subsetting: Try to use specific class -- if feasible :
as_dense <- function(x) {
    as(x, paste(.M.kind(x), .dense.prefixes[.M.shape(x)], "Matrix", sep=''))
}

as_Csparse <- function(x) {
    as(x, paste(.M.kind(x), .sparse.prefixes[.M.shape(x)], "CMatrix", sep=''))
}
as_Rsparse <- function(x) {
    as(x, paste(.M.kind(x), .sparse.prefixes[.M.shape(x)], "RMatrix", sep=''))
}
as_Tsparse <- function(x) {
    as(x, paste(.M.kind(x), .sparse.prefixes[.M.shape(x)], "TMatrix", sep=''))
}

as_geClass <- function(x, cl) {
    if	   (extends(cl, "diagonalMatrix")  && isDiagonal(x))
	as(x, cl)
    else if(extends(cl, "symmetricMatrix") &&  isSymmetric(x))
	as(x, class2(cl, kind, do.sub= kind != "d"))
    else if(extends(cl, "triangularMatrix") && isTriangular(x))
	as(x, cl)
    else
	as(x, paste(.M.kind(x), "geMatrix", sep=''))
}

as_CspClass <- function(x, cl) {
    if ((extends(cl, "diagonalMatrix")	&& isDiagonal(x)) ||
	(extends(cl, "symmetricMatrix") &&  isSymmetric(x)) ||
	(extends(cl, "triangularMatrix")&& isTriangular(x)))
	as(x, cl)
    else as(x, paste(.M.kind(x), "gCMatrix", sep=''))
}


## -> ./ddenseMatrix.R :
d2l_Matrix <- function(from) {
    stopifnot(is(from, "dMatrix"))
    fixupDense(new(sub("^d", "l", class(from)), # no need for dClass2 here
		   Dim = from@Dim, Dimnames = from@Dimnames),
	       from)
    ## FIXME: treat 'factors' smartly {not for triangular!}
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


## For *dense* matrices
isTriMat <- function(object, upper = NA) {
    ## pretest: is it square?
    d <- dim(object)
    if(d[1] != d[2]) return(FALSE)
    ## else slower test
    if(!is.matrix(object))
	object <- as(object,"matrix")
    ## == 0 even works for logical & complex:
    if(is.na(upper)) {
	if(all(object[lower.tri(object)] == 0))
	    structure(TRUE, kind = "U")
	else if(all(object[upper.tri(object)] == 0))
	    structure(TRUE, kind = "L")
	else FALSE
    } else if(upper)
	all(object[lower.tri(object)] == 0)
    else ## upper is FALSE
	all(object[upper.tri(object)] == 0)
}

## For Csparse matrices
isTriC <- function(x, upper = NA) {
    ## pretest: is it square?
    d <- dim(x)
    if(d[1] != d[2]) return(FALSE)
    ## else
    if(d[1] == 0) return(TRUE)
    ni <- 1:d[2]
    ## the row indices split according to column:
    ilist <- split(x@i, factor(rep.int(ni, diff(x@p)), levels= ni))
    lil <- unlist(lapply(ilist, length), use.names = FALSE)
    if(any(lil == 0)) {
	pos <- lil > 0
	if(!any(pos)) ## matrix of all 0's
	    return(TRUE)
	ilist <- ilist[pos]
	ni <- ni[pos]
    }
    ni0 <- ni - 1:1 # '0-based ni'
    if(is.na(upper)) {
	if(all(sapply(ilist, max, USE.NAMES = FALSE) <= ni0))
	    structure(TRUE, kind = "U")
	else if(all(sapply(ilist, min, USE.NAMES = FALSE) >= ni0))
	    structure(TRUE, kind = "L")
	else FALSE
    } else if(upper) {
	all(sapply(ilist, max, USE.NAMES = FALSE) <= ni0)
    } else { ## 'lower'
	all(sapply(ilist, min, USE.NAMES = FALSE) >= ni0)
    }
}

.is.diagonal <- function(object) {
    d <- dim(object)
    if(d[1] != (n <- d[2])) FALSE
    else all(object[rep(c(FALSE, rep.int(TRUE,n)), length = n^2)] == 0)
}

diagU2N <- function(x)
{
    ## Purpose: Transform a *unit diagonal* triangular matrix
    ##	into one with explicit diagonal entries '1'
    xT <- as(x, "dgTMatrix")
    ## leave it as  T* - the caller can always coerce to C* if needed:
    new("dtTMatrix", x = xT@x, i = xT@i, j = xT@j, Dim = x@Dim,
	Dimnames = x@Dimnames, uplo = x@uplo, diag = "N")
}

## FIXME: this should probably be dropped / replaced by as_Csparse
.as.dgC.Fun <- function(x, na.rm = FALSE, dims = 1) {
    x <- as(x, "dgCMatrix")
    callGeneric()
}

.as.dgT.Fun <- function(x, na.rm = FALSE, dims = 1) {
    ## used e.g. inside colSums() etc methods
    x <- as(x, "dgTMatrix")
    callGeneric()
}


### Fast much simplified version of tapply()
tapply1 <- function (X, INDEX, FUN = NULL, ..., simplify = TRUE) {
    sapply(split(X, INDEX), FUN, ..., simplify = simplify, USE.NAMES = FALSE)
}

## tapply.x <- function (X, n, INDEX, FUN = NULL, ..., simplify = TRUE) {
##     tapply1(X, factor(INDEX, 0:(n-1)), FUN = FUN, ..., simplify = simplify)
## }

