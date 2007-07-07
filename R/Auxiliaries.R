#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)


if(R.version$`svn rev` < 41863) ## use the fixed one
    ## will be hidden in namespace and can be removed when DEPENDS: R >= 2.5.1
    callGeneric <- function(...)
{
    frame <- sys.parent()
    envir <- parent.frame()
    call <- sys.call(frame)

    ## the  lines below this comment do what the previous version
    ## did in the expression fdef <- sys.function(frame)
    if(exists(".Generic", envir = envir, inherits = FALSE))
	fname <- get(".Generic", envir = envir)
    else { # in a local method (special arguments), or	an error
	## FIXME:  this depends on the .local mechanism, which should change
	if(identical(as.character(call[[1]]), ".local"))
	    call <- sys.call(sys.parent(2))
	fname <- as.character(call[[1]])
    }
    fdef <- get(fname, envir = envir)

    if(is.primitive(fdef)) {
        if(nargs() == 0)
            stop("'callGeneric' with a primitive needs explicit arguments (no formal args defined)")
        else {
            fname <- as.name(fname)
            call <- substitute(fname(...))
        }
    }
    else {
        env <- environment(fdef)
        if(!exists(".Generic", env, inherits = FALSE))
            stop("'callGeneric' must be called from a generic function or method")
        f <- get(".Generic", env, inherits = FALSE)
        fname <- as.name(f)
        if(nargs() == 0) {
            call[[1]] <- as.name(fname) # in case called from .local
            call <- match.call(fdef, call)
            anames <- names(call)
            matched <- !is.na(match(anames, names(formals(fdef))))
            for(i in seq_along(anames))
                if(matched[[i]])
                    call[[i]] <- as.name(anames[[i]])
        }
        else {
            call <- substitute(fname(...))
        }
    }
    eval(call, sys.frame(sys.parent()))
}



## Need to consider NAs ;  "== 0" even works for logical & complex:
is0  <- function(x) !is.na(x) & x == 0
isN0 <- function(x)  is.na(x) | x != 0
all0 <- function(x) !any(is.na(x)) && all(x == 0)

allTrue  <- function(x) !any(is.na(x)) && all(x)
allFalse <- function(x) !any(is.na(x)) && !any(x)

## maybe we should have this in base, maybe via an .Internal(paste0(.)) -> do_paste(.. op=2)
paste0 <- function(...) paste(..., sep = '')

## For %*% (M = Matrix; v = vector (double or integer {complex maybe?}):
.M.v <- function(x, y) callGeneric(x, as.matrix(y))
.v.M <- function(x, y) callGeneric(rbind(x), y)

.M.DN <- function(x) if(!is.null(dn <- dimnames(x))) dn else list(NULL,NULL)

.if.NULL <- function(x, orElse) if(!is.null(x)) x else orElse

.has.DN <- ## has non-trivial Dimnames slot?
    function(x) !identical(list(NULL,NULL), x@Dimnames)

.bail.out.1 <- function(fun, cl) {
    stop(gettextf('not-yet-implemented method for %s(<%s>).\n ->>  Ask the package authors to implement the missing feature.', fun, cl),
	 call. = FALSE)
}
.bail.out.2 <- function(fun, cl1, cl2) {
    stop(gettextf('not-yet-implemented method for %s(<%s>, <%s>).\n ->>  Ask the package authors to implement the missing feature.',
		  fun, cl1, cl2), call. = FALSE)
}

## This should be done in C and be exported by 'methods':  [FIXME - ask JMC ]
copyClass <- function(x, newCl, sNames =
		      intersect(slotNames(newCl), slotNames(x))) {
    r <- new(newCl)
    for(n in sNames)
	slot(r, n) <- slot(x, n)
    r
}

## chol() via "dpoMatrix"
cholMat <- function(x, pivot, ...) {
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
    is(x, "denseMatrix") &&
    ## unneeded(!): any("x" == slotNames(x)) &&
    length(x@x) < prod(dim(x))
}

emptyColnames <- function(x, msg.if.not.empty = FALSE)
{
    ## Useful for compact printing of (parts) of sparse matrices
    ## possibly	 dimnames(x) "==" NULL :
    dn <- dimnames(x)
    d <- dim(x)
    if(msg.if.not.empty && is.list(dn) && length(dn) >= 2 &&
       is.character(cn <- dn[[2]]) && any(cn != "")) {
	lc <- length(cn)
	message(sprintf("   [[ suppressing %d column names %s%s ]]", d[2],
			paste(sQuote(cn[1:min(3, lc)]), collapse = ", "),
			if(lc > 3) " ..." else ""))
    }
    dimnames(x) <- list(dn[[1]], rep("", d[2]))
    x
}

### TODO:  write in C and port to base (or 'utils') R
indTri <- function(n, upper = TRUE) {
    ## == which(upper.tri(diag(n)) or
    ##	  which(lower.tri(diag(n)) -- but much more efficiently for largish 'n'
    stopifnot(length(n) == 1, n == (n. <- as.integer(n)), (n <- n.) >= 0)
    if(n <= 2)
	return(if(n == 2) as.integer(if(upper) n+1 else n) else integer(0))
    ## First, compute the 'diff(.)'  fast.  Use integers
    one <- 1:1 ; two <- 2:2
    n1 <- n - one
    n2 <- n1 - one
    r <- rep.int(one, n*n1/two - one)
    r[cumsum(if(upper) 1:n2 else c(n1, if(n >= 4) n2:two))] <- if(upper) n:3 else 3:n
    ## now have "dliu" difference; revert to "liu":
    cumsum(c(if(upper) n+one else two, r))
}


prTriang <- function(x, digits = getOption("digits"),
                     maxp = getOption("max.print"),
		     justify = "none", right = TRUE)
{
    ## modeled along stats:::print.dist
    upper <- x@uplo == "U"

    m <- as(x, "matrix")
    cf <- format(m, digits = digits, justify = justify)
    if(upper)
	cf[row(cf) > col(cf)] <- "."
    else
	cf[row(cf) < col(cf)] <- "."
    print(cf, quote = FALSE, right = right, max = maxp)
    invisible(x)
}

prMatrix <- function(x, digits = getOption("digits"),
                     maxp = getOption("max.print")) {
    d <- dim(x)
    cl <- class(x)
    cat(sprintf('%d x %d Matrix of class "%s"\n', d[1], d[2], cl))
    if(prod(d) <= maxp) {
	if(is(x, "triangularMatrix"))
	    prTriang(x, digits = digits, maxp = maxp)
	else {
            print(as(x, "matrix"), digits = digits, max = maxp)
	}
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

nonFALSE <- function(x) {
    ## typically used for lMatrices:  (TRUE,NA,FALSE) |-> (TRUE,FALSE)
    if(any(ix <- is.na(x))) x[ix] <- TRUE
    x
}

nz.NA <- function(x, na.value) {
    ## Non-Zeros of x
    ## na.value: TRUE: NA's give TRUE, they are not 0
    ##             NA: NA's are not known ==> result := NA
    ##          FALSE: NA's give FALSE, could be 0
    stopifnot(is.logical(na.value) && length(na.value) == 1)
    if(is.na(na.value)) x != 0
    else  if(na.value)	isN0(x)
    else		x != 0 & !is.na(x)
}

## Number of "structural" non-zeros --- this is  nnzmax() in Matlab
##        of effectively  non-zero values =      nnz()     "   "

## Our nnzero() is like Matlab's nnz() -- but more sophisticated because of NAs
## This is now exported!
nnzero <- function(x, na.counted = NA) {
    ## na.counted: TRUE: NA's are counted, they are not 0
    ##		     NA: NA's are not known (0 or not) ==>  result := NA
    ##		  FALSE: NA's are omitted before counting
    cl <- class(x)
    ## speedup:
    cld <- getClassDef(cl)
    if(!extends(cld, "Matrix"))
	sum(nz.NA(x, na.counted))
    else { ## Matrix
       iSym <- extends(cld, "symmetricMatrix")
       if(extends(cld, "pMatrix")) # is "sparse" too
	   nrow(x)
       else if(extends(cld, "sparseMatrix")) {
	   nn <-
	       if(extends(cld, "nMatrix")) # <==> no 'x' slot
		   switch(.sp.class(cl),
			  "CsparseMatrix" = length(x@i),
			  "TsparseMatrix" = length(x@i),
			  "RsparseMatrix" = length(x@j))
	       else ## consider NAs in 'x' slot:
		   sum(nz.NA(x@x, na.counted))
	   if(iSym) (nn+nn - sum(nz.NA(diag(x), na.counted))) else nn
       }
       else if(extends(cld, "diagonalMatrix"))
	    sum(nz.NA(diag(x), na.counted))
       else { ## dense, not diagonal: Can use 'x' slot;
           nn <- sum(nz.NA(as_geClass(x, cl)@x, na.counted))
           if(iSym && length(x@x) < prod(dim(x))) ## packed symmetric
               ## n(n+1)/2  |--> n^2
               nn <- (nn + nn) - as.integer(sqrt(2*nn))
           nn
       }
    }
}

## For sparseness handling, return a
## 2-column (i,j) matrix of 0-based indices of non-zero entries:
non0ind <- function(x, classDef.x = getClassDef(class(x)))
{
    if(is.numeric(x))
	return(if((n <- length(x))) (0:(n-1))[isN0(x)] else integer(0))
    ## else

    stopifnot(extends(classDef.x, "sparseMatrix"))

    non0.i <- function(M, cM = class(M)) {
	if(extends(cM, "TsparseMatrix"))
	    return(unique(cbind(M@i,M@j)))
	if(extends(cM, "pMatrix"))
	    return(cbind(seq_len(nrow(M)), M@perm) - 1:1)
	## else: C* or R*
	isC <- extends(cM, "CsparseMatrix")
	.Call(compressed_non_0_ij, M, isC)
    }

    if(extends(classDef.x, "symmetricMatrix")) { # also get "other" triangle
	ij <- non0.i(x, classDef.x)
	notdiag <- ij[,1] != ij[,2]# but not the diagonals again
	rbind(ij, ij[notdiag, 2:1])
    }
    else if(extends(classDef.x, "triangularMatrix")) { # check for "U" diag
	if(x@diag == "U") {
	    i <- seq_len(dim(x)[1]) - 1:1
	    rbind(non0.i(x, classDef.x), cbind(i,i))
	} else non0.i(x, classDef.x)
    }
    else
	non0.i(x, classDef.x)
}

## nr= nrow: since  i in {0,1,.., nrow-1}  these are 1:1 "decimal" encodings:
## Further, these map to and from the usual "Fortran-indexing" (but 0-based)
encodeInd  <- function(ij,  nr) ij[,1] + ij[,2] * nr
encodeInd2 <- function(i,j, nr) i      +  j     * nr

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
	   ## do we need this for "logical" ones, there's no sum() there!
	   "ngTMatrix" = as(as(x, "ngCMatrix"), "ngTMatrix"),
	   "nsTMatrix" = as(as(x, "nsCMatrix"), "nsTMatrix"),
	   "ntTMatrix" = as(as(x, "ntCMatrix"), "ntTMatrix"),
	   ## otherwise:
	   stop("not yet implemented for class ", class.x))
}

## Note: maybe, using
## ----    xj <- .Call(Matrix_expand_pointers, x@p)
## would be slightly more efficient than as( <dgC> , "dgTMatrix")
## but really efficient would be to use only one .Call(.) for uniq(.) !

drop0 <- function(x, clx = c(class(x))) {
    if(!extends(clx, "CsparseMatrix"))
        clx <- sub(".Matrix$", "CMatrix", clx)
    ## FIXME: Csparse_drop should do this, but it
    ##	      drops triangularity and symmetry :
    ## .Call(Csparse_drop, as_CspClass(x, clx), 0)
    as_CspClass(.Call(Csparse_drop, as_CspClass(x, clx), 0.),
		clx)
}

uniq <- function(x) {
    if(is(x, "TsparseMatrix")) uniqTsparse(x) else
    if(is(x, "sparseMatrix")) drop0(x) else x
}

asTuniq <- function(x) {
    if(is(x, "TsparseMatrix")) uniqTsparse(x) else as(x,"TsparseMatrix")
}

## is 'x' a uniq Tsparse Matrix ?
is_not_uniqT <- function(x, nr = nrow(x))
    is.unsorted(x@j) || any(duplicated(encodeInd2(x@i, x@j, nr)))

## is 'x' a TsparseMatrix with no duplicated entries (to be *added* for uniq):
is_duplicatedT <- function(x, nr = nrow(x))
    any(duplicated(encodeInd2(x@i, x@j, nr)))


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

fixupDense <- function(m, from, cldm = getClassDef(class(m))) {
    if(extends(cldm, "triangularMatrix")) {
	m@uplo <- from@uplo
	m@diag <- from@diag
    } else if(extends(cldm, "symmetricMatrix")) {
	m@uplo <- from@uplo
    }
    m
}

## -> ./ldenseMatrix.R :
l2d_Matrix <- function(from, cl = class(from), cld = getClassDef(cl)) {
    ## stopifnot(is(from, "lMatrix"))
    fixupDense(new(sub("^l", "d", cl),
		   x = as.double(from@x),
		   Dim = from@Dim, Dimnames = from@Dimnames),
	       from, cld)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}

## -> ./ndenseMatrix.R :
n2d_Matrix <- function(from, cl = class(from), cld = getClassDef(cl)) {
    ## stopifnot(is(from, "nMatrix"))
    fixupDense(new(sub("^n", "d", cl), x = as.double(from@x),
		   Dim = from@Dim, Dimnames = from@Dimnames),
	       from, cld)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}
n2l_Matrix <- function(from, cl = class(from), cld = getClassDef(cl)) {
    fixupDense(new(sub("^n", "l", cl),
		   x = from@x, Dim = from@Dim, Dimnames = from@Dimnames),
	       from, cld)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}
## -> ./ddenseMatrix.R :
d2l_Matrix <- function(from, cl = class(from), cld = getClassDef(cl)) {
    fixupDense(new(sub("^d", "l", cl), x = as.logical(from@x),
                   Dim = from@Dim, Dimnames = from@Dimnames),
	       from, cld)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}

n2l_spMatrix <- function(from) {
    ## stopifnot(is(from, "nMatrix"))
    new(sub("^n", "l", class(from)),
        ##x = as.double(from@x),
        Dim = from@Dim, Dimnames = from@Dimnames)
}

tT2gT <- function(x, cl = class(x),
                  toClass = paste(substr(cl,1,1), "tTMatrix", sep=''),# "d" | "l"|"i"|"z"
                  cld = getClassDef(cl)) {
    ## coerce *tTMatrix to *gTMatrix {triangular -> general}
    d <- x@Dim
    if(uDiag <- x@diag == "U")	     # unit diagonal, need to add '1's
	uDiag <- (n <- d[1]) > 0
    if(extends(cld, "nMatrix")) # no 'x' slot
	new("ngTMatrix", Dim = d, Dimnames = x@Dimnames,
	    i = c(x@i, if(uDiag) 0:(n-1)),
	    j = c(x@j, if(uDiag) 0:(n-1)))
    else
	new(toClass, Dim = d, Dimnames = x@Dimnames,
	    i = c(x@i, if(uDiag) 0:(n-1)),
	    j = c(x@j, if(uDiag) 0:(n-1)),
	    x = c(x@x, if(uDiag) rep.int(1,n)))
}

gT2tT <- function(x, uplo, diag,
		  cl = class(x),
		  toClass = paste(substr(cl,1,1), "tTMatrix", sep=''),# d,l,i,z
		  cld = getClassDef(cl)) {
    ## coerce *gTMatrix to *tTMatrix {general -> triangular}
    i <- x@i
    j <- x@j
    sel <-
	if(uplo == "U") {
	    if(diag == "U") i < j else i <= j
	} else {
	    if(diag == "U") i > j else i >= j
	}
    i <- i[sel]
    j <- j[sel]
    if(extends(cld, "nMatrix")) # no 'x' slot
	new("ntTMatrix", i = i, j = j, uplo = uplo, diag = diag,
	    Dim = x@Dim, Dimnames = x@Dimnames)
    else
	new(toClass, i = i, j = j, uplo = uplo, diag = diag,
	    x = x@x[sel], Dim = x@Dim, Dimnames = x@Dimnames)
}

check.gT2tT <- function(from, cl = class(from),
			toClass = paste(substr(cl,1,1), "tTMatrix", sep=''),# d,l,i,z
			cld = getClassDef(cl)) {
    if(isTr <- isTriangular(from)) {
        force(cl)
	gT2tT(from, uplo = .if.NULL(attr(isTr, "kind"), "U"),
	      diag = "N", ## improve: also test for unit diagonal
	      cl = cl, toClass = toClass, cld = cld)
    } else stop("not a triangular matrix")
}

if(FALSE)# unused
l2d_meth <- function(x) {
    cl <- class(x)
    as(callGeneric(as(x, sub("^l", "d", cl))), cl)
}

## return "d" or "l" or "n" or "z"
.M.kind <- function(x, clx = class(x)) {
    ## 'clx': class() *or* class definition of x
    if(is.matrix(x) || is.atomic(x)) { ## 'old style' matrix or vector
	if     (is.numeric(x)) "d"
	else if(is.logical(x)) "l" ## FIXME ? "n" if no NA ??
	else if(is.complex(x)) "z"
	else stop("not yet implemented for matrix w/ typeof ", typeof(x))
    }
    else .M.kindC(clx)
}

## the same as .M.kind, but also knows "i"
.V.kind <- function(x, clx = class(x)) {
    ## 'clx': class() *or* class definition of x
    if(is.matrix(x) || is.atomic(x)) { ## 'old style' matrix or vector
	if     (is.integer(x)) "i"
	else if (is.numeric(x)) "d"
	else if(is.logical(x)) "l" ## FIXME ? "n" if no NA ??
	else if(is.complex(x)) "z"
	else stop("not yet implemented for matrix w/ typeof ", typeof(x))
    }
    else .M.kindC(clx)
}

.M.kindC <- function(clx) { ## 'clx': class() *or* classdefinition
    if(is.character(clx))		# < speedup: get it once
        clx <- getClassDef(clx)
    if(extends(clx, "sparseVector")) ## shortcut
	substr(as.character(clx@className), 1,1)
    else if(extends(clx, "dMatrix")) "d"
    else if(extends(clx, "nMatrix")) "n"
    else if(extends(clx, "lMatrix")) "l"
    else if(extends(clx, "zMatrix")) "z"
    else if(extends(clx, "pMatrix")) "n" # permutation -> pattern
    else stop(" not yet be implemented for ", clx@className)
}


## typically used as .type.kind[.M.kind(x)]:
.type.kind <- c("d" = "double",
		"i" = "integer",
		"l" = "logical",
		"n" = "logical",
		"z" = "complex")

.M.shape <- function(x, clx = class(x)) {
    ## 'clx': class() *or* class definition of x
    if(is.matrix(x)) { ## 'old style matrix'
	if     (isDiagonal  (x)) "d"
	else if(isTriangular(x)) "t"
	else if(isSymmetric (x)) "s"
	else "g" # general
    }
    else {
	if(is.character(clx)) # < speedup: get it once
	    clx <- getClassDef(clx)
	if(extends(clx, "diagonalMatrix"))  "d"
	else if(extends(clx, "triangularMatrix"))"t"
	else if(extends(clx, "symmetricMatrix")) "s"
	else "g"
    }
}


class2 <- function(cl, kind = "l", do.sub = TRUE) {
    ## Find "corresponding" class; since pos.def. matrices have no pendant:
    if	   (cl == "dpoMatrix") paste(kind, "syMatrix", sep='')
    else if(cl == "dppMatrix") paste(kind, "spMatrix", sep='')
    else if(do.sub) sub("^[a-z]", kind, cl)
    else cl
}

## see also as_geClass() below
geClass <- function(x) {
    if     (is(x, "dMatrix")) "dgeMatrix"
    else if(is(x, "lMatrix")) "lgeMatrix"
    else if(is(x, "nMatrix")) "ngeMatrix"
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
as_dense <- function(x, cld = if(isS4(x)) getClassDef(class(x))) {
    as(x, paste(.M.kind(x, cld), .dense.prefixes[.M.shape(x, cld)], "Matrix", sep=''))
}

.sp.class <- function(x) { ## find and return the "sparseness class"
    if(!is.character(x)) x <- class(x)
    for(cl in paste(c("C","T","R"), "sparseMatrix", sep=''))
	if(extends(x, cl))
	    return(cl)
    ## else (should rarely happen)
    as.character(NA)
}


### Goal: Eventually get rid of these --- want to foster coercions
### ----  *to* virtual classes whenever possible, i.e.
##  as(*, "CsparseMatrix"),  etc

## Here, getting the class definition and passing it, should be faster
as_Csparse <- function(x, cld = if(isS4(x)) getClassDef(class(x))) {
    as(x, paste(.M.kind(x, cld),
                .sparse.prefixes[.M.shape(x, cld)], "CMatrix", sep=''))
}

as_Csparse2 <- function(x, cld = if(isS4(x)) getClassDef(class(x))) {
    ## Csparse + U2N when needed
    sh <- .M.shape(x, cld)
    x <- as(x, paste(.M.kind(x, cld), .sparse.prefixes[sh], "CMatrix", sep=''))
    if(sh == "t") .Call(Csparse_diagU2N, x) else x
}


## 'cl'   : class() *or* class definition of from
as_gCsimpl2 <- function(from, cl = class(from))
    as(from, paste(.M.kind(from, cl), "gCMatrix", sep=''))
## to be used directly in setAs(.) needs one-argument-only  (from) :
as_gCsimpl <- function(from) as(from, paste(.M.kind(from), "gCMatrix", sep=''))

## slightly smarter:
as_Sp <- function(from, shape, cl = class(from)) {
    if(is.character(cl)) cl <- getClassDef(cl)
    as(from, paste(.M.kind(from, cl),
		   shape,
		   if(extends(cl, "TsparseMatrix")) "TMatrix" else "CMatrix",
		   sep=''))
}
as_gSparse <- function(from) as_Sp(from, "g", getClassDef(class(from)))
as_sSparse <- function(from) as_Sp(from, "s", getClassDef(class(from)))
as_tSparse <- function(from) as_Sp(from, "t", getClassDef(class(from)))

as_geSimpl2 <- function(from, cl = class(from))
    as(from, paste(.M.kind(from, cl), "geMatrix", sep=''))
## to be used directly in setAs(.) needs one-argument-only  (from) :
as_geSimpl <- function(from) as(from, paste(.M.kind(from), "geMatrix", sep=''))

## smarter, (but sometimes too smart!) compared to geClass() above:
as_geClass <- function(x, cl) {
    if(missing(cl)) return(as_geSimpl(x))
    ## else
    cld <- getClassDef(cl)
    if(extends(cld, "diagonalMatrix")  && isDiagonal(x))
	as(x, cl)
    else if(extends(cld, "symmetricMatrix") &&  isSymmetric(x)) {
        kind <- .M.kind(x, cld)
	as(x, class2(cl, kind, do.sub= kind != "d"))
    } else if(extends(cld, "triangularMatrix") && isTriangular(x))
	as(x, cl)
    else ## revert to
	as_geSimpl2(x, cld)
}

as_CspClass <- function(x, cl) {
    if (## diagonal is *not* sparse:
	##(extends(cl, "diagonalMatrix") && isDiagonal(x)) ||
	(extends(cl, "symmetricMatrix") && isSymmetric(x)) ||
	(extends(cl, "triangularMatrix")&& isTriangular(x)))
	as(x, cl)
    else if(is(x, "CsparseMatrix")) x
    else as(x, paste(.M.kind(x, cl), "gCMatrix", sep=''))
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
    if(is.na(upper)) {
	if(all0(object[lower.tri(object)]))
	    structure(TRUE, kind = "U")
	else if(all0(object[upper.tri(object)]))
	    structure(TRUE, kind = "L")
	else FALSE
    } else if(upper)
	all0(object[lower.tri(object)])
    else ## upper is FALSE
	all0(object[upper.tri(object)])
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
    ## "matrix" or "denseMatrix" (but not "diagonalMatrix")
    d <- dim(object)
    if(d[1] != (n <- d[2])) FALSE
    else if(is.matrix(object))
        ## requires that "vector-indexing" works for 'object' :
        all0(object[rep(c(FALSE, rep.int(TRUE,n)), length = n^2)])
    else ## "denseMatrix" -- packed or unpacked
        if(is(object, "generalMatrix")) # "dge", "lge", ...
            all0(object@x[rep(c(FALSE, rep.int(TRUE,n)), length = n^2)])
        else { ## "dense" but not {diag, general}, i.e. triangular or symmetric:
            ## -> has 'uplo'  differentiate between packed and unpacked

### .......... FIXME ...............

            packed <- isPacked(object)
            if(object@uplo == "U") {
            } else { ## uplo == "L"
            }

### very cheap workaround
	    all0(as.matrix(object)[rep(c(FALSE, rep.int(TRUE,n)), length = n^2)])
        }
}


## Purpose: Transform a *unit diagonal* sparse triangular matrix
##	into one with explicit diagonal entries '1'

## fast no-test version:
.diagU2N <- function(x, cl)
{
    if(extends(cl, "CsparseMatrix"))
	.Call(Csparse_diagU2N, x)
    else {
	kind <- .M.kind(x, cl)
	xT <- as(x, paste(kind, "gTMatrix", sep=''))
	## leave it as	T* - the caller can always coerce to C* if needed:
	new(paste(kind, "tTMatrix", sep=''), x = xT@x, i = xT@i, j = xT@j,
	    Dim = x@Dim, Dimnames = x@Dimnames, uplo = x@uplo, diag = "N")
    }
}

diagU2N <- function(x, cl = getClassDef(class(x)))
{
    if(extends(cl, "triangularMatrix") && x@diag == "U")
	.diagU2N(x, cl)
    else x
}



.as.dgC.0.factors <- function(x) {
    if(!is(x, "dgCMatrix"))
	as(x, "dgCMatrix") # will not have 'factors'
    else ## dgCMatrix
	if(!length(x@factors)) x else { x@factors <- list() ; x }
}


### Fast, much simplified version of tapply()
tapply1 <- function (X, INDEX, FUN = NULL, ..., simplify = TRUE) {
    sapply(unname(split(X, INDEX)), FUN, ...,
	   simplify = simplify, USE.NAMES = FALSE)
}

## tapply.x <- function (X, n, INDEX, FUN = NULL, ..., simplify = TRUE) {
##     tapply1(X, factor(INDEX, 0:(n-1)), FUN = FUN, ..., simplify = simplify)
## }

### MM: Unfortunately, these are still pretty slow for large sparse ...

sparsapply <- function(x, MARGIN, FUN, sparseResult = TRUE, ...)
{
    ## Purpose: "Sparse Apply": better utility than tapply1() for colSums() etc :
    ##    NOTE: Only applicable sum()-like where the "zeros do not count"
    ## ----------------------------------------------------------------------
    ## Arguments: x: sparseMatrix;  others as in *apply()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 16 May 2007
    stopifnot(MARGIN %in% 1:2)
    xi <- if(MARGIN == 1) x@i else x@j
    ui <- unique(xi)
    n <- x@Dim[MARGIN]
    ## FIXME: Here we assume 'FUN' to return  'numeric' !
    r <- if(sparseResult) new("dsparseVector", length = n) else numeric(n)
    r[ui + 1L] <- sapply(ui, function(i) FUN(x@x[xi == i], ...))
    r
}

sp.colMeans <- function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
{
    nr <- nrow(x)
    if(na.rm) ## use less than nrow(.) in case of NAs
	nr <- nr - sparsapply(x, 2, function(u) sum(is.na(u)),
			      sparseResult=sparseResult)
    sparsapply(x, 2, sum, sparseResult=sparseResult, na.rm=na.rm) / nr
}

sp.rowMeans <- function(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
{
    nc <- ncol(x)
    if(na.rm) ## use less than ncol(.) in case of NAs
	nc <- nc - sparsapply(x, 1, function(u) sum(is.na(u)),
			      sparseResult=sparseResult)
    sparsapply(x, 1, sum, sparseResult=sparseResult, na.rm=na.rm) / nc
}
