#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)

## Always FALSE when building for CRAN ... documented since 2015
## NB: keep ../NAMESPACE synchronized
.Matrix.avoiding.as.matrix <- FALSE

## Some reverse dependencies built with Matrix <= 1.4.1 cache methods
## referring to objects that we no longer strictly need in the namespace.
## TRUE ensures that these "unused" objects continue to exist, so that
## the cached methods continue to work as before if called ...
## NB: keep Matrix_SupportingCachedMethods in ../src/Mutils.h synchronized
.Matrix.supporting.cached.methods <- TRUE

## These would be faster by a factor ~2 if done in C:
if(FALSE) {
## Need to consider NAs ;  "== 0" even works for logical & complex:
## Note that "!x" is faster than "x == 0", but does not (yet!) work for complex
is0  <- function(x) !is.na(x) & x == 0
isN0 <- function(x)  is.na(x) | x != 0
is1  <- function(x) !is.na(x) & x
} else {
## MJ: These seem more natural ... ?
is0  <- function(x) !(is.na(x) | x)
isN0 <- function(x)   is.na(x) | x
is1  <- function(x)  !is.na(x) & x == 1
isN1 <- function(x)   is.na(x) | x != 1
isT  <- function(x)  !is.na(x) &  x
isNT <- function(x)   is.na(x) | !x
}

##allFalse <- function(x) !any(x) && !any(is.na(x))## ~= all0, but allFalse(NULL) = TRUE w/warning
##all0 <- function(x) !any(is.na(x)) && all(!x) ## ~= allFalse
allFalse <- function(x) if(is.atomic(x)) .Call(R_all0, x) else !any(x) && !any(is.na(x))
all0     <- function(x) if(is.atomic(x)) .Call(R_all0, x) else all(!x) && !any(is.na(x))

##anyFalse <- function(x) isTRUE(any(!x))		 ## ~= any0
## any0 <- function(x) isTRUE(any(x == 0))	      ## ~= anyFalse
anyFalse <-
any0 <- function(x) if(is.atomic(x)) .Call(R_any0, x) else isTRUE(any(!x))

## These work "identically" for  1 ('==' TRUE)  and 0 ('==' FALSE)
##	(but give a warning for "double"  1 or 0)
## TODO: C versions of these would be faster
allTrue  <- function(x) all(x) && !anyNA(x)


## Note that mode(<integer>) = "numeric" -- as0(), as1() return "double"
## which is good *AS LONG AS* we do not really have i..Matrix integer matrices
as1 <- function(x, mod=mode(x))
    switch(mod, "integer"= 1L, "double"=, "numeric"= 1, "logical"= TRUE,
	   "complex"= 1+0i, stop(gettextf("invalid 'mod': %s", mod), domain = NA))
as0 <- function(x, mod=mode(x))
    switch(mod, "integer"= 0L, "double"=, "numeric"= 0, "logical"= FALSE,
	   "complex"= 0+0i, stop(gettextf("invalid 'mod': %s", mod), domain = NA))

##' equivalent to   extends(cl, classes[1]) || extends(cl, classes[2]) || ....
extends1of <- function(class, classes, ...) {
    if(is.character(class))
	class <- getClassDef(class[[1L]])
    for(c2 in classes)
	if(extends(class, c2, ...)) return(TRUE)
    ## otherwise return
    FALSE
}

## MJ: This .M.kind() is 3-6 times faster than the variants used previously
##     (see farther below).  It seems to be _slower_ than .M.kindC(cld, ext)
##     but only with 'cld' _and_ 'ext' already computed.  That would not be
##     the case if R_check_class_etc() were optimized as described here:
##     https://stat.ethz.ch/pipermail/r-devel/2022-August/081952.html

## "[nlidz]" for Matrix, sparseVector, logical, integer, double, complex 'x'
.M.kind  <- function(x) .Call(R_Matrix_kind, x, TRUE)  # integer -> "d"
.V.kind  <- function(x) .Call(R_Matrix_kind, x, FALSE) # integer -> "i"
## "[gstd]" for Matrix, sparseVector 'x'
.M.shape <- function(x) .Call(R_Matrix_shape, x)
## "[CRT]" for [CRT]sparseMatrix 'x', "" for other S4 'x'
.M.repr  <- function(x) .Call(R_Matrix_repr, x)
.isCRT   <- function(x) nzchar(.M.repr(x))

##' Should the matrix/Matrix  x  or a combination of x and y   be treated as  'sparse' ?
## sparseDefault <- function(x, y=NULL) {
##     if(is.null(y))
##         prod(dim(x)) > 2*sum(isN0(as(x, "matrix")))
##     else ## nrow / ncol ... differentiate  this would be for  rbind / cbind --> ./bind2.R
##         (nnzero(x) + nnzero(y)) * 2 < (nrow(x)+nrow(y)) * nc
## }

if(FALSE) {
sparseDefault <- function(x) prod(dim(x)) > 2 * sum(isN0(as(x, "matrix")))
} else {
## MJ: This version uses 'nnzero' which can be rather more efficient than
##     the above ... it also works for vectors without a 'dim' attribute
sparseDefault <- function(x) length(x) > 2 * nnzero(x, na.counted = TRUE)
}

##' return 'x' unless it is NULL where you'd use 'orElse'
`%||%` <- function(x, orElse) if(!is.null(x)) x else orElse

##'  not %in%  :
`%nin%` <- function (x, table) is.na(match(x, table))

nonTRUEoption <- function(ch) is.null(v <- getOption(ch)) || !isTRUE(v)


##' @title Check identical(i, 0:n) {or identical(i, 1:n) when Ostart is false}
##' @param i an integer vector, to be compared with 0:n or 1:n
##' @param n an integer number
##' @param Ostart logical indicating if comparison with 0:n or 1:n should be made
##' @return TRUE or FALSE
##' @author Martin Maechler
isSeq <- function(i, n, Ostart = TRUE) {
    ## FIXME: Port to C, use simple .Call() which is much faster notably in FALSE cases
    ##        and then *export* (and hence document)
    identical(i, if(Ostart) 0L:n else seq_len(n))
}


.bail.out.1 <- function(fun, cl) {
    stop(gettextf(
     'not-yet-implemented method for %s(<%s>).\n ->>  Ask the package authors to implement the missing feature.',
		  fun, cl[1L]), call. = FALSE, domain=NA)
}
.bail.out.2 <- function(fun, cl1, cl2) {
    stop(gettextf(
     'not-yet-implemented method for %s(<%s>, <%s>).\n ->>  Ask the package authors to implement the missing feature.',
		  fun, cl1[1L], cl2[1L]), call. = FALSE, domain=NA)
}

Matrix.verbose <- function()
    getOption("Matrix.verbose", .MatrixEnv[["verbose"]])

Matrix.warn <- function()
    getOption("Matrix.warn", .MatrixEnv[["warn"]])

## MJ: or maybe a general one?
if(FALSE) {
Matrix.getOption <- function(x, default = .Matrix.Env[[x]])
    getOption(paste0("Matrix.", x), default)
} ## MJ

if(FALSE) {
Matrix.msg <- function(..., .M.level = 1) {
    if(!is.null(v <- getOption("Matrix.verbose")) && v >= .M.level)
        message(...)
}
} else {
## MJ: a bit smarter
Matrix.msg <- function(..., .M.level = 1, call. = FALSE, domain = NULL) {
    if(Matrix.verbose() >= .M.level) {
        m <-
            if((w <- Matrix.warn()) < 1)
                function(..., call., domain) message(..., domain = domain)
            else if(w < 2)
                warning
            else stop
        m(..., call. = call., domain = domain)
    }
}
}

## currently unused; see also .solve.dsC.status() in ./solve.R
if(FALSE) {
Matrix.msg12 <- function(m1, m2, ...) {
    if(!is.null(v <- getOption("Matrix.verbose")) && v >= 1)
        message(if(v >= 2) m2 else m1, ...)
}
} ## unused

## we can set this to FALSE and possibly measure speedup:
.copyClass.check <- TRUE

## This should be done in C and be exported by 'methods':  [FIXME - ask JMC ]
copyClass <- function(x, newCl, sNames =
		      intersect(slotNames(newCl), slotNames(x)),
                      check = .copyClass.check)
{
    r <- new(newCl)
    ## Equivalent of
    ##   for(n in sNames) slot(r, n, check=check) <- slot(x, n)  :
    if(check) for(n in sNames) slot(r, n) <- slot(x, n)
    else for(n in sNames) # don't check, be fast
	## .Call("R_set_slot", r, n, slot(x,n), PACKAGE = "methods")
	## "ugly", but not using .Call(*, "methods")
	attr(r, n) <- attr(x, n)
    r
}

##' Return the (maybe super-)class of class 'cl'   from "Matrix", returning  character(0) if there is none.
##'
##' @title The Matrix (Super-) Class of a Class
##' @param cl string, class name
##' @param cld its class definition
##' @param ...Matrix if TRUE, the result must be of pattern "[dlniz]..Matrix"
##'     where the first letter "[dlniz]" denotes the content kind.
##' @param dropVirtual
##' @param ... other arguments are passed to .selectSuperClasses()
##' @return a character string
##' @author Martin Maechler, Date: 24 Mar 2009
MatrixClass <- function(cl, cld = getClassDef(cl),
			...Matrix = TRUE, dropVirtual = TRUE, ...)
{
    ## stopifnot(is.character(cl))
    ## Hmm, packageSlot(cl)  *can* be misleading --> use  cld@package  first:
    if(is.null(pkg <- cld@package)) {
	if(is.null(pkg <- packageSlot(cl))) return(character())
	## else we use 'pkg'
    }
    if(identical(pkg, "Matrix") &&
       (!...Matrix || (cl != "indMatrix" && identical(1L, grep("^[dlniz]..Matrix$", cl)))))
	cl
    else { ## possibly recursively
	r <- .selectSuperClasses(cld@contains, dropVirtual = dropVirtual,
				 namesOnly = TRUE, ...)
	if(length(r)) {
	    while(!length(r1 <- Recall(r[1], ...Matrix = ...Matrix, dropVirtual = dropVirtual))
		  && length(r) > 1) r <- r[-1]
	    r1
	} else r
    }
}

attrSlotNames <- function(m, factors = TRUE) {
    ## slotnames of Matrix objects which are *not* directly content related
    sn <- slotNames(m)
    sn[sn %nin% c("x","i","j","p", if(!factors) "factors")]
}

##' @param m
##' @return the slots of 'm' which are "attributes" of some kind.
attrSlots <- function(m, factors = TRUE)
    sapply(attrSlotNames(m, factors=factors),
           function(sn) slot(m, sn), simplify = FALSE)

##' @return { NULL | TRUE | character | list(.) }
attr.all_Mat <- function(target, current,
			 check.attributes = TRUE, factorsCheck = FALSE, ...) {
    msg <- if(check.attributes)
	all.equal(attrSlots(target,  factors=factorsCheck),
		  attrSlots(current, factors=factorsCheck),
		  check.attributes = TRUE, ...) ## else NULL
    if(!identical((c1 <- class(target)), (c2 <- class(current))))
	## list(): so we can easily check for this
	list(c(if(!isTRUE(msg)) msg, paste0("class(target) is ", c1, ", current is ", c2)))
    else msg
}

##' @return combination for  all.equal() functions in ./Matrix.R & ./sparseMatrix.R
.a.e.comb <- function(msg, r) {
    if((is.null(msg) || isTRUE(msg)) & (r.ok <- isTRUE(r))) TRUE
    else c(if(!isTRUE(msg)) msg, if(!r.ok) r)
}

identicalSlots <- function(x, y, slots, ...) {
    for (name in slots)
        if (!identical(slot(x, name), slot(y, name), ...))
            return(FALSE)
    TRUE
}

invPerm.R <- function(p) { p[p] <- seq_along(p) ; p }
## how much faster would this be in C? -- less than a factor of two?
invPerm <- function(p, zero.p = FALSE, zero.res = FALSE)
    .Call(inv_permutation, p, zero.p, zero.res)

##  sign( <permutation> ) == determinant( <pMatrix>)
signPerm <- function(p)
{
    ## Purpose: sign(<permutation>) via the cycles
    ## ----------------------------------------------------------------------
    ## Arguments: a permutation of 1:n
    ## ----------------------------------------------------------------------
    ## Author: Peter Dalgaard, 14 Apr 2008 // speedup: Martin Maechler 2008-04-16

    n <- length(p)
    x <- integer(n)
    ii <- seq_len(n)
    for (i in ii) {
	z <- ii[!x][1]             # index of first unmarked x[] entry
	if (is.na(z)) break
	repeat { ## mark x[] <- i  for those in i-th cycle
	    x[z] <- i
	    z <- p[z]
	    if (x[z]) break
	}
    }
    ## Now, table(x) gives the cycle lengths,
    ## where  split(seq_len(n), x)  would give the cycles list
    ## tabulate(x, i - 1L) is quite a bit faster than the equivalent
    ## table(x)
    clen <- tabulate(x, i - 1L)
    ## The sign is -1 (<==>  permutation is odd)  iff
    ## the cycle factorization contains an odd number of even-length cycles:
    1L - (sum(clen %% 2 == 0) %% 2L)*2L
}

checkDim <- function(da, db) {
    if(any(da != db))
	stop(gettextf("non-conformable matrix dimensions in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE, domain = NA)
    da
}

mmultDim <- function(d.a, d.b, type = 1L) {
    ## Return the 'dim' of the product indicated by 'type':
    ##     type 1:    a  %*%   b
    ##          2:  t(a) %*%   b    {crossprod}
    ##          3:    a  %*% t(b)  {tcrossprod}
    ## after asserting that ncol(<left operand>) == nrow(<right operand>)
    i.a <- 1L + (type != 2L)
    i.b <- 1L + (type == 3L)
    if(d.a[i.a] != d.b[i.b])
	stop(gettextf("non-conformable matrix dimensions in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE, domain = NA)
    c(d.a[-i.a], d.b[-i.b])
}

mmultDimnames <- function(dn.a, dn.b, type = 1L) {
    ## Return the 'dimnames' of the product indicated by 'type':
    ##     type 1:    a  %*%   b
    ##          2:  t(a) %*%   b    {crossprod}
    ##          3:    a  %*% t(b)  {tcrossprod}
    c(if(is.null(dn.a)) list(NULL) else dn.a[2L - (type != 2L)],
      if(is.null(dn.b)) list(NULL) else dn.b[2L - (type == 3L)])
}

## Still used in many places (for now):
dimCheck <- function(a, b) checkDim(dim(a), dim(b))
mmultCheck <- function(a, b, kind = 1L) mmultDim(dim(a), dim(b), type = kind)

##' Constructs "sensical" dimnames for something like  a + b ;
##' assume dimCheck() has happened before
##'
##' NOTA BENE:   R's  ?Arithmetic  says
##' ---------
##'>  For arrays (and an array result) the dimensions and dimnames are taken from
##'>  first argument if it is an array, otherwise the second.
##' but that's not quite correct:
##' The dimnames are taken from second *if* the first are NULL.
##'
##' @title Construct dimnames for  a  o  b
##' @param a matrix
##' @param b matrix
##' @param useFirst logical indicating if dimnames(a), the first, is taken, unless NULL
##' @param check logical indicating if a warning should be signalled for mismatches
##' @return a \code{\link{list}} of length two with dimnames
##' @author Martin Maechler
dimNamesCheck <- function(a, b, useFirst = TRUE, check = FALSE) {
    nullDN <- list(NULL,NULL)
    h.a <- !identical(nullDN, dna <- dimnames(a))
    h.b <- !identical(nullDN, dnb <- dimnames(b))
    if(h.a || h.b) {
        if(useFirst) {
            if(!h.a) dnb else dna
        } else {
            if (!h.b) dna
            else if(!h.a) dnb
            else { ## both have non-trivial dimnames
                r <- dna # "default" result
                for(j in 1:2) if(!is.null(dn <- dnb[[j]])) {
                    if(is.null(r[[j]]))
                        r[[j]] <- dn
                    else if(check && !identical(r[[j]], dn))
                        warning(gettextf("dimnames [%d] mismatch in %s", j,
                                         deparse(sys.call(sys.parent()))),
                                call. = FALSE, domain=NA)
                }
                r
            }
        }
    }
    else
	nullDN
}

##' valid Matrix-class @Dimnames slot {assuming only NULL needs to be transformed}
.M.DN <- function(x)
    dimnames(x) %||% list(NULL, NULL)

## NB: Now exported and documented in ../man/is.null.DN.Rd:
is.null.DN <- function(dn) {
    if(is.null(dn))
        return(TRUE)
    if(!is.null(names(dn)))
        names(dn) <- NULL
    identical(dn, list(NULL, NULL)) ||
        identical(dn, list(ch0 <- character(0L), NULL)) ||
        identical(dn, list(NULL, ch0)) ||
        identical(dn, list(ch0, ch0))
}

## Is 'dn' valid in the sense of 'validObject(<Matrix>)'?
## It is assumed that 'dim' is a length-2 non-negative integer vector.
validDN <- function(dn, dim)
    .Call(R_DimNames_validate, dn, dim)

validDim <- function(dim)
    .Call(R_Dim_validate, dim)

## Is 'dn' symmetric?
## This allows, e.g., list(NULL, nms), _unlike_ identical(dn[1], dn[2]),
## the definition used by base::isSymmetric.matrix ...
isSymmetricDN <- function(dn)
    .Call(R_DimNames_is_symmetric, dn)

symmDN <- function(dn)
    .Call(R_symmDN, dn)

##' @title Symmetrize dimnames
##' @param x A square matrix.
##' @return
##' \code{y} identical to \code{x} except with \code{dny <- dimnames(y)}
##' given by \code{rep(dimnames(x)[J], 2)} rather than \code{dimnames(x)}
##' (where \code{J} is 1 if \code{x} has row names but not column names,
##' and 2 otherwise) and thus satisfying \code{identical(dny[1], dny[2])}.
##' @author Martin Maechler and Mikael Jagan
symmetrizeDimnames <- function(x) {
    if(isS4(x)) # assuming is(x, "Matrix")
        `dimnames<-`(x, symmDN(x@Dimnames))
    else if(!is.null(dn <- dimnames(x))) # assuming list of length 2
        `dimnames<-`(x, symmDN(dn))
    else x
}

## MJ: no longer ... see above
if(FALSE) {
symmDN <- function(dn, col=TRUE, names=TRUE) {
    if(is.null(dn) || identical(dn[1L], dn[2L]))
	return(dn)
    J <-
        if(col) {
            if(is.null(dn[[2L]])) 1L else 2L
        } else { ## !col : row
            if(is.null(dn[[1L]])) 2L else 1L
        }

    if(!is.null(n <- names(dn))) {
	if(length(n) != 2)
	    stop("names(dimnames(<matrix>)) must be NULL or of length two")
	if(n[1L] != n[2L])
	    names(dn) <- if(names) n[c(J,J)] # else NULL
    }
    dn[c(J,J)]
}

symmetrizeDimnames <- function(x, col=TRUE, names=TRUE) {
    dimnames(x) <- symmDN(dimnames(x), col=col, names=names)
    x
}
} ## MJ

.M2symm <- function(from, ...) {
    if(isSymmetric(from, ...))
        forceSymmetric(from)
    else
        stop("matrix is not symmetric; consider forceSymmetric() or symmpart()")
}
..M2symm <- function(from) { # for setAs()
    if(isSymmetric(from))
        forceSymmetric(from)
    else
        stop("matrix is not symmetric; consider forceSymmetric() or symmpart()")
}

.M2tri <- function(from, ...) {
    if(!(it <- isTriangular(from, ...)))
        stop("matrix is not triangular; consider triu() or tril()")
    else if(attr(it, "kind") == "U")
        triu(from)
    else
        tril(from)
}
..M2tri <- function(from) { # for setAs()
    if(!(it <- isTriangular(from)))
        stop("matrix is not triangular; consider triu() or tril()")
    else if(attr(it, "kind") == "U")
        triu(from)
    else
        tril(from)
}

.M2diag <- function(from, check = TRUE) {
    if (check && !isDiagonal(from))
        stop("matrix is not diagonal; consider Diagonal(x=diag(.))")
    x <- diag(from, names = FALSE)
    cl <- switch(typeof(x),
                 double = {
                     unit <- allTrue(x == 1)
                     "ddiMatrix" },
                 integer = {
                     unit <- allTrue(x == 1L)
                     storage.mode(x) <- "double"
                     "ddiMatrix" },
                 ## integer = {
                 ##     unit <- allTrue(x == 1L)
                 ##     "idiMatrix" },
                 logical = {
                     unit <- allTrue(x)
                     "ldiMatrix" },
                 complex =
                     stop("complex \"diagonalMatrix\" not yet implemented"),
                 ## complex = {
                 ##     unit <- allTrue(x == 1+0i)
                 ##     "zdiMatrix" },
                 stop(gettextf("cannot coerce matrix of type \"%s\" to \"diagonalMatrix\"",
                               typeof(x)),
                      domain = NA))
    new(cl, Dim = dim(from), Dimnames = .M.DN(from),
        diag = if(unit) "U" else "N", x = if(unit) x[FALSE] else x)
}
..M2diag <- function(from) # for setAs()
    .M2diag(from, check = TRUE)

.dense2g <- function(from, kind = ".")
    .Call(R_dense_as_general, from, kind)

.sparse2g <- function(from)
    .Call(R_sparse_as_general, from)

.dense2m <- function(from, ndense = FALSE)
    .Call(R_dense_as_matrix, from, ndense)

.sparse2m <- function(from)
    .Call(R_sparse_as_matrix, from)

.diag2m <- function(from) {
    D <- base::diag(if(from@diag == "N") from@x else as1(from@x),
                    nrow = from@Dim[1L])
    if(!identical(dn <- from@Dimnames, list(NULL, NULL)))
        dimnames(D) <- dn
    D
}

.ind2m <- function(from) {
    P <- array(FALSE, d <- from@Dim)
    if((n <- d[1L]) > 0L)
        P[seq_len(n) + (from@perm - 1L) * as.double(n)] <- TRUE
    if(!identical(dn <- from@Dimnames, list(NULL, NULL)))
        dimnames(P) <- dn
    P
}

.ge2m <- function(from, ndense = FALSE)
    .Call(R_geMatrix_as_matrix, from, ndense)

.dense2v <- function(from, ndense = FALSE)
    .Call(R_dense_as_vector, from, ndense)

.sparse2v <- function(from)
    .Call(R_sparse_as_vector, from)

.diag2v <- function(from) {
    n <- from@Dim[1L]
    x <- from@x
    m <- mode(x)
    r <- vector(m, length = n^2)
    if(n > 0L)
        r[1 + 0:(n - 1L) * (n + 1)] <- if(from@diag == "N") x else as1(mod = m)
    r
}

.ind2v <- function(from) {
    x <- logical(prod(d <- from@Dim))
    if((n <- d[1L]) > 0L)
        x[seq_len(n) + (from@perm - 1L) * as.double(n)] <- TRUE
    x
}

.ge2v <- function(from, ndense = FALSE)
    .Call(R_geMatrix_as_vector, from, ndense)

.dense2kind <- function(from, kind)
    .Call(R_dense_as_kind, from, kind)

..dense2d <- function(from) # for setAs() but used widely:
    .Call(R_dense_as_kind, from, "d")

..dense2l <- function(from)
    .Call(R_dense_as_kind, from, "l")

..dense2n <- function(from)
    .Call(R_dense_as_kind, from, "n")

.sparse2kind <- function(from, kind, drop0 = FALSE)
    .Call(R_sparse_as_kind, from, kind, drop0)

..sparse2d <- function(from) # for setAs() but used widely:
    .Call(R_sparse_as_kind, from, "d", FALSE)

..sparse2l <- function(from)
    .Call(R_sparse_as_kind, from, "l", FALSE)

..sparse2n <- function(from)
    .Call(R_sparse_as_kind, from, "n", FALSE)

.diag2kind <- function(from, kind)
    .Call(R_diagonal_as_kind, from, kind)

..diag2d <- function(from)
    .Call(R_diagonal_as_kind, from, "d")

..diag2l <- function(from)
    .Call(R_diagonal_as_kind, from, "l")

.dense2sparse <- function(from, code)
    .Call(R_dense_as_sparse, from, code, NULL, NULL)

.diag2sparse <- function(from, code, uplo = "U", drop0 = TRUE)
    .Call(R_diagonal_as_sparse, from, code, uplo, drop0)

.sparse2dense <- function(from, packed = FALSE)
    .Call(R_sparse_as_dense, from, packed)

.diag2dense <- function(from, code, uplo = "U")
    .Call(R_diagonal_as_dense, from, code, uplo)

.m2ge <- function(from, kind = ".")
    .Call(R_matrix_as_dense, from, `substr<-`(".ge", 1L, 1L, kind), NULL, NULL)

.m2dense <- function(from, code, uplo = "U", diag = "N")
    .Call(R_matrix_as_dense, from, code, uplo, diag)

.m2dense.checking <- function(from, kind = ".", ...) {
    switch(typeof(from), logical =, integer =, double = NULL,
           stop(gettextf("matrix of invalid type \"%s\" to .m2dense.checking()",
                         typeof(from)),
                domain = NA))
    if(kind != ".") {
        ## These must happen before isSymmetric() call
        storage.mode(from) <-
            switch(kind, n =, l = "logical", d = "double",
                   stop(gettextf("invalid kind \"%s\" to .m2dense.checking()",
                                 kind),
                        domain = NA))
        if(kind == "n" && anyNA(from))
            from[is.na(from)] <- TRUE
    }
    if(isSymmetric(from, ...))
        .m2dense(from, paste0(kind, "sy"), "U", NULL)
    else if(it <- isTriangular(from))
        .m2dense(from, paste0(kind, "tr"), attr(it, "kind"), "N")
    else
        .m2dense(from, paste0(kind, "ge"), NULL, NULL)
}

.m2sparse <- function(from, code, uplo = "U", diag = "N")
    .Call(R_dense_as_sparse, from, code, uplo, diag)

.m2sparse.checking <- function(from, kind = ".", repr = "C", ...) {
    switch(typeof(from), logical =, integer =, double = NULL,
           stop(gettextf("matrix of invalid type \"%s\" to .m2sparse.checking()",
                         typeof(from)),
                domain = NA))
    if(kind != ".") {
        ## These must happen before isSymmetric() call
        storage.mode(from) <-
            switch(kind, n =, l = "logical", d = "double",
                   stop(gettextf("invalid kind \"%s\" to .m2sparse.checking()",
                                 kind),
                        domain = NA))
        if(kind == "n" && anyNA(from))
            from[is.na(from)] <- TRUE
    }
    if(isSymmetric(from, ...))
        .m2sparse(from, paste0(kind, "s", repr), "U", NULL)
    else if(it <- isTriangular(from))
        .m2sparse(from, paste0(kind, "t", repr), attr(it, "kind"), "N")
    else
        .m2sparse(from, paste0(kind, "g", repr), NULL, NULL)
}

.CR2T <- function(from) .Call(CRsparse_as_Tsparse, from)
.T2CR <- function(from, Csparse) .Call(Tsparse_as_CRsparse, from, Csparse)

.T2C  <- function(from) .Call(Tsparse_as_CRsparse, from,  TRUE)
.T2R  <- function(from) .Call(Tsparse_as_CRsparse, from, FALSE)

.tCR2RC <- function(from) .Call(tCRsparse_as_RCsparse, from)
.CR2RC <- function(from) {
    to <- .tCR2RC(.Call(R_sparse_transpose, from))
    if(.hasSlot(from, "factors"))
        to@factors <- from@factors
    to
}

drop0.notol <- function(x)
    .Call(R_sparse_drop0, x)

if(.Matrix.supporting.cached.methods) {
.C.2.R <- .CR2RC
.R.2.C <- .CR2RC
.R.2.T <- .CR2T
.T.2.C <- .T2C
.tC.2.R <- function(m, cl, clx) .tCR2RC(m)
.tR.2.C <- .tCR2RC
}

rowCheck <- function(a, b) {
    da <- dim(a)
    db <- dim(b)
    if(da[1] != db[1])
	stop(gettextf("Matrices must have same number of rows in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE, domain=NA)
    ## return the common nrow()
    da[1]
}

colCheck <- function(a, b) {
    da <- dim(a)
    db <- dim(b)
    if(da[2] != db[2])
	stop(gettextf("Matrices must have same number of columns in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE, domain=NA)
    ## return the common ncol()
    da[2]
}

## MJ: no longer needed ... can now do is(x, "packedMatrix")
if(FALSE) {
## Note: !isPacked(.)  i.e. `full' still contains
## ----  "*sy" and "*tr" which have "undefined" lower or upper part
isPacked <- function(x) {
    ## Is 'x' a packed (dense) matrix ?
    is(x, "denseMatrix") &&
    ## unneeded(!): any("x" == slotNames(x)) &&
    length(x@x) < prod(x@Dim)
}
} ## MJ

## Is 'x' a packed, dense matrix?
## MJ: Fast (not checking class) but "wrong" for n-by-n packedMatrix if n < 2,
##     e.g., .isPacked(new("dtpMatrix")) == FALSE ... FIXME ??
.isPacked <- function(x) length(x@x) < prod(x@Dim)

emptyColnames <- function(x, msg.if.not.empty = FALSE) {
    ## Useful for compact printing of (parts) of sparse matrices
    ## possibly	 dimnames(x) "==" NULL :
    if((nd <- length(d <- dim(x))) < 2L)
        return(x)
    nc <- d[2L]
    if(is.null(dn <- dimnames(x)))
        dn <- vector("list", nd)
    else if(msg.if.not.empty &&
            is.character(cn <- dn[[2L]]) &&
            any(nzchar(cn)))
        message(gettextf("  [[ suppressing %d column name%s %s ... ]]",
                         nc,
                         if(nc == 1L) "" else "s",
                         paste0(sQuote(if(nc <= 3L) cn else cn[1:3]),
                                collapse = ", ")),
                domain = NA)
    dn[[2L]] <- character(nc)
    dimnames(x) <- dn
    x
}

## MJ: unused
if(FALSE) {
## The i-th unit vector  e[1:n] with e[j] = \delta_{i,j}
## .E.i.log <- function(i,n)  i == (1:n)
## .E.i <- function(i,n)
##     r <- numeric(n)
##     r[i] <- 1.
##     r
## }
idiag <- function(n, p=n)
{
    ## Purpose: diag() returning  *integer*
    ## --------------------------------------------------------
    ## Author: Martin Maechler, Date:  8 Dec 2007, 23:13
    r <- matrix(0L, n,p)
    if ((m <- min(n, p)) > 0)
	r[1 + 0:(m - 1) * (n + 1)] <- 1L
    r
}

ldiag <- function(n, p=n)
{
    ## Purpose: diag() returning  *logical*
    r <- matrix(FALSE, n,p)
    if ((m <- min(n, p)) > 0)
	r[1 + 0:(m - 1) * (n + 1)] <- TRUE
    r
}
} ## MJ

indDiag <- function(n, upper = TRUE, packed = FALSE)
    .Call(R_index_diagonal, n, upper, packed)

indTri <- function(n, upper = TRUE, diag = FALSE, packed = FALSE)
    .Call(R_index_triangle, n, upper, diag, packed)

## MJ: now done in C, above, with options for packedMatrix
if(FALSE) {
## The indices of the diagonal entries of an  n x n matrix,  n >= 1
## i.e. indDiag(n) === which(diag(n) == 1)
indDiag <- function(n) cumsum(c(1L, rep.int(n+1L, n-1)))

### TODO:  write in C and port to base (or 'utils') R
### -----
### "Theory" behind this: /u/maechler/R/MM/MISC/lower-tri-w.o-matrix.R
## NB: also have "abIndex" version:  abIindTri() --> ./abIndex.R
## Size problem:  indTri(n) is of size ~ n^2/2  which may be too large!
indTri <- function(n, upper = TRUE, diag = FALSE) {
    ## Indices of (strict) upper/lower triangular part
    ## == which(upper.tri(diag(n), diag=diag) or
    ##	  which(lower.tri(diag(n), diag=diag) -- but
    ## much more efficiently for largish 'n'
    stopifnot(length(n) == 1, n == (n. <- as.integer(n)), (n <- n.) >= 0)
    if(n <= 2) {
        if(n == 0) return(integer(0))
        if(n == 1) return(if(diag) 1L else integer(0))
        ## else n == 2
        v <- if(upper) 3L else 2L
	return(if(diag) c(1L, v, 4L) else v)
    }

    ## n >= 3 [also for n == 2 && diag (==TRUE)] :

    ## First, compute the 'diff(.)' of the result [fast, using integers]
    n. <- if(diag) n else n - 1L
    n1 <- n. - 1L
    ## all '1' but a few
    r <- rep.int(1L, choose(n.+1, 2) - 1)
    tt <- if(diag) 2L else 3L
    r[cumsum(if(upper) 1:n1 else n.:2)] <- if(upper) n:tt else tt:n
    ## now have differences; revert to "original":
    cumsum(c(if(diag) 1L else if(upper) n+1L else 2L, r))
}
} ## MJ

prTriang <- function(x, digits = getOption("digits"),
                     maxp = getOption("max.print"),
		     justify = "none", right = TRUE) {
    ## modeled along stats:::print.dist
    upper <- x@uplo == "U"
    m <- as(x, "matrix")
    cf <- format(m, digits = digits, justify = justify)
    cf[if(upper) row(cf) > col(cf)
	else	 row(cf) < col(cf)] <- "."
    print(cf, quote = FALSE, right = right, max = maxp)
    invisible(x)
}

prMatrix <- function(x, digits = getOption("digits"),
                     maxp = getOption("max.print")) {
    d <- dim(x)
    cl <- class(x) ## cld <- getClassDef(cl)
    tri <- extends(cl, "triangularMatrix")
    xtra <- if(tri && x@diag == "U") " (unitriangular)" else ""
    cat(sprintf('%d x %d Matrix of class "%s"%s\n',
		d[1], d[2], cl, xtra))
    if(prod(d) <= maxp) {
	if(tri)
	    prTriang(x, digits = digits, maxp = maxp)
	else
	    print(as(x, "matrix"), digits = digits, max = maxp)
    }
    else { ## d[1] > maxp / d[2] >= nr :
	m <- as(x, "matrix")
	nr <- maxp %/% d[2]
	n2 <- ceiling(nr / 2)
	print(head(m, max(1, n2)))
	cat("\n ..........\n\n")
	print(tail(m, max(1, nr - n2)))
	cat("\n ..........\n\n")

    }
    ## DEBUG: cat("str(.):\n") ; str(x)
    invisible(x)# as print() S3 methods do
}

## MJ: no longer used
if(FALSE) {
nonFALSE <- function(x) {
    ## typically used for lMatrices:  (TRUE,NA,FALSE) |-> (TRUE,TRUE,FALSE)
    if(any(ix <- is.na(x))) x[ix] <- TRUE
    x
}

nz.NA <- function(x, na.value) {
    ## Non-Zeros of x
    ## na.value: TRUE: NA's give TRUE, they are not 0
    ##             NA: NA's are not known ==> result := NA
    ##          FALSE: NA's give FALSE, could be 0
    stopifnot(is.logical(na.value), length(na.value) == 1)
    if(is.na(na.value)) x != 0
    else  if(na.value)	isN0(x)
    else		x != 0 & !is.na(x)
}

### This assumes that e.g. the i-slot in Csparse is *not* over-allocated:
nnzSparse <- function(x, cl = class(x), cld = getClassDef(cl))
{
    ## Purpose: number of *stored* / structural non-zeros {NA's counted too}
    ## ----------------------------------------------------------------------
    ## Arguments: x sparseMatrix
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, 18 Apr 2008
    if(extends1of(cld, c("CsparseMatrix", "TsparseMatrix")))
	length(x@i)
    else if(extends(cld, "RsparseMatrix"))
	length(x@j)
    else if(extends(cld, "indMatrix"))	# is "sparse" too
	x@Dim[1]
    else stop(gettext("'x' must be \"sparseMatrix\""), domain=NA)
}
} ## MJ

## For sparseness handling, return a
## 2-column (i,j) matrix of 0-based indices of non-zero entries:

##' the workhorse for non0ind(.), but occasionally used directly
non0.i <- function(M, cM = class(M), uniqT = TRUE) {
    cld <- getClassDef(cM)
    if(extends(cld, "CsparseMatrix"))
        .Call(compressed_non_0_ij, M, TRUE)
    else if(extends(cld, "RsparseMatrix"))
        .Call(compressed_non_0_ij, M, FALSE)
    else if(extends(cld, "TsparseMatrix")) {
        if(uniqT && is_not_uniqT(M))
	    .Call(compressed_non_0_ij, .T2C(M), TRUE)
	else cbind(M@i, M@j, deparse.level = 0L)
    } else if(extends(cld, "diagonalMatrix")) {
        i <- seq.int(from = 0L, length.out = M@Dim[1L])
        if(M@diag == "N")
	    i <- i[isN0(M@x)]
	cbind(i, i, deparse.level = 0L)
    } else if(extends(cld, "indMatrix")) {
        i <- seq.int(from = 0L, length.out = M@Dim[1L])
        cbind(i, M@perm - 1L, deparse.level = 0L)
    } else
        stop(gettextf("non0.i() not yet implemented for class %s", dQuote(cM)),
             domain = NA)
}

##' the "more versatile / user" function (still not exported):
non0ind <- function(x, cld = getClassDef(class(x)),
		    uniqT = TRUE, xtendSymm = TRUE, check.Udiag = TRUE)
{
    if(is.numeric(x)) {
        n <- length(x)
        return(if(n == 0L)
                   integer(0L)
               else if(is.matrix(x))
                   arrayInd(seq_len(n)[isN0(x)], dim(x)) - 1L
               else (0:(n-1L))[isN0(x)])
    }
    stopifnot(extends(cld, "sparseMatrix"))

    ij <- non0.i(x, cld, uniqT = uniqT)
    if(xtendSymm && extends(cld, "symmetricMatrix")) {
        ## also get "other" triangle, but not the diagonal again
	notdiag <- ij[, 1L] != ij[, 2L]
	rbind(ij, ij[notdiag, 2:1], deparse.level = 0L)
    } else if(check.Udiag && extends(cld, "triangularMatrix") &&
              x@diag == "U") {
        i <- seq.int(from = 0L, length.out = x@Dim[1L])
        rbind(ij, cbind(i, i, deparse.level = 0L), deparse.level = 0L)
    } else ij
}

if(FALSE) { ## -- now have  .Call(m_encodeInd, ...) etc :

## nr= nrow: since  i in {0,1,.., nrow-1}  these are 1L "decimal" encodings:
## Further, these map to and from the usual "Fortran-indexing" (but 0-based)
encodeInd  <- function(ij, di) {
    stopifnot(length(di) == 2)
    nr <- di[1L]
    ## __ care against integer overflow __
    if(prod(di) >= .Machine$integer.max) nr <- as.double(nr)
    ij[,1] + ij[,2] * nr
}
encodeInd2 <- function(i,j, di) {
    stopifnot(length(di) == 2)
    nr <- di[1L]
    ## __ care against integer overflow __
    if(prod(di) >= .Machine$integer.max) nr <- as.double(nr)
    i +  j * nr
}

} else {

##' Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
##'
##' @param ij 2-column integer matrix
##' @param dim dim(.), i.e. length 2 integer vector
##' @param checkBnds logical indicating  0 <= ij[,k] < di[k]  need to be checked.
##'
##' @return encoded index; integer if prod(dim) is small; double otherwise
encodeInd <- function(ij, dim, orig1=FALSE, checkBnds=TRUE)
    .Call(m_encodeInd, ij, dim, orig1, checkBnds)
## --> in ../src/Mutils.c :  m_encodeInd(ij, di, orig_1, chk_bnds)
##                           ~~~~~~~~~~~

##' Here, 1-based indices (i,j) are default:
encodeInd2 <- function(i, j, dim, orig1=TRUE, checkBnds=TRUE)
    .Call(m_encodeInd2, i,j, dim, orig1, checkBnds)

##' Decode "encoded" (i,j) indices back to  cbind(i,j)
##' This is the inverse of encodeInd(.)
##'
##' @title Decode "Encoded" (i,j) Indices
##' @param code integer in 0:((n x m - 1)  <==> encodeInd(.) result
##' @param nr the number of rows
##' @return
##' @author Martin Maechler
decodeInd <- function(code, nr)
    cbind(as.integer(code %% nr), as.integer(code %/% nr),
          deparse.level = 0L)

complementInd <- function(ij, dim, orig1=FALSE, checkBnds=FALSE) {
    ## Purpose: Compute the complement of the 2-column 0-based ij-matrix
    ##		but as 1-based indices
    n <- prod(dim)
    if(n == 0L)
        return(integer(0L))
    seq_len(n)[-(1L + .Call(m_encodeInd, ij, dim, orig1, checkBnds))]
}

unionInd <- function(ij1, ij2)
    unique(rbind(ij1, ij2))

intersectInd <- function(ij1, ij2, di, orig1=FALSE, checkBnds=FALSE) {
    ## from 2-column (i,j) matrices where i in {0,.., nrow-1},
    ## return only the *common* entries
    decodeInd(intersect(.Call(m_encodeInd, ij1, di, orig1, checkBnds),
			.Call(m_encodeInd, ij2, di, orig1, checkBnds)),
              nr = di[1L])
}

WhichintersectInd <- function(ij1, ij2, di, orig1=FALSE, checkBnds=FALSE) {
    ## from 2-column (i,j) matrices where i \in {0,.., nrow-1},
    ## find *where*  common entries are in ij1 & ij2
    m1 <- match(.Call(m_encodeInd, ij1, di, orig1, checkBnds),
                .Call(m_encodeInd, ij2, di, orig1, checkBnds))
    ni <- !is.na(m1)
    list(which(ni), m1[ni])
}

### There is a test on this in ../tests/dgTMatrix.R !

uniqTsparse <- function(x, class.x = c(class(x))) {
    ## Purpose: produce a *unique* triplet representation:
    ##		by having (i,j) sorted and unique
    ## -----------------------------------------------------------
    ## The following is not quite efficient, but easy to program,
    ## and much based on C code
    ##
    ## TODO: faster for the case where 'x' is already 'uniq'?  if(anyDuplicatedT(.))
    if(!extends(class.x, "TsparseMatrix"))
	stop(gettextf("not yet implemented for class \"%s\"", dQuote(class.x)),
	     domain = NA)
    .CR2T(.T2C(x))
}

##' non-exported version with*OUT* check -- called often only  if(anyDuplicatedT(.))
.uniqTsparse <- function(x) .CR2T(.T2C(x))

## Note: maybe, using
## ----    xj <- .Call(Matrix_expand_pointers, x@p)
## would be slightly more efficient than as( <dgC> , "dgTMatrix")
## but really efficient would be to use only one .Call(.) for uniq(.) !

drop0 <- function(x, tol = 0, is.Csparse = NA) {
    ## ## MJ: My R_sparse_drop0() handles all [CRT]sparseMatrix without
    ## ##     coercion, but we should assess how many packages depend on
    ## ##     > is(drop0(...), "CsparseMatrix")
    ## ##     before going ahead and using it ...
    ## if(tol <= 0)
    ##     return(.Call(R_sparse_drop0, x))
    if(is.na(is.Csparse))
        is.Csparse <- is(x, "CsparseMatrix")
    .Call(Csparse_drop, if(is.Csparse) x else as(x, "CsparseMatrix"), tol)
}

uniq <- function(x) {
    cld <- getClassDef(class(x))
    if(extends(cld, "TsparseMatrix"))
        uniqTsparse(x)
    else if(extends(cld, "sparseMatrix"))
        drop0(x)
    else x
}

asTuniq <- function(x) {
    if(is(x, "TsparseMatrix"))
        uniqTsparse(x)
    else as(x, "TsparseMatrix")
}

## is 'x' a uniq Tsparse Matrix ?
is_not_uniqT <- function(x, di = dim(x))
    is.unsorted(x@j) || anyDuplicatedT(x, di)

## is 'x' a TsparseMatrix with duplicated entries (to be *added* for uniq):
anyDuplicatedT <- function(x, di = dim(x))
    anyDuplicated(.Call(m_encodeInd2, x@i, x@j, di, FALSE, FALSE))

}

## MJ: no longer needed ... replacement in ./unpackedMatrix.R
if(FALSE) {
t_geMatrix <- function(x) {
    x@x <- as.vector(t(array(x@x, dim = x@Dim))) # no dimnames here
    x@Dim <- x@Dim[2:1]
    x@Dimnames <- x@Dimnames[2:1]
    x@factors <- list() ## FIXME -- do better, e.g., for "LU"?
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
} ## MJ

## MJ: no longer needed ... replaced by .dense2kind() above
if(FALSE) {
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
l2d_Matrix <- function(from, cl = MatrixClass(class(from)), cld = getClassDef(cl)) {
    ## stopifnot(is(from, "lMatrix"))
    fixupDense(new(sub("^l", "d", cl),
		   x = as.double(from@x),
		   Dim = from@Dim, Dimnames = from@Dimnames),
	       from, cld)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}

## -> ./ndenseMatrix.R :
n2d_Matrix <- function(from, cl = MatrixClass(class(from)), cld = getClassDef(cl)) {
    ## stopifnot(is(from, "nMatrix"))
    fixupDense(new(sub("^n", "d", cl), x = as.double(from@x),
		   Dim = from@Dim, Dimnames = from@Dimnames),
	       from, cld)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}
n2l_Matrix <- function(from, cl = MatrixClass(class(from)), cld = getClassDef(cl)) {
    fixupDense(new(sub("^n", "l", cl),
		   x = from@x, Dim = from@Dim, Dimnames = from@Dimnames),
	       from, cld)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}
## -> ./ddenseMatrix.R :
d2l_Matrix <- function(from, cl = MatrixClass(class(from)), cld = getClassDef(cl)) {
    fixupDense(new(sub("^d", "l", cl), x = as.logical(from@x),
                   Dim = from@Dim, Dimnames = from@Dimnames),
	       from, cld)
    ## FIXME: treat 'factors' smartly {not for triangular!}
}

n2l_spMatrix <- function(from) {
    ## stopifnot(is(from, "nMatrix"))
    new(sub("^n", "l", MatrixClass(class(from))),
        ##x = as.double(from@x),
        Dim = from@Dim, Dimnames = from@Dimnames)
}
} ## MJ

## MJ: no longer needed
if(FALSE) {
.R.2.C  <- function(from) .Call(R_to_CMatrix, from)
.C.2.R  <- function(from)
    .tC.2.R(.Call(Csparse_transpose, from, is(from, "triangularMatrix")))
## slightly less efficient than above, but preserves symmetry correctly
.viaC.2.R <- function(from) .tC.2.R(as(t(from), "CsparseMatrix"))
## .R.2.T() fails on 32bit--enable-R-shlib with segfault {Kurt}
.R.2.T  <- function(from) .Call(compressed_to_TMatrix, from, FALSE)

## in ../src/Tsparse.c :  |-> cholmod_T -> cholmod_C -> chm_sparse_to_SEXP
## adjusted for triangular matrices not represented in cholmod
.T.2.C <- function(from) {
    to <- .Call(Tsparse_to_Csparse, from, is(from, "triangularMatrix"))
    if(.hasSlot(from, "factors"))
        to@factors <- from@factors
    to
}

.T.2.R <- function(from) {
    to <- .tC.2.R(.Call(Tsparse_to_Csparse, t(from), is(from, "triangularMatrix")))
    if(.hasSlot(from, "factors"))
        to@factors <- from@factors
    to
}

.T2Cmat <- function(from, isTri = is(from, "triangularMatrix"))
    .Call(Tsparse_to_Csparse, from, isTri)

tT2gT <- function(x, cl = class(x), toClass, cld = getClassDef(cl)) {
    ## coerce *tTMatrix to *gTMatrix {triangular -> general}
    d <- x@Dim
    if(uDiag <- x@diag == "U")	     # unit diagonal, need to add '1's
        uDiag <- (n <- d[1]) > 0
    if(missing(toClass)) {
        do.n <- extends(cld, "nMatrix")
        toKind <- if(do.n) "n" else substr(MatrixClass(cl), 1,1) # "d" | "l"|"i"|"z"
        toClass <- paste0(toKind, "gTMatrix")
    } else {
        do.n <- extends(toClass, "nMatrix")
        toKind <- if(do.n) "n" else substr(toClass, 1,1)
    }

    if(do.n) ## no 'x' slot
	new(toClass, # == "ngTMatrix"
            Dim = d, Dimnames = x@Dimnames,
	    i = c(x@i, if(uDiag) 0:(n-1)),
	    j = c(x@j, if(uDiag) 0:(n-1)))
    else
	new(toClass, Dim = d, Dimnames = x@Dimnames,
	    i = c(x@i, if(uDiag) 0:(n-1)),
	    j = c(x@j, if(uDiag) 0:(n-1)),
	    x = c(x@x, if(uDiag) rep.int(if(toKind == "l") TRUE else 1, n)))
}
## __TODO__
## Hack for the above, possibly considerably faster:
## Just *modify* the 'x' object , using attr(x, "class') <- toClass

## Fast very special one ../src/Tsparse.c -- as_cholmod_triplet()
## in ../src/chm_common.c
## 'x' *must* inherit from TsparseMatrix!
.gT2tC <- function(x, uplo, diag="N") .Call(Tsparse_to_tCsparse, x, uplo, diag)
## Ditto in ../src/Csparse.c :
.gC2tT <- function(x, uplo, diag="N") .Call(Csparse_to_tTsparse, x, uplo, diag)
.gC2tC <- function(x, uplo, diag="N") .Call(Csparse_to_tCsparse, x, uplo, diag)

.gC2sC <- function(x, uplo) .Call(Csparse_general_to_symmetric, x, uplo, TRUE)

gT2tT <- function(x, uplo, diag, toClass,
		  do.n = extends(toClass, "nMatrix")) {
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
    if(do.n) ## no 'x' slot
	new("ntTMatrix", i = i, j = j, uplo = uplo, diag = diag,
	    Dim = x@Dim, Dimnames = x@Dimnames)
    else
	new(toClass, i = i, j = j, uplo = uplo, diag = diag,
	    x = x@x[sel], Dim = x@Dim, Dimnames = x@Dimnames)
}

check.gT2tT <- function(from, toClass, do.n = extends(toClass, "nMatrix")) {
    if(isTr <- isTriangular(from)) {
	gT2tT(from, uplo = attr(isTr, "kind") %||% "U",
	      diag = "N", ## improve: also test for unit diagonal
	      toClass = toClass, do.n = do.n)
    } else stop("not a triangular matrix")
}

gT2sT <- function(x, toClass, do.n = extends(toClass, "nMatrix")) {
    upper <- x@i <= x@j
    i <- x@i[upper]
    j <- x@j[upper]
    if(do.n) ## no 'x' slot
	new("nsTMatrix", Dim = x@Dim, Dimnames = x@Dimnames,
	    i = i, j = j, uplo = "U")
    else
	new(toClass, Dim = x@Dim, Dimnames = x@Dimnames,
	    i = i, j = j, x = x@x[upper], uplo = "U")
}

check.gT2sT <- function(x, toClass, do.n = extends(toClass, "nMatrix"))
{
    if(isSymmetric(x))
        gT2sT(x, toClass, do.n)
    else
	stop("not a symmetric matrix; consider forceSymmetric() or symmpart()")
}
} ## MJ

## MJ: no longer needed ... replaced above
if(FALSE) {
## return "d" or "l" or "n" or "z"
.M.kind <- function(x, clx = class(x)) {
    ## 'clx': class() *or* class definition of x
    if(is.matrix(x) || is.atomic(x)) { ## 'old style' matrix or vector
	if     (is.numeric(x)) "d" ## also for integer: see .V.kind(), .M.kindC()
	else if(is.logical(x)) "l" ## FIXME ? "n" if no NA ??
	else if(is.complex(x)) "z"
	else stop(gettextf("not yet implemented for matrix with typeof %s",
			   typeof(x)), domain = NA)
    }
    else .M.kindC(clx)
}

##' *V*ector kind (as .M.kind, but also knows "i")
.V.kind <- function(x, clx = class(x)) {
    ## 'clx': class() *or* class definition of x
    if(is.matrix(x) || is.atomic(x)) { ## 'old style' matrix or vector
	if     (is.integer(x)) "i"
	else if(is.numeric(x)) "d"
	else if(is.logical(x)) "l" ## FIXME ? "n" if no NA ??
	else if(is.complex(x)) "z"
	else stop(gettextf("not yet implemented for matrix with typeof %s",
			   typeof(x)), domain = NA)
    }
    else .M.kindC(clx)
}

.M.kindC <- function(clx, ex = extends(clx)) { ## 'clx': class() *or* classdefinition
    if(is.character(clx))		# < speedup: get it once
	clx <- getClassDef(clx)
    if(any(ex == "sparseVector")) {
	## must work for class *extending* "dsparseVector" ==> cannot use  (clx@className) !
	if     (any(ex == "dsparseVector")) "d"
	else if(any(ex == "nsparseVector")) "n"
	else if(any(ex == "lsparseVector")) "l"
	else if(any(ex == "zsparseVector")) "z"
	else if(any(ex == "isparseVector")) "i"
	else stop(gettextf(" not yet implemented for %s", clx@className),
		  domain = NA)
    }
    else if(any(ex == "dMatrix")) "d"
    else if(any(ex == "nMatrix")) "n"
    else if(any(ex == "lMatrix")) "l"
    else if(any(ex == "indMatrix")) "n" # permutation -> pattern
    else if(any(ex == "zMatrix")) "z"
    else if(any(ex == "iMatrix")) "i"
    else stop(gettextf(" not yet implemented for %s", clx@className),
	      domain = NA)
}

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
        ex <- extends(clx)
	if(     any(ex == "diagonalMatrix"))  "d"
	else if(any(ex == "triangularMatrix"))"t"
	else if(any(ex == "symmetricMatrix")) "s"
	else "g"
    }
}

## a faster simpler version [for sparse matrices, i.e., never diagonal]
.M.shapeC <- function(x, clx = class(x)) {
    if(is.character(clx)) # < speedup: get it once
	clx <- getClassDef(clx)
    if	   (extends(clx, "triangularMatrix")) "t"
    else if(extends(clx, "symmetricMatrix"))  "s" else "g"
}
} ## MJ

## MJ: unused
if(FALSE) {
## TODO: faster via C, either R's  R_data_class() [which needs to become API !]
##       or even direct  getAttrib(x, R_ClassSymbol); ..
##' class - single string, no "package" attribute,..
.class0 <- function(x)  as.vector(class(x))
} ## MJ

class2 <- function(cl, kind = "l", do.sub = TRUE) {
    ## Find "corresponding" class; since pos.def. matrices have no pendant:
    cl <- MatrixClass(cl)
    if(cl %in% c("dpoMatrix", "corMatrix"))
	paste0(kind, "syMatrix")
    else if(cl == "dppMatrix")
	paste0(kind, "spMatrix")
    else if(do.sub) sub("^[a-z]", kind, cl)
    else cl
}

## MJ: no longer used
if(FALSE) {
geClass <- function(x) {
    if     (is(x, "dMatrix")) "dgeMatrix"
    else if(is(x, "lMatrix")) "lgeMatrix"
    else if(is(x, "nMatrix") || is(x, "indMatrix")) "ngeMatrix"
    else if(is(x, "zMatrix")) "zgeMatrix"
    else stop(gettextf("general Matrix class not yet implemented for %s",
		       dQuote(class(x))), domain = NA)
}
} ## MJ

## typically used as .type.kind[.M.kind(x)]:
.type.kind <- c("d" = "double",
		"i" = "integer",
		"l" = "logical",
		"n" = "logical",
		"z" = "complex")

## the reverse, a "version of" .M.kind(.):
.kind.type <- setNames(names(.type.kind), as.vector(.type.kind))

.dense.prefixes <- c("d" = "tr", ## map diagonal to triangular
                     "t" = "tr",
                     "s" = "sy",
                     "g" = "ge")

.sparse.prefixes <- c("d" = "t", ## map diagonal to triangular
                      "t" = "t",
                      "s" = "s",
                      "g" = "g")

## MJ: no longer needed
if(FALSE) {
as_M.kind <- function(x, clx) {
    if(is.character(clx)) # < speedup: get it once
	clx <- getClassDef(clx)
    if(is(x, clx)) x else as(x, paste0(.M.kindC(clx), "Matrix"))
}

## Used, e.g. after subsetting: Try to use specific class -- if feasible :
as_dense <- function(x, cld = if(isS4(x)) getClassDef(class(x))) {
    as(x, paste0(.M.kind(x, cld), .dense.prefixes[.M.shape(x, cld)], "Matrix"))
}

if(FALSE) {
## This is "general" but slower than the next definition
.sp.class <- function(x) {
    if(!is.character(x)) x <- class(x)
    for(cl in paste0(c("C", "T", "R"), "sparseMatrix"))
	if(extends(x, cl))
	    return(cl)
    NA_character_
}
} else {
## find and return the "sparseness class" (aka "representation")
.sp.class <- function(x) {
    cl <- MatrixClass(if(is.character(x)) x else class(x))
    if(match(repr <- substr(cl, 3L, 3L), c("C", "T", "R"), 0L))
        return(paste0(repr, "sparseMatrix"))
    NA_character_
}
}

### Goal: Eventually get rid of these --- want to foster coercions
### ----  *to* virtual classes whenever possible, e.g. as(*, "CsparseMatrix")
## 2007-12: better goal: use them only for "matrix" [maybe speed them up later]

## Here, getting the class definition and passing it, should be faster
as_Csparse <- function(x, cld = if(isS4(x)) getClassDef(class(x))) {
    as(x, paste0(.M.kind(x, cld),
		 .sparse.prefixes[.M.shape(x, cld)], "CMatrix"))
}

if(FALSE) # replaced by .Call(dense_to_Csparse, *) which is perfect for "matrix"
as_Csparse2 <- function(x, cld = if(isS4(x)) getClassDef(class(x))) {
    ## Csparse + U2N when needed
    sh <- .M.shape(x, cld)
    x <- as(x, paste0(.M.kind(x, cld), .sparse.prefixes[sh], "CMatrix"))
    if(sh == "t") .Call(Csparse_diagU2N, x) else x
}

## 'cl'   : class() *or* class definition of from
as_gCsimpl2 <- function(from, cl = class(from))
    as(from, paste0(.M.kind(from, cl), "gCMatrix"))
## to be used directly in setAs(.) needs one-argument-only  (from) :
as_gCsimpl <- function(from) as(from, paste0(.M.kind(from), "gCMatrix"))

## slightly smarter:
as_Sp <- function(from, shape, cl = class(from)) {
    if(is.character(cl)) cl <- getClassDef(cl)
    as(from, paste0(.M.kind(from, cl),
		    shape,
		    if(extends(cl, "TsparseMatrix")) "TMatrix" else "CMatrix"))
}
## These are used in ./sparseMatrix.R:
as_gSparse <- function(from) as_Sp(from, "g", getClassDef(class(from)))
as_tSparse <- function(from) as_Sp(from, "t", getClassDef(class(from)))
as_sSparse <- function(from) as_Sp(from, "s", getClassDef(class(from)))

as_geSimpl2 <- function(from, cl = class(from))
    as(from, paste0(.M.kind(from, cl), "geMatrix"))
## to be used directly in setAs(.) needs one-argument-only  (from) :
as_geSimpl <- function(from) as(from, paste0(.M.kind(from), "geMatrix"))
} ## MJ

if(.Matrix.supporting.cached.methods) {
as_gCsimpl <- function(from) as(as(from, "CsparseMatrix"), "generalMatrix")
}

## (matrix|denseMatrix)->denseMatrix as similar as possible to "target"
as_denseClass <- function(x, cl, cld = getClassDef(cl)) {
    kind <- .M.kind(x)
    symmetric <- extends(cld, "symmetricMatrix") && isSymmetric(x)
    triangular <- !symmetric &&
        (extends(cld, "triangularMatrix") && (it <- isTriangular(x)))
    if(!(symmetric || triangular))
        return(.dense2g(x, kind))
    y <- if(symmetric)
             forceSymmetric(x)
         else if (attr(it, "kind") == "U")
             triu(x)
         else tril(x)
    if(extends(cld, "packedMatrix"))
        y <- pack(y)
    .dense2kind(y, kind)
}

## (matrix|sparseMatrix)->CsparseMatrix as similar as possible to "target"
as_CspClass <- function(x, cl, cld = getClassDef(cl)) {
    x <- as(x, "CsparseMatrix")
    x <- if(extends(cld, "symmetricMatrix") && isSymmetric(x))
             forceSymmetric(x)
         else if(!(extends(cld, "triangularMatrix") && (it <- isTriangular(x))))
             .sparse2g(x)
         else if(attr(it, "kind") == "U")
             triu(x)
         else tril(x)
    .sparse2kind(x, .M.kind(x))
}

## as(<Matrix>, <non-unit diagonal CsparseMatrix>)
asCspN <- function(x, cl = class(x), cld = getClassDef(cl)) {
    if(!extends(cld, "CsparseMatrix"))
        cld <- getClassDef(class(x <- as(x, "CsparseMatrix")))
    .Call(R_sparse_diag_U2N, x)
}

## MJ: no longer used
if(FALSE) {
asTri <- function(from, newclass) {
    ## TODO: also check for unit-diagonal: 'diag = "U"'
    isTri <- isTriangular(from)
    if(isTri)
 	new(newclass, x = from@x, Dim = from@Dim,
	    Dimnames = from@Dimnames, uplo = attr(isTri, "kind"))
    else stop("not a triangular matrix")
}

mat2tri <- function(from, sparse=NA) {
    isTri <- isTriangular(from)
    if(isTri) {
	d <- dim(from)
	if(is.na(sparse))
	    sparse <- prod(d) > 2*sum(isN0(from)) ## <==> sparseDefault() above
	if(sparse)
	    as(as(from, "sparseMatrix"), "triangularMatrix")
	else
	    new(paste0(.M.kind(from),"trMatrix"), x = base::as.vector(from),
                Dim = d, Dimnames = .M.DN(from), uplo = attr(isTri, "kind"))
    }
    else stop("not a triangular matrix")
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
} ## MJ

## MJ: no longer needed ... replacement in ./(un)?packedMatrix.R
if(FALSE) {
## For *dense* matrices
isTriMat <- function(object, upper = NA, ...) {
    ## pretest: is it square?
    d <- dim(object)
    if(d[1] != d[2]) return(FALSE)
    TRUE.U <- structure(TRUE, kind = "U")
    if(d[1] == 0) return(TRUE.U)
    ## else slower test
    TRUE.L <- structure(TRUE, kind = "L")
    if(!is.matrix(object))
	object <- as(object,"matrix")
    if(is.na(upper)) {
	if(all0(object[lower.tri(object)]))
	    TRUE.U
	else if(all0(object[upper.tri(object)]))
	    TRUE.L
	else FALSE
    } else if(upper)
	if(all0(object[lower.tri(object)])) TRUE.U else FALSE
    else ## upper is FALSE
	if(all0(object[upper.tri(object)])) TRUE.L else FALSE
}
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
## For Tsparse matrices:
isTriT <- function(object, upper = NA, ...) {
    ## pretest: is it square?
    d <- dim(object)
    if(d[1] != d[2]) return(FALSE)
    ## else
    TRUE.U <- structure(TRUE, kind = "U")
    if(d[1] == 0) return(TRUE.U)
    TRUE.L <- structure(TRUE, kind = "L")
    if(is.na(upper)) {
	if(all(object@i <= object@j))
	    TRUE.U
	else if(all(object@i >= object@j))
	    TRUE.L
	else FALSE
    } else if(upper) {
	if(all(object@i <= object@j)) TRUE.U else FALSE
    } else { ## 'lower'
	if(all(object@i >= object@j)) TRUE.L else FALSE
    }
}

## For Csparse matrices
isTriC <- function(object, upper = NA, ...) {
    ## pretest: is it square?
    d <- dim(object)
    if(d[1] != d[2]) return(FALSE)
    ## else
    TRUE.U <- structure(TRUE, kind = "U")
    if((n <- d[1]) == 0) return(TRUE.U)
    TRUE.L <- structure(TRUE, kind = "L")
    ## Need this, since 'i' slot of symmetric looks like triangular :
    if(is(object, "symmetricMatrix")) # triangular only iff diagonal :
        return(if(length(oi <- object@i) == n && isSeq(oi, n-1L)
                  && isSeq(object@p, n))
               structure(TRUE, kind = object@uplo) else FALSE)
    ## else
    ni <- 1:n
    ## the row indices split according to column:
    ilist <- split(object@i, factor(rep.int(ni, diff(object@p)), levels= ni))
    lil <- lengths(ilist, use.names = FALSE)
    if(any(lil == 0)) {
	pos <- lil > 0
	if(!any(pos)) ## matrix of all 0's
	    return(TRUE.U)
	ilist <- ilist[pos]
	ni <- ni[pos]
    }
    ni0 <- ni - 1L # '0-based ni'
    if(is.na(upper)) {
	if(all(sapply(ilist, max, USE.NAMES = FALSE) <= ni0))
	    TRUE.U
	else if(all(sapply(ilist, min, USE.NAMES = FALSE) >= ni0))
	    TRUE.L
	else FALSE
    } else if(upper) {
	if(all(sapply(ilist, max, USE.NAMES = FALSE) <= ni0))
	    TRUE.U else FALSE
    } else { ## 'lower'
	if(all(sapply(ilist, min, USE.NAMES = FALSE) >= ni0))
	    TRUE.L else FALSE
    }
}
} ## MJ

## MJ: no longer needed ... replacement in ./(un)?packedMatrix.R
if(FALSE) {
## When the matrix is known to be [n x n] aka "square"
## (need "vector-indexing" work for 'M'):
.is.diagonal.sq.matrix <- function(M, n = dim(M)[1L])
    all0(M[rep_len(c(FALSE, rep.int(TRUE,n)), n^2)])

.is.diagonal <- function(object) {
    ## "matrix" or "denseMatrix" (but not "diagonalMatrix")
    d <- dim(object)
    if(d[1L] != (n <- d[2L])) FALSE
    else if(is.matrix(object)) .is.diagonal.sq.matrix(object, n)
    else ## "denseMatrix" -- packed or unpacked
        if(is(object, "generalMatrix")) # "dge", "lge", ...
	    .is.diagonal.sq.matrix(object@x, n)
        else { ## "dense" but not {diag, general}, i.e. triangular or symmetric:
            ## -> has 'uplo'  differentiate between packed and unpacked

### .......... FIXME ...............
	    ## packed <- isPacked(object)
	    ## if(object@uplo == "U") {
	    ## } else { ## uplo == "L"
	    ## }
### very cheap workaround
	    all0(as(object,"matrix")[rep_len(c(FALSE, rep.int(TRUE,n)), n^2)])
        }
}
} ## MJ

## MJ: no longer used
if(FALSE) {
## for "dtC*", "ltC* ..: directly
xtC.diagU2N <- function(x) if(x@diag == "U") .Call(Csparse_diagU2N, x) else x

##' @title uni-diagonal to "regular" triangular Matrix
##'
##' NOTE:   class is *not* checked here! {speed}
##' @param x a dense unidiagonal (x@diag == "U") triangular Matrix
##'     ("ltrMatrix", "dtpMatrix", ...).
##' @param kind character indicating content kind: "d","l",..
##' @param isPacked logical indicating if 'x' is packed
##' @return Matrix "like" x, but with x@diag == "N" (and 1 or TRUE values "filled" in .@x)
##' @author Martin Maechler
.dense.diagU2N <- function(x, kind = .M.kind(x), isPacked = length(x@x) < n^2) {
    ## FIXME: Move this to C ----- (possibly with an option of *not* copying)
    ## For denseMatrix, .@diag = "U"  means the 'x' slot can have wrong values
    ## which are documented to never be accessed
    n <- x@Dim[1]
    if(n > 0) {
	one <- if(kind == "d") 1. else TRUE
	if(isPacked) { ## { == isPacked(x)) } : dtp-, ltp-, or "ntpMatrix":
	    ## x@x is of length	 n*(n+1)/2
	    if(n == 1)
		x@x <- one
	    else {
		di <- if(x@uplo == "U") seq_len(n) else c(1L,n:2L)
		x@x[cumsum(di)] <- one
	    }
	} else {
	    ## okay: now have  'x' slot of length n x n
	    x@x[1L+ (0:(n-1L))*(n+1L)] <- one # even for "n..Matrix"
	}
    }
    x@diag <- "N"
    x
}

##' @title coerce triangular Matrix to uni-diagonal
##'
##' NOTE: class is *not* checked here! {speed}
##' @param x a dense triangular Matrix ("ltrMatrix", "dtpMatrix", ...).
##' @return Matrix "like" x, but with x@diag == "U"
.dense.diagN2U <- function(x) {
    ## as we promise that the diagonal entries are not accessed when
    ##	diag = "U",   we don't even need to set them to one !!
    ## and *contrary* to the sparseMatrix case, we keep the diagonal entries in @x !
    x@diag <- "U"
    x
}

## This one and the following are fast "no-test" versions,
## which _know_ that 'x' is formally a triangularMatrix
## having a unit or non-unit diagonal
.diagU2N <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(extends(cl, "CsparseMatrix"))
	.Call(Csparse_diagU2N, x)
    if(extends(cl, "RsparseMatrix"))
	.tCR2RC(.Call(Csparse_diagU2N, .tCR2RC(x)))
    else if(extends(cl, "TsparseMatrix"))
	.Call(Tsparse_diagU2N, x)
    else if(checkDense && extends(cl, "denseMatrix"))
        .dense.diagU2N(x)
    else # still possibly dense
        .Call(Csparse_diagU2N, as(x, "CsparseMatrix"))
    ## ^leaving as CsparseMatrix ... caller can coerce as necessary
}

.diagN2U <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(extends(cl, "CsparseMatrix"))
        .Call(Csparse_diagN2U, x)
    else if(extends(cl, "RsparseMatrix"))
        .tCR2RC(.Call(Csparse_diagN2U, .tCR2RC(x)))
    else if(extends(cl, "TsparseMatrix"))
        .CR2T(.Call(Csparse_diagN2U, .T2C(x)))
    else if(checkDense && extends(cl, "denseMatrix"))
	.dense.diagN2U(x)
    else # still possibly dense
	.Call(Csparse_diagN2U, as(x, "CsparseMatrix"))
    ## ^leaving as CsparseMatrix ... caller can coerce as necessary
}
} ## MJ

diagU2N <- function (x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(extends(cl, "triangularMatrix") && x@diag == "U")
        .diagU2N(x, cl = cl, checkDense = checkDense)
    else x
}

.diagU2N <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(!checkDense && extends(cl, "denseMatrix"))
        x <- as(x, "CsparseMatrix")
    ..diagU2N(x)
}

..diagU2N <- function(x) {
    diag(x) <- TRUE
    x
}

diagN2U <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(extends(cl, "triangularMatrix") && x@diag == "N")
	.diagN2U(x, cl = cl, checkDense = checkDense)
    else x
}

.diagN2U <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(!checkDense & (isDense <- extends(cl, "denseMatrix")))
        ..diagN2U(as(x, "CsparseMatrix"), sparse = TRUE)
    else ..diagN2U(x, sparse = !isDense)
}

..diagN2U <- function(x, sparse) {
    if(sparse && x@Dim[1L] > 0L)
        x <- switch(x@uplo,
                    U = .Call(R_sparse_band, x, 1L, NULL),
                    L = .Call(R_sparse_band, x, NULL, -1L),
                    stop("invalid 'uplo'"))
    x@diag <- "U"
    x
}

if(.Matrix.supporting.cached.methods) {
.dgC.0.factors <- function(x) {
    if(length(x@factors))
        x@factors <- list()
    x
}
}

# MJ: no longer used
if(FALSE) {
.as.dgC.0.factors <- function(x) {
    if(is(x, "dgCMatrix"))
        .dgC.0.factors(x)
    else as(x, "dgCMatrix") # will not have 'factors'
}
} ## MJ

## Caches 'value' in the 'factors' slot of 'x', i.e. modifies 'x', and returns 'value'
## WARNING:: for updating the '@ factors' slot of a function *argument* [CARE!]
.set.factors <- function(x, name, value, warn.no.slot=FALSE)
    .Call(R_set_factor, x, name, value, warn.no.slot)

##' Change function *argument* 'x', emptying its 'factors' slot; USE with CARE! __ DANGER ! __
##' @return TRUE iff 'x' is modified, FALSE if not.
.empty.factors <- function(x, warn.no.slot=FALSE)
    .Call(R_empty_factors, x, warn.no.slot)

##' The *SAFE* regular function version:  empty the factor slot
.drop.factors <- function(x, check=FALSE)
   `slot<-`(x, "factors", check=check, value=list())

## MJ: not used
if(FALSE) {
### Fast, much simplified version of tapply()
tapply1 <- function (X, INDEX, FUN = NULL, ..., simplify = TRUE) {
    sapply(unname(split(X, INDEX)), FUN, ...,
	   simplify = simplify, USE.NAMES = FALSE)
}
tapply.x <- function (X, n, INDEX, FUN = NULL, ..., simplify = TRUE) {
    tapply1(X, factor(INDEX, 0:(n-1)), FUN = FUN, ..., simplify = simplify)
}

### MM: Unfortunately, these are still pretty slow for large sparse ...
sparsapply <- function(x, MARGIN, FUN, sparseResult = TRUE, ...) {
    ## Purpose: "Sparse Apply": better than tapply1() for colSums(), etc.:
    ##    NOTE: Only applicable sum()-like where the "zeros do not count"
    ## ----------------------------------------------------------------------
    ## Arguments: x: sparseMatrix;  others as in *apply()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 16 May 2007
    stopifnot(MARGIN %in% 1:2)
    xi <- if(MARGIN == 1) x@i else x@j
    ui <- unique.default(xi)
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

all0Matrix <- function(n,m) {
    ## an  all-0 matrix	 -- chose what Matrix() also gives -- "most efficiently"
    n <- as.integer(n)
    m <- as.integer(m)
    new(if(n == m) "dsCMatrix" else "dgCMatrix",
	Dim = c(n,m),
	p = rep.int(0L, m+1L))
}
} ## MJ

## all-0 matrix from 'x' which must inherit from Matrix
.setZero <- function(x, kind = .M.kind(x)) {
    cl <- if(.hasSlot(x, "diag"))
              ".tCMatrix"
          else if(.hasSlot(x, "uplo"))
              ".sCMatrix"
          else ".gCMatrix"
    substr(cl, 1L, 1L) <- kind
    d <- x@Dim
    new(cl, Dim = d, Dimnames = x@Dimnames, p = rep.int(0L, d[2L] + 1))
}

##' Subsetting matrix/vector in "vector style", e.g. M[0], M[TRUE], M[1:2], M[-7]
##' @param x any matrix/Matrix/(sparse)vector, to be subset
##' @param i integer (incl negative!) or logical 'index'.
##' @param allowSparse logical indicating if the result may be a
##'   \code{"sparseVector"}; the default is false for reasons of back
##'   compatibility  (against efficiency here).
##' @note 2018-03: Now partially based on \code{as(x, "sparseVector")[i]}
##'   which has been improved itself.
.M.vectorSub <- function(x, i, allowSparse=FALSE) {
    if(prod(dim(x)) == 0)
	as(x, "matrix")[i]
    else if(any(as.logical(i))) {
	if(inherits(x, "denseMatrix"))
	    as(x, "matrix")[i]
	else { ## sparse ...
	    ## if(is.logical(i)) # unfortunately, this is not-yet-implemented!
	    ##     x[as(i, "sparseVector")]
	    ## else if(all(i >= 0))
	    if(is.numeric(i) && all(i >= 0))
		subset.ij(x, ij = arrayInd(i, .dim=dim(x), useNames=FALSE))
	    else if(allowSparse) # more efficient here
		as(x, "sparseVector")[i]
	    else # sparse result not allowed
		sp2vec(as(x, "sparseVector")[i])
	}
    } else ## save memory (for large sparse M):
	as.vector(x[1,1])[FALSE]
}

## MJ: not used
if(FALSE) {
##' Compute the three "parts" of two sets:
##' @param x arbitrary vector; possibly with duplicated values,
##' @param y (ditto)
##' @param uniqueCheck
##' @param check
##'
##' @return list(x.only = setdiff(x,y),
##'              y.only = setdiff(y,x),
##'	         int = intersect(x,y))
setparts <- function(x,y, uniqueCheck = TRUE, check = TRUE) {
    if(check) {
        x <- as.vector(x)
        y <- as.vector(y)
    }
    if(uniqueCheck) {
        x <- unique.default(x)
        y <- unique.default(y)
    }
    .setparts(x,y)
}
} ## MJ

.setparts <- function(x, y) {
    n1 <- length(m1 <- match(x, y, 0L))
    n2 <- length(m2 <- match(y, x, 0L))
    ix <- which(m1 == 0L)
    iy <- which(m2 == 0L)
    list(x.only = x[ix], ix.only = ix, mx = m1,
         y.only = y[iy], iy.only = iy, my = m2,
         int = if(n1 < n2) y[m1] else x[m2])
}

## MJ: no longer needed ... now using base::chkDots()
if(FALSE) {
##' @title Warn about extraneous arguments in the "..."  (of its caller).
##' A merger of my approach and the one in seq.default() -- FIXME: now have base::chkDots()
##' @author Martin Maechler, June 2012, May 2014
##' @param ...
##' @param which.call passed to sys.call().  A caller may use -2 if the message should
##' mention *its* caller
chk.s <- function(..., which.call = -1,
		  depCtrl = if(exists("..deparseOpts")) "niceNames")
{
    if(nx <- length(list(...)))
	warning(sprintf(ngettext(nx,
                                 "extra argument %s will be disregarded in\n %s",
                                 "extra arguments %s will be disregarded in\n %s"),
                        sub(")$", '', sub("^list\\(", '',
                                          deparse1(list(...), control=depCtrl))),
                        deparse1(sys.call(which.call), control=depCtrl)),
                call. = FALSE, domain=NA)
}
} ## MJ

##' *Only* to be used as function in
##'    setMethod.("Compare", ...., .Cmp.swap)  -->  ./Ops.R  & ./diagMatrix.R
.Cmp.swap <- function(e1, e2) {
    ## "swap RHS and LHS" and use the method below:
    switch(.Generic,
	   "==" =,
           "!=" = callGeneric(e2, e1),
	   "<"	= e2 >	e1,
	   "<=" = e2 >= e1,
	   ">"	= e2 <	e1,
	   ">=" = e2 <= e1)
}


### These two are very similar, the first one has the advantage
### to be applicable to 'Chx' directly:

## FIXME:  kind = "diagBack" is not yet implemented
##	would be much more efficient, but there's no CHOLMOD UI (?)

## "used" currently only in ../tests/factorizing.R
.diag.dsC <- function(x, Chx = Cholesky(x, LDL=TRUE), res.kind = "diag") {
    force(Chx)
    if(!missing(Chx)) stopifnot(.isLDL(Chx), is.integer(Chx@p), is.double(Chx@x))
    .Call(diag_tC, Chx, res.kind)
    ##    ^^^^^^^ from ../src/Csparse.c
    ## => res.kind in ("trace", "sumLog", "prod", "min", "max", "range", "diag", "diagBack")
}

## MJ: unused
if(FALSE) {
## here, we  *could* allow a 'mult = 0' factor :
.CHM.LDL.D <- function(x, perm = TRUE, res.kind = "diag") {
    .Call(dsCMatrix_LDL_D, x, perm, res.kind)
    ##    ^^^^^^^^^^^^^^^^ from ../src/dsCMatrix.c
}
} ## MJ

dimScale <- function(x, d1 = sqrt(1/diag(x, names = FALSE)), d2 = d1) {
    dim.x <- dim(x)
    D1 <- Diagonal(n = dim.x[1L], x = d1)
    D2 <- if(missing(d2)) D1 else Diagonal(n = dim.x[2L], x = d2)
    y <- D1 %*% x %*% D2 # inefficient for symmetricMatrix 'x', but "general"
    if(isS4(x) && is(x, "symmetricMatrix") && identical(d1, d2))
        y <- forceSymmetric(y, x@uplo)
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}

rowScale <- function(x, d) {
    y <- Diagonal(n = nrow(x), x = d) %*% x
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}

colScale <- function(x, d) {
    y <- x %*% Diagonal(n = ncol(x), x = d)
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}
