#### Will be sourced by several R scripts in ../tests/

### ------- Part I --  unrelated to "Matrix" classes ---------------

paste0 <- function(...) paste(..., sep = '')

identical3 <- function(x,y,z)	identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d) identical(a,b) && identical3(b,c,d)

## Make sure errors are signaled
assertError <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- try(expr, silent = TRUE)
    if(!inherits(t.res, "try-error"))
	stop(d.expr, "\n\t did not give an error", call. = FALSE)
    invisible(t.res)
}

isValid <- function(x, class) validObject(x, test=TRUE) && is(x, class)

is.all.equal3 <- function(x,y,z, tol = .Machine$double.eps^0.5)
    isTRUE(all.equal(x,y, tol=tol)) && isTRUE(all.equal(y,z, tol=tol))

is.all.equal4 <- function(x,y,z,u, tol = .Machine$double.eps^0.5)
    is.all.equal3(x,y,z, tol=tol) && isTRUE(all.equal(z,u, tol=tol))

## A version of all.equal() for the slots
all.slot.equal <- function(x,y, ...) {
    slts <- slotNames(x)
    for(sl in slts) {
        aeq <- all.equal(slot(x,sl), slot(y,sl), ...)
        if(!identical(TRUE, aeq))
            return(paste("slot '",sl,"': ", aeq, sep=''))
    }
    TRUE
}

## The relative error typically returned by all.equal:
relErr <- function(target, current) { ## make this work for 'Matrix'
    ## ==> no mean() ..
    n <- length(current)
    if(length(target) < n)
        target <- rep(target, length.out = n)
    sum(abs(target - current)) / sum(abs(target))
}

## is.R22 <- (paste(R.version$major, R.version$minor, sep=".") >= "2.2")

pkgRversion <- function(pkgname)
    substring(packageDescription(pkgname)[["Built"]], 3,5)


### ------- Part II  -- related to matrices, but *not* "Matrix" -----------

add.simpleDimnames <- function(m) {
    stopifnot(length(d <- dim(m)) == 2)
    dimnames(m) <- list(paste0("r", seq_len(d[1])),
                        paste0("c", seq_len(d[2])))
    m
}

as.mat <- function(m) {
    ## as(., "matrix")	but with no extraneous empty dimnames
    m <- as(m, "matrix")
    if(identical(dimnames(m), list(NULL,NULL)))
	dimnames(m) <- NULL
    m
}

assert.EQ.mat <- function(M, m, tol = if(show) 0 else 1e-15, show=FALSE) {
    ## Purpose: check equality of  'Matrix' M with  'matrix' m
    ## ----------------------------------------------------------------------
    ## Arguments: M: is(., "Matrix") typically {but just needs working as(., "matrix")}
    ##            m: is(., "matrix")
    ##            show: if TRUE, return (and hence typically print) all.equal(...)
    validObject(M)
    MM <- as.mat(M)                     # as(M, "matrix")
    if(is.logical(MM) && is.numeric(m))
	storage.mode(MM) <- "integer"
    attr(MM, "dimnames") <- attr(m, "dimnames") <- NULL
    if(show) all.equal(MM, m, tol = tol)
    else if(!isTRUE(r <- all.equal(MM, m, tol = tol)))
	stop("all.equal() |->  ", r)
}

chk.matrix <- function(M) {
    ## check object; including coercion to "matrix" :
    cl <- class(M)
    cat("class ", dQuote(cl), " [",nrow(M)," x ",ncol(M),"]; slots (",
	paste(slotNames(M), collapse=","), ")\n", sep='')
    stopifnot(validObject(M),
	      dim(M) == c(nrow(M), ncol(M)),
	      identical(dim(m <- as(M, "matrix")), dim(M))
	      )
}

isOrthogonal <- function(x, tol = 1e-15) {
    all.equal(diag(as(zapsmall(crossprod(x)), "diagonalMatrix")),
              rep(1, ncol(x)), tol = tol)
}


### ------- Part III --  "Matrix" (classes) specific ----------------------

lsM <- function(...) {
    for(n in ls(..., envir=parent.frame()))
        if(is((. <- get(n)),"Matrix"))
            cat(sprintf("%5s: '%s' [%d x %d]\n",n,class(.), nrow(.),ncol(.)))
}

asD <- function(m) { ## as "Dense"
    if(canCoerce(m, "denseMatrix")) as(m, "denseMatrix")
    else if(canCoerce(m, (cl <- paste(.M.kind(m), "denseMatrix", sep=''))))
        as(m, cl)
    else if(canCoerce(m, "dgeMatrix")) as(m, "dgeMatrix")
    else stop("cannot coerce to a typical dense Matrix")
}

## "normal" sparse Matrix: Csparse, no diag="U"
asCsp <- function(x) Matrix:::diagU2N(as(x, "CsparseMatrix"))

Qidentical <- function(x,y) {
    ## quasi-identical - for 'Matrix' matrices
    if(class(x) != class(y)) return(FALSE)
    slts <- slotNames(x)
    if("factors" %in% slts) { ## allow one empty and one non-empty 'factors'
        slts <- slts[slts != "factors"]
        ## if both are not empty, they must be the same:
        if(length(xf <- x@factors) && length(yf <- y@factors))
            if(!identical(xf, yf)) return(FALSE)
    }
    for(sl in slts)
        if(!identical(slot(x,sl), slot(y,sl)))
            return(FALSE)
    TRUE
}

Q.C.identical <- function(x,y, sparse = is(x,"sparseMatrix")) {
    if(sparse) Qidentical(as(x,"CsparseMatrix"), as(y,"CsparseMatrix"))
    else Qidentical(x,y)
}

## Useful Matrix constructors for testing:

rspMat <- function(n, m = n, density = 1/4, nnz = round(density * n*m))
{
    ## Purpose: random sparse Matrix
    ## ----------------------------------------------------------------------
    ## Arguments: (n,m) : dimension [default m=n ==> *square* matrix}
    ##           density: the desired sparseness density:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 5 Mar 2008, 11:07
    stopifnot(length(n) == 1, n == as.integer(n),
              length(m) == 1, m == as.integer(m),
              0 <= density, density <= 1,
              0 <= nnz, nnz <= n*m)
    x <- numeric(n*m)
    ## entries 2 : (nnz+1) {so we can have '1' as 'special'}
    x[sample(n*m, nnz, replace=FALSE)] <- as.numeric(1L + seq_len(nnz))
    Matrix(x, n,m, sparse=TRUE)
}



rUnitTri <- function(n, upper = TRUE, ...)
{
    ## Purpose: random unit-triangular sparse Matrix .. built from rspMat()
    ## ----------------------------------------------------------------------
    ## Arguments:  n: matrix dimension
    ##         upper: logical indicating if upper or lower triangular
    ##         ...  : further arguments passed to rspMat(), eg. 'density'
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  5 Mar 2008, 11:35

    r <- (if(upper) triu else tril)(rspMat(n, ...))
    ## make sure the diagonal is empty
    diag(r) <- 0
    r <- drop0(r)
    r@diag <- "U"
    r
}

## This is related to rUnitTri(), ver
mkLDL <- function(n, density = 1/3) {
    ## Purpose: make nice artifical   A = L D L'  (with exact numbers) decomp
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 15 Mar 2008
    stopifnot(n == round(n))
    n <- as.integer(n)
    L <- Matrix(0, n,n)

    nnz <- round(n*n * density)
    L[sample(n*n, nnz)] <- seq_len(nnz)
    L <- tril(L,-1)
    diag(L) <- 1
    d.half <- sample(10*(n:1))# random permutation ; use '10*' to be "different" from L entries
    D <- Diagonal(x = d.half * d.half)
    A <- tcrossprod(L * rep(d.half, each=n))
    ## = as(L %*% D %*% t(L), "symmetricMatrix")
    list(A = A, L = L, d.half = d.half, D = D)
}

eqDeterminant <- function(m1, m2, ...) {
    d1 <- determinant(m1)
    d2 <- determinant(m2)
    (d1$modulus == -Inf && d2$modulus == -Inf) || all.equal(d1, d2, ...)
}

##--- Compatibility tests "Matrix" =!= "traditional Matrix" ---
checkMatrix <- function(m, m.m = if(do.matrix) as(m, "matrix"),
			do.matrix = !(isSparse || isDiag) || prod(dim(m)) < 1e6,
			do.t = TRUE, doNorm = TRUE, doOps = TRUE, doSummary = TRUE,
			doCoerce = TRUE,
			doCoerce2 = doCoerce && !extends(cld, "RsparseMatrix"),
			verbose = TRUE, catFUN = cat)
{
    ## Purpose: Compatibility tests "Matrix" <-> "traditional matrix"
    ## ----------------------------------------------------------------------
    ## Arguments: m: is(., "Matrix")
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 Apr 2008; building on tests originally
    ##	       in dotestMat()  ../tests/Class+Meth.R
    stopifnot(is(m, "Matrix"))
    validObject(m)
    clNam <- class(m)
    cld <- getClassDef(clNam) ## extends(cld, FOO) is faster than is(m, FOO)
    isCor    <- extends(cld, "corMatrix")
    isSparse <- extends(cld, "sparseMatrix")
    isSym    <- extends(cld, "symmetricMatrix")
    isDiag   <- extends(cld, "diagonalMatrix")

    Cat	 <- function(...) if(verbose) cat(...)
    CatF <- function(...) if(verbose) catFUN(...)
    warnNow <- function(...) warning(..., call. = FALSE, immediate. = TRUE)
    if(getRversion() < "2.7.0") {
	## Equivalent of fixing base::determinant.matrix() :
	determinant.matrix <- function(x, logarithm = TRUE, ...)
	{
	    if ((n <- ncol(x)) != nrow(x))
		stop("'x' must be a square matrix")
	    if (n < 1)
		return(structure(list(modulus =
				      structure(if(logarithm) 0 else 1,
						logarithm = logarithm),
				      sign = 1L),
				 class = "det"))
	    if (is.complex(x))
		stop("determinant not currently defined for complex matrices")
	    storage.mode(x) <- "double"
	    .Call("det_ge_real", x, logarithm, PACKAGE = "base")
	}
	setMethod(determinant, c("matrix","ANY"), determinant.matrix)
	setMethod(determinant, c("matrix","logical"), determinant.matrix)
    }

    if(getRversion() < "2.7.1") {
	## "Fixup base::diag()", i.e. the default diag() method, locally:
	b.diag <- function(x = 1, nrow, ncol)
	{
	    if (is.matrix(x) && nargs() == 1) {
		if((m <- min(dim(x))) == 0)
		    return(vector(typeof(x), 0)) # logical, integer, also list ..

		y <- c(x)[1 + 0:(m - 1) * (dim(x)[1] + 1)]
		nms <- dimnames(x)
		if (is.list(nms) && !any(sapply(nms, is.null)) &&
		    identical((nm <- nms[[1]][1:m]), nms[[2]][1:m]))
		    names(y) <- nm
		return(y)
	    }
	    if(is.array(x) && length(dim(x)) != 1)
		stop("first argument is array, but not matrix.")

	    if(missing(x))
		n <- nrow
	    else if(length(x) == 1 && missing(nrow) && missing(ncol)) {
		n <- as.integer(x)
		x <- 1
	    }
	    else n <- length(x)
	    if(!missing(nrow))
		n <- nrow
	    if(missing(ncol))
		ncol <- n
	    p <- ncol
	    y <- array(0, c(n, p))
	    if((m <- min(n, p)) > 0) y[1 + 0:(m - 1) * (n + 1)] <- x
	    y
	}
	## .old.base.diag <- base::diag
	on.exit(setMethod(diag, "ANY", base::diag))
	setMethod(diag, "ANY", b.diag)
    }## end{Fix base::diag() for R <= 2.7.0}

    DO.m <- function(expr) if(do.matrix) eval(expr) else TRUE
    ina <- is.na(m)
    if(do.matrix) {
	stopifnot(all(ina == is.na(m.m)),
		  all(m == m | ina)) ## check all() , "==" [Compare], "|" [Logic]
	if(any(m != m & !ina)) stop(" any (m != m) should not be true")
    } else {
	if(any(m != m)) stop(" any (m != m) should not be true")
    }
    if(do.t) {
	tm <- t(m)
	if(isSym) ## check that t() swaps 'uplo'  L <--> U :
	    stopifnot(c("L","U") == sort(c(m@uplo, tm@uplo)))
	ttm <- t(tm)
	if(extends(cld, "CsparseMatrix") ||
	   extends(cld, "generalMatrix") || isDiag)
            stopifnot(Qidentical(m, ttm))
	else if(do.matrix)
	    stopifnot(class(ttm) == clNam, all(m == ttm | ina))
        ## else : not testing
    }
    if(!do.matrix) {
	CatF(" will *not* coerce to 'matrix' since do.matrix is FALSE\n")
    } else if(doNorm) {
	CatF(sprintf(" norm(m [%d x %d]) :", nrow(m), ncol(m)))
	for(typ in c("1","I","F","M")) {
	    Cat('', typ, '')
	    stopifnot(all.equal(norm(m,typ), norm(m.m,typ)))
	}
	Cat(" ok\n")
    }
    if(do.matrix && doSummary) {
	summList <- lapply(getGroupMembers("Summary"), get,
			   envir = asNamespace("Matrix"))
	CatF(" Summary: ")
	for(f in summList) {
	    ## suppressWarnings():  e.g. any(<double>)	would warn here:
	    r <- suppressWarnings(if(isCor) all.equal(f(m), f(m.m)) else
				  identical(f(m), f(m.m)))
	    if(!isTRUE(r)) {
		f.nam <- sub("..$", '', sub("^\\.Primitive..", '', format(f)))
		(if(f.nam == "prod") warnNow else stop)(
		    sprintf("%s(m) differs from %s(m.m)", f.nam, f.nam))
	    }
	}
	if(verbose) cat(" ok\n")
    }

    ## and test 'dim()' as well:
    d <- dim(m)
    isSqr <- d[1] == d[2]
    isPerm <- extends(cld, "pMatrix")
    if(do.t) stopifnot(identical(diag(m), diag(t(m))))
    ## TODO: also === diag(band(m,0,0))
    if(do.matrix)
    stopifnot(identical(dim(m.m), dim(m)),
	      ## base::diag() keeps names [Matrix FIXME]
	      if(isPerm) {
		  identical(as.integer(unname(diag(m))), unname(diag(m.m)))
	      } else
	      identical(unname(diag(m)),
			unname(diag(m.m))),## not for NA: diag(m) == diag(m.m),
	      identical(nnzero(m), sum(m.m != 0)),
	      identical(nnzero(m, na.= FALSE), sum(m.m != 0, na.rm = TRUE)),
	      identical(nnzero(m, na.= TRUE),  sum(m.m != 0 | is.na(m.m)))
	      )

    if(isSparse) {
	n0m <- drop0(m) #==> n0m is Csparse
	has0 <- !Qidentical(n0m, as(m,"CsparseMatrix"))
    }
    ## use non-square matrix when "allowed":

    ## m12: sparse and may have 0s even if this is not: if(isSparse && has0)
    m12 <- as(as(  m, "lMatrix"),"CsparseMatrix")
    m12 <- drop0(m12)
    if(do.matrix) {
	## "!" should work (via as(*, "l...")) :
	m11 <- as(as(!!m,"CsparseMatrix"), "lMatrix")
	if(!Qidentical(m11, m12))
	    stopifnot(Qidentical(as(m11, "generalMatrix"),
				 as(m12, "generalMatrix")))
    }
    if(isSparse && !extends(cld,"nMatrix")) {
	## ensure that as(., "nMatrix") gives nz-pattern
	CatF("as(., \"nMatrix\") giving full nonzero-pattern: ")
	n1 <- as(m, "nMatrix")
	ns <- as(m, "nsparseMatrix")
	stopifnot(identical(n1,ns),
		  isDiag || ((if(isSym) Matrix:::nnzSparse else sum)(n1) ==
			     length(if(isPerm) m@perm else Matrix:::diagU2N(m)@x)))
        Cat("ok\n")
    }

    if(doOps) {
	## makes sense with non-trivial m (!)
	CatF("2*m =?= m+m: ")
	if(identical(2*m, m+m)) Cat("identical\n")
	else if(do.matrix) {
	    stopifnot(as(2*m,"matrix") == as(m+m, "matrix"))
	    Cat("ok\n")
	} else {# !do.matrix
	    identical(as(2*m, "CsparseMatrix"),
		      as(m+m, "CsparseMatrix"))
	    Cat("ok\n")
	}
	if(do.matrix) {
	    ## m == m etc, now for all, see above
	    CatF("m >= m for all: "); stopifnot(all(m >= m | ina)); Cat("ok\n")
	}
	if(prod(d) > 0) {
	    CatF("m < m for none: ")
	    mlm <- m < m
	    if(!any(ina)) stopifnot(!any(mlm))
	    else if(do.matrix) stopifnot(!any(mlm & !ina))
	    else { ## !do.matrix & any(ina) :  !ina can *not* be used
		mlm[ina] <- FALSE
		stopifnot(!any(mlm))
	    }
	    Cat("ok\n")
	}

	if(isSqr) {
	    if(do.matrix) {
		## determinant(<dense>) "fails" for triangular with NA such as
		## (m <- matrix(c(1:0,NA,1), 2))
		CatF("det...(): ")
		if(any(is.na(m.m)) && extends(cld, "triangularMatrix"))
		    Cat(" skipped: is triang. and has NA")
		else
		    stopifnot(eqDeterminant(m, m.m))
		Cat("ok\n")
	    }
	} else assertError(determinant(m))
    }

    if(doCoerce && do.matrix && canCoerce("matrix", clNam)) {
	CatF("as(<matrix>, ",clNam,"): ", sep='')
	m3 <- as(m.m, clNam)
	Cat("valid:", validObject(m3), "\n")
	## m3 should ``ideally'' be identical to 'm'
    }

    if(doCoerce2 && do.matrix) { ## not for large m:  !m will be dense

	if(extends(cld, "nMatrix")) {
	    stopifnot(identical(m, as(as(m, "dMatrix"),"nMatrix")),
		      identical(m, as(as(m, "lMatrix"),"nMatrix")))
	}
	else if(extends(cld, "lMatrix")) { ## should fulfill even with NA:
	    stopifnot(all(m | !m | ina), !any(!m & m & !ina))
	    if(extends(cld, "TsparseMatrix")) # allow modify, since at end here
		m <- Matrix:::uniqTsparse(m, clNam)
	    stopifnot(identical(m, m & TRUE),
		      identical(m, FALSE | m))
	    ## also check the  coercions to [dln]Matrix
	    m. <- if(isSparse && has0) n0m else m
	    m1. <- m. # replace NA by 1 in m1. , carefully not changing class:
	    if(any(ina)) m1.@x[is.na(m1.@x)] <- TRUE
	    stopifnot(identical(m. , as(as(m. , "dMatrix"),"lMatrix")),
		      clNam == "ldiMatrix" || # <- there's no "ndiMatrix"
		      ## coercion to n* and back: only identical when no extra 0s:
		      identical(m1., as(as(m1., "nMatrix"),"lMatrix")))
	}
	else if(extends(cld, "dMatrix")) {
	    m. <- if(isSparse && has0) n0m else m
	    m1 <- (m. != 0)*1
	    if(!isSparse && substr(clNam,1,3) == "dpp")
		## no "nppMatrix" possible
		m1 <- unpack(m1)

	    m1. <- m1 # replace NA by 1 in m1. , carefully not changing class:
	    if(any(ina)) m1.@x[is.na(m1.@x)] <- 1
	    ## coercion to n* (nz-pattern!) and back: only identical when no extra 0s and no NAs:
	    stopifnot(Q.C.identical(m1., as(as(m., "nMatrix"),"dMatrix"), isSparse),
		      Q.C.identical(m1 , as(as(m., "lMatrix"),"dMatrix"), isSparse))
	}

	if(extends(cld, "triangularMatrix")) {
	    mm. <- m
	    i0 <- if(m@uplo == "L")
		upper.tri(mm.) else lower.tri(mm.)
	    mm.[i0] <- 0 # ideally, mm. remained triangular, but can be dge*
	    CatF("as(<triangular (ge)matrix>, ",clNam,"): ", sep='')
	    tm <- as(as(mm., "triangularMatrix"), clNam)
	    Cat("valid:", validObject(tm), "\n")
	    if(m@uplo == tm@uplo) ## otherwise, the matrix effectively was *diagonal*
		## note that diagU2N(<dtr>) |-> dtC :
		stopifnot(Qidentical(tm, as(Matrix:::diagU2N(m), clNam)))
	}
	else if(extends(cld, "diagonalMatrix")) {

	    ## TODO

	} else {

	    ## TODO
	}
    }
    invisible(TRUE)
}
