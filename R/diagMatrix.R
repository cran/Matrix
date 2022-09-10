## Methods for "diagonalMatrix" and its subclasses, currently "[dl]diMatrix"

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("Matrix", "diagonalMatrix", ..M2diag)
setAs("matrix", "diagonalMatrix", ..M2diag)

## MJ: no longer needed ... replacement above
if(FALSE) {
setAs("matrix", "diagonalMatrix",
      function(from) {
	  d <- dim(from)
	  if(d[1] != (n <- d[2])) stop("non-square matrix")
	  if(any(from[row(from) != col(from)] != 0))
	      stop("matrix with non-zero off-diagonals cannot be coerced to \"diagonalMatrix\"")
	  x <- diag(from); names(x) <- NULL # don't want them in 'x' slot
	  if(is.logical(x)) {
	      cl <- "ldiMatrix"
	      uni <- allTrue(x) ## uni := {is it unit-diagonal ?}
	  } else {
	      cl <- "ddiMatrix"
	      uni <- allTrue(x == 1)
	      storage.mode(x) <- "double"
	  } ## TODO: complex
	  new(cl, Dim = c(n,n), diag = if(uni) "U" else "N",
	      x = if(uni) x[FALSE] else x, Dimnames = .M.DN(from))
      })

## ``generic'' coercion to  diagonalMatrix : build on  isDiagonal() and diag()
setAs("Matrix", "diagonalMatrix",
      function(from) {
	  d <- dim(from)
	  if(d[1] != (n <- d[2])) stop("non-square matrix")
	  if(!isDiagonal(from)) stop("matrix is not diagonal")
	  ## else:
	  x <- diag(from); names(x) <- NULL # don't want them in 'x' slot
	  if(is.logical(x)) {
	      cl <- "ldiMatrix"
	      uni <- allTrue(x)
	  } else {
	      cl <- "ddiMatrix"
	      uni <- allTrue(x == 1)
	      storage.mode(x) <- "double"
	  } ## TODO: complex
	  new(cl, Dim = c(n,n), diag = if(uni) "U" else "N",
	      x = if(uni) x[FALSE] else x, Dimnames = from@Dimnames)
      })
} ## MJ


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

..diag2dsparse <- function(from)
    .Call(R_diagonal_as_sparse, from, "dtT", "U", TRUE)
..diag2lsparse <- function(from)
    .Call(R_diagonal_as_sparse, from, "ltT", "U", TRUE)
..diag2nsparse <- function(from)
    .Call(R_diagonal_as_sparse, from, "ntT", "U", TRUE)

..diag2tC <- function(from)
    .Call(R_diagonal_as_sparse, from, ".tC", "U", TRUE)
..diag2tR <- function(from)
    .Call(R_diagonal_as_sparse, from, ".tR", "U", TRUE)
..diag2tT <- function(from)
    .Call(R_diagonal_as_sparse, from, ".tT", "U", TRUE)
..diag2sC <- function(from)
    .Call(R_diagonal_as_sparse, from, ".sC", "U", TRUE)
..diag2sR <- function(from)
    .Call(R_diagonal_as_sparse, from, ".sR", "U", TRUE)
..diag2sT <- function(from)
    .Call(R_diagonal_as_sparse, from, ".sT", "U", TRUE)
..diag2gC <- function(from)
    .Call(R_diagonal_as_sparse, from, ".gC", "U", TRUE)
..diag2gR <- function(from)
    .Call(R_diagonal_as_sparse, from, ".gR", "U", TRUE)
..diag2gT <- function(from)
    .Call(R_diagonal_as_sparse, from, ".gT", "U", TRUE)

## .diag2[ts]T() are exported ...
.diag2tT <- function(from, uplo = "U", kind = ".", drop0 = TRUE)
    .Call(R_diagonal_as_sparse, from,
          `substr<-`(".tT", 1L, 1L, kind), uplo, drop0)
.diag2sT <- function(from, uplo = "U", kind = ".", drop0 = TRUE)
    .Call(R_diagonal_as_sparse, from,
          `substr<-`(".sT", 1L, 1L, kind), uplo, drop0)

## For group methods
.diag2tT.smart <- function(from, x, uplo = "U", kind = ".", drop0 = TRUE) {
    .Call(R_diagonal_as_sparse, from,
          `substr<-`(".sT", 1L, 1L, kind),
          if(is(x, "triangularMatrix")) x@uplo else uplo,
          drop0)
}
.diag2T.smart <- function(from, x, uplo = "U", kind = ".", drop0 = TRUE) {
    symmetric <- extends(cld <- getClassDef(class(x)), "symmetricMatrix")
    .Call(R_diagonal_as_sparse, from,
          `substr<-`(if(symmetric) ".sT" else ".tT", 1L, 1L, kind),
          if(symmetric || extends(cld, "triangularMatrix")) x@uplo else uplo,
          drop0)
}

..diag2tp  <- function(from) .diag2dense(from, ".tp", "U")
..diag2tr  <- function(from) .diag2dense(from, ".tr", "U")
..diag2dtr <- function(from) .diag2dense(from, "dtr", "U")
..diag2ltr <- function(from) .diag2dense(from, "ltr", "U")
..diag2ntr <- function(from) .diag2dense(from, "ntr", "U")
..diag2dsy <- function(from) .diag2dense(from, "dsy", "U")
..diag2lsy <- function(from) .diag2dense(from, "lsy", "U")
..diag2nsy <- function(from) .diag2dense(from, "nsy", "U")
..diag2dge <- function(from) .diag2dense(from, "dge", NULL)
..diag2lge <- function(from) .diag2dense(from, "lge", NULL)
..diag2nge <- function(from) .diag2dense(from, "nge", NULL)

setAs("diagonalMatrix",          "dMatrix", ..diag2d)
setAs("diagonalMatrix",          "lMatrix", ..diag2l)
setAs("diagonalMatrix",          "nMatrix", ..diag2nsparse)

setAs("diagonalMatrix",    "dsparseMatrix", ..diag2dsparse)
setAs("diagonalMatrix",    "lsparseMatrix", ..diag2lsparse)
setAs("diagonalMatrix",    "nsparseMatrix", ..diag2nsparse)

setAs("diagonalMatrix",    "CsparseMatrix", ..diag2tC)
setAs("diagonalMatrix",    "RsparseMatrix", ..diag2tR)
setAs("diagonalMatrix",    "TsparseMatrix", ..diag2tT)

setAs("diagonalMatrix", "triangularMatrix", ..diag2tC)
setAs("diagonalMatrix",  "symmetricMatrix", ..diag2sC)
setAs("diagonalMatrix",    "generalMatrix", ..diag2gC)

setAs("diagonalMatrix",      "denseMatrix", ..diag2tr)
setAs("diagonalMatrix",   "unpackedMatrix", ..diag2tr)
setAs("diagonalMatrix",     "packedMatrix", ..diag2tp)
setAs("diagonalMatrix",     "ddenseMatrix", ..diag2dtr)
setAs("diagonalMatrix",     "ldenseMatrix", ..diag2ltr)
setAs("diagonalMatrix",     "ndenseMatrix", ..diag2ntr)
setAs("diagonalMatrix",           "matrix", .diag2m)
setAs("diagonalMatrix",           "vector", .diag2v)

setMethod("as.vector", signature(x = "diagonalMatrix"),
          function(x, mode) as.vector(.diag2v(x), mode))

setMethod("as.numeric", signature(x = "diagonalMatrix"),
          function(x, ...) as.double(.diag2v(x)))
setMethod("as.numeric", signature(x = "ddiMatrix"),
          function(x, ...) .diag2v(x))

setMethod("as.logical", signature(x = "diagonalMatrix"),
          function(x, ...) as.logical(.diag2v(x)))
setMethod("as.logical", signature(x = "ldiMatrix"),
          function(x, ...) .diag2v(x))

## DEPRECATED IN 1.4-2; see ./zzz.R
if(FALSE) {
.kinds <- c("d", "l")
for (.kind in .kinds) {
    ## ddi->[^d]di and similar
    for (.otherkind in .kinds[.kinds != .kind])
        setAs(paste0(     .kind, "diMatrix"),
              paste0(.otherkind, "diMatrix"),
              get(paste0("..diag2", .otherkind),
                  mode = "function", inherits = FALSE))
    ## ddi->d[tsg][CRT] and similar
    for (.x in c("t", "s", "g"))
        for (.y in c("C", "R", "T"))
            setAs(paste0(.kind, "diMatrix"),
                  paste0(.kind, .x, .y, "Matrix"),
                  get(paste0("..diag2", .x, .y),
                      mode = "function", inherits = FALSE))
    ## ddi->d(tr|sy|ge) and similar
    for (.xy in c("tr", "sy", "ge"))
        setAs(paste0(.kind, "diMatrix"),
              paste0(.kind, .xy, "Matrix"),
              get(paste0("..diag2", .kind, .xy),
                  mode = "function", inherits = FALSE))
}
rm(.kinds, .kind, .otherkind, .x, .y, .xy)
} ## DEPRECATED IN 1.4-2; see ./zzz.R

rm(..diag2dsparse, ..diag2lsparse, ..diag2nsparse,
   ..diag2tC, ..diag2tR, ..diag2tT,
   ..diag2sC, ..diag2sR, ..diag2sT,
   ..diag2gC, ..diag2gR, ..diag2gT,
   ..diag2tp, ..diag2tr,
   ..diag2dtr, ..diag2ltr, ..diag2ntr,
   ..diag2dsy, ..diag2lsy, ..diag2nsy,
   ..diag2dge, ..diag2lge, ..diag2nge)

## MJ: no longer needed ... replacement above
if(FALSE) {
.diag2tT <- function(from, uplo = "U", kind = .M.kind(from), drop0 = TRUE) {
    ## to triangular Tsparse
    x <- from@x
    i <- if(from@diag == "U")
	      integer(0L)
	  else if(drop0 & any0(x)) {
	      ii <- which(isN0(x))
	      x <- x[ii]
	      ii - 1L
	  }
	  else
	      seq_len(from@Dim[1]) - 1L
    new(paste0(kind, "tTMatrix"),
	diag = from@diag, Dim = from@Dim, Dimnames = from@Dimnames,
	uplo = uplo,
	x = x, # <- ok for diag = "U" and "N" (!)
	i = i, j = i)
}

.diag2sT <- function(from, uplo = "U", kind = .M.kind(from)) {
    ## to symmetric Tsparse
    n <- from@Dim[1]
    i <- seq_len(n) - 1L
    new(paste0(kind, "sTMatrix"),
	Dim = from@Dim, Dimnames = from@Dimnames,
	i = i, j = i, uplo = uplo,
	x = if(from@diag == "N") from@x else ## "U"-diag
	rep.int(switch(kind,
		       "d" = 1.,
		       "l" =,
		       "n" = TRUE,
		       ## otherwise
		       stop(gettextf("%s kind not yet implemented",
				     sQuote(kind)), domain=NA)),
		n))
}

## diagonal -> triangular,  upper / lower depending on "partner" 'x':
diag2tT.u <- function(d, x, kind = .M.kind(d), drop0 = TRUE)
    .diag2tT(d, uplo = if(is(x,"triangularMatrix")) x@uplo else "U", kind, drop0)

## diagonal -> sparse {triangular OR symmetric} (upper / lower) depending on "partner":
diag2Tsmart <- function(d, x, kind = .M.kind(d)) {
    clx <- getClassDef(class(x))
    if(extends(clx, "symmetricMatrix"))
	.diag2sT(d, uplo = x@uplo, kind)
    else
	.diag2tT(d, uplo = if(extends(clx,"triangularMatrix")) x@uplo else "U", kind)
}

## In order to evade method dispatch ambiguity warnings,
## and because we can save a .M.kind() call, we use this explicit
## "hack"  instead of signature  x = "diagonalMatrix" :
##
## ddi*:
di2tT <- function(from) .diag2tT(from, "U", "d")
setAs("ddiMatrix", "triangularMatrix", di2tT)
##_no_longer_ setAs("ddiMatrix", "sparseMatrix", di2tT)
## needed too (otherwise <dense> -> Tsparse is taken):
setAs("ddiMatrix", "TsparseMatrix", di2tT)
setAs("ddiMatrix", "dsparseMatrix", di2tT)
ddi2Csp <- function(from) .T2Cmat(.diag2tT(from, "U", "d"), isTri=TRUE) #-> dtC*
setAs("ddiMatrix", "dtCMatrix",     ddi2Csp)
setAs("ddiMatrix", "CsparseMatrix", ddi2Csp)
## Such that  as(Matrix(0, d,d), "dgCMatrix")  continues working:
setAs("ddiMatrix", "dgCMatrix", function(from) .dtC2g(ddi2Csp(from)))

setAs("ddiMatrix", "symmetricMatrix", function(from) .diag2sT(from, "U", "d"))
##
## ldi*:
ldi2tT <- function(from) .diag2tT(from, "U", "l")
setAs("ldiMatrix", "triangularMatrix", ldi2tT)
##_no_longer_ setAs("ldiMatrix", "sparseMatrix", di2tT)
## needed too (otherwise <dense> -> Tsparse is taken):
setAs("ldiMatrix", "TsparseMatrix", ldi2tT)
setAs("ldiMatrix", "lsparseMatrix", ldi2tT)
setAs("ldiMatrix", "CsparseMatrix",
      function(from) .T2Cmat(.diag2tT(from, "U", "l"), isTri=TRUE))
setAs("ldiMatrix", "symmetricMatrix", function(from) .diag2sT(from, "U", "l"))
rm(ldi2tT)

setAs("diagonalMatrix", "nMatrix",
      di2nMat <- function(from) {
	  i <- if(from@diag == "U") integer(0) else which(isN0(from@x)) - 1L
	  new("ntTMatrix", i = i, j = i, diag = from@diag,
	      Dim = from@Dim, Dimnames = from@Dimnames)
      })
setAs("diagonalMatrix", "nsparseMatrix", function(from) as(from, "nMatrix"))

##' A version of diag(x,n) which *does* preserve the mode of x, where diag() "fails"
mkDiag <- function(x, n) {
    y <- matrix(as0(mod=mode(x)), n,n)
    if (n > 0) y[1L + 0:(n - 1L) * (n + 1)] <- x
    y
}
## NB: diag(x,n) is really faster for n >= 20, and even more for large n
## --> using diag() where possible, ==> .ddi2mat()

.diag2mat <- function(from)
    ## want "ldiMatrix" -> <logical> "matrix"  (but integer -> <double> for now)
    mkDiag(if(from@diag == "U") as1(from@x) else from@x, n = from@Dim[1])

.ddi2mat <- function(from)
    `dimnames<-`(base::diag(if(from@diag == "U") as1(from@x) else from@x, nrow = from@Dim[1]),
                 from@Dimnames)

setAs("ddiMatrix", "matrix", .ddi2mat)
## the non-ddi diagonalMatrix -- only "ldiMatrix" currently:
setAs("diagonalMatrix", "matrix", .diag2mat)

setAs("diagonalMatrix", "generalMatrix", # prefer sparse:
      function(from) as(as(from, "CsparseMatrix"), "generalMatrix"))

setAs("diagonalMatrix", "denseMatrix",
      function(from) as(as(from, "CsparseMatrix"), "denseMatrix"))

setAs("ddiMatrix", "dgeMatrix", function(from) ..2dge(from))

setAs("ddiMatrix", "ddenseMatrix", #-> "dtr"
      function(from) as(as(from, "triangularMatrix"),"denseMatrix"))
setAs("ldiMatrix", "ldenseMatrix", #-> "ltr"
      function(from) as(as(from, "triangularMatrix"),"denseMatrix"))
} ## MJ


## ~~~~ CONSTRUCTORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Purpose: Constructor of diagonal matrices -- ~= diag() ,
##          but *not* diag() extractor!
Diagonal <- function(n, x = NULL)
{
    ## Allow  Diagonal(4), Diagonal(x=1:5), and  Diagonal(4, TRUE)
    n <- if(missing(n)) length(x) else {
	stopifnot(length(n) == 1, n == as.integer(n), n >= 0)
	as.integer(n)
    }

    if(missing(x)) ## unit diagonal matrix
	new("ddiMatrix", Dim = c(n,n), diag = "U")
    else {
	lx <- length(x)
	lx.1 <- lx == 1L
	stopifnot(lx.1 || lx == n) # but keep 'x' short for now
	if(is.logical(x))
	    cl <- "ldiMatrix"
	else if(is.numeric(x)) {
	    cl <- "ddiMatrix"
	    x <- as.numeric(x)
	}
	else if(is.complex(x)) {
	    cl <- "zdiMatrix"  # will not yet work
	} else stop("'x' has invalid data type")
	if(lx.1 && !is.na(x) && x == 1) # cheap check for uni-diagonal..
	    new(cl, Dim = c(n,n), diag = "U")
	else
	    new(cl, Dim = c(n,n), diag = "N",
		x = if(lx.1) rep.int(x,n) else x)
    }
}

.sparseDiagonal <- function(n, x = 1, uplo = "U",
			    shape = if(missing(cols)) "t" else "g",
			    unitri, kind,
			    cols = if(n) 0:(n - 1L) else integer(0))
{
    n <- if (missing(n)) length(x) else {
	stopifnot(length(n) == 1, n == as.integer(n), n >= 0)
	as.integer(n)
    }
    if(!(mcols <- missing(cols)))
	stopifnot(0 <= (cols <- as.integer(cols)), cols < n)
    m <- length(cols)
    if(missing(kind))
	kind <-
	    if(is.double(x)) "d"
	    else if(is.logical(x)) "l"
	    else { ## for now
		storage.mode(x) <- "double"
		"d"
	    }
    else stopifnot(any(kind == c("d","l","n")))
    stopifnot(is.character(shape), nchar(shape) == 1,
	      any(shape == c("t","s","g"))) # triangular / symmetric / general
    if((missing(unitri) || unitri) && shape == "t" &&
       (mcols || cols == 0:(n-1L)) &&
       ((any(kind == c("l", "n")) && allTrue(x)) ||
	(    kind == "d"	  && allTrue(x == 1)))) { ## uni-triangular
	new(paste0(kind,"tCMatrix"), Dim = c(n,n),
		   uplo = uplo, diag = "U", p = rep.int(0L, n+1L))
    }
    else if(kind == "n") {
	if(shape == "g")
	    new("ngCMatrix", Dim = c(n,m), i = cols, p = 0:m)
	else new(paste0("n", shape, "CMatrix"), Dim = c(n,m), uplo = uplo,
		 i = cols, p = 0:m)
    }
    else { ## kind != "n" -- have x slot :
	if((lx <- length(x)) == 1) x <- rep.int(x, m)
	else if(lx != m) stop("length(x) must be either 1 or #{cols}")
	if(shape == "g")
	    new(paste0(kind, "gCMatrix"), Dim = c(n,m),
		x = x, i = cols, p = 0:m)
	else new(paste0(kind, shape, "CMatrix"), Dim = c(n,m), uplo = uplo,
		 x = x, i = cols, p = 0:m)
    }
}

## Pkg 'spdep' had (relatively slow) versions of this as_dsCMatrix_I()
.symDiagonal <- function(n, x = rep.int(1,n), uplo = "U", kind)
    .sparseDiagonal(n, x, uplo, shape = "s", kind = kind)

## NOTA BENE: .triDiagonal() would be misleading (<=> banded tri-diagonal matrix !)
# instead of   diagU2N(as(Diagonal(n), "CsparseMatrix")), diag = "N" in any case:
.trDiagonal <- function(n, x = 1, uplo = "U", unitri = TRUE, kind)
    .sparseDiagonal(n, x, uplo, shape = "t", unitri=unitri, kind=kind)


## This is modified from a post of Bert Gunter to R-help on  1 Sep 2005.
## Bert's code built on a post by Andy Liaw who most probably was influenced
## by earlier posts, notably one by Scott Chasalow on S-news, 16 Jan 2002
## who posted his bdiag() function written in December 1995.
if(FALSE)##--- no longer used:
.bdiag <- function(lst) {
    ## block-diagonal matrix [a dgTMatrix] from list of matrices
    stopifnot(is.list(lst), length(lst) >= 1)
    dims <- vapply(lst, dim, 1L, USE.NAMES=FALSE)
    ## make sure we had all matrices:
    if(!(is.matrix(dims) && nrow(dims) == 2))
	stop("some arguments are not matrices")
    csdim <- rbind(rep.int(0L, 2),
                   apply(dims, 1, cumsum))
    r <- new("dgTMatrix")
    r@Dim <- as.integer(csdim[nrow(csdim),])
    add1 <- matrix(1:0, 2,2)
    for(i in seq_along(lst)) {
	indx <- apply(csdim[i:(i+1),] + add1, 2, function(n) n[1]:n[2])
	if(is.null(dim(indx))) ## non-square matrix
	    r[indx[[1]],indx[[2]]] <- lst[[i]]
	else ## square matrix
	    r[indx[,1], indx[,2]] <- lst[[i]]
    }
    r
}
## expand(<mer>) needed something like bdiag() for lower-triangular
## (Tsparse) Matrices; hence Doug Bates provided a much more efficient
##  implementation for those; now extended and generalized:
.bdiag <- function(lst) {
    ## block-diagonal matrix [a dgTMatrix] from list of matrices
    stopifnot(is.list(lst), (nl <- length(lst)) >= 1L)

### FIXME: next line is *slow* when lst = list of 75'000  dense 3x3 matrices
    Tlst <- lapply(unname(lst), function(x) .CR2T(asCspN(x)))
    if(nl == 1L)
        return(Tlst[[1L]])
    ## else
    i_off <- c(0L, cumsum(vapply(Tlst, function(x) x@Dim[1L], 1L)))
    j_off <- c(0L, cumsum(vapply(Tlst, function(x) x@Dim[2L], 1L)))

    clss <- vapply(Tlst, class, "")
    ## NB ("FIXME"): this requires the component classes to be *called*
    ## -- "dgTMatrix" | "dnTMatrix" etc (and not just *extend* those)!
    typ <- substr(clss, 2L, 2L)
    knd <- substr(clss, 1L, 1L)
    sym <- typ == "s" # symmetric ones
    tri <- typ == "t" # triangular ones
    use.n <- any(is.n <- knd == "n")
    if(use.n && !(use.n <- all(is.n))) {
	Tlst[is.n] <- lapply(Tlst[is.n], ..sparse2l)
	knd [is.n] <- "l"
    }
    use.l <- !use.n && all(knd == "l")
    if(all(sym)) { ## result should be *symmetric*
	uplos <- vapply(Tlst, slot, "", "uplo") ## either "U" or "L"
	tLU <- table(uplos)# of length 1 or 2 ..
	if(length(tLU) == 1L) { ## all "U" or all "L"
	    useU <- uplos[1L] == "U"
	} else { ## length(tLU) == 2, counting "L" and "U"
	    useU <- diff(tLU) >= 0L
	    if(useU && (hasL <- tLU[1L] > 0L))
		Tlst[hasL] <- lapply(Tlst[hasL], t)
	    else if(!useU && (hasU <- tLU[2L] > 0L))
		Tlst[hasU] <- lapply(Tlst[hasU], t)
	}
	if(use.n) { ## return nsparseMatrix :
	    r <- new("nsTMatrix")
	} else {
	    r <- new(paste0(if(use.l) "l" else "d", "sTMatrix"))
	    r@x <- unlist(lapply(Tlst, slot, "x"), FALSE, FALSE)
	}
	r@uplo <- if(useU) "U" else "L"
    }
    else if(all(tri) && { ULs <- vapply(Tlst, slot, "", "uplo")##  "U" or "L"
			  all(ULs[1L] == ULs[-1L]) } ## all upper or all lower
       ){ ## *triangular* result

	if(use.n) { ## return nsparseMatrix :
	    r <- new("ntTMatrix")
	} else {
	    r <- new(paste0(if(use.l) "l" else "d", "tTMatrix"))
	    r@x <- unlist(lapply(Tlst, slot, "x"), FALSE, FALSE)
	}
	r@uplo <- ULs[1L]
    }
    else {
	if(any(sym))
	    Tlst[sym] <- lapply(Tlst[sym], .sparse2g)
	if(use.n) { ## return nsparseMatrix :
	    r <- new("ngTMatrix")
	} else {
	    r <- new(paste0(if(use.l) "l" else "d", "gTMatrix"))
	    r@x <- unlist(lapply(Tlst, slot, "x"), FALSE, FALSE)
	}
    }
    r@Dim <- c(i_off[nl+1], j_off[nl + 1])
    r@i <- unlist(lapply(1:nl, function(k) Tlst[[k]]@i + i_off[k]),
                  FALSE, FALSE)
    r@j <- unlist(lapply(1:nl, function(k) Tlst[[k]]@j + j_off[k]),
                  FALSE, FALSE)
    r
}

bdiag <- function(...) {
    if((nA <- nargs()) == 0L) return(new("dgCMatrix"))
    if(nA == 1L && !is.list(...))
	return(as(..., "CsparseMatrix"))
    alis <- if(nA == 1L && is.list(..1)) ..1 else list(...)
    if(length(alis) == 1L)
	return(as(alis[[1L]], "CsparseMatrix"))
    ## else : two or more arguments
    .T2C(.bdiag(alis))
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("band", signature(x = "diagonalMatrix"),
          function(x, k1, k2, ...)
              if(k1 <= 0L && k2 >= 0L) x else .setZero(x))

setMethod("tril", signature(x = "diagonalMatrix"),
          function(x, k = 0, ...)
              if(k >= 0L) x else .setZero(x))

setMethod("triu", signature(x = "diagonalMatrix"),
          function(x, k = 0, ...)
              if(k <= 0L) x else .setZero(x))

setMethod("forceSymmetric", signature(x = "diagonalMatrix", uplo = "character"),
          function(x, uplo) .diag2sparse(x, ".sC", uplo = uplo))

setMethod("forceSymmetric", signature(x = "diagonalMatrix", uplo = "missing"),
          function(x, uplo) .diag2sparse(x, ".sC", uplo = "U"))

..diag.x <- function(m)                   rep.int(as1(m@x), m@Dim[1])
.diag.x  <- function(m) if(m@diag == "U") rep.int(as1(m@x), m@Dim[1]) else m@x

.diag.2N <- function(m) {
    if(m@diag == "U") m@diag <- "N"
    m
}

diag.x <- function(x, nrow, ncol, names=TRUE) {
    y <- .diag.x(x)
    if(names) {
        nms <- dimnames(x)
        if(is.list(nms) && !any(vapply(nms, is.null, NA)) &&
           identical((nm <- nms[[1L]][im <- seq_len(min(dim(x)))]), nms[[2L]][im]))
            names(y) <- nm
    }
    y
}
setMethod("diag", signature(x = "diagonalMatrix"), diag.x)

subDiag <- function(x, i, j, ..., drop) {
    x <- .diag2sparse(x, ".gC") ## was ->TsparseMatrix but C* is faster now
    x <- if(missing(i))
	x[, j, drop=drop]
    else if(missing(j))
	if(nargs() == 4L) x[i, , drop=drop] else x[i, drop=drop]
    else
	x[i,j, drop=drop]
    if(isS4(x) && isDiagonal(x)) as(x, "diagonalMatrix") else x
}

setMethod("[", signature(x = "diagonalMatrix", i = "index",
			 j = "index", drop = "logical"), subDiag)
setMethod("[", signature(x = "diagonalMatrix", i = "index",
			 j = "missing", drop = "logical"),
	  function(x, i, j, ..., drop) {
	      na <- nargs()
	      Matrix.msg("diag[i,m,l] : nargs()=", na, .M.level = 2)
	      if(na == 4L)
		   subDiag(x, i=i, , drop=drop)
	      else subDiag(x, i=i,   drop=drop)
	  })
setMethod("[", signature(x = "diagonalMatrix", i = "missing",
			 j = "index", drop = "logical"),
	  function(x, i, j, ..., drop) subDiag(x, j=j, drop=drop))

## When you assign to a diagonalMatrix, the result should be
## diagonal or sparse ---
replDiag <- function(x, i, j, ..., value) {
## FIXME: if   (i == j)  &&  isSymmetric(value) then -- want symmetricMatrix result! -- or diagMatrix
    x <- .diag2sparse(x, ".gC") # was ->TsparseMatrix till 2012-07
    if(missing(i))
	x[, j] <- value
    else if(missing(j)) { ##  x[i , ] <- v  *OR*   x[i] <- v
        na <- nargs()
##         message("diagnosing replDiag() -- nargs()= ", na)
	if(na == 4L)
            x[i, ] <- value
	else if(na == 3L)
            x[i] <- value
	else stop(gettextf("Internal bug: nargs()=%d; please report",
			   na), domain=NA)
    } else
	x[i,j] <- value
    ## TODO: the following is a bit expensive; have cases above e.g. [i,] where
    ## ----- we could check *much* faster :
    if(isDiagonal(x))
        .M2diag(x, check = FALSE)
    else if(isSymmetric(x))
        forceSymmetric(x)
    else if(!(it <- isTriangular(x)))
        x
    else if(attr(it, "kind") == "U")
        triu(x)
    else tril(x)
}

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index",
				j = "index", value = "replValue"), replDiag)

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index",
				j = "missing", value = "replValue"),
		 function(x,i,j, ..., value) {
                     ## message("before replDiag() -- nargs()= ", nargs())
                     if(nargs() == 3L)
                         replDiag(x, i=i, value=value)
                     else ## nargs() == 4 :
                         replDiag(x, i=i, , value=value)
                 })

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing",
				j = "index", value = "replValue"),
		 function(x,i,j, ..., value) replDiag(x, j=j, value=value))

## x[] <- value :
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing",
				j = "missing", value = "ANY"),
		 function(x,i,j, ..., value)
	     {
	      if(all0(value)) { # be faster
		  r <- new(paste0(.M.kind(x), "tTMatrix")) # of all "0"
		  r@Dim <- x@Dim
		  r@Dimnames <- x@Dimnames
		  r
	      } else { ## typically non-sense: assigning to full sparseMatrix
		  x[TRUE] <- value
		  x
	      }
	  })


setReplaceMethod("[", signature(x = "diagonalMatrix",
                                i = "matrix", # 2-col.matrix
				j = "missing", value = "replValue"),
		 function(x,i,j, ..., value) {
		     if(ncol(i) == 2L) {
			 if(all((ii <- i[,1L]) == i[,2L])) { # replace in diagonal only
			     if(x@diag == "U") {
				 one <- as1(x@x)
				 if(any(value != one | is.na(value))) {
				     x@diag <- "N"
				     x@x <- rep.int(one, x@Dim[1L])
				 } else return(x)
			     }
			     x@x[ii] <- value
			     x
			 } else { ## no longer diagonal, but remain sparse:
###  FIXME:  use  uplo="U" or uplo="L"  (or *not* "triangularMatrix") depending on LE <- i <= j
###          all(LE) //  all(!LE) // remaining cases
			     x <- .diag2sparse(x, ".tC") # was ->TsparseMatrix
			     x[i] <- value
			     x
			 }
		     }
		     else { # behave as "base R": use as if vector
			 x <- as(x, "matrix")
			 x[i] <- value
			 Matrix(x)
		     }
		 })


## value = "sparseMatrix":
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing", j = "index",
				value = "sparseMatrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, , j=j, value = as(value, "sparseVector")))
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "missing",
				value = "sparseMatrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, i=i, , value = as(value, "sparseVector")))
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "index",
				value = "sparseMatrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, i=i, j=j, value = as(value, "sparseVector")))

## value = "sparseVector":
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "missing", j = "index",
				value = "sparseVector"),
		 replDiag)
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "missing",
				value = "sparseVector"),
		 replDiag)
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "index",
				value = "sparseVector"),
		 replDiag)


setMethod("t", signature(x = "diagonalMatrix"),
          function(x) { x@Dimnames <- x@Dimnames[2:1]; x })

setMethod("isSymmetric", signature(object = "diagonalMatrix"),
          function(object, checkDN = TRUE, ...) {
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...) check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              TRUE
          })

setMethod("isTriangular", signature(object = "diagonalMatrix"),
          function(object, upper = NA, ...)
              if(is.na(upper)) `attr<-`(TRUE, "kind", "U") else TRUE)

setMethod("isDiagonal", signature(object = "diagonalMatrix"),
          function(object) TRUE)

setMethod("symmpart", signature(x = "diagonalMatrix"),
          function(x) forceSymmetric(..diag2d(x)))

setMethod("skewpart", signature(x = "diagonalMatrix"),
          function(x) symmetrizeDimnames(.setZero(x, "d")))

## FIXME: Many of these products are not handling 'Dimnames' appropriately ...

## Basic Matrix Multiplication {many more to add}
##       ---------------------
## Note that "ldi" logical are treated as numeric
diagdiagprod <- function(x, y) {
    dimCheck(x,y)
    if(x@diag != "U") {
	if(y@diag != "U") {
	    nx <- x@x * y@x
	    if(is.numeric(nx) && !is.numeric(x@x))
		x <- ..diag2d(x)
	    x@x <- as.double(nx)
	}
	x
    } else ## x is unit diagonal
	y
}
setMethod("%*%", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
	  diagdiagprod)

##' Boolean Algebra/Arithmetic Product of Diagonal Matrices
##'  %&%
diagdiagprodBool <- function(x, y) {
    dimCheck(x,y)
    if(x@diag != "U") {
	if(!is.logical(x@x))
            x <- ..diag2l(x)
	if(y@diag != "U") {
	    nx <- x@x & y@x
	    x@x <- as.logical(nx)
	}
	x
    } else { ## x is unit diagonal: return y
	if(!is.logical(y@x))
            y <- ..diag2l(y)
	y
    }
}
setMethod("%&%", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
	  diagdiagprodBool) # giving ldiMatrix as do *not* have "ndiMatrix" !

##' Both Numeric or Boolean Algebra/Arithmetic Product of Diagonal Matrices
diagdiagprodFlexi <- function(x, y=NULL, boolArith = NA, ...)
{
    dimCheck(x,y)
    bool <- isTRUE(boolArith)
    if(x@diag != "U") {
	if(bool && !is.logical(x@x))
            x <- ..diag2l(x)
	if(y@diag != "U") {
	    if(bool) {
		nx <- x@x & y@x
		x@x <- as.logical(nx)
	    } else { ## boolArith is NA or FALSE: ==> numeric, as have *no* "diagMatrix" patter[n]:
		nx <- x@x * y@x
		if(is.numeric(nx) && !is.numeric(x@x))
		    x <- ..diag2d(x)
		x@x <- as.double(nx)
	    }
	}
	x
    } else { ## x is unit diagonal: return y
	if(bool && !is.logical(y@x))
            y <- ..diag2l(y)
	y
    }
}
setMethod("crossprod", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
	  diagdiagprodFlexi)
setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
	  diagdiagprodFlexi)

##' crossprod(x) := x'x
diagprod <- function(x, y = NULL, boolArith = NA, ...) {
    bool <- isTRUE(boolArith)
    if(bool && !is.logical(x@x))
        x <- ..diag2l(x)
    if(x@diag != "U") {
        if(bool) {
            nx <- x@x & y@x
            x@x <- as.logical(nx)
        } else { ## boolArith is NA or FALSE: ==> numeric, as have *no* "diagMatrix" patter[n]:
            nx <- x@x * x@x
            if(is.numeric(nx) && !is.numeric(x@x))
                x <- ..diag2d(x)
            x@x <- as.double(nx)
        }
    }
    x
}
setMethod( "crossprod", signature(x = "diagonalMatrix", y = "missing"),
          diagprod)
setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "missing"),
          diagprod)

## analogous to matdiagprod() below:
diagmatprod <- function(x, y) {
    ## x is diagonalMatrix
    dy <- dim(y)
    if(x@Dim[2L] != dy[1L]) stop("non-matching dimensions")
    if(prod(dy))
	Matrix(if(x@diag == "U") y else x@x * y)
    else
	Matrix(if(x@diag == "U") y else x@x * y, nrow=dy[1L], ncol=dy[2L])
}
setMethod("%*%", signature(x = "diagonalMatrix", y = "matrix"), diagmatprod)

##formals(diagmatprod) <- alist(x=, y=NULL, boolArith = NA, ...=) ## FIXME boolArith
diagmatprod2 <- function(x, y=NULL, boolArith = NA, ...) {
    ## x is diagonalMatrix
    dy <- dim(y)
    if(x@Dim[2L] != dy[1L]) stop("non-matching dimensions")
    if(prod(dy))
	Matrix(if(x@diag == "U") y else x@x * y)
    else
	Matrix(if(x@diag == "U") y else x@x * y, nrow=dy[1L], ncol=dy[2L])
}
setMethod("crossprod",  signature(x = "diagonalMatrix", y = "matrix"), diagmatprod2)

diagGeprod <- function(x, y) {
    if(x@Dim[2L] != y@Dim[1L]) stop("non-matching dimensions")
    if(x@diag != "U") {
        if(!is.numeric(y@x))
            y <- ..dense2d(y)
        y@x <- x@x * y@x
    }
    y
}
setMethod("%*%", signature(x= "diagonalMatrix", y= "dgeMatrix"), diagGeprod)
setMethod("%*%", signature(x= "diagonalMatrix", y= "lgeMatrix"), diagGeprod)

diagGeprodBool <- function(x, y) {
    if(x@Dim[2L] != y@Dim[1L]) stop("non-matching dimensions")
    if(!is.logical(y@x)) # MJ: hmm ... what about nge* with NA in 'x' slot?
        y <- ..dense2l(y)
    if(x@diag != "U")
        y@x <- x@x & y@x
    y
}
setMethod("%&%", signature(x= "diagonalMatrix", y= "geMatrix"), diagGeprodBool)

diagGeprod2 <- function(x, y=NULL, boolArith = NA, ...) {
    if(x@Dim[2L] != y@Dim[1L]) stop("non-matching dimensions")
    bool <- isTRUE(boolArith)
    if(bool && !is.logical(y@x))
        y <- ..dense2l(y)
    else if(!bool && !is.numeric(y@x))
        y <- ..dense2d(y)
    if(x@diag != "U")
        y@x <- if(bool) x@x & y@x else x@x * y@x
    y
}
setMethod("crossprod", signature(x = "diagonalMatrix", y = "dgeMatrix"),
          diagGeprod2)
setMethod("crossprod", signature(x = "diagonalMatrix", y = "lgeMatrix"),
          diagGeprod2)

## analogous to diagmatprod() above:
matdiagprod <- function(x, y) {
    dx <- dim(x)
    if(dx[2L] != y@Dim[1L]) stop("non-matching dimensions")
    Matrix(if(y@diag == "U") x else x * rep(y@x, each = dx[1L]))
}
setMethod("%*%", signature(x = "matrix", y = "diagonalMatrix"), matdiagprod)

gediagprod <- function(x, y) {
    dx <- dim(x)
    if(dx[2L] != y@Dim[1L]) stop("non-matching dimensions")
    if(y@diag == "N") {
        if(!is.numeric(x@x)) x <- as(x, "dMatrix")
	x@x <- x@x * rep(y@x, each = dx[1L])
    }
    x
}
setMethod("%*%", signature(x= "dgeMatrix", y= "diagonalMatrix"), gediagprod)
setMethod("%*%", signature(x= "lgeMatrix", y= "diagonalMatrix"), gediagprod)

gediagprodBool <- function(x, y) {
    dx <- dim(x)
    if(dx[2L] != y@Dim[1L]) stop("non-matching dimensions")
    if(!is.logical(x@x))
        x <- ..dense2l(x)
    if(y@diag == "N")
	x@x <- x@x & rep(y@x, each = dx[1L])
    x
}
setMethod("%&%", signature(x= "geMatrix", y= "diagonalMatrix"), gediagprodBool)

setMethod("tcrossprod",signature(x = "matrix", y = "diagonalMatrix"),
          function(x, y=NULL, boolArith = NA, ...) {
              dx <- dim(x)
              if(dx[2L] != y@Dim[1L]) stop("non-matching dimensions")
              bool <- isTRUE(boolArith)
              if(bool && !is.logical(y@x))
                  y <- ..diag2l(y)
              Matrix(if(y@diag == "U") x else
                     if(bool) x & rep(y@x, each = dx[1L])
                     else     x * rep(y@x, each = dx[1L]))
          })

setMethod("crossprod", signature(x = "matrix", y = "diagonalMatrix"),
	  function(x, y=NULL, boolArith = NA, ...) {
	      dx <- dim(x)
	      if(dx[1L] != y@Dim[1L]) stop("non-matching dimensions")
              bool <- isTRUE(boolArith)
              if(bool && !is.logical(y@x))
                  y <- ..diag2l(y)
	      Matrix(if(y@diag == "U") t(x) else
		     if(bool) t(rep.int(y@x, dx[2L]) & x)
		     else     t(rep.int(y@x, dx[2L]) * x))
	  })


gediagprod2 <- function(x, y=NULL, boolArith = NA, ...) {
    dx <- dim(x)
    if(dx[2L] != y@Dim[1L]) stop("non-matching dimensions")
    bool <- isTRUE(boolArith)
    if(bool && !is.logical(x@x))
        x <- ..dense2l(x)
    else if(!bool && !is.numeric(x@x))
        x <- ..dense2d(x)
    if(y@diag == "N")
	x@x <- if(bool) x@x & rep(y@x, each = dx[1L])
	       else     x@x * rep(y@x, each = dx[1L])
    x
}
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "diagonalMatrix"), gediagprod2)
setMethod("tcrossprod", signature(x = "lgeMatrix", y = "diagonalMatrix"), gediagprod2)


## crossprod {more of these}

## tcrossprod --- all are not yet there: do the dense ones here:

setMethod("%*%", signature(x = "diagonalMatrix", y = "denseMatrix"),
	  function(x, y) if(x@diag == "U") y else x %*% .dense2g(y))
setMethod("%*%", signature(x = "denseMatrix", y = "diagonalMatrix"),
	  function(x, y) if(y@diag == "U") x else .dense2g(x) %*% y)


## FIXME:
## setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "denseMatrix"),
## 	  function(x, y = NULL) {
##           })

##' @param x CsparseMatrix
##' @param y diagonalMatrix
##' @return x %*% y
Cspdiagprod <- function(x, y, boolArith = NA, ...) {
    if((m <- ncol(x)) != y@Dim[1L]) stop("non-matching dimensions")
    if(y@diag == "N") { ## otherwise: y == Diagonal(n) : multiplication is identity
	x <- .Call(Csparse_diagU2N, x)
	cx <- getClass(class(x))
	if(!all(y@x[1L] == y@x[-1L]) && extends(cx, "symmetricMatrix"))
	    x <- as(x, "generalMatrix")
	ind <- rep.int(seq_len(m), x@p[-1] - x@p[-m-1L])
	if(isTRUE(boolArith)) {
	    if(extends(cx, "nMatrix")) x <- as(x, "lMatrix") # so, has y@x
	    x@x <- r <- x@x & y@x[x@i + 1L]
	    if(!anyNA(r) && !extends(cx, "diagonalMatrix")) x <- as(drop0(x), "nMatrix")
	} else {
	    if(!extends(cx, "dMatrix")) x <- as(x, "dMatrix") # <- FIXME if we have zMatrix
	    x@x <- x@x * y@x[ind]
	}
        if(.hasSlot(x, "factors") && length(x@factors)) {# drop cashed ones
	    ## instead of dropping all factors, be smart about some
	    ## TODO ......
	    x@factors <- list()
	}
        x
    } else { #	y is unit-diagonal ==> "return x"
	cx <- getClass(class(x))
	if(isTRUE(boolArith)) {
	    is.l <- if(extends(cx, "dMatrix")) { ## <- FIXME: extend once we have iMatrix, zMatrix
		x <- as(x, "lMatrix"); TRUE } else extends(cx, "lMatrix")
	    if(is.l && !anyNA(x@x)) as(drop0(x), "nMatrix")
	    else if(is.l) x else # defensive:
	    as(x, "lMatrix")
	} else {
	    ## else boolArith is  NA or FALSE {which are equivalent here, das diagonal = "numLike"}
	    if(extends1of(cx, c("nMatrix", "lMatrix")))
		as(x, "dMatrix") else x
	}
    }
}

##' @param x diagonalMatrix
##' @param y CsparseMatrix
##' @return x %*% y
diagCspprod <- function(x, y, boolArith = NA, ...) {
    if(x@Dim[2L] != y@Dim[1L]) stop("non-matching dimensions")
    if(x@diag == "N") {
	y <- .Call(Csparse_diagU2N, y)
	cy <- getClass(class(y))
	if(!all(x@x[1L] == x@x[-1L]) && extends(cy, "symmetricMatrix"))
	    y <- .sparse2g(y)
	if(isTRUE(boolArith)) {
	    if(extends(cy, "nMatrix"))
                y <- ..sparse2l(y) # so, has y@x
	    y@x <- r <- y@x & x@x[y@i + 1L]
	    if(!anyNA(r) && !extends(cy, "diagonalMatrix"))
                y <- ..sparse2n(drop0(y))
	} else {
	    if(!extends(cy, "dMatrix"))
                y <- ..sparse2d(y) # <- FIXME if we have zMatrix
	    y@x <- y@x * x@x[y@i + 1L]
	}
	if(.hasSlot(y, "factors") && length(y@factors)) {
            ## if(.hasSlot(y, "factors") && length(yf <- y@factors)) { ## -- TODO? --
	    ## instead of dropping all factors, be smart about some
	    ## keep <- character()
	    ## if(any(names(yf) == "LU")) { ## <- not easy: y = P'LUQ,  x y = xP'LUQ => LU ???
	    ##     keep <- "LU"
	    ## }
	    ## y@factors <- yf[keep]
	    y@factors <- list()
        }
        y
    } else { ## x @ diag  == "U"
	cy <- getClass(class(y))
	if(isTRUE(boolArith)) {
	    is.l <-
                if(extends(cy, "dMatrix")) { ## <- FIXME: extend once we have iMatrix, zMatrix
                    y <- ..sparse2l(y)
                    TRUE
                }      else extends(cy, "lMatrix")
            if(is.l && !anyNA(y@x))
                ..sparse2n(drop0(y))
            else if(is.l)
                y
            else ..sparse2l(y) # defensive
        } else {
            ## else boolArith is  NA or FALSE {which are equivalent here, das diagonal = "numLike"}
            if(extends1of(cy, c("nMatrix", "lMatrix")))
                ..sparse2d(y)
            else y
        }
    }
}

## + 'boolArith' argument  { ==> .local() is used in any case; keep formals simple :}
setMethod("crossprod", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
	  function(x, y = NULL, boolArith = NA, ...)
              diagCspprod(x, y, boolArith = boolArith))

setMethod("crossprod", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y = NULL, boolArith = NA, ...)
	      diagCspprod(x, as(y, "CsparseMatrix"), boolArith = boolArith))

## Prefer calling diagCspprod to Cspdiagprod if going to transpose anyway
##  x'y == (y'x)'
setMethod("crossprod", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
	  function(x, y = NULL, boolArith = NA, ...)
              t(diagCspprod(y, x, boolArith = boolArith)))

setMethod("crossprod", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y = NULL, boolArith = NA, ...)
              t(diagCspprod(y, as(x, "Csparsematrix"), boolArith = boolArith)))

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
	  function(x, y = NULL, boolArith = NA, ...)
              diagCspprod(x, t(y), boolArith = boolArith))

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y=NULL, boolArith=NA, ...) diagCspprod(x, t(as(y, "CsparseMatrix")), boolArith=boolArith))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
	  function(x, y = NULL, boolArith = NA, ...)
              Cspdiagprod(x, y, boolArith = boolArith))

setMethod("tcrossprod", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y = NULL, boolArith = NA, ...)
              Cspdiagprod(as(x, "CsparseMatrix"), y, boolArith = boolArith))

setMethod("%*%", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
	  function(x, y) diagCspprod(x, y, boolArith=NA))
setMethod("%&%", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
	  function(x, y) diagCspprod(x, y, boolArith=TRUE))

## instead of "sparseMatrix", use: [RT]sparse.. ("closer" in method dispatch)
for(cl in c("TsparseMatrix", "RsparseMatrix")) {

setMethod("%*%", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y) diagCspprod(as(x, "CsparseMatrix"), y, boolArith=NA))

setMethod("%*%", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y) Cspdiagprod(as(x, "CsparseMatrix"), y, boolArith=NA))

setMethod("%&%", signature(x = "diagonalMatrix", y = "sparseMatrix"),
	  function(x, y) diagCspprod(as(x, "CsparseMatrix"), y, boolArith=TRUE))

setMethod("%&%", signature(x = "sparseMatrix", y = "diagonalMatrix"),
	  function(x, y) Cspdiagprod(as(x, "CsparseMatrix"), y, boolArith=TRUE))
}
rm(cl)

setMethod("%*%", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
	  function(x, y) Cspdiagprod(x, y, boolArith=NA))
setMethod("%&%", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
	  function(x, y) Cspdiagprod(x, y, boolArith=TRUE))

## TODO: Write tests in ./tests/ which ensure that many "ops" with diagonal*
##       do indeed work by going through sparse (and *not* ddense)!


###---------------- <Ops> (<Arith>, <Logic>, <Compare> ) ----------------------

## Use as S4 method for several signatures ==>  using callGeneric()
diagOdiag <- function(e1,e2) {
    ## result should also be diagonal _ if possible _
    r <- callGeneric(.diag.x(e1), .diag.x(e2)) # error if not "compatible"
    ## Check what happens with non-diagonals, i.e. (0 o 0), (FALSE o 0), ...:
    r00 <- callGeneric(if(is.numeric(e1@x)) 0 else FALSE,
		       if(is.numeric(e2@x)) 0 else FALSE)
    if(is0(r00)) { ##  r00 == 0 or FALSE --- result *is* diagonal
	if(is.numeric(r)) { # "double" *or* "integer"
            if(!is.double(r))
                r <- as.double(r)
            if(is.double(e2@x)) {
		e2@x <- r
                return(.diag.2N(e2))
            }
	    if(!is.double(e1@x))
		## e.g. e1, e2 are logical;
		e1 <- ..diag2d(e1)
	}
	else if(is.logical(r))
	    e1 <- ..diag2l(e1)
	else stop(gettextf("intermediate 'r' is of type %s",
			   typeof(r)), domain=NA)
	e1@x <- r
	.diag.2N(e1)
    }
    else { ## result not diagonal, but at least symmetric:
        ## e.g., m == m
	isNum <- (is.numeric(r) || is.numeric(r00))
	isLog <- (is.logical(r) || is.logical(r00))
        Matrix.msg("exploding <diag> o <diag> into dense matrix", .M.level = 2)
	d <- e1@Dim
	n <- d[1L]
	stopifnot(length(r) == n)
	if(isNum && !is.double(r))
            r <- as.double(r)
	## faster (?) than  m <- matrix(r00,n,n); diag(m) <- r ; as.vector(m)
        xx <- rbind(r, matrix(r00,n,n), deparse.level=0L)[seq_len(n*n)]
	newcl <-
	    paste0(if(isNum) "d" else if(isLog) {
		if(!anyNA(r) && !anyNA(r00)) "n" else "l"
	    } else stop("not yet implemented .. please report"), "syMatrix")

	new(newcl, Dim = e1@Dim, Dimnames = e1@Dimnames, x = xx)
    }
}

### This would be *the* way, but we get tons of "ambiguous method dispatch"
## we use this hack instead of signature  x = "diagonalMatrix" :
diCls <- names(getClass("diagonalMatrix")@subclasses)
if(FALSE) {
setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "diagonalMatrix"),
          diagOdiag)
} else { ## These are just for method disambiguation:
    for(c1 in diCls)
	for(c2 in diCls)
	    setMethod("Ops", signature(e1 = c1, e2 = c2), diagOdiag)
}

## diagonal  o  triangular  |-->  triangular
## diagonal  o  symmetric   |-->  symmetric
##    {also when other is sparse: do these "here" --
##     before conversion to sparse, since that loses "diagonality"}
diagOtri <- function(e1,e2) {
    ## result must be triangular
    r <- callGeneric(d1 <- .diag.x(e1), diag(e2)) # error if not "compatible"
    ## Check what happens with non-diagonals, i.e. (0 o 0), (FALSE o 0), ...:
    e1.0 <- if(is.numeric(d1)) 0 else FALSE
    r00 <- callGeneric(e1.0, if(.n2 <- is.numeric(e2[0L])) 0 else FALSE)
    if(is0(r00)) { ##  r00 == 0 or FALSE --- result *is* triangular
        diag(e2) <- r
        ## check what happens "in the triangle"
        e2.2 <- if(.n2) 2 else TRUE
        if(!callGeneric(e1.0, e2.2) == e2.2) { # values "in triangle" can change:
            n <- dim(e2)[1L]
            it <- indTri(n, upper = (e2@uplo == "U"))
            e2[it] <- callGeneric(e1.0, e2[it])
        }
        e2
    }
    else { ## result not triangular ---> general
        rr <- as(e2, "generalMatrix")
        diag(rr) <- r
        rr
    }
}


setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "triangularMatrix"),
          diagOtri)
## For the reverse,  Ops == "Arith" | "Compare" | "Logic"
##   'Arith'  :=  '"+"', '"-"', '"*"', '"^"', '"%%"', '"%/%"', '"/"'
setMethod("Arith", signature(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          function(e1,e2)
      { ## this must only trigger for *dense* e1
	  switch(.Generic,
		 "+" = .Call(dtrMatrix_addDiag,
                             unpack(..dense2d(e1)),   .diag.x(e2)),
		 "-" = .Call(dtrMatrix_addDiag,
                             unpack(..dense2d(e1)), - .diag.x(e2)),
		 "*" = {
		     n <- e2@Dim[1L]
		     d2 <- if(e2@diag == "U") { # unit-diagonal
			 d <- rep.int(as1(e2@x), n)
			 e2@x <- d
			 e2@diag <- "N"
			 d
		     } else e2@x
		     e2@x <- diag(e1) * d2
		     e2
		 },
		 "^" = { ## will be dense ( as  <ANY> ^ 0 == 1 ):
		     e1 ^ .diag2dense(e2, ".ge")
		 },
		 ## otherwise:
		 callGeneric(e1, .diag2T.smart(e2, e1)))
})

## Compare --> 'swap' (e.g.   e1 < e2   <==>  e2 > e1 ):
setMethod("Compare", signature(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
	  .Cmp.swap)
## '&' and "|'  are commutative:
setMethod("Logic", signature(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          function(e1, e2) callGeneric(e2, e1))

## For almost everything else, diag* shall be treated "as sparse" :
## These are cheap implementations via coercion

## For disambiguation --- define this for "sparseMatrix" , then for "ANY";
## and because we can save an .M.kind() call, we use this explicit
## "hack" for all diagonalMatrix *subclasses* instead of just "diagonalMatrix" :
##
## ddi*:
setMethod("Ops", signature(e1 = "ddiMatrix", e2 = "sparseMatrix"),
	  function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "d"), e2))
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "ddiMatrix"),
	  function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind = "d")))
## ldi*
setMethod("Ops", signature(e1 = "ldiMatrix", e2 = "sparseMatrix"),
	  function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "l"), e2))
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "ldiMatrix"),
	  function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind = "l")))

## Ops:	 Arith	--> numeric : "dMatrix"
##	 Compare --> logical
##	 Logic	 --> logical: "lMatrix"

## Other = "numeric" : stay diagonal if possible
## ddi*: Arith: result numeric, potentially ddiMatrix
for(arg2 in c("numeric","logical"))
setMethod("Arith", signature(e1 = "ddiMatrix", e2 = arg2),
	  function(e1,e2) {
	      n <- e1@Dim[1L]
	      if(length(e2) == 0L)
		  return(if(n) numeric() else e1)
	      f0 <- callGeneric(0, e2)
	      if(all0(f0)) { # remain diagonal
		  if(e1@diag == "U") {
		      if(any((r <- callGeneric(1, e2)) != 1)) {
			  e1@diag <- "N"
			  e1@x[seq_len(n)] <- r # possibly recycling r
		      } ## else: result = e1  (is "U" diag)
		  } else if(n) {
		      L1 <- (le <- length(e2)) == 1L
		      r <- callGeneric(e1@x, e2)
		      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
		      e1@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
		  }
		  e1
	      } else
		  callGeneric(.diag2tT.smart(e1, e2, kind = "d"), e2)
	  })

for(arg1 in c("numeric","logical"))
setMethod("Arith", signature(e1 = arg1, e2 = "ddiMatrix"),
	  function(e1,e2) {
	      n <- e2@Dim[1L]
	      if(length(e1) == 0L)
		  return(if(n) numeric() else e2)
	      f0 <- callGeneric(e1, 0)
	      if(all0(f0)) { # remain diagonal
		  if(e2@diag == "U") {
		      if(any((r <- callGeneric(e1, 1)) != 1)) {
			  e2@diag <- "N"
			  e2@x[seq_len(n)] <- r # possibly recycling r
		      } ## else: result = e2  (is "U" diag)
		  } else {
		      L1 <- (le <- length(e1)) == 1L
		      r <- callGeneric(e1, e2@x)
		      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
		      e2@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
		  }
		  e2
	      } else
		  callGeneric(e1, .diag2tT.smart(e2, e1, kind = "d"))
	  })

## ldi* Arith --> result numeric, potentially ddiMatrix
for(arg2 in c("numeric","logical"))
setMethod("Arith", signature(e1 = "ldiMatrix", e2 = arg2),
	  function(e1,e2) {
	      n <- e1@Dim[1L]
	      if(length(e2) == 0L)
		  return(if(n) numeric()
			 else copyClass(e1, "ddiMatrix", c("diag", "Dim", "Dimnames"), check=FALSE))
	      f0 <- callGeneric(0, e2)
	      if(all0(f0)) { # remain diagonal
		  E <- copyClass(e1, "ddiMatrix", c("diag", "Dim", "Dimnames"), check=FALSE)
		  ## storage.mode(E@x) <- "double"
		  if(e1@diag == "U") {
		      if(any((r <- callGeneric(1, e2)) != 1)) {
			  E@diag <- "N"
			  E@x[seq_len(n)] <- r # possibly recycling r
		      } ## else: result = E  (is "U" diag)
		  } else if(n) {
		      L1 <- (le <- length(e2)) == 1L
		      r <- callGeneric(e1@x, e2)
		      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
		      E@x[seq_len(n)] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
		  }
		  E
	      } else
		  callGeneric(.diag2tT.smart(e1, e2, kind = "l"), e2)
	  })

for(arg1 in c("numeric","logical"))
setMethod("Arith", signature(e1 = arg1, e2 = "ldiMatrix"),
	  function(e1,e2) {
	      n <- e2@Dim[1L]
	      if(length(e1) == 0L)
		  return(if(n) numeric()
			 else copyClass(e2, "ddiMatrix", c("diag", "Dim", "Dimnames"), check=FALSE))
	      f0 <- callGeneric(e1, 0)
	      if(all0(f0)) { # remain diagonal
		  E <- copyClass(e2, "ddiMatrix", c("diag", "Dim", "Dimnames"), check=FALSE)
		  ## storage.mode(E@x) <- "double"
		  if(e2@diag == "U") {
		      if(any((r <- callGeneric(e1, 1)) != 1)) {
			  E@diag <- "N"
			  E@x[seq_len(n)] <- r # possibly recycling r
		      } ## else: result = E  (is "U" diag)
		  } else if(n) {
		      L1 <- (le <- length(e1)) == 1L
		      r <- callGeneric(e1, e2@x)
		      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
		      E@x[seq_len(n)] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
		  }
		  E
	      } else
		  callGeneric(e1, .diag2tT.smart(e2, e1, kind = "l"))
	  })

## ddi*: for "Ops" without "Arith": <Compare> or <Logic> --> result logical, potentially ldi
##
## Note that  ("numeric", "ddiMatrix")  is simply swapped, e.g.,
if(FALSE) {
    selectMethod("<", c("numeric","lMatrix"))# Compare
    selectMethod("&", c("numeric","lMatrix"))# Logic
} ## so we don't need to define a method here :

for(arg2 in c("numeric","logical"))
setMethod("Ops", signature(e1 = "ddiMatrix", e2 = arg2),
	  function(e1,e2) {
	      n <- e1@Dim[1L]
	      if(length(e2) == 0L)
		  return(if(n) logical()
			 else copyClass(e1, "ldiMatrix", c("diag", "Dim", "Dimnames"), check=FALSE))
	      f0 <- callGeneric(0, e2)
	      if(all0(f0)) { # remain diagonal
		  E <- copyClass(e1, "ldiMatrix", c("diag", "Dim", "Dimnames"), check=FALSE)
		  ## storage.mode(E@x) <- "logical"
		  if(e1@diag == "U") {
		      if(any((r <- callGeneric(1, e2)) != 1)) {
			  E@diag <- "N"
			  E@x[seq_len(n)] <- r # possibly recycling r
		      } ## else: result = E  (is "U" diag)
		  } else if(n) {
		      L1 <- (le <- length(e2)) == 1L
		      r <- callGeneric(e1@x, e2)
		      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
		      E@x[seq_len(n)] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
		  }
		  E
	      } else
		  callGeneric(.diag2tT.smart(e1, e2, kind = "d"), e2)
	  })

## ldi*: for "Ops" without "Arith": <Compare> or <Logic> --> result logical, potentially ldi
for(arg2 in c("numeric","logical"))
setMethod("Ops", signature(e1 = "ldiMatrix", e2 = arg2),
	  function(e1,e2) {
	      n <- e1@Dim[1L]
	      if(length(e2) == 0L)
                  return(if(n) logical() else e1)
	      f0 <- callGeneric(FALSE, e2)
	      if(all0(f0)) { # remain diagonal
		  if(e1@diag == "U") {
		      if(any((r <- callGeneric(TRUE, e2)) != 1)) {
			  e1@diag <- "N"
			  e1@x[seq_len(n)] <- r # possibly recycling r
		      } ## else: result = e1  (is "U" diag)
		  } else if(n) {
		      L1 <- (le <- length(e2)) == 1L
		      r <- callGeneric(e1@x, e2)
		      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
		      e1@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
		  }
		  e1
	      } else
		  callGeneric(.diag2tT.smart(e1, e2, kind = "l"), e2)
	  })


## Not {"sparseMatrix", "numeric} :  {"denseMatrix", "matrix", ... }
for(other in c("ANY", "Matrix", "dMatrix")) {
    ## ddi*:
    setMethod("Ops", signature(e1 = "ddiMatrix", e2 = other),
	      function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind="d"), e2))
    setMethod("Ops", signature(e1 = other, e2 = "ddiMatrix"),
	      function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind="d")))
    ## ldi*:
    setMethod("Ops", signature(e1 = "ldiMatrix", e2 = other),
	      function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind="l"), e2))
    setMethod("Ops", signature(e1 = other, e2 = "ldiMatrix"),
	      function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind="l")))
}

## Direct subclasses of "denseMatrix": currently ddenseMatrix, ldense... :
if(FALSE) # now also contains "geMatrix"
dense.subCl <- local({ dM.scl <- getClass("denseMatrix")@subclasses
		       names(dM.scl)[vapply(dM.scl, slot, 0, "distance") == 1] })
dense.subCl <- paste0(c("d","l","n"), "denseMatrix")
for(DI in diCls) {
    dMeth <- if(extends(DI, "dMatrix"))
	function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "d"), e2)
    else # "lMatrix", the only other kind for now
	function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "l"), e2)
    for(c2 in c(dense.subCl, "Matrix")) {
	for(Fun in c("*", "&")) {
	    setMethod(Fun, signature(e1 = DI, e2 = c2),
		      function(e1,e2) callGeneric(e1, Diagonal(x = diag(e2))))
	    setMethod(Fun, signature(e1 = c2, e2 = DI),
		      function(e1,e2) callGeneric(Diagonal(x = diag(e1)), e2))
	}
	setMethod("^", signature(e1 = c2, e2 = DI),
		  function(e1,e2) callGeneric(Diagonal(x = diag(e1)), e2))
	for(Fun in c("%%", "%/%", "/")) ## 0 <op> 0 |--> NaN  for these.
	    setMethod(Fun, signature(e1 = DI, e2 = c2), dMeth)
    }
}

## Group methods "Math", "Math2" in			--> ./Math.R

### "Summary" : "max"   "min"   "range" "prod"  "sum"   "any"   "all"
### ----------   the last 4: separately here
for(cl in diCls) {
setMethod("any", cl,
	  function (x, ..., na.rm) {
	      if(any(x@Dim == 0)) FALSE
	      else if(x@diag == "U") TRUE else any(x@x, ..., na.rm = na.rm)
	  })
setMethod("all",  cl, function (x, ..., na.rm) {
    n <- x@Dim[1L]
    if(n >= 2) FALSE
    else if(n == 0 || x@diag == "U") TRUE
    else all(x@x, ..., na.rm = na.rm)
})
setMethod("prod", cl, function (x, ..., na.rm) {
    n <- x@Dim[1L]
    if(n >= 2) 0
    else if(n == 0 || x@diag == "U") 1
    else ## n == 1, diag = "N" :
	prod(x@x, ..., na.rm = na.rm)
})

setMethod("sum", cl,
	  function(x, ..., na.rm) {
	      r <- sum(x@x, ..., na.rm = na.rm)# double or integer, correctly
	      if(x@diag == "U" && !is.na(r)) r + x@Dim[1L] else r
	  })
}

## The remaining ones are  max, min, range :

setMethod("Summary", "ddiMatrix",
	  function(x, ..., na.rm) {
	      if(any(x@Dim == 0)) callGeneric(numeric(0), ..., na.rm=na.rm)
	      else if(x@diag == "U")
		  callGeneric(x@x, 0, 1, ..., na.rm=na.rm)
	      else callGeneric(x@x, 0, ..., na.rm=na.rm)
	  })
setMethod("Summary", "ldiMatrix",
	  function(x, ..., na.rm) {
	      if(any(x@Dim == 0)) callGeneric(logical(0), ..., na.rm=na.rm)
	      else if(x@diag == "U")
		  callGeneric(x@x, FALSE, TRUE, ..., na.rm=na.rm)
	      else callGeneric(x@x, FALSE, ..., na.rm=na.rm)
	  })



## similar to prTriang() in ./Auxiliaries.R :
prDiag <-
    function(x, digits = getOption("digits"), justify = "none", right = TRUE)
{
    cf <- array(".", dim = x@Dim, dimnames = x@Dimnames)
    cf[row(cf) == col(cf)] <-
	vapply(diag(x), format, "", digits = digits, justify = justify)
    print(cf, quote = FALSE, right = right)
    invisible(x)
}

## somewhat consistent with "print" for sparseMatrix :
setMethod("print", signature(x = "diagonalMatrix"), prDiag)

setMethod("show", signature(object = "diagonalMatrix"),
	  function(object) {
	      d <- dim(object)
	      cl <- class(object)
	      cat(sprintf('%d x %d diagonal matrix of class "%s"',
			  d[1L], d[2L], cl))
	      if(d[1L] < 50) {
		  cat("\n")
		  prDiag(object)
	      } else {
		  cat(", with diagonal entries\n")
		  show(diag(object))
		  invisible(object)
	      }
	  })

rm(arg1, arg2, other, DI, Fun, cl, c1, c2,
   dense.subCl, diCls)# not used elsewhere

setMethod("summary", signature(object = "diagonalMatrix"),
	  function(object, ...) {
	      d <- dim(object)
	      r <- summary(object@x, ...)
	      attr(r, "header") <-
		  sprintf('%d x %d diagonal Matrix of class "%s"',
			  d[1L], d[2L], class(object))
	      ## use ole' S3 technology for such a simple case
	      class(r) <- c("diagSummary", class(r))
	      r
	  })

print.diagSummary <- function (x, ...) {
    cat(attr(x, "header"),"\n")
    class(x) <- class(x)[-1]
    print(x, ...)
    invisible(x)
}
