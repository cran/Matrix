#### Sparse Matrices in Compressed row-oriented format
####                               --- "R"

### ``mainly for completeness'' --- we *do* favour Csparse
##    - - - - - - - - - - - -   hence only "minimal" methods here !
##  see also ./SparseM-conv.R

### contains = "dMatrix"

## MJ: all in C now, and moved to ./Auxiliaries.R
if(FALSE) {
.R.2.C <- function(from)
{
    cl <- class(from)
    valid <- c("dgRMatrix", "dsRMatrix", "dtRMatrix",
               "lgRMatrix", "lsRMatrix", "ltRMatrix",
               "ngRMatrix", "nsRMatrix", "ntRMatrix",
               "zgRMatrix", "zsRMatrix", "ztRMatrix")
    icl <- match(cl, valid) - 1L
    if(is.na(icl)) stop(gettextf("invalid class: %s", dQuote(cl)), domain=NA)
    Ccl <- sub("^(..)R","\\1C", cl)  # corresponding Csparse class name
    r <- new(Ccl)
    r@Dim <- from@Dim[2:1]
    if(icl %/% 3 != 2) ## not "n..Matrix" --> has 'x' slot
        r@x <- from@x
    if(icl %% 3 != 0) {                 # symmetric or triangular
        r@uplo <- from@uplo
        if(icl %% 3 == 2)               # triangular
            r@diag <- from@diag
    }
    r@i <- from@j
    r@p <- from@p
    r <- t(r)
    r@Dimnames <- from@Dimnames
    r
}

## However, a quick way to "treat a t(<R..>) as corresponding <C..> " :
.tR.2.C <- function(from) {
    cl <- class(from)
    valid <- c("dgRMatrix", "dsRMatrix", "dtRMatrix",
               "lgRMatrix", "lsRMatrix", "ltRMatrix",
               "ngRMatrix", "nsRMatrix", "ntRMatrix",
               "zgRMatrix", "zsRMatrix", "ztRMatrix")
    icl <- match(cl, valid) - 1L
    if(is.na(icl)) stop(gettextf("invalid class: %s", dQuote(cl)), domain=NA)
    Ccl <- sub("^(..)R","\\1C", cl)  # corresponding Csparse class name
    r <- new(Ccl)
    r@i <- from@j
    ##-         -
    r@p <- from@p
    r@Dim      <- from@Dim[2:1]
    r@Dimnames <- from@Dimnames[2:1]

    if(icl %/% 3 != 2) ## not "n..Matrix" --> has 'x' slot
        r@x <- from@x
    if(icl %% 3 != 0) {                 # symmetric or triangular
        r@uplo <- if(from@uplo == "U") "L" else "U"
        if(icl %% 3 == 2)               # triangular
            r@diag <- from@diag
    }
    r
}

.tC.2.R <- function(m, cl = class(m), clx = getClassDef(cl)) {
    has.x <- !extends(clx, "nsparseMatrix")## <==> has 'x' slot
    sh <- .M.shapeC(m,clx)
    r <- new(paste0(.M.kindC(clx), sh, "RMatrix"))
    r@Dim      <- m@Dim[2:1]
    r@Dimnames <- m@Dimnames[2:1]
    r@p <- m@p
    r@j <- m@i
    if(has.x)
	r@x <- m@x
    if(sh != "g") {
	r@uplo <- if(m@uplo != "U") "U" else "L"
	if(sh == "t")
	    r@diag <- m@diag
    }
    r
}
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
## coercion to other virtual classes --- the functionality we want to encourage

setAs("RsparseMatrix", "TsparseMatrix", .R.2.T)
setAs("RsparseMatrix", "CsparseMatrix", .R.2.C)

setAs("RsparseMatrix", "denseMatrix",
      function(from) as(.R.2.C(from), "denseMatrix"))

setAs("RsparseMatrix", "dsparseMatrix",
      function(from) as(.R.2.C(from), "dsparseMatrix"))
setAs("RsparseMatrix", "lsparseMatrix",
      function(from) as(.R.2.C(from), "lsparseMatrix"))
setAs("RsparseMatrix", "nsparseMatrix",
      function(from) as(.R.2.C(from), "nsparseMatrix"))

setAs("RsparseMatrix", "dMatrix",
      function(from) as(.R.2.C(from), "dMatrix"))
setAs("RsparseMatrix", "lMatrix",
      function(from) as(.R.2.C(from), "lMatrix"))
setAs("RsparseMatrix", "nMatrix",
      function(from) as(.R.2.C(from), "nMatrix"))

setAs("RsparseMatrix", "generalMatrix",
      function(from) as(.R.2.C(from), "generalMatrix"))

## for printing etc:
setAs("RsparseMatrix", "dgeMatrix",
      function(from) as(.R.2.C(from), "dgeMatrix"))
setAs("RsparseMatrix", "matrix",
      function(from) as(.R.2.C(from), "matrix"))

setAs("sparseMatrix", "RsparseMatrix", .viaC.2.R)
setAs("CsparseMatrix", "RsparseMatrix", .C.2.R)

##setAs("dgRMatrix", "dgeMatrix",
##      function(from) .Call(csc_to_dgeMatrix, from))
} ## MJ

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
## **VERY** cheap substitute:  work via dgC and t(.)
.viaC.2.dgR <- function(from) {
    m <- as(t(from), "dgCMatrix")
    new("dgRMatrix", Dim = dim(from), Dimnames = .M.DN(from),
	p = m@p, j = m@i, x = m@x)
}

## one of the few coercions "to <specific>" {tested in ../tests/Class+Meth.R}
setAs("matrix", "dgRMatrix", .viaC.2.dgR)
## setAs("dtCMatrix", "dtRMatrix", .viaC.to.dgR) # should work; can NOT use 'p'

setAs("matrix",      "RsparseMatrix", .viaC.2.R)
setAs("denseMatrix", "RsparseMatrix", .viaC.2.R)

##setAs("matrix", "dgRMatrix",
##      function(from) {
##          storage.mode(from) <- "double"
##          .Call(matrix_to_csc, from)
##      })
} ## MJ

## MJ: "fixed" in ./sparseMatrix.R
if(FALSE) {
## symmetric: can use same 'p' slot
setAs("dsCMatrix", "dsRMatrix",
      function(from) new("dsRMatrix", Dim = dim(from), Dimnames = .M.DN(from),
	      p = from@p, j = from@i, x = from@x,
	      uplo = if (from@uplo == "U") "L" else "U"))
## FIXME: if this makes sense, do it for "l" and "n" as well as "d"
} ## MJ

## MJ: no longer needed ... now inherited from Matrix
if(FALSE) {
setMethod("image", "dgRMatrix", function(x, ...) image(.CR2T(x), ...))
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
setMethod("t", "RsparseMatrix", function(x) .C.2.R(.tR.2.C(x)))

## Want tril(), triu(), band() --- just as "indexing" ---
## return a "close" class:
setMethod("tril", "RsparseMatrix",
	  function(x, k = 0, ...)
	  as(tril(.R.2.C(x), k = k, ...), "RsparseMatrix"))
setMethod("triu", "RsparseMatrix",
	  function(x, k = 0, ...)
	  as(triu(.R.2.C(x), k = k, ...), "RsparseMatrix"))
setMethod("band", "RsparseMatrix",
	  function(x, k1, k2, ...)
	  as(band(.R.2.C(x), k1 = k1, k2 = k2, ...), "RsparseMatrix"))

##setMethod("diag", signature(x = "dgRMatrix"),
##          function(x = 1, nrow, ncol = n) .Call(csc_getDiag, x))
} ## MJ

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "index", j = "missing",
				value = "replValue"),
		 function (x, i, j, ..., value)
		 replTmat(.CR2T(x), i=i, , value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "missing", j = "index",
				value = "replValue"),
		 function (x, i, j, ..., value)# extra " , ": want nargs() == 4
		 replTmat(.CR2T(x), , j=j, value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "index", j = "index",
				value = "replValue"),
		 function (x, i, j, ..., value)
		 replTmat(.CR2T(x), i=i, j=j, value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "index", j = "missing",
				value = "sparseVector"),
		 function (x, i, j, ..., value)
		 replTmat(.CR2T(x), i=i, value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "missing", j = "index",
				value = "sparseVector"),
		 function (x, i, j, ..., value)# extra " , ": want nargs() == 4
		 replTmat(.CR2T(x), , j=j, value=value))

setReplaceMethod("[", signature(x = "RsparseMatrix", i = "index", j = "index",
				value = "sparseVector"),
		 function (x, i, j, ..., value)
		 replTmat(.CR2T(x), i=i, j=j, value=value))


setReplaceMethod("[", signature(x = "RsparseMatrix", i = "matrix", j = "missing",
				value = "replValue"),
		 function (x, i, j, ..., value)
		 .TM.repl.i.mat(.CR2T(x), i=i, value=value))
