## METHODS FOR CLASS: diagonalMatrix (virtual)
## diagonal matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("Matrix", "diagonalMatrix", .M2diag)
setAs("matrix", "diagonalMatrix", .M2diag)


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
..diag2gC <- function(from)
    .Call(R_diagonal_as_sparse, from, ".gC", "U", TRUE)

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
setAs("diagonalMatrix",           "matrix",  .diag2m)
setAs("diagonalMatrix",           "vector",  .diag2v)

setMethod("as.vector", signature(x = "diagonalMatrix"),
          function(x, mode = "any") as.vector(.diag2v(x), mode))

setMethod("as.numeric", signature(x = "diagonalMatrix"),
          function(x, ...) as.double(.diag2v(x)))
setMethod("as.numeric", signature(x = "ddiMatrix"),
          function(x, ...) .diag2v(x))

setMethod("as.logical", signature(x = "diagonalMatrix"),
          function(x, ...) as.logical(.diag2v(x)))
setMethod("as.logical", signature(x = "ldiMatrix"),
          function(x, ...) .diag2v(x))

rm(..diag2dsparse, ..diag2lsparse, ..diag2nsparse,
   ..diag2tC, ..diag2tR, ..diag2tT, ..diag2sC, ..diag2gC,
   ..diag2tp, ..diag2tr, ..diag2dtr, ..diag2ltr, ..diag2ntr)


## ~~~~ CONSTRUCTORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## diagonalMatrix constructor, allowing either or both of 'n' and 'x' to be
## missing ... like base::diag() but _not_ also extracting diagonal entries
Diagonal <- function(n, x = NULL, names = FALSE) {
    nx <- length(x)
    if(missing(n))
        n <- nx
    else if(!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0L)
        stop("'n' must be a non-negative integer")
    if(is.double(n) && n >= .Machine$integer.max + 1)
        stop("dimensions cannot exceed 2^31-1")
    n <- as.integer(n) # discarding attributes
    if(is.null(x)) {
        r <- new("ddiMatrix")
        r@diag <- "U"
        if(n > 0L) {
            r@Dim <- c(n, n)
            if(is.character(names) && length(names) == n)
                r@Dimnames <- list(names, names)
        }
        return(r)
    }
    if(is.object(x))
        stop(gettextf("'x' has unsupported class \"%s\"", class(x)[1L]),
             domain = NA)
    names.x <- names(x) # keeping for later
    r <- new(switch(typeof(x),
                    ## discarding attributes, incl. 'dim' and 'names'
                    logical = { x <- as.logical(x); "ldiMatrix" },
                    integer =,
                    double = { x <- as.double(x); "ddiMatrix" },
                    stop(gettextf("'x' has unsupported type \"%s\"", typeof(x)),
                         domain = NA)))
    if(n == 0L)
        return(r)
    if(nx != 1L)
        r@x <-
            if(nx == n)
                x
            else if(nx > 0L)
                rep_len(x, n)
            else stop("attempt to recycle 'x' of length 0 to length 'n' (n > 0)")
    else if(is.na(x) || x != 1)
        r@x <- rep.int(x, n)
    else r@diag <- "U"
    r@Dim <- c(n, n)
    if(is.character(names)) {
        if(length(names) == n)
            r@Dimnames <- list(names, names)
    } else if(isTRUE(names) && !is.null(names.x)) {
        names.x <- rep_len(names.x, n) # we know length(names.x) > 0L
        r@Dimnames <- list(names.x, names.x)
    }
    r
}

.sparseDiagonal <- function(n, x = NULL, uplo = "U", shape = "t", unitri = TRUE,
                            kind, cols) {
    if(missing(n))
        n <- length(x)
    else if(!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0L)
        stop("'n' must be a non-negative integer")
    if(is.double(n) && n >= .Machine$integer.max + 1)
        stop("dimensions cannot exceed 2^31-1")
    n <- nj <- as.integer(n) # stripping attributes

    if(!(missing(shape) ||
         (is.character(shape) && length(shape) == 1L && !is.na(shape) &&
          any(shape == c("g", "t", "s")))))
        stop("'shape' must be one of \"g\", \"t\", \"s\"")

    if(!((m.kind <- missing(kind)) ||
         (is.character(kind) && length(kind) == 1L && !is.na(kind) &&
          any(kind == c("d", "l", "n")))))
        stop("'kind' must be one of \"d\", \"l\", \"n\"")

    if(m.kind || kind != "n") {
        if(is.null(x))
           x <- if(m.kind) { kind <- "d"; 1 } else switch(kind, d = 1, l = TRUE)
        else if(is.object(x))
            stop(gettextf("'x' has unsupported class \"%s\"",
                          class(x)[1L]),
                 domain = NA)
        else {
            kind. <- switch(typeof(x),
                            ## discarding attributes, incl. 'dim' in array case
                            logical = { x <- as.logical(x); "l" },
                            integer =,
                            double = { x <- as.double(x); "d" },
                            stop(gettextf("'x' has unsupported type \"%s\"",
                                          typeof(x)),
                                 domain = NA))
            if(m.kind)
                kind <- kind.
            else if(kind != kind.) {
                warning(gettextf("mismatch between typeof(x)=\"%s\" and kind=\"%s\"; using kind=\"%s\"",
                                 typeof(x), kind, kind.),
                        domain = NA)
                kind <- kind.
            }
        }
    }

    if(!(m.cols <- missing(cols))) {
        if(!is.numeric(cols))
            stop("'cols' must be numeric")
        else if((nj <- length(cols)) > 0L &&
                (n == 0L || anyNA(rj <- range(cols)) ||
                 rj[1L] < 0L || rj[2L] >= n))
            stop("'cols' has elements not in seq(0, length.out = n)")
        else {
            cols <- as.integer(cols)
            shape <- "g"
        }
    }

    r <- new(paste0(kind, shape, "CMatrix"))
    r@Dim <- c(n, nj)
    if(shape != "g") {
        if(!missing(uplo)) {
            if(is.character(uplo) && length(uplo) == 1L && !is.na(uplo) &&
               any(uplo == c("U", "L")))
                r@uplo <- uplo
            else stop("'uplo' must be \"U\" or \"L\"")
        }
        if(shape == "t" && unitri &&
           (kind == "n" || (!anyNA(x) && all(if(kind == "l") x else x == 1)))) {
            r@diag <- "U"
            r@p <- integer(nj + 1)
            return(r)
        }
    }
    if(nj > 0L) {
        r@p <- 0:nj
        r@i <- if(m.cols) 0:(nj - 1L) else cols
        if(kind != "n") {
            x <-
                if((nx <- length(x)) == n)
                    x
                else if(nx > 0L)
                    rep_len(x, n)
                else stop("attempt to recycle 'x' of length 0 to length 'n' (n > 0)")
            r@x <- if(m.cols) x else x[1L + cols]
        }
    }
    r
}

.trDiagonal <- function(n, x = NULL, uplo = "U", unitri = TRUE, kind)
    .sparseDiagonal(n, x, uplo, shape = "t", unitri = unitri, kind = kind)

.symDiagonal <- function(n, x = NULL, uplo = "U", kind)
    .sparseDiagonal(n, x, uplo, shape = "s", kind = kind)

.bdiag <- function(lst) {
    if(!is.list(lst))
        stop("'lst' must be a list")
    if((n <- length(lst)) == 0L)
        return(new("dgTMatrix"))
    if(n == 1L)
        return(.CR2T(asCspN(lst[[1L]])))

### FIXME? this is _slow_ when 'lst' is list of 75000 3-by-3 dense matrices
    lst <- unname(lapply(lst, function(x) .CR2T(asCspN(x))))

    ## NB: class(.CR2T(.)) is always "[dln][gts]TMatrix"
    cl <- vapply(lst, class, "")
    kind  <- substr(cl, 1L, 1L) # "d", "l", or "n"
    shape <- substr(cl, 2L, 2L) # "g", "t", or "s"

    if(!(any(kind == (kind. <- "d")) || any(kind == (kind. <- "l"))))
        kind. <- "n"
    else if(any(z <- kind == "n"))
        lst[z] <- lapply(lst[z], .sparse2kind, kind.)

    shape. <-
        if(all(symmetric <- shape == "s"))
            "s"
        else if(all(shape == "t"))
            "t"
        else "g"

    if(shape. != "g") {
        uplo <- vapply(lst, slot, "", "uplo") # "U" or "L"
        if(shape. == "s")
            uplo. <-
                if(all(z <- uplo == "U"))
                    "U"
                else if(!any(z))
                    "L"
                else {
                    uplo.. <- if(2 * sum(z) >= n) { z <- !z; "U" } else "L"
                    lst[z] <- lapply(lst[z],
                                     function(x) .Call(R_sparse_transpose, x))
                    uplo..
                }
        else if(any(uplo != (uplo. <- uplo[1L])))
            shape. <- "g"
    }

    i_off <- c(0L, cumsum(vapply(lst, function(x) x@Dim[1L], 0L)))
    j_off <- c(0L, cumsum(vapply(lst, function(x) x@Dim[2L], 0L)))

    r <- new(paste0(kind., shape., "TMatrix"))
    r@Dim <- r@Dim <- c(i_off[n + 1L], j_off[n + 1L])
    if(shape. == "g")
        lst[symmetric] <- lapply(lst[symmetric], .sparse2g)
    else r@uplo <- uplo.
    r@i <- unlist(lapply(seq_len(n), function(k) i_off[k] + lst[[k]]@i),
                  FALSE, FALSE)
    r@j <- unlist(lapply(seq_len(n), function(k) j_off[k] + lst[[k]]@j),
                  FALSE, FALSE)
    if(kind. != "n")
        r@x <- unlist(lapply(lst, slot, "x"), FALSE, FALSE)
    r
}

bdiag <- function(...) {
    if((n <- ...length()) == 0L)
        new("dgCMatrix")
    else if(n > 1L)
        .T2C(.bdiag(list(...)))
    else if(!is.list(x <- ..1))
        as(x, "CsparseMatrix")
    else if(length(x) == 1L)
        as(x[[1L]], "CsparseMatrix")
    else .T2C(.bdiag(x))
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 .diag.x <- function(m) if(m@diag != "N") rep.int(as1(m@x), m@Dim[1L]) else m@x
..diag.x <- function(m)                   rep.int(as1(m@x), m@Dim[1L])

setMethod("diag", signature(x = "diagonalMatrix"),
          function(x, nrow, ncol, names = TRUE) {
              r <- .diag.x(x)
              if(names &&
                 !any(vapply(dn <- x@Dimnames, is.null, NA)) &&
                 {
                     i <- seq_len(min(x@Dim))
                     identical(nms <- dn[[1L]][i], dn[[2L]][i])
                 })
                  names(r) <- nms
              r
          })

setMethod("diag<-", signature(x = "diagonalMatrix"),
          function(x, value) {
              n <- x@Dim[1L]
              nv <- length(value)
              if(nv != 1L && nv != n)
                  stop("replacement diagonal has wrong length")
              x@x <-
                  if(is.logical(x@x))
                      switch(typeof(value),
                             logical = rep_len(value, n),
                             integer =,
                             double =
                                 {
                                     x <- ..diag2d(x)
                                     rep_len(as.double(x), n)
                                 },
                             stop(gettextf("replacement diagonal has incompatible type \"%s\"", typeof(value)),
                                  domain = NA))
                  else
                      switch(typeof(value),
                             logical =,
                             integer =,
                             double = rep_len(as.double(value), n),
                             stop(gettextf("replacement diagonal has incompatible type \"%s\"", typeof(value)),
                                  domain = NA))
              x@diag <- "N"
              x
          })

setMethod("t", signature(x = "diagonalMatrix"),
          function(x) { x@Dimnames <- x@Dimnames[2:1]; x })

setMethod("band", signature(x = "diagonalMatrix"),
          function(x, k1, k2, ...)
              if(k1 <= 0L && k2 >= 0L) x else .setZero(x))

setMethod("triu", signature(x = "diagonalMatrix"),
          function(x, k = 0L, ...)
              if(k <= 0L) x else .setZero(x))

setMethod("tril", signature(x = "diagonalMatrix"),
          function(x, k = 0L, ...)
              if(k >= 0L) x else .setZero(x))

setMethod("forceSymmetric", signature(x = "diagonalMatrix", uplo = "character"),
          function(x, uplo) .diag2sparse(x, ".sC", uplo = uplo))

setMethod("forceSymmetric", signature(x = "diagonalMatrix", uplo = "missing"),
          function(x, uplo) .diag2sparse(x, ".sC", uplo = "U"))

setMethod("symmpart", signature(x = "diagonalMatrix"),
          function(x) forceSymmetric(..diag2d(x)))

setMethod("skewpart", signature(x = "diagonalMatrix"),
          function(x) symmetrizeDimnames(.setZero(x, "d")))

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
        forceDiagonal(x)
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
                 function(x,i,j, ..., value) {
                     if(all0(value)) { # be faster
                         r <- new(paste0(.M.kind(x), "tTMatrix")) # of all "0"
                         r@Dim <- x@Dim
                         r@Dimnames <- x@Dimnames
                         r
                     } else {
                         ## typically non-sense: assigning to full sparseMatrix
                         x[TRUE] <- value
                         x
                     }
                 })


setReplaceMethod("[", signature(x = "diagonalMatrix",
                                i = "matrix", # 2-col.matrix
                                j = "missing", value = "replValue"),
                 function(x,i,j, ..., value) {
                     if(ncol(i) == 2L) {
                         if(all((ii <- i[,1L]) == i[,2L])) {
                             ## replace in diagonal only
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
### FIXME:  use  uplo="U" or uplo="L"  (or *not* "triangularMatrix")
### depending on LE <- i <= j
### all(LE) //  all(!LE) // remaining cases
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
                     callGeneric(x=x, , j=j, value=as(value, "sparseVector")))

setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "missing",
                                value = "sparseMatrix"),
                 function (x, i, j, ..., value)
                     callGeneric(x=x, i=i, , value=as(value, "sparseVector")))
setReplaceMethod("[", signature(x = "diagonalMatrix", i = "index", j = "index",
                                value = "sparseMatrix"),
                 function (x, i, j, ..., value)
                     callGeneric(x=x, i=i, j=j, value=as(value, "sparseVector")))

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

## FIXME: Many of these products are not handling 'Dimnames' appropriately ...

.prod.diag.missing <- function(x, boolArith) {
    if(boolArith) {
        if(!is.logical(x@x))
            x <- ..diag2l(x)
    } else {
        if(!is.double(x@x))
            x <- ..diag2d(x)
        if(x@diag == "N")
            x@x <- x@x * x@x
    }
    x
}

setMethod( "crossprod", signature(x = "diagonalMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <- .prod.diag.missing(x, boolArith = isTRUE(boolArith))
              r@Dimnames <- r@Dimnames[c(2L, 2L)]
              r
          })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <- .prod.diag.missing(x, boolArith = isTRUE(boolArith))
              r@Dimnames <- r@Dimnames[c(1L, 1L)]
              r
          })

.prod.diag.diag <- function(x, y, boolArith) {
    if(boolArith) {
        if(x@diag == "N") {
            if(!is.logical(x@x))
                x <- ..diag2l(x)
            if(y@diag == "N")
                x@x <- x@x & y@x
            x
        } else if(is.logical(y@x))
            y
        else ..diag2l(y)
    } else {
        if(x@diag == "N") {
            if(!is.double(x@x))
                x <- ..diag2d(x)
            if(y@diag == "N")
                x@x <- x@x * y@x
            x
        } else if(is.double(y@x))
            y
        else ..diag2d(y)
    }
}

setMethod("%*%", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.diag(x, y, boolArith = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.diag(x, y, boolArith = TRUE)
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.diag.diag(x, y, boolArith = isTRUE(boolArith))
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.diag.diag(x, y, boolArith = isTRUE(boolArith))
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 3L)
              r
          })

.prod.diag.m <- function(x, y, boolArith, trans) {
    ## MJ: .m2ge() avoids a copy when argument is unreferenced,
    ##     so it is more efficient than Matrix() here
    if(boolArith) {
        kind <- "n"
        op <- `&`
    } else {
        kind <- "d"
        op <- `*`
    }
    .m2ge(if(x@diag == "N")
              op(x@x, if(trans) t(y) else y)
          else if(trans)
              t(y)
          else y,
          kind)
}

setMethod("%*%", signature(x = "diagonalMatrix", y = "matrix"),
          function(x, y) {
              mmultDim(x@Dim, dim(y), type = 1L)
              r <- .prod.diag.m(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "diagonalMatrix", y = "matrix"),
          function(x, y) {
              mmultDim(x@Dim, dim(y), type = 1L)
              r <- .prod.diag.m(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "diagonalMatrix", y = "matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, dim(y), type = 2L)
              r <- .prod.diag.m(x, y, boolArith = isTRUE(boolArith),
                                trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, dim(y), type = 3L)
              r <- .prod.diag.m(x, y, boolArith = isTRUE(boolArith),
                                trans = TRUE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 3L)
              r
          })

.prod.m.diag <- function(x, y, boolArith, trans) {
    ## MJ: .m2ge() avoids a copy when argument is unreferenced,
    ##     so it is more efficient than Matrix() here
    if(boolArith) {
        kind <- "n"
        op <- `&`
    } else {
        kind <- "d"
        op <- `*`
    }
    .m2ge(if(y@diag == "N")
              op(if(trans) t(x) else x,
                 rep(y@x, each = dim(x)[1L + trans]))
          else if(trans)
              t(x)
          else x,
          kind)
}

setMethod("%*%", signature(x = "matrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(dim(x), y@Dim, type = 1L)
              r <- .prod.m.diag(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "matrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(dim(x), y@Dim, type = 1L)
              r <- .prod.m.diag(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "matrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(dim(x), y@Dim, type = 2L)
              r <- .prod.m.diag(x, y, boolArith = isTRUE(boolArith),
                                trans = TRUE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "matrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(dim(x), y@Dim, type = 3L)
              r <- .prod.m.diag(x, y, boolArith = isTRUE(boolArith),
                                trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 3L)
              r
          })

## FIXME: <unit diagonalMatrix> %*% <symmetricMatrix> should not be symmetric,
##        because the inherited 'rownames' and 'colnames' could differ ...

.prod.diag.dense <- function(x, y, boolArith, trans) {
    if(boolArith) {
        y <- ..dense2n(y)
        op <- `&`
        one <- TRUE
    } else {
        if(!is.double(y@x))
            y <- ..dense2d(y)
        op <- `*`
        one <- 1
    }
    if(x@diag == "N") {
        y@x <-
            if(!.hasSlot(y, "uplo")) {
                ## y=[nd]geMatrix
                if(trans)
                    y <- t(y)
                y@factors <- list()
                op(x@x, y@x)
            } else if(.hasSlot(y, "diag")) {
                ## y=[nd]t[rp]Matrix
                if(trans)
                    y <- t(y)
                if(y@diag != "N")
                    diag(y) <- one
                if(length(y@x) == (n <- y@Dim[1L])^2)
                    op(x@x, y@x)
                else if(y@uplo == "U")
                    op(x@x[sequence.default(1:n, rep.int(1L, n))], y@x)
                else
                    op(x@x[sequence.default(n:1,            1:n)], y@x)
            } else {
                ## y=[nd]s[yp]Matrix
                y <- .dense2g(y)
                y@factors <- list()
                op(x@x, y@x)
            }
        y
    } else if(trans)
        t(y)
    else y
}

setMethod("%*%", signature(x = "diagonalMatrix", y = "denseMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.dense(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "diagonalMatrix", y = "denseMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.dense(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "diagonalMatrix", y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.diag.dense(x, y, boolArith = isTRUE(boolArith),
                                    trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.diag.dense(x, y, boolArith = isTRUE(boolArith),
                                    trans = TRUE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 3L)
              r
          })

.prod.dense.diag <- function(x, y, boolArith, trans) {
    if(boolArith) {
        x <- ..dense2n(x)
        op <- `&`
        one <- TRUE
    } else {
        if(!is.double(x@x))
            x <- ..dense2d(x)
        op <- `*`
        one <- 1
    }
    if(y@diag == "N") {
        x@x <-
            if(!.hasSlot(x, "uplo")) {
                ## x=[nd]geMatrix
                if(trans)
                    x <- t(x)
                x@factors <- list()
                op(x@x, rep(y@x, each = x@Dim[1L]))
            } else if(.hasSlot(x, "diag")) {
                ## x=[nd]t[rp]Matrix
                if(trans)
                    x <- t(x)
                if(x@diag != "N")
                    diag(x) <- one
                if(length(x@x) == (n <- x@Dim[1L])^2)
                    op(x@x, rep(y@x, each = x@Dim[1L]))
                else if(x@uplo == "U")
                    op(x@x, rep.int(y@x, 1:n))
                else
                    op(x@x, rep.int(y@x, n:1))
            } else {
                ## x=[nd]s[yp]Matrix
                x <- .dense2g(x)
                x@factors <- list()
                op(x@x, rep(y@x, each = x@Dim[1L]))
            }
        x
    } else if(trans)
        t(x)
    else x
}

setMethod("%*%", signature(x = "denseMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.dense.diag(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "denseMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.dense.diag(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "denseMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.dense.diag(x, y, boolArith = isTRUE(boolArith),
                                    trans = TRUE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "denseMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.dense.diag(x, y, boolArith = isTRUE(boolArith),
                                    trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 3L)
              r
          })

.prod.diag.Csparse <- function(x, y, boolArith, trans) {
    if(x@diag == "N") {
        y <- .sparse2kind(y, if(boolArith) "l" else "d", drop0 = FALSE)
        if(!.hasSlot(y, "uplo")) {
            ## y=[ld]gCMatrix
            if(trans)
                y <- t(y)
            y@factors <- list()
        } else if(.hasSlot(y, "diag")) {
            ## y=[ld]tCMatrix
            if(trans)
                y <- t(y)
            if(y@diag != "N")
                y <- ..diagU2N(y)
        } else {
            ## y=[ld]sCMatrix
            y <- .sparse2g(y)
            y@factors <- list()
        }
        op <- if(boolArith) `&` else `*`
        y@x <- op(x@x[y@i + 1L], y@x)
        if(boolArith) .sparse2kind(y, "n", drop0 = TRUE) else y
    } else
        (if(trans) t else identity)(
            if(boolArith)
                .sparse2kind(y, "n", drop0 = TRUE)
            else .sparse2kind(y, "d", drop0 = FALSE))
}

setMethod("%*%", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.Csparse(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.Csparse(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.diag.Csparse(x, y, boolArith = isTRUE(boolArith),
                                      trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.diag.Csparse(x, y, boolArith = isTRUE(boolArith),
                                      trans = TRUE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 3L)
              r
          })

.prod.Csparse.diag <- function(x, y, boolArith, trans) {
    if(y@diag == "N") {
        x <- .sparse2kind(x, if(boolArith) "l" else "d", drop0 = FALSE)
        if(!.hasSlot(x, "uplo")) {
            ## x=[ld]gCMatrix
            if(trans)
                x <- t(x)
            x@factors <- list()
        } else if(.hasSlot(x, "diag")) {
            ## x=[ld]tCMatrix
            if(trans)
                x <- t(x)
            if(x@diag != "N")
                x <- ..diagU2N(x)
        } else {
            ## x=[ld]sCMatrix
            x <- .sparse2g(x)
            x@factors <- list()
        }
        dp <- if((n <- length(p <- x@p)) > 1L) p[-1L] - p[-n] else integer(0L)
        x@x <- (if(boolArith) `&` else `*`)(x@x, rep.int(y@x, dp))
        if(boolArith) .sparse2kind(x, "n", drop0 = TRUE) else x
    } else
        (if(trans) t else identity)(
            if(boolArith)
                .sparse2kind(x, "n", drop0 = TRUE)
            else .sparse2kind(x, "d", drop0 = FALSE))
}

setMethod("%*%", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.Csparse.diag(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.Csparse.diag(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.Csparse.diag(x, y, boolArith = isTRUE(boolArith),
                                      trans = TRUE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.Csparse.diag(x, y, boolArith = isTRUE(boolArith),
                                      trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 3L)
              r
          })

.prod.diag.Rsparse <- function(x, y, boolArith, trans) {
    if(x@diag == "N") {
        y <- .sparse2kind(y, if(boolArith) "l" else "d", drop0 = FALSE)
        if(!.hasSlot(y, "uplo")) {
            ## y=[ld]gRMatrix
            if(trans)
                y <- t(y)
            y@factors <- list()
        } else if(.hasSlot(y, "diag")) {
            ## y=[ld]tRMatrix
            if(trans)
                y <- t(y)
            if(y@diag != "N")
                y <- ..diagU2N(y)
        } else {
            ## y=[ld]sRMatrix
            y <- .sparse2g(y)
            y@factors <- list()
        }
        dp <- if((n <- length(p <- x@p)) > 1L) p[-1L] - p[-n] else integer(0L)
        y@x <- (if(boolArith) `&` else `*`)(rep.int(x@x, dp), y@x)
        if(boolArith) .sparse2kind(y, "n", drop0 = TRUE) else y
    } else
        (if(trans) t else identity)(
            if(boolArith)
                .sparse2kind(y, "n", drop0 = TRUE)
            else .sparse2kind(y, "d", drop0 = FALSE))
}

setMethod("%*%", signature(x = "diagonalMatrix", y = "RsparseMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.Rsparse(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "diagonalMatrix", y = "RsparseMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.Rsparse(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "diagonalMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.diag.Rsparse(x, y, boolArith = isTRUE(boolArith),
                                      trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.diag.Rsparse(x, y, boolArith = isTRUE(boolArith),
                                      trans = TRUE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 3L)
              r
          })

.prod.Rsparse.diag <- function(x, y, boolArith, trans) {
    if(y@diag == "N") {
        x <- .sparse2kind(x, if(boolArith) "l" else "d", drop0 = FALSE)
        if(!.hasSlot(x, "uplo")) {
            ## x=[ld]gRMatrix
            if(trans)
                x <- t(x)
            x@factors <- list()
        } else if(.hasSlot(x, "diag")) {
            ## x=[ld]tRMatrix
            if(trans)
                x <- t(x)
            if(x@diag != "N")
                x <- ..diagU2N(x)
        } else {
            ## x=[ld]sRMatrix
            x <- .sparse2g(x)
            x@factors <- list()
        }
        op <- if(boolArith) `&` else `*`
        x@x <- op(x@x, y@x[x@j + 1L])
        if(boolArith) .sparse2kind(x, "n", drop0 = TRUE) else x
    } else
        (if(trans) t else identity)(
            if(boolArith)
                .sparse2kind(x, "n", drop0 = TRUE)
            else .sparse2kind(x, "d", drop0 = FALSE))
}

setMethod("%*%", signature(x = "RsparseMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.Rsparse.diag(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "RsparseMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.Rsparse.diag(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "RsparseMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.Rsparse.diag(x, y, boolArith = isTRUE(boolArith),
                                      trans = TRUE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.Rsparse.diag(x, y, boolArith = isTRUE(boolArith),
                                      trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 3L)
              r
          })

.prod.diag.Tsparse <- function(x, y, boolArith, trans) {
    if(x@diag == "N") {
        y <- .sparse2kind(y, if(boolArith) "l" else "d", drop0 = FALSE)
        if(!.hasSlot(y, "uplo")) {
            ## y=[ld]gTMatrix
            if(trans)
                y <- t(y)
            y@factors <- list()
        } else if(.hasSlot(y, "diag")) {
            ## y=[ld]tTMatrix
            if(trans)
                y <- t(y)
            if(y@diag != "N")
                y <- ..diagU2N(y)
        } else {
            ## y=[ld]sTMatrix
            y <- .sparse2g(y)
            y@factors <- list()
        }
        op <- if(boolArith) `&` else `*`
        y@x <- op(x@x[y@i + 1L], y@x)
        if(boolArith) .sparse2kind(y, "n", drop0 = TRUE) else y
    } else
        (if(trans) t else identity)(
            if(boolArith)
                .sparse2kind(y, "n", drop0 = TRUE)
            else .sparse2kind(y, "d", drop0 = FALSE))
}

setMethod("%*%", signature(x = "diagonalMatrix", y = "TsparseMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.Tsparse(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "diagonalMatrix", y = "TsparseMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.diag.Tsparse(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "diagonalMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.diag.Tsparse(x, y, boolArith = isTRUE(boolArith),
                                      trans = FALSE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.diag.Tsparse(x, y, boolArith = isTRUE(boolArith),
                                      trans = TRUE)
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 3L)
              r
          })

.prod.Tsparse.diag <- function(x, y, boolArith, trans) {
    if(y@diag == "N") {
        x <- .sparse2kind(x, if(boolArith) "l" else "d", drop0 = FALSE)
        if(!.hasSlot(x, "uplo")) {
            ## x=[ld]gTMatrix
            if(trans)
                x <- t(x)
            x@factors <- list()
        } else if(.hasSlot(x, "diag")) {
            ## x=[ld]tTMatrix
            if(trans)
                x <- t(x)
            if(x@diag != "N")
                x <- ..diagU2N(x)
        } else {
            ## x=[ld]sTMatrix
            x <- .sparse2g(x)
            x@factors <- list()
        }
        op <- if(boolArith) `&` else `*`
        x@x <- op(x@x, y@x[x@j + 1L])
        if(boolArith) .sparse2kind(x, "n", drop0 = TRUE) else x
    } else
        (if(trans) t else identity)(
            if(boolArith)
                .sparse2kind(x, "n", drop0 = TRUE)
            else .sparse2kind(x, "d", drop0 = FALSE))
}

setMethod("%*%", signature(x = "TsparseMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.Tsparse.diag(x, y, boolArith = FALSE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "TsparseMatrix", y = "diagonalMatrix"),
          function(x, y) {
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- .prod.Tsparse.diag(x, y, boolArith = TRUE, trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod( "crossprod", signature(x = "TsparseMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              r <- .prod.Tsparse.diag(x, y, boolArith = isTRUE(boolArith),
                                      trans = TRUE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              r <- .prod.Tsparse.diag(x, y, boolArith = isTRUE(boolArith),
                                      trans = FALSE)
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 3L)
              r
          })



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
                e2@diag <- "N"
                return(e2)
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
        e1@diag <- "N"
        e1
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
            paste0(if(isNum) "d"
                   else if(isLog) {
                       if(!anyNA(r) && !anyNA(r00)) "n" else "l"
                   } else stop("not yet implemented .. please report"), "syMatrix")

        new(newcl, Dim = e1@Dim, Dimnames = e1@Dimnames, x = xx)
    }
}

### This would be *the* way, but we get tons of "ambiguous method dispatch"
## we use this hack instead of signature  x = "diagonalMatrix" :
diCls <- names(getClassDef("diagonalMatrix")@subclasses)
if(FALSE) {
setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "diagonalMatrix"),
          diagOdiag)
} else { ## These are just for method disambiguation:
    for(c1 in diCls)
        for(c2 in diCls)
            setMethod("Ops", signature(e1 = c1, e2 = c2), diagOdiag)
    rm(c1, c2)
}
rm(diagOdiag)

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
rm(diagOtri)

## For the reverse,  Ops == "Arith" | "Compare" | "Logic"
##   'Arith'  :=  '"+"', '"-"', '"*"', '"^"', '"%%"', '"%/%"', '"/"'
setMethod("Arith", signature(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          function(e1, e2) { ## this must only trigger for *dense* e1
              switch(.Generic,
                     "+" = `diag<-`(e1, as.double(diag(e1, names=FALSE) + .diag.x(e2))),
                     "-" = `diag<-`(e1, as.double(diag(e1, names=FALSE) - .diag.x(e2))),
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
rm(arg2)

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
rm(arg1)

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
rm(arg2)

for(arg1 in c("numeric","logical"))
setMethod("Arith", signature(e1 = arg1, e2 = "ldiMatrix"),
          function(e1,e2) {
              n <- e2@Dim[1L]
              if(length(e1) == 0L)
                  return(if(n) numeric()
                         else copyClass(e2, "ddiMatrix",
                                        c("diag", "Dim", "Dimnames"),
                                        check=FALSE))
              f0 <- callGeneric(e1, 0)
              if(all0(f0)) { # remain diagonal
                  E <- copyClass(e2, "ddiMatrix",
                                 c("diag", "Dim", "Dimnames"),
                                 check=FALSE)
                  ## storage.mode(E@x) <- "double"
                  if(e2@diag == "U") {
                      if(any((r <- callGeneric(e1, 1)) != 1)) {
                          E@diag <- "N"
                          E@x[seq_len(n)] <- r # possibly recycling r
                      } ## else: result = E  (is "U" diag)
                  } else if(n) {
                      L1 <- (le <- length(e1)) == 1L
                      r <- callGeneric(e1, e2@x)
                      ## "future fixme": if we have idiMatrix,
                      ## and r is 'integer', use idiMatrix
                      E@x[seq_len(n)] <-
                          if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  E
              } else
                  callGeneric(e1, .diag2tT.smart(e2, e1, kind = "l"))
          })
rm(arg1)

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
                         else copyClass(e1, "ldiMatrix",
                                        c("diag", "Dim", "Dimnames"),
                                        check=FALSE))
              f0 <- callGeneric(0, e2)
              if(all0(f0)) { # remain diagonal
                  E <- copyClass(e1, "ldiMatrix",
                                 c("diag", "Dim", "Dimnames"),
                                 check=FALSE)
                  ## storage.mode(E@x) <- "logical"
                  if(e1@diag == "U") {
                      if(any((r <- callGeneric(1, e2)) != 1)) {
                          E@diag <- "N"
                          E@x[seq_len(n)] <- r # possibly recycling r
                      } ## else: result = E  (is "U" diag)
                  } else if(n) {
                      L1 <- (le <- length(e2)) == 1L
                      r <- callGeneric(e1@x, e2)
                      ## "future fixme": if we have idiMatrix,
                      ### and r is 'integer', use idiMatrix
                      E@x[seq_len(n)] <-
                          if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  E
              } else
                  callGeneric(.diag2tT.smart(e1, e2, kind = "d"), e2)
          })
rm(arg2)

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
                      ## "future fixme": if we have idiMatrix,
                      ## and r is 'integer', use idiMatrix
                      e1@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  e1
              } else
                  callGeneric(.diag2tT.smart(e1, e2, kind = "l"), e2)
          })
rm(arg2)

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
rm(other)

## Direct subclasses of "denseMatrix": currently ddenseMatrix, ldense... :
if(FALSE) # now also contains "geMatrix"
dense.subCl <- local({ dM.scl <- getClassDef("denseMatrix")@subclasses
    names(dM.scl)[vapply(dM.scl, slot, 0, "distance") == 1] })
dense.subCl <- paste0(c("d","l","n"), "denseMatrix")
for(DI in diCls) {
    dMeth <-
        if(extends(DI, "dMatrix"))
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
rm(dense.subCl, DI, dMeth, c2, Fun)

## Group methods "Math", "Math2" in			--> ./Math.R

### "Summary" : "max"   "min"   "range" "prod"  "sum"   "any"   "all"
### ----------   the last 4: separately here
for(cl in diCls) {
setMethod("any", cl,
          function (x, ..., na.rm) {
              if(any(x@Dim == 0)) FALSE
              else if(x@diag == "U") TRUE else any(x@x, ..., na.rm = na.rm)
          })
setMethod("all",  cl,
          function (x, ..., na.rm) {
              n <- x@Dim[1L]
              if(n >= 2) FALSE
              else if(n == 0 || x@diag == "U") TRUE
              else all(x@x, ..., na.rm = na.rm)
          })
setMethod("prod", cl,
          function (x, ..., na.rm) {
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
rm(cl, diCls)

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
function(x, digits = getOption("digits"), justify = "none", right = TRUE) {
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
