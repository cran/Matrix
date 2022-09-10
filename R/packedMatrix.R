## METHODS FOR CLASS: packedMatrix (virtual)
## dense triangular or symmetric matrices with packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.pM.subclasses <- names(getClassDef("packedMatrix")@subclasses)


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## as(<packedMatrix>,             "matrix") inherited from denseMatrix
## as(<packedMatrix>,        "[dln]Matrix") inherited from denseMatrix
## as(<packedMatrix>,   "[dln]denseMatrix") inherited from denseMatrix
## as(<packedMatrix>,  "[dln]sparseMatrix") inherited from denseMatrix
## as(<packedMatrix>,      "generalMatrix") inherited from denseMatrix
## as(<packedMatrix>,   "triangularMatrix") inherited from      Matrix
## as(<packedMatrix>,    "symmetricMatrix") inherited from      Matrix
## as(<packedMatrix>,     "unpackedMatrix") inherited from denseMatrix
## as(<packedMatrix>, "[CRT]?sparseMatrix") inherited from denseMatrix
## as(<packedMatrix>,     "diagonalMatrix") inherited from      Matrix
## as(<packedMatrix>,          "indMatrix") inherited from      Matrix


## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## as(           <matrix>, "packedMatrix") in ./denseMatrix.R
## as(      <denseMatrix>, "packedMatrix") in ./denseMatrix.R
## as(<[CRT]sparseMatrix>, "packedMatrix") in ./sparseMatrix.R
## as(   <diagonalMatrix>, "packedMatrix") in ./diagMatrix.R
## as(        <indMatrix>, "packedMatrix") in ./indMatrix.R


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("pack", signature(x = "packedMatrix"),
          function(x, ...) x)
setMethod("unpack", signature(x = "packedMatrix"),
          function(x, ...) .Call(packedMatrix_unpack, x, TRUE))

setMethod("forceSymmetric", signature(x = "packedMatrix", uplo = "missing"),
          function(x, uplo) .Call(packedMatrix_force_symmetric, x, NULL))
setMethod("forceSymmetric", signature(x = "packedMatrix", uplo = "character"),
          function(x, uplo) .Call(packedMatrix_force_symmetric, x, uplo))

## Not all of these .pM.is.* are used, because all packedMatrix inherit
## from symmetricMatrix or triangularMatrix, and those classes have
## their own methods.  They are retained here somewhat for completeness ...

.pM.is.sy <- function(object, checkDN = TRUE, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    ## requiring exact symmetry (fast):
    .Call(packedMatrix_is_symmetric, object, checkDN)
}

.pM.is.sy.dz <- function(object, tol = 100 * .Machine$double.eps,
                         tol1 = 8 * tol, checkDN = TRUE, ...) {
    if (tol <= 0)
        .
    else {
        ## going via all.equal (slow):
        isSymmetric(unpack(object), tol = tol, tol1 = tol1,
                    checkDN = checkDN, ...)
    }
}
body(.pM.is.sy.dz) <-
    do.call(substitute, list(body(.pM.is.sy.dz), list(. = body(.pM.is.sy))))

.pM.is.tr <- function(object, upper = NA, ...)
    .Call(packedMatrix_is_triangular, object, upper)

.pM.is.di <- function(object) .Call(packedMatrix_is_diagonal, object)

## method for     .spMatrix in ./symmetricMatrix.R
## method for [lni]tpMatrix in ./triangularMatrix.R
for (.cl in grep("^[dz]tpMatrix$", .pM.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl), .pM.is.sy.dz)

setMethod("isDiagonal", signature(object = "packedMatrix"), .pM.is.di)

rm(.pM.is.sy, .pM.is.sy.dz, .pM.is.tr, .pM.is.di, .cl)

setMethod("t", signature(x = "packedMatrix"),
          function(x) .Call(packedMatrix_transpose, x))
setMethod("diag", signature(x = "packedMatrix"),
          function(x, nrow, ncol, names) .Call(packedMatrix_diag_get, x, names))
setMethod("diag<-", signature(x = "packedMatrix"),
          function(x, value) .Call(packedMatrix_diag_set, x, value))

setMethod("symmpart", signature(x = "packedMatrix"),
          function(x) .Call(packedMatrix_symmpart, x))
setMethod("skewpart", signature(x = "packedMatrix"),
          function(x) .Call(packedMatrix_skewpart, x))


## ~~~~ UTILITIES FOR INDEXING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.pM.error.oob <- function() {
    stop("subscript out of bounds")
}
.pM.error.ist <- function(i) {
    if (isS4(i)) {
        stop("invalid subscript (S4) class ", sQuote(class(i)))
    } else {
        stop("invalid subscript type ", sQuote(typeof(i)))
    }
}
.pM.error.dim <- function() {
    stop("incorrect number of dimensions")
}
.pM.error.neg <- function() {
    stop("negative values are not allowed in a matrix subscript")
}
.pM.error.toolong <- function() {
    stop("logical subscript too long")
}

## NOTATION:
## utility .pM.sub<arity>.<kind> subsets "packedMatrix"
## arity=
## 1: x[i]
## 2: x[i, , drop=], x[, j, drop=], x[i, j, drop=]
## kind=
## vec: indexing by logical, numeric, or character vector
##      _or_ array to be treated as dimensionless ...
##      <vector>[<array>] is equivalent to <vector>[as.vector(<array>)]
##      while avoiding a copy
## mat: (only when arity=1)
##      indexing by logical, numeric, or character matrix; "lMatrix" ...
##      dispatches to kind=vec when not numeric or character matrix
##      with 2 columns

.pM.sub1.vec <- function(x, i) {
    n <- x@Dim[1L]
    if (n <= 1L) {
        ## in this case x@x and as(x, "[dln]geMatrix")@x are the same
        return(x@x[i])
    }
    ni <- length(i)
    if (ni == 0L) {
        return(x@x[0L])
    }
    ## obvious inefficiencies here that we could avoid at C level,
    ## but not without rewriting parts of src/main/subscript.c and
    ## src/main/subset.c ...
    switch(mode(i),
           numeric =
               {
                   if (is.object(i)) {
                       ## use numeric representation of factors, etc.
                       class(i) <- NULL
                   }
                   if (any(i <= -1L, na.rm = TRUE)) {
                       i <- seq_len(as.double(n) * n)[i]
                   } else {
                       i <- i[i >= 1L]
                   }
                   .Call(packedMatrix_sub1, x, i)
               },
           logical =
               {
                   i <- seq_len(as.double(n) * n)[i]
                   .Call(packedMatrix_sub1, x, i)
               },
           character =
               {
                   ## emulating 'stringSubscript' in src/main/subscript.c,
                   ## but this case is quite pathological and need (should?)
                   ## not be maintained ...
                   rep.int(x@x[1L][NA], ni)
               },
           .pM.error.ist(i))
}

.pM.sub2.vec <- function(x, i, j, drop) {
    n <- x@Dim[1L]
    index <- list(if (missing(i)) NULL else if (is.null(i)) integer(0L) else i,
                  if (missing(j)) NULL else if (is.null(j)) integer(0L) else j)
    for (pos in 1:2) {
        if (is.null(k <- index[[pos]])) {
            next
        }
        if ((nk <- length(k)) == 0L) {
            index[[pos]] <- integer(0L)
            next
        }
        index[[pos]] <-
            switch(mode(k),
                   numeric =
                       {
                           if (is.object(k)) {
                               ## use numeric representation of factors, etc.
                               class(k) <- NULL
                           }
                           if (any(k >= n + 1, na.rm = TRUE)) {
                               .pM.error.oob()
                           }
                           seq_len(n)[k]
                       },
                   logical =
                       {
                           if (nk > n) {
                               .pM.error.toolong()
                           }
                           seq_len(n)[k]
                       },
                   character =
                       {
                           nms <- dimnames(x)[[pos]]
                           if (is.null(nms) || anyNA(k <- match(k, nms))) {
                               .pM.error.oob()
                           }
                           k
                       },
                   .pM.error.ist(k))
    }
    .Call(packedMatrix_sub2, x, index[[1L]], index[[2L]],
          if (missing(drop)) TRUE else drop)
}

## could easily support indexing by "[dn]Matrix" and "array",
## but leaving out for now
.pM.sub1.mat <- function(x, i) {
    if (is(i, "lMatrix")) {
        return(.pM.sub1.vec(x, as.vector(i)))
    }
    if (is.logical(i) || length(di <- dim(i)) != 2L || di[2L] != 2L) {
        return(.pM.sub1.vec(x, i))
    }
    if (is.numeric(i)) {
        if (is.double(i)) {
            i <- as.integer(i)
            dim(i) <- di
        }
        ## rows containing 0 are deleted, rows containing NA result in NA,
        ## rows containing both are handled according to the first column
        i <- i[i[, 1L] != 0L, , drop = FALSE] # [NA,j] -> [NA,NA]
        i <- i[i[, 2L] != 0L, , drop = FALSE]
        if (dim(i)[1L] == 0L) {
            return(x@x[0L])
        }
        if (any(i < 1L, na.rm = TRUE)) {
            .pM.error.neg()
        }
        if (any(i > x@Dim[1L], na.rm = TRUE)) {
            .pM.error.oob()
        }
        .Call(packedMatrix_sub1_mat, x, i)
    } else if (is.character(i)) {
        if (di[1L] == 0L) {
            return(x@x[0L])
        }
        dn <- dimnames(x)
        m <- c(match(i[, 1L], dn[[1L]]), match(i[, 2L], dn[[2L]]))
        dim(m) <- di
        if (any(rowSums(is.na(i)) == 0L & rowSums(is.na(m)) > 0L)) {
            ## error if character row contains zero NA but integer row
            ## contains at least one NA, indicating nonmatch that cannot
            ## be ignored
            .pM.error.oob()
        }
        .Call(packedMatrix_sub1_mat, x, m)
    } else {
        .pM.error.ist(i)
    }
}


## ~~~~ METHODS FOR INDEXING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("[", signature(x = "packedMatrix", i = "missing", j = "missing", drop = "missing"),
          function(x, i, j, ..., drop) {
              na <- nargs()
              Matrix.msg(sprintf("pM[%s, %s, %s] : nargs() = %d",
                                 "missing", "missing", "missing", na),
                         .M.level = 2)
	      if (na == 2L) {
                  ## x[]
                  x
              } else if (na == 3L) {
                  ## x[, ]
                  if (x@Dim[1L] == 1L) x@x else x # drop=TRUE implicit
              } else {
                  ## x[, , ], etc.
                  .pM.error.dim()
              }
          })

setMethod("[", signature(x = "packedMatrix", i = "missing", j = "missing", drop = "logical"),
          function(x, i, j, ..., drop) {
              na <- nargs()
              Matrix.msg(sprintf("pM[%s, %s, %s] : nargs() = %d",
                                 "missing", "missing", drop, na),
                         .M.level = 2)
              if (na < 4L) {
                  ## x[drop=], x[, drop=], x[drop=, ]
                  x
              } else if (na == 4L) {
                  ## x[, , drop=], x[, drop=, ], x[drop=, , ]
                  if (!isFALSE(drop[1L]) && x@Dim[1L] == 1L) x@x else x
              } else {
                  ## x[, , , drop=], etc.
                  .pM.error.dim()
              }
          })

.cl <- expand.grid(x = "packedMatrix",
                   i = c("missing", "NULL", "index", "matrix", "lMatrix"),
                   j = c("missing", "NULL", "index"),
                   drop = c("missing", "logical"),
                   stringsAsFactors = FALSE)
.ms <- lapply(.cl[c("i", "j", "drop")], `==`, "missing")

## some abstraction here to avoid repetition ...
for (.k in seq_len(nrow(.cl))) {
    if (.ms$i[.k] && .ms$j[.k]) {
        ## both 'i' and 'j' are missing ... methods set outside of loop
        next
    } else if (.ms$i[.k] || .ms$j[.k]) {
        ## exactly one of 'i' and 'j' is missing
        .i1 <- if (.ms$i[.k]) quote(j) else quote(i)
        .f1 <- if (grepl("[mM]atrix", .cl[.k, as.character(.i1)])) quote(.pM.sub1.mat) else quote(.pM.sub1.vec)
        .definition <- eval(bquote({
            function(x, i, j, ..., drop) {
                na <- nargs()
                Matrix.msg(sprintf("pM[%s, %s, %s] : nargs() = %d",
                                   .(.cl$i[.k]), .(.cl$j[.k]), .(.cl$drop[.k]),
                                   na),
                           .M.level = 2)
                if (na == .(2L + !.ms$drop[.k])) {
                    ## x[i],  x[i,  drop=],
                    ## x[j=], x[j=, drop=]
                    .(.f1)(x, .(.i1))
                } else if (na == .(3L + !.ms$drop[.k])) {
                    ## x[i, ], x[i, , drop=],
                    ## x[, j], x[, j, drop=]
                    .pM.sub2.vec(x, i, j, drop)
                } else {
                    ## x[i, , ], x[i, , , drop=],
                    ## x[, j, ], x[, j, , drop=], etc.
                    .pM.error.dim()
                }
            }
        }))
    } else {
        ## neither 'i' nor 'j' is missing
        .definition <- eval(bquote({
            function(x, i, j, ..., drop) {
                na <- nargs()
                Matrix.msg(sprintf("pM[%s, %s, %s] : nargs() = %d",
                                   .(.cl$i[.k]), .(.cl$j[.k]), .(.cl$drop[.k]),
                                   na),
                           .M.level = 2)
                if (na == .(3L + !.ms$drop[.k])) {
                    ## x[i, j], x[i, j, drop=]
                    .pM.sub2.vec(x, i, j, drop)
                } else {
                    ## x[i, j, , ], x[i, j, , drop=], etc.
                    .pM.error.dim()
                }
            }
        }))
    }
    setMethod("[", do.call(signature, .cl[.k, ]), .definition)
}

rm(.pM.subclasses, .cl, .ms, .k, .i1, .f1, .definition)
