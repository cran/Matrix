## METHODS FOR GENERIC: cbind2, rbind2
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## GOAL: write methods that preallocate, i.e., do _not_ use [cr]bind2,
##       and maybe implement c.Matrix and *.sparseVector similarly ... ?
if(TRUE) {
.cbind <- function(..., deparse.level = 1)
    .External(R_bind, deparse.level, 1L, substitute(list(...)), ...)
.rbind <- function(..., deparse.level = 1)
    .External(R_bind, deparse.level, 0L, substitute(list(...)), ...)
}

bindDim <- function(d.x, d.y, margin) {
    r <- d.x
    if(d.x[margin] != d.y[margin]) {
        if(margin == 1L)
            stop("number of rows of matrices must match")
        else stop("number of columns of matrices must match")
    }
    n.x <- d.x[-margin]
    n.y <- d.y[-margin]
    if(n.y > .Machine$integer.max - n.x)
        stop("dimensions cannot exceed 2^31-1")
    r[-margin] <- n.x + n.y
    r
}

bindDimnames <- function(dn.x, dn.y, d.x, d.y, margin) {
    r <- list(NULL, NULL)
    if(!(is.null(tmp <- dn.x[[margin]]) && is.null(tmp <- dn.y[[margin]])))
        r[[margin]] <- tmp
    nms.x <- dn.x[[-margin]]
    nms.y <- dn.y[[-margin]]
    if(!(is.null(nms.x) && is.null(nms.y)))
        r[[-margin]] <-
            c(if(is.null(nms.x)) character(d.x[-margin]) else nms.x,
              if(is.null(nms.y)) character(d.y[-margin]) else nms.y)
    r
}


## ==== Trivial special cases ==========================================

setMethod("cbind2", signature(x = "Matrix", y = "missing"),
          function(x, y, ...) x)
setMethod("cbind2", signature(x = "Matrix", y = "NULL"),
          function(x, y, ...) x)
setMethod("cbind2", signature(x = "NULL", y = "Matrix"),
          function(x, y, ...) y)
if(FALSE) {
## Correct, but breaks evclust ... leaving for 1.6-2
setMethod("cbind2", signature(x = "Matrix", y = "vector"),
          function(x, y, ...) cbind2(x, matrix(y, x@Dim[1L], 1L)))
setMethod("cbind2", signature(x = "vector", y = "Matrix"),
          function(x, y, ...) cbind2(matrix(x, y@Dim[1L], 1L), y))
} else {
setMethod("cbind2", signature(x = "Matrix", y = "vector"),
          function(x, y, ...) cbind2(x, matrix(y, nrow = x@Dim[1L])))
setMethod("cbind2", signature(x = "vector", y = "Matrix"),
          function(x, y, ...) cbind2(matrix(x, nrow = y@Dim[1L]), y))
}

setMethod("rbind2", signature(x = "Matrix", y = "missing"),
          function(x, y, ...) x)
setMethod("rbind2", signature(x = "Matrix", y = "NULL"),
          function(x, y, ...) x)
setMethod("rbind2", signature(x = "NULL", y = "Matrix"),
          function(x, y, ...) y)
if(FALSE) {
## Correct, but breaks evclust ... leaving for 1.6-2
setMethod("rbind2", signature(x = "Matrix", y = "vector"),
          function(x, y, ...) rbind2(x, matrix(y, 1L, x@Dim[2L])))
setMethod("rbind2", signature(x = "vector", y = "Matrix"),
          function(x, y, ...) rbind2(matrix(x, 1L, y@Dim[2L]), y))
} else {
setMethod("rbind2", signature(x = "Matrix", y = "vector"),
          function(x, y, ...) rbind2(x, matrix(y, ncol = x@Dim[2L])))
setMethod("rbind2", signature(x = "vector", y = "Matrix"),
          function(x, y, ...) rbind2(matrix(x, ncol = y@Dim[2L]), y))
}


###-- General -----------------------------------------------------------

###-- Dense, incl Diagonal ----------------------------------------------

###-- Sparse ------------------------------------------------------------

setMethod("cbind2", signature(x = "sparseMatrix", y = "matrix"),
          function(x, y, ...) cbind2(x, .m2sparse(y, ".gC")))
setMethod("cbind2", signature(x = "matrix", y = "sparseMatrix"),
          function(x, y, ...) cbind2(.m2sparse(x, ".gC"), y))
setMethod("rbind2", signature(x = "sparseMatrix", y = "matrix"),
          function(x, y, ...) rbind2(x, .m2sparse(y, ".gC")))
setMethod("rbind2", signature(x = "matrix", y = "sparseMatrix"),
          function(x, y, ...) rbind2(.m2sparse(x, ".gC"), y))

## originally from ./Matrix.R : -------------------------------

## Makes sure one gets x decent error message for the unimplemented cases:
setMethod("cbind2", signature(x = "Matrix", y = "Matrix"),
          function(x, y, ...) {
              bindDim(x@Dim, y@Dim, 1L)
              .bail.out.2("cbind2", class(x), class(y))
          })

## Use a working fall back {particularly useful for sparse}:
## FIXME: implement rbind2 via "cholmod" for C* and Tsparse ones
setMethod("rbind2", signature(x = "Matrix", y = "Matrix"),
          function(x, y, ...) {
              bindDim(x@Dim, y@Dim, 2L)
              t(cbind2(t(x), t(y)))
          })

## originally from ./denseMatrix.R : -------------------------------

### cbind2
setMethod("cbind2", signature(x = "denseMatrix", y = "numeric"),
          function(x, y, ...) {
              d <- dim(x); nr <- d[1]; nc <- d[2]
              y <- rep_len(y, nr) # 'silent procrustes'
              ## beware of (packed) triangular, symmetric, ...
              x <- .M2gen(x)
              x@x <- c(x@x, as.double(y))
              x@Dim[2] <- nc + 1L
              if(is.character(dn <- x@Dimnames[[2]]))
                  x@Dimnames[[2]] <- c(dn, "")
              x
          })
## the same, (x,y) <-> (y,x):
setMethod("cbind2", signature(x = "numeric", y = "denseMatrix"),
          function(x, y, ...) {
              d <- dim(y); nr <- d[1]; nc <- d[2]
              x <- rep_len(x, nr)
              y <- .M2gen(y)
              y@x <- c(as.double(x), y@x)
              y@Dim[2] <- nc + 1L
              if(is.character(dn <- y@Dimnames[[2]]))
                  y@Dimnames[[2]] <- c("", dn)
              y
          })


setMethod("cbind2", signature(x = "denseMatrix", y = "matrix"),
          function(x, y, ...) cbind2(x, .m2dense(y, ".ge")))
setMethod("cbind2", signature(x = "matrix", y = "denseMatrix"),
          function(x, y, ...) cbind2(.m2dense(x, ".ge"), y))

setMethod("cbind2", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y, ...) {
              d.x <- x@Dim
              d.y <- y@Dim
              d.r <- bindDim(d.x, d.y, 1L)
              ## beware of (packed) triangular, symmetric, ...
              x <- .M2gen(x)
              y <- .M2gen(y)
              xx <- c(x@x, y@x)
              ## be careful, e.g., if we have an 'n' and 'd'
              if(identical((tr <- typeof(xx)), typeof(x@x))) {
                  x@x <- xx
                  x@Dim <- d.r
                  x@Dimnames <- bindDimnames(dimnames(x), dimnames(y), d.x, d.y, 1L)
                  x
              } else if(identical(tr, typeof(y@x))) {
                  y@x <- xx
                  y@Dim <- d.r
                  y@Dimnames <- bindDimnames(dimnames(x), dimnames(y), d.x, d.y, 1L)
                  y
              } else stop("resulting x-slot has different type than x's or y's")
          })

### rbind2 -- analogous to cbind2 --- more to do for @x though:

setMethod("rbind2", signature(x = "denseMatrix", y = "numeric"),
          function(x, y, ...) {
              if(is.character(dn <- x@Dimnames[[1]]))
                  dn <- c(dn, "")
              y <- rbind2(as(x,"matrix"), y)
              new(paste0(.M.kind(y), "geMatrix"), x = c(y),
                  Dim = x@Dim + 1:0, Dimnames = list(dn, x@Dimnames[[2]]))
          })
## the same, (x,y) <-> (y,x):
setMethod("rbind2", signature(x = "numeric", y = "denseMatrix"),
          function(x, y, ...) {
              if(is.character(dn <- y@Dimnames[[1]]))
                  dn <- c("", dn)
              x <- rbind2(x, as(y,"matrix"))
              new(paste0(.M.kind(x), "geMatrix"), x = c(x),
                  Dim = y@Dim + 1:0, Dimnames = list(dn, y@Dimnames[[2]]))
          })

setMethod("rbind2", signature(x = "denseMatrix", y = "matrix"),
          function(x, y, ...) rbind2(x, .m2dense(y, ".ge")))
setMethod("rbind2", signature(x = "matrix", y = "denseMatrix"),
          function(x, y, ...) rbind2(.m2dense(x, ".ge"), y))

setMethod("rbind2", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y, ...) {
              d.x <- x@Dim
              d.y <- y@Dim
              d.r <- bindDim(d.x, d.y, 2L)
              ## beware of (packed) triangular, symmetric, ...
              x <- .M2gen(x)
              y <- .M2gen(y)
              xx <- .Call(R_rbind2_vector, x, y)
              ## be careful, e.g., if we have an 'n' and 'd'
              if(identical((tr <- typeof(xx)), typeof(x@x))) {
                  x@x <- xx
                  x@Dim <- d.r
                  x@Dimnames <- bindDimnames(dimnames(x), dimnames(y), d.x, d.y, 2L)
                  x
              } else if(identical(tr, typeof(y@x))) {
                  y@x <- xx
                  y@Dim <- d.r
                  y@Dimnames <- bindDimnames(dimnames(x), dimnames(y), d.x, d.y, 2L)
                  y
              } else stop("resulting x-slot has different type than x's or y's")
          })

## originally from ./diagMatrix.R : --------------------------------------

## For diagonalMatrix:  preserve sparseness {not always optimal, but "the law"}

setMethod("cbind2", signature(x = "diagonalMatrix", y = "sparseMatrix"),
          function(x, y, ...)
              cbind2(.diag2sparse(x, "g", "C"), .M2C(y)))
setMethod("cbind2", signature(x = "sparseMatrix", y = "diagonalMatrix"),
          function(x, y, ...)
              cbind2(.M2C(x), .diag2sparse(y, "g", "C")))
setMethod("rbind2", signature(x = "diagonalMatrix", y = "sparseMatrix"),
          function(x, y, ...)
              rbind2(.diag2sparse(x, "g", "C"), .M2C(y)))
setMethod("rbind2", signature(x = "sparseMatrix", y = "diagonalMatrix"),
          function(x, y, ...)
              rbind2(.M2C(x), .diag2sparse(y, "g", "C")))

## in order to evade method dispatch ambiguity, but still remain "general"
## we use this hack instead of signature  x = "diagonalMatrix"
for(cls in names(getClassDef("diagonalMatrix")@subclasses)) {

setMethod("cbind2", signature(x = cls, y = "matrix"),
          function(x, y, ...)
              cbind2(.diag2sparse(x, "g", "C"), .m2sparse(y, ".gC")))
setMethod("cbind2", signature(x = "matrix", y = cls),
          function(x, y, ...)
              cbind2(.m2sparse(x, ".gC"), .diag2sparse(y, "g", "C")))
setMethod("rbind2", signature(x = cls, y = "matrix"),
          function(x, y, ...)
              rbind2(.diag2sparse(x, "g", "C"), .m2sparse(y, ".gC")))
setMethod("rbind2", signature(x = "matrix", y = cls),
          function(x, y, ...)
              rbind2(.m2sparse(x, ".gC"), .diag2sparse(y, "g", "C")))

 ## These are already defined for "Matrix"
 ## -- repeated here for method dispatch disambiguation {"design-FIXME" ?}
setMethod("cbind2", signature(x = cls, y = "vector"),
          function(x, y, ...) cbind2(x, matrix(y, nrow = nrow(x))))
setMethod("cbind2", signature(x = "vector", y = cls),
          function(x, y, ...) cbind2(matrix(x, nrow = nrow(y)), y))
setMethod("rbind2", signature(x = cls, y = "vector"),
          function(x, y, ...) rbind2(x, matrix(y, ncol = ncol(x))))
setMethod("rbind2", signature(x = "vector", y = cls),
          function(x, y, ...) rbind2(matrix(x, ncol = ncol(y)), y))

}
rm(cls)

## originally from ./dsparseMatrix.R : --------------------------------

## FIXME: dimnames() handling should happen in C code
## ------> ../src/Csparse.c

## Fast - almost non-checking methods
.cbind2Csp <- function(x, y) .Call(Csparse_horzcat, asCspN(x), asCspN(y))
.rbind2Csp <- function(x, y) .Call(Csparse_vertcat, asCspN(x), asCspN(y))

cbind2sparse <- function(x, y) {
    ## beware of (packed) triangular, symmetric, ...
    if(identical(c(dnx <- dimnames(x),
                   dny <- dimnames(y)),
                 list(NULL, NULL, NULL, NULL)))
        ## keep empty dimnames
        .cbind2Csp(x, y)
    else {
        ## R and S+ are different in which names they take
        ## if they differ -- but there's no warning in any case
        rn <-
            if(!is.null(dnx[[1]]))
                dnx[[1]]
            else dny[[1]]
        cx <- dnx[[2]]
        cy <- dny[[2]]
        cn <-
            if(is.null(cx) && is.null(cy))
                NULL
            else c(if(!is.null(cx)) cx else character(ncol(x)),
                   if(!is.null(cy)) cy else character(ncol(y)))
        ans <- .cbind2Csp(x, y)
        ans@Dimnames <- list(rn, cn)
        ans
    }
}
setMethod("cbind2", signature(x = "sparseMatrix", y = "sparseMatrix"),
          function(x, y, ...) {
              bindDim(x@Dim, y@Dim, 1L)
              cbind2sparse(x, y)
          })

rbind2sparse <- function(x, y) {
    ## beware of (packed) triangular, symmetric, ...
    if(identical(c(dnx <- dimnames(x),
                   dny <- dimnames(y)),
                 list(NULL, NULL, NULL, NULL)))
        ## keep empty dimnames
        .rbind2Csp(x, y)
    else {
        ## R and S+ are different in which names they take
        ## if they differ -- but there's no warning in any case
        cn <-
            if(!is.null(dnx[[2]]))
                dnx[[2]]
            else dny[[2]]
        rx <- dnx[[1]] ; ry <- dny[[1]]
        rn <- if(is.null(rx) && is.null(ry))
                  NULL
              else c(if(!is.null(rx)) rx else character(nrow(x)),
                     if(!is.null(ry)) ry else character(nrow(y)))
        ans <- .rbind2Csp(x, y)
        ans@Dimnames <- list(rn, cn)
        ans
    }
}
setMethod("rbind2", signature(x = "sparseMatrix", y = "sparseMatrix"),
          function(x, y, ...) {
              bindDim(x@Dim, y@Dim, 2L)
              rbind2sparse(x, y)
          })

setMethod("cbind2", signature(x = "sparseMatrix", y = "denseMatrix"),
          function(x, y, sparse = NA, ...) {
              d.r <- bindDim(x@Dim, y@Dim, 1L)
              if(is.na(sparse))
                  sparse <-
                      2 * (nnzero(x, na.counted = TRUE) +
                           nnzero(y, na.counted = TRUE)) <
                      as.double(d.r[1L]) * (ncol(x) + ncol(y))
              if(sparse)
                  cbind2sparse(x, y)
              else cbind2(as(x, "denseMatrix"), y)
          })
setMethod("cbind2", signature(x = "denseMatrix", y = "sparseMatrix"),
          function(x, y, sparse = NA, ...) {
              d.r <- bindDim(x@Dim, y@Dim, 1L)
              if(is.na(sparse))
                  sparse <-
                      2 * (nnzero(x, na.counted = TRUE) +
                           nnzero(y, na.counted = TRUE)) <
                      as.double(d.r[1L]) * (ncol(x) + ncol(y))
              if(sparse)
                  cbind2sparse(x, y)
              else cbind2(x, as(y, "denseMatrix"))
          })
setMethod("rbind2", signature(x = "sparseMatrix", y = "denseMatrix"),
          function(x, y, sparse = NA, ...) {
              d.r <- bindDim(x@Dim, y@Dim, 2L)
              if(is.na(sparse))
                  sparse <-
                      2 * (nnzero(x, na.counted = TRUE) +
                           nnzero(y, na.counted = TRUE)) <
                      (nrow(x) + nrow(y)) * as.double(d.r[2L])
              if(sparse)
                  rbind2sparse(x, y)
              else rbind2(as(x, "denseMatrix"), y)
          })
setMethod("rbind2", signature(x = "denseMatrix", y = "sparseMatrix"),
          function(x, y, sparse = NA, ...) {
              d.r <- bindDim(x@Dim, y@Dim, 2L)
              if(is.na(sparse))
                  sparse <-
                      2 * (nnzero(x, na.counted = TRUE) +
                           nnzero(y, na.counted = TRUE)) <
                      (nrow(x) + nrow(y)) * as.double(d.r[2L])
              if(sparse)
                  rbind2sparse(x, y)
              else rbind2(x, as(y, "denseMatrix"))
          })

if(FALSE) {
## FIXME
##------------- maybe a bit faster --- but too much to maintain
## would have to be done for "rbind2" as well ...
setMethod("cbind2", signature(x = "sparseMatrix", y = "numeric"),
          function(x, y, ...) {
              d <- dim(x); nr <- d[1]; nc <- d[2]; cl <- class(x)
              x <- as(x, "CsparseMatrix")
              if(nr > 0) {
                  y <- rep_len(y, nr) # 'silent procrustes'
                  n0y <- y != 0
                  n.e <- length(x@i)
                  x@i <- c(x@i, (0:(nr-1))[n0y])
                  x@p <- c(x@p, n.e + sum(n0y))
                  x@x <- c(x@x, y[n0y])
              }
              x@Dim[2] <- nc + 1L
              if(is.character(dn <- x@Dimnames[[2]]))
                  x@Dimnames[[2]] <- c(dn, "")
              x
          })
## the same, (x,y) <-> (y,x):
setMethod("cbind2", signature(x = "numeric", y = "sparseMatrix"),
          function(x, y, ...) {
              d <- dim(y); nr <- d[1]; nc <- d[2]; cl <- class(y)
              y <-  as(y, "CsparseMatrix")
              if(nr > 0) {
                  x <- rep_len(x, nr) # 'silent procrustes'
                  n0x <- x != 0
                  y@i <- c((0:(nr-1))[n0x], y@i)
                  y@p <- c(0L, sum(n0x) + y@p)
                  y@x <- c(x[n0x], y@x)
              }
              y@Dim[2] <- nc + 1L
              if(is.character(dn <- y@Dimnames[[2]]))
                  y@Dimnames[[2]] <- c(dn, "")
              y
          })
}## -- no longer

setMethod("rbind2", signature(x = "indMatrix", y = "indMatrix"),
          function(x, y, ...) {
              if(x@margin != 1L || y@margin != 1L)
                  return(rbind2(as(x, "RsparseMatrix"), as(y, "RsparseMatrix")))
              d.x <- x@Dim
              d.y <- y@Dim
              r <- new("indMatrix")
              r@Dim <- bindDim(d.x, d.y, 2L)
              r@Dimnames <- bindDimnames(x@Dimnames, y@Dimnames, d.x, d.y, 2L)
              r@perm <- c(x@perm, y@perm)
              r
          })

setMethod("cbind2", signature(x = "indMatrix", y = "indMatrix"),
          function(x, y, ...) {
              if(x@margin == 1L || y@margin == 1L)
                  return(cbind2(as(x, "CsparseMatrix"), as(y, "CsparseMatrix")))
              d.x <- x@Dim
              d.y <- y@Dim
              r <- new("indMatrix")
              r@Dim <- bindDim(d.x, d.y, 1L)
              r@Dimnames <- bindDimnames(x@Dimnames, y@Dimnames, d.x, d.y, 1L)
              r@perm <- c(x@perm, y@perm)
              r@margin <- 2L
              r
          })
