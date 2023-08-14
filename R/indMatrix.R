## METHODS FOR CLASS: indMatrix
## index matrices, i.e., matrices with standard unit vectors
## for all rows _or_ all columns
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.perm2ind <- function(perm, n, margin = 1L, check.p = FALSE) {
    perm.i <- perm
    if(!is.numeric(perm))
        stop("'perm' must be numeric")
    else if(anyNA(r <- range(perm)) || r[1L] < 1L ||
            (is.double(perm) && any(perm != (perm.i <- as.integer(perm)))))
        stop("elements of 'perm' must be positive integers")
    else if((m.i <- length(perm)) > (M <- .Machine$integer.max) || r[2L] > M)
        stop("dimensions cannot exceed 2^31-1")

    if(missing(n))
        n.i <- as.integer(r[2L])
    else {
        n.i <- n
        if(!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0L ||
           (is.double(n) && n != (n.i <- as.integer(n))))
            stop("'n' must be a non-negative integer")
        else if(n > M)
            stop("dimensions cannot exceed 2^31-1")
        else if(r[2L] > n)
            stop("elements of 'perm' cannot exceed 'n'")
    }

    if(!is.numeric(margin) || length(margin) != 1L || is.na(margin) ||
       (margin != 1L && margin != 2L))
        stop("'margin' must be 1 or 2")
    margin.i <- as.integer(margin)

    give.p <- check.p && m.i == n.i &&
        (m.i == 0L || (all(r == c(1L, m.i)) && !anyDuplicated.default(perm.i)))

    J <- new(if(give.p) "pMatrix" else "indMatrix")
    nms <- names(perm)
    if(margin.i == 1L) {
        J@Dim <- c(m.i, n.i)
        J@Dimnames = list(nms, if(give.p) nms)
    } else {
        J@Dim <- c(n.i, m.i)
        J@Dimnames = list(if(give.p) nms, nms)
        J@margin <- 2L
    }
    J@perm <- perm.i
    J
}

setAs("numeric", "indMatrix",
      function(from) {
          J <- new("indMatrix")
          if((m <- length(from)) == 0L)
              return(J)
          from.i <- from
          if(anyNA(r <- range(from)) || r[1L] < 1L ||
             (is.double(from) && any(from != (from.i <- as.integer(from)))))
              stop("elements of 'perm' slot must be positive integers")
          if(m > (M <- .Machine$integer.max) || r[2L] > M)
              stop("dimensions cannot exceed 2^31-1")
          J@Dim <- c(m, as.integer(r[2L]))
          J@Dimnames <- list(names(from), NULL)
          J@perm <- from.i
          J
      })

## FIXME: deprecate this method and export more general function .perm2ind
setAs("list", "indMatrix",
      function(from) do.call(.perm2ind, unname(from)))

setAs("nsparseMatrix", "indMatrix",
      function(from) {
          from <- .M2gen(from)
          J <- new("indMatrix")
          J@Dim <- from@Dim
          J@Dimnames <- from@Dimnames
          from. <- as(from, "RsparseMatrix")
          p <- from.@p
          m <- length(p) - 1L
          if(all(p == 0:m)) {
              J@perm <- from.@j + 1L
              return(J)
          }
          from. <- as(from, "CsparseMatrix")
          p <- from.@p
          n <- length(p) - 1L
          if(all(p == 0:n)) {
              J@perm <- from.@i + 1L
              J@margin <- 2L
              return(J)
          }
          stop("matrix must have exactly one nonzero element in each row or column")
      })

setMethod("isSymmetric", signature(object = "indMatrix"),
          function(object, checkDN = TRUE, ...) {
              d <- object@Dim
              if((n <- d[1L]) != d[2L])
                  return(FALSE)
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...) check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              perm <- object@perm
              all(perm[perm] == seq_len(n))
          })

setMethod("isTriangular", signature(object = "indMatrix"),
          function(object, upper = NA, ...) {
              d <- object@Dim
              if((n <- d[1L]) != d[2L])
                  return(FALSE)
              if(object@margin == 1L) {
                  i <- seq_len(n)
                  j <- object@perm
              } else {
                  i <- object@perm
                  j <- seq_len(n)
              }
              if(is.na(upper)) {
                  if(all(j >= i))
                      return(`attr<-`(TRUE, "kind", "U"))
                  if(all(i <= j))
                      return(`attr<-`(TRUE, "kind", "L"))
                  FALSE
              } else if(upper) {
                  all(j >= i)
              } else {
                  all(i <= j)
              }
          })

setMethod("isDiagonal", signature(object = "indMatrix"),
          function(object) {
              d <- object@Dim
              if((n <- d[1L]) != d[2L])
                  return(FALSE)
              all(object@perm == seq_len(n))
          })

setMethod("t", signature(x = "indMatrix"),
          function(x) {
              r <- new("indMatrix")
              r@Dim <- x@Dim[2:1]
              r@Dimnames = x@Dimnames[2:1]
              r@perm <- x@perm
              if(x@margin == 1L)
                  r@margin <- 2L
              r
          })

setMethod("diag", signature(x = "indMatrix"),
          function(x, nrow, ncol, names = TRUE) {
              if((m <- min(x@Dim)) == 0L)
                  return(logical(0L))
              i <- seq_len(m)
              r <- x@perm[i] == i
              if(names &&
                 !any(vapply(dn <- x@Dimnames, is.null, NA)) &&
                 identical(nms <- dn[[1L]][i], dn[[2L]][i]))
                  names(r) <- nms
              r
          })

setMethod("diag<-", signature(x = "indMatrix"),
          function(x, value) `diag<-`(as(x, "nsparseMatrix"), value))

setMethod("band", signature(x = "indMatrix"),
          function(x, k1, k2, ...) band(as(x, "nsparseMatrix"), k1, k2))

setMethod("triu", signature(x = "indMatrix"),
          function(x, k = 0L, ...) triu(as(x, "nsparseMatrix")))

setMethod("tril", signature(x = "indMatrix"),
          function(x, k = 0L, ...) tril(as(x, "nsparseMatrix")))

setMethod("forceSymmetric", signature(x = "indMatrix", uplo = "missing"),
          function(x, uplo) forceSymmetric(as(x, "nsparseMatrix")))

setMethod("forceSymmetric", signature(x = "indMatrix", uplo = "character"),
          function(x, uplo) forceSymmetric(as(x, "nsparseMatrix"), uplo))

setMethod("symmpart", signature(x = "indMatrix"),
          function(x) symmpart(as(x, "dsparseMatrix")))

setMethod("skewpart", signature(x = "indMatrix"),
          function(x) skewpart(as(x, "dsparseMatrix")))

setMethod("%*%", signature(x = "indMatrix", y = "indMatrix"),
          function(x, y) {
              mx <- x@margin
              my <- y@margin
              px <- x@perm
              py <- y@perm
              r <- new(if(mx == my)
                           "indMatrix"
                       else if(mx == 1L)
                           "dgeMatrix"
                       else "dgTMatrix")
              r@Dim <- mmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              if(mx == my)
                  r@perm <- if(mx == 1L) py[px] else { r@margin <- 2L; px[py] }
              else if(mx == 1L)
                  r@x <- as.double(px == rep(py, each = length(px)))
              else {
                  r@i <- px - 1L
                  r@j <- py - 1L
                  r@x <- rep.int(1, length(px))
              }
              r
          })

setMethod("%*%", signature(x = "indMatrix", y = "matrix"),
          function(x, y) {
              if(x@margin != 1L)
                  return(as(x, "dsparseMatrix") %*% y)
              mmultDim(x@Dim, dim(y), type = 1L)
              r <- .m2dense(y[x@perm, , drop = FALSE], "dge")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%*%", signature(x = "matrix", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %*% as(y, "dsparseMatrix"))
              mmultDim(dim(x), y@Dim, type = 1L)
              r <- .m2dense(x[, y@perm, drop = FALSE], "dge")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%*%", signature(x = "indMatrix", y = "Matrix"),
          function(x, y) {
              if(x@margin != 1L)
                  return(as(x, "dsparseMatrix") %*% y)
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- as(y[x@perm, , drop = FALSE], "dMatrix")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%*%", signature(x = "Matrix", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %*% as(y, "dsparseMatrix"))
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- as(x[, y@perm, drop = FALSE], "dMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "indMatrix", y = "indMatrix"),
          function(x, y) {
              mx <- x@margin
              my <- y@margin
              px <- x@perm
              py <- y@perm
              r <- new(if(mx == my)
                           "indMatrix"
                       else if(mx == 1L)
                           "ngeMatrix"
                       else "ngTMatrix")
              r@Dim <- mmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- mmultDimnames(x@Dimnames, y@Dimnames, type = 1L)
              if(mx == my)
                  r@perm <- if(mx == 1L) py[px] else { r@margin <- 2L; px[py] }
              else if(mx == 1L)
                  r@x <- px == rep(py, each = length(px))
              else {
                  r@i <- px - 1L
                  r@j <- py - 1L
              }
              r
          })

setMethod("%&%", signature(x = "indMatrix", y = "matrix"),
          function(x, y) {
              if(x@margin != 1L)
                  return(as(x, "nsparseMatrix") %&% y)
              mmultDim(x@Dim, dim(y), type = 1L)
              r <- .m2dense(y[x@perm, , drop = FALSE], "nge")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "matrix", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %&% as(y, "nsparseMatrix"))
              mmultDim(dim(x), y@Dim, type = 1L)
              r <- .m2dense(x[, y@perm, drop = FALSE], "nge")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "indMatrix", y = "Matrix"),
          function(x, y) {
              if(x@margin != 1L)
                  return(as(x, "nsparseMatrix") %&% y)
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- as(y[x@perm, , drop = FALSE], "nMatrix")
              r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "Matrix", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %&% as(y, "nsparseMatrix"))
              mmultDim(x@Dim, y@Dim, type = 1L)
              r <- as(x[, y@perm, drop = FALSE], "nMatrix")
              r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("crossprod", signature(x = "indMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              if(x@margin != 1L)
                  return(tcrossprod(t(x), boolArith = boolArith, ...))
              n <- x@Dim[2L]
              tt <- tabulate(x@perm, n)
              if(isTRUE(boolArith)) {
                  r <- new("ldiMatrix")
                  r@x <- as.logical(tt)
              } else {
                  r <- new("ddiMatrix")
                  r@x <- as.double(tt)
              }
              r@Dim <- c(n, n)
              r@Dimnames <- x@Dimnames[c(2L, 2L)]
              r
          })

setMethod("crossprod", signature(x = "indMatrix", y = "matrix"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(isTRUE(boolArith)) `%&%` else `%*%`)(t(x), y))

setMethod("crossprod", signature(x = "matrix", y = "indMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(dim(x), y@Dim, type = 2L)
              boolArith <- isTRUE(boolArith)
              if(y@margin == 1L) {
                  cl <- if(boolArith) "nsparseMatrix" else "dsparseMatrix"
                  r <- crossprod(x, as(y, cl), boolArith = boolArith, ...)
              } else {
                  r <- .m2dense(t(x)[, y@perm, drop = FALSE],
                                if(boolArith) "nge" else "dge")
                  r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames,
                                              type = 2L)
              }
              r
          })

setMethod("crossprod", signature(x = "indMatrix", y = "Matrix"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(isTRUE(boolArith)) `%&%` else `%*%`)(t(x), y))

setMethod("crossprod", signature(x = "Matrix", y = "indMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 2L)
              boolArith <- isTRUE(boolArith)
              if(y@margin == 1L) {
                  cl <- if(boolArith) "nsparseMatrix" else "dsparseMatrix"
                  r <- crossprod(x, as(y, cl), boolArith = boolArith, ...)
              } else {
                  cl <- if(boolArith) "nMatrix" else "dMatrix"
                  r <- as(t(x)[, y@perm, drop = FALSE], cl)
                  r@Dimnames <- mmultDimnames(dimnames(x), y@Dimnames,
                                              type = 2L)
              }
              r
          })

setMethod("tcrossprod", signature(x = "indMatrix", y = "missing"),
          function(x, y = NULL, boolArith = TRUE, ...) {
              if(x@margin != 1L)
                  return(crossprod(t(x), boolArith = boolArith, ...))
              if(isTRUE(boolArith)) {
                  r <- new("ngeMatrix")
                  r@x <- as.vector(
                      `storage.mode<-`(as(x, "matrix"), "logical")[, x@perm])
              } else {
                  r <- new("dgeMatrix")
                  r@x <- as.vector(
                      `storage.mode<-`(as(x, "matrix"),  "double")[, x@perm])
              }
              r@Dim <- x@Dim[c(1L, 1L)]
              r@Dimnames <- x@Dimnames[c(1L, 1L)]
              r
          })

setMethod("tcrossprod", signature(x = "indMatrix", y = "matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, dim(y), type = 3L)
              boolArith <- isTRUE(boolArith)
              if(y@margin == 1L) {
                  r <- .m2dense(t(y)[x@perm, , drop = FALSE],
                                if(boolArith) "nge" else "dge")
                  r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y),
                                              type = 3L)
              } else {
                  cl <- if(boolArith) "nsparseMatrix" else "dsparseMatrix"
                  r <- tcrossprod(as(x, cl), y, boolArith = boolArith, ...)
              }
              r
          })

setMethod("tcrossprod", signature(x = "matrix", y = "indMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(isTRUE(boolArith)) `%&%` else `%*%`)(x, t(y)))

setMethod("tcrossprod", signature(x = "indMatrix", y = "Matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              mmultDim(x@Dim, y@Dim, type = 3L)
              boolArith <- isTRUE(boolArith)
              if(y@margin == 1L) {
                  cl <- if(boolArith) "nMatrix" else "dMatrix"
                  r <- as(t(y)[x@perm, , drop = FALSE], cl)
                  r@Dimnames <- mmultDimnames(x@Dimnames, dimnames(y),
                                              type = 3L)
              } else {
                  cl <- if(boolArith) "nsparseMatrix" else "dsparseMatrix"
                  r <- tcrossprod(as(x, cl), y, boolArith = boolArith, ...)
              }
              r
          })

setMethod("tcrossprod", signature(x = "Matrix", y = "indMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(isTRUE(boolArith)) `%&%` else `%*%`)(x, t(y)))


.indMatrix.sub <- function(x, i, j, ..., value) {
    x <- as(x, "TsparseMatrix")
    callGeneric()
}
for (.i in c("missing", "index"))
for (.j in c("missing", "index"))
setReplaceMethod("[", signature(x = "indMatrix", i = .i, j = .j, value = "ANY"),
                 .indMatrix.sub)
rm(.indMatrix.sub, .i, .j)
