## METHODS FOR GENERIC: [
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.subscript.invalid <- function(i) {
    if(is.object(i))
        gettextf("invalid subscript class \"%s\"", class(i)[1L])
    else gettextf("invalid subscript type \"%s\"", typeof(i))
}

.subscript.recycle <- function(i, n, pattern) {
    ## Return integer or double vector corresponding
    ## to [nl]sparseVector 'i' recycled to length 'n' :
    if(length(i.i <- i@i) == 0L)
        integer(0L)
    else if((i.length <- length(i)) >= n) {
        if(i.length > n) {
            if(n < 0x1p+53) {
                if(i.i[length(i.i)] >= n + 1)
                    i.i[i.i >= n + 1] <- NA
            } else {
                if(i.i[length(i.i)] > n)
                    i.i[i.i > n] <- NA
            }
        }
        if(pattern) i.i else i.i[i@x]
    } else {
        r <- ceiling(n / i.length)
        n. <- r * i.length
        i.i <-
            if(n. <= .Machine$integer.max)
                rep.int(as.integer(i.i), r) +
                    rep(seq.int(from = 0L,
                                by = as.integer(i.length),
                                length.out = r),
                        each = length(i.i))
            else if(i.i[length(i.i)] + (r - 1) * i.length <= 0x1p+53)
                rep.int(as.double(i.i), r) +
                    rep(seq.int(from = 0,
                                by = as.double(i.length),
                                length.out = r),
                        each = length(i.i))
            else stop(gettextf("recycled %s would have maximal index exceeding %s",
                               "[nl]sparseVector", "2^53"),
                      domain = NA)
        if(pattern) {
            if(n. > n) i.i[      i.i <= n] else i.i
        } else {
            if(n. > n) i.i[i@x & i.i <= n] else i.i[i@x]
        }
    }
}

## x[i] where 'i' is NULL or any vector or sparseVector
.subscript.1ary <- function(x, i) {
    x.length <- prod(x@Dim)
    if(is.null(i))
        i <- integer(0L)
    else if(isS4(i)) {
        if(!.isVector(i))
            stop(.subscript.invalid(i), domain = NA)
        kind <- .M.kind(i)
        if((pattern <- kind == "n") || kind == "l") {
            ## [nl]sparseVector
            i <- .subscript.recycle(i, x.length, pattern)
            return(..subscript.1ary(x, i, unsorted = !pattern && anyNA(i)))
        }
        i <- i@x
    }
    switch(typeof(i),
           double =
               {
                   r <- min(1, i, na.rm = TRUE)
                   if(r < 1)
                       i <- if(r <= -1)
                                seq_len(x.length)[i] # FIXME
                            else i[i >= 1]
                   ..subscript.1ary(x, i)
               },
           integer =
               {
                   r <- min(1L, i, na.rm = TRUE)
                   if(r < 1L)
                       i <- if(r <= -1L)
                                seq_len(x.length)[i] # FIXME
                            else i[i >= 1L]
                   ..subscript.1ary(x, i)
               },
           logical =
               {
                   if((i.length <- length(i)) && !is.na(a <- all(i)) && a) {
                       if(i.length <= x.length)
                           as.vector(x)
                       else c(as.vector(x), rep.int(NA, i.length - x.length))
                   } else .subscript.1ary(x, .m2V(i)) # recursively
               },
           character =
               {
                   rep.int(if(.hasSlot(x, "x")) x@x[NA_integer_] else NA,
                           length(i))
               },
           stop(.subscript.invalid(i), domain = NA))
}

## x[i] where 'i' is vector of type "integer" or "double"
## with elements greater than or equal to 1 (or NA)
..subscript.1ary <- function(x, i,
                             shape = .M.shape(x),
                             repr = .M.repr(x),
                             unsorted = is.unsorted(i)) {
    if(!any(repr == c("C", "R", "T")))
        return(.Call(R_subscript_1ary, x, i))
    if(shape == "t" && x@diag != "N")
        x <- ..diagU2N(x)
    if(shape == "s" || repr == "R") {
        x.length <- prod(d <- x@Dim)
        if(x.length < 0x1p+53) {
            r <- max(x.length, i, na.rm = TRUE)
            if(r >= x.length + 1)
                i[i >= x.length + 1] <- NA
        } else if(is.double(i)) {
            r <- max(0x1p+53, i, na.rm = TRUE)
            if(r > 0x1p+53) {
                if(any(i > 0x1p+53 && i <= x.length, na.rm = TRUE))
                    ## could be avoided in C, which has 64-bit integers :
                    warning(gettextf("subscripts exceeding %s replaced with NA", "2^53"),
                            domain = NA)
                i[i > 0x1p+53] <- NA
            }
        }
        i1s <- i - 1L
        m <- d[1L]
        i. <- as.integer(i1s %% m)
        if(shape == "s") {
            j. <- as.integer(i1s %/% m)
            op <- if(x@uplo == "U") `>` else `<`
            if(length(w <- which(op(i., j.)))) {
                i.. <- i.[w]
                j.. <- j.[w]
                if(repr == "R")
                    i.[w] <- j..
                if(x.length > .Machine$integer.max)
                    m <- as.double(m)
                i[w] <- m * i.. + j.. + 1L
            }
        }
    }
    o <-
        if(repr == "R")
            order(i., i)
        else if(is.na(unsorted) || unsorted)
            sort.list(i)
        else return(.Call(R_subscript_1ary, x, i))
    if(is.unsorted(o)) {
        s <- .Call(R_subscript_1ary, x, i[o])
        s[o] <- s
        s
    } else .Call(R_subscript_1ary, x, i)
}

## x[i] where 'i' is any array or Matrix
.subscript.1ary.mat <- function(x, i) {
    if(isS4(i)) {
        if(!.isMatrix(i))
            stop(.subscript.invalid(i), domain = NA)
        if((logic <- any(.M.kind(i) == c("n", "l"))) || i@Dim[2L] != 2L) {
            if(logic && all(i@Dim) && !is.na(a <- all(i)) && a) {
                x <- as.vector(x)
                if((i.length <- length(i)) <= (x.length <- length(x)))
                    return(x)
                else return(c(x, rep.int(NA, i.length - x.length)))
            }
            i <- if(.isDense(i)) .M2v(i) else .M2V(i)
            return(.subscript.1ary(x, i))
        }
        i <- as(i, "matrix")
    } else if(is.logical(i) || length(di <- dim(i)) != 2L || di[2L] != 2L)
        return(.subscript.1ary(x, i))
    switch(typeof(i),
           double =,
           integer =
               {
                   storage.mode(i) <- "integer"
                   if(min(1L, i, na.rm = TRUE) < 1L)
                       stop("negative values are not allowed in a matrix subscript")
                   dx <- x@Dim
                   m <- dx[1L]
                   n <- dx[2L]
                   if(m == n) {
                       if(max(n, i, na.rm = TRUE) > n)
                           stop("subscript out of bounds")
                   } else {
                       if(max(m, i[, 1L], na.rm = TRUE) > m ||
                          max(n, i[, 2L], na.rm = TRUE) > n)
                           stop("subscript out of bounds")
                   }
                   ## * rows containing 0 are deleted
                   ## * rows containing NA are kept
                   ## * rows containing both 0 and NA are handled
                   ##   according to value in first column
                   if(is.na(a <- all(i. <- i[, 1L])) || !a)
                       i <- i[i. > 0L, , drop = FALSE]
                   if(!all(j. <- i[, 2L], na.rm = TRUE))
                       i <- i[j. > 0L, , drop = FALSE]
                   ..subscript.1ary.mat(x, i)
               },
           character =
               {
                   dnx <- dimnames(x)
                   m <- c(match(i[, 1L], dnx[[1L]]), match(i[, 2L], dnx[[2L]]))
                   dim(m) <- di
                   if(any(!rowSums(is.na(i)) & rowSums(is.na(m))))
                       ## error if character row contains zero NA and
                       ## integer row contains at least one NA,
                       ## indicating non-match that should not be ignored
                       stop("subscript out of bounds")
                   ..subscript.1ary.mat(x, m)
               },
           stop(.subscript.invalid(i), domain = NA))
}

## x[i] where 'i' is a 2-column matrix of type "integer"
## with i[, 1L] in 1:m (or NA) and i[, 2L] in 1:n (or NA)
..subscript.1ary.mat <- function(x, i,
                                 shape = .M.shape(x),
                                 repr = .M.repr(x)) {
    if(!any(repr == c("C", "R", "T")))
        return(.Call(R_subscript_1ary_mat, x, i))
    if(shape == "t" && x@diag != "N")
        x <- ..diagU2N(x)
    i. <- i[, 1L]
    j. <- i[, 2L]
    if(shape == "s") {
        op <- if(x@uplo == "U") `>` else `<`
        if(length(w <- which(op(i., j.)))) {
            i[w, ] <- i[w, 2:1]
            i. <- i[, 1L]
            j. <- i[, 2L]
        }
    }
    o <-
        if(repr == "R") {
            if(anyNA(j.))
                i.[is.na(j.)] <- NA
            order(i., j.)
        } else {
            if(anyNA(i.))
                j.[is.na(i.)] <- NA
            order(j., i.)
        }
    if(is.unsorted(o)) {
        s <- .Call(R_subscript_1ary_mat, x, i[o, , drop = FALSE])
        s[o] <- s
        s
    } else .Call(R_subscript_1ary_mat, x, i)
}

## x[i, j, drop] where 'i' and 'j' are NULL or any vector
.subscript.2ary <- function(x, i, j, drop) {
    d <- x@Dim
    l <- list(if(missing(i)) NULL else if(is.null(i)) integer(0L) else i,
              if(missing(j)) NULL else if(is.null(j)) integer(0L) else j)
    for(pos in 1:2) {
        if(!is.null(k <- l[[pos]])) {
            l[pos] <- list(
                switch(typeof(k),
                       double =
                           {
                               r <- d[pos]
                               if(max(r, k, na.rm = TRUE) >= r + 1)
                                   stop("subscript out of bounds")
                               if(min(1, k, na.rm = TRUE) < 1)
                                   seq_len(r)[k]
                               else as.integer(k)
                           },
                       integer =
                           {
                               r <- d[pos]
                               if(max(r, k, na.rm = TRUE) > r)
                                   stop("subscript out of bounds")
                               if(min(1L, k, na.rm = TRUE) < 1L)
                                   seq_len(r)[k]
                               else k
                           },
                       logical =
                           {
                               r <- d[pos]
                               if(length(k) > r)
                                   stop("logical subscript too long")
                               if(length(k) && !is.na(a <- all(k)) && a)
                                   NULL
                               else seq_len(r)[k]
                           },
                       character =
                           {
                               if(length(k) == 0L)
                                   integer(0L)
                               else if(is.null(nms <- dimnames(x)[[pos]]) ||
                                       anyNA(k <- match(k, nms)))
                                   stop("subscript out of bounds")
                               else k
                           },
                       stop(.subscript.invalid(k), domain = NA)))
        }
    }
    if(is.double(lengths(l, use.names = FALSE)))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"), domain = NA)
    ..subscript.2ary(x, l[[1L]], l[[2L]], drop = drop[1L])
}

## x[i, j, drop] where 'i' and 'j' are vectors of type "integer"
## of length not exceeding 2^31-1 with 'i' in 1:m (or NA) and 'j'
## in 1:n (or NA) ... NULL => missing
..subscript.2ary <- function(x, i, j, drop) {
    if(is.null(i) && is.null(j))
        r <- x
    else {
        r <- .Call(R_subscript_2ary, x, i, j)
        dn <- dimnames(x)
        if(!(is.null(i) || is.null(rn <- dn[[1L]])))
            dn[1L] <- list(if(length(i)) rn[i] else NULL)
        if(!(is.null(j) || is.null(cn <- dn[[2L]])))
            dn[2L] <- list(if(length(j)) cn[j] else NULL)
        r@Dimnames <- dn
    }
    if((is.na(drop) || drop) && any(r@Dim == 1L)) drop(as(r, "matrix")) else r
}

setMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na == 2L)
                  ## x[]
                  x
              else if(na == 3L)
                  ## x[, ]
                  drop(x)
              else
                  ## x[, , ], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                         drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na < 4L)
                  ## x[drop=], x[, drop=], x[drop=, ]
                  x
              else if(na == 4L)
                  ## x[, , drop=], x[, drop=, ], x[drop=, , ]
                  if(is.na(drop <- drop[1L]) || drop) drop(x) else x
              else
                  ## x[, , , drop=], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na == 2L)
                  ## x[i=]
                  .subscript.1ary(x, i)
              else if(na == 3L)
                  ## x[i=, ], x[, i=]
                  .subscript.2ary(x, i, , drop = TRUE)
              else
                  ## x[i=, , ], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "missing",
                         drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na == 3L)
                  ## x[i=, drop=]
                  .subscript.1ary(x, i)
              else if(na == 4L)
                  ## x[i=, , drop=], x[, i=, drop=]
                  .subscript.2ary(x, i, , drop = drop)
              else
                  ## x[i=, , , drop=], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "index",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na == 2L)
                  ## x[j=]
                  .subscript.1ary(x, j)
              else if(na == 3L)
                  ## x[j=, ], x[, j=]
                  .subscript.2ary(x, , j, drop = TRUE)
              else
                  ## x[, j=, ], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "index",
                         drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na == 3L)
                  ## x[j=, drop=]
                  .subscript.1ary(x, j)
              else if(na == 4L)
                  ## x[j=, , drop=], x[, j=, drop=]
                  .subscript.2ary(x, , j, drop = drop)
              else
                  ## x[, j=, , drop=], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "index",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na == 3L)
                  ## x[i=, j=], x[j=, i=]
                  .subscript.2ary(x, i, j, drop = TRUE)
              else
                  ## x[i=, j=, ], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "index",
                         drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na == 4L)
                  ## x[i=, j=, drop=], x[j=, i=, drop=]
                  .subscript.2ary(x, i, j, drop = drop)
              else
                  ## x[i=, j=, , drop=], etc.
                  stop("incorrect number of dimensions")
          })

for(.cl in c("matrix", "nMatrix", "lMatrix"))
setMethod("[", signature(x = "Matrix", i = .cl, j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if(na == 2L)
                  ## x[i=]
                  .subscript.1ary.mat(x, i)
              else if(na == 3L)
                  ## x[i=, ], x[, i=]
                  .subscript.2ary(x, i, , drop = TRUE)
              else
                  ## x[i=, , ], etc.
                  stop("incorrect number of dimensions")
          })
rm(.cl)

setMethod("[", signature(x = "Matrix", i = "NULL", j = "ANY",
                         drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              i <- integer(0L)
              callGeneric()
          })

setMethod("[", signature(x = "Matrix", i = "ANY", j = "NULL",
                         drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              j <- integer(0L)
              callGeneric()
          })

setMethod("[", signature(x = "Matrix", i = "NULL", j = "NULL",
                         drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              i <- integer(0L)
              j <- integer(0L)
              callGeneric()
          })

setMethod("[", signature(x = "sparseVector", i = "missing", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              if(nargs() != 2L)
                  stop("incorrect number of dimensions")
              x
          })

setMethod("[", signature(x = "sparseVector", i = "index", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              if(nargs() != 2L)
                  stop("incorrect number of dimensions")
              x.length <- length(x)
              pattern <- .M.kind(x) == "n"
              switch(typeof(i),
                     double =
                         {
                             r <- min(1, i, na.rm = TRUE)
                             if(r <= -1) {
                                 if(r <= -x.length - 1)
                                     i <- i[i > -x.length - 1]
                                 r <- max(-1, i)
                                 if(is.na(r) || r >= 1)
                                     stop("only zeros may be mixed with negative subscripts")
                                 if(r > -1)
                                     i <- i[i <= -1]
                                 d <- unique.default(sort.int(-trunc(i)))
                                 k <- match(x@i, d, 0L) == 0L
                                 x@length <- length(x) - length(d)
                                 x@i <-
                                     {
                                         tmp <- x@i[k]
                                         tmp - findInterval(tmp, d) # !!
                                     }
                                 if(!pattern)
                                     x@x <- x@x[k]
                             } else {
                                 if(r < 1)
                                     i <- i[i >= 1]
                                 if(max(0, i, na.rm = TRUE) >= x.length + 1)
                                     i[i >= x.length + 1] <- NA
                                 if((a <- anyNA(i)) && pattern) {
                                     x <- .V2kind(x, "l")
                                     pattern <- FALSE
                                 }
                                 j <- match(trunc(i), x@i, 0L)
                                 x@length <- length(i)
                                 x@i <-
                                     if(!a)
                                         which(j != 0L)
                                     else {
                                         i. <- is.na(i)
                                         j[i.] <- NA
                                         which(j != 0L | i.)
                                     }
                                 if(!pattern)
                                     x@x <- x@x[j]
                             }
                             x
                         },
                     integer =
                         {
                             r <- min(1L, i, na.rm = TRUE)
                             if(r <= -1L) {
                                 if(r < -x.length)
                                     i <- i[i >= -x.length]
                                 r <- max(-1L, i)
                                 if(is.na(r) || r >= 1L)
                                     stop("only zeros may be mixed with negative subscripts")
                                 if(r > -1L)
                                     i <- i[i <= -1L]
                                 d <- unique.default(sort.int(-i))
                                 k <- is.na(match(x@i, d))
                                 x@length <- length(x) - length(d)
                                 x@i <-
                                     {
                                         tmp <- x@i[k]
                                         tmp - findInterval(tmp, d) # !!
                                     }
                                 if(!pattern)
                                     x@x <- x@x[k]
                             } else {
                                 if(r < 1L)
                                     i <- i[i >= 1L]
                                 if(max(0L, i, na.rm = TRUE) > x.length)
                                     i[i > x.length] <- NA
                                 if((a <- anyNA(i)) && pattern) {
                                     x <- .V2kind(x, "l")
                                     pattern <- FALSE
                                 }
                                 j <- match(i, x@i, 0L)
                                 x@length <- length(i)
                                 x@i <-
                                     if(!a)
                                         which(j != 0L)
                                     else {
                                         i. <- is.na(i)
                                         j[i.] <- NA
                                         which(j != 0L | i.)
                                     }
                                 if(!pattern)
                                     x@x <- x@x[j]
                             }
                             x
                         },
                     logical =
                         {
                             if((i.length <- length(i)) && !is.na(a <- all(i)) && a) {
                                 if(i.length > x.length) {
                                     if(pattern)
                                         x <- .V2kind(x, "l")
                                     x@length <- i.length
                                     x@i <- c(x@i, (x.length + 1):i.length)
                                     x@x <- c(x@x, rep.int(NA, i.length - x.length))
                                 }
                                 x
                             } else x[.m2V(i)] # recursively
                         },
                     stop(.subscript.invalid(i), domain = NA))
          })

setMethod("[", signature(x = "sparseVector", i = "nsparseVector", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              if(nargs() != 2L)
                  stop("incorrect number of dimensions")
              x[.subscript.recycle(i, length(x), TRUE)]
          })

setMethod("[", signature(x = "sparseVector", i = "lsparseVector", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              if(nargs() != 2L)
                  stop("incorrect number of dimensions")
              x[.subscript.recycle(i, length(x), FALSE)]
          })

setMethod("[", signature(x = "sparseVector", i = "NULL", j = "ANY",
                         drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              i <- integer(0L)
              callGeneric()
          })


## METHODS FOR GENERIC: head
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("head", signature(x = "Matrix"),
          head.matrix)

setMethod("head", signature(x = "sparseVector"),
          function(x, n = 6L, ...) {
              stopifnot(is.numeric(n), length(n) == 1L, !is.na(n))
              len <- length(x)
              n <- if(n < 0L) max(len + n, 0L) else min(n, len)
              if(n >= len)
                  return(x)
              nnz <- length(i <- x@i)
              x@length <- n <- if(is.integer(i)) as.integer(n) else trunc(n)
              if(nnz > 0L && i[nnz] > n) {
                  pattern <- .M.kind(x) == "n"
                  if(i[1L] > n) {
                      x@i <- integer(0L)
                      if(!pattern)
                          x@x <- x@x[0L]
                  } else {
                      ii <- 1L:(which.max(i > n) - 1L)
                      x@i <- i[ii]
                      if(!pattern)
                          x@x <- x@x[ii]
                  }
              }
              x
          })


## METHODS FOR GENERIC: tail
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("tail", signature(x = "Matrix"),
          tail.matrix)

setMethod("tail", signature(x = "sparseVector"),
          function(x, n = 6L, ...) {
              stopifnot(is.numeric(n), length(n) == 1L, !is.na(n))
              len <- length(x)
              n <- if(n < 0L) max(len + n, 0L) else min(n, len)
              if(n >= len)
                  return(x)
              nnz <- length(i <- x@i)
              x@length <- n <- if(is.integer(i)) as.integer(n) else trunc(n)
              if(nnz > 0L && i[1L] <= (k <- len - n)) {
                  pattern <- .M.kind(x) == "n"
                  if(i[nnz] <= k) {
                      x@i <- integer(0L)
                      if(!pattern)
                          x@x <- x@x[0L]
                  } else {
                      ii <- which.min(i <= k):nnz
                      x@i <- i[ii] - k
                      if(!pattern)
                          x@x <- x@x[ii]
                  }
              }
              x
          })
