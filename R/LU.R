## METHODS FOR GENERIC: lu
## pivoted LU factorization, returning denseLU or sparseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lu", signature(x = "matrix"),
          function(x, ...) lu(.m2dense(x, ",ge"), ...))

setMethod("lu", signature(x = "denseMatrix"),
          function(x, ...) lu(.M2kind(x, ","), ...))

setMethod("lu", signature(x = "dgeMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(dgeMatrix_trf, x, as.logical(warnSing)))

setMethod("lu", signature(x = "dsyMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["denseLU"]]))
                  return(ch)
              r <- lu(.M2gen(x), ...)
              if(cache) .set.factor(x, "denseLU", r) else r
          })

setMethod("lu", signature(x = "dspMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["denseLU"]]))
                  return(ch)
              r <- lu(.M2gen(x), ...)
              if(cache) .set.factor(x, "denseLU", r) else r
          })

for(.cl in c("dtrMatrix", "dtpMatrix"))
setMethod("lu", signature(x = .cl),
          function(x, ...) {
              if(x@uplo == "U" || x@diag == "U") {
                  r <- new("denseLU")
                  r@Dim <- d <- x@Dim
                  r@Dimnames <- x@Dimnames
                  r@perm <- seq_len(d[1L])
                  r@x <- .M2gen(x)@x
                  r
              } else lu(.M2gen(x), ...)
          })
rm(.cl)

setMethod("lu", signature(x = "sparseMatrix"),
          function(x, ...)
              lu(.M2kind(.M2C(x), ","), ...))

setMethod("lu", signature(x = "dgCMatrix"),
          function(x, errSing = TRUE, order = NA_integer_, tol = 1, ...)
              .Call(dgCMatrix_trf, x, order, tol, errSing))

setMethod("lu", signature(x = "dsCMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2gen(x), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", "dtCMatrix",
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@Dim <- d
                  r@Dimnames <- x@Dimnames
                  r@L <- if(upper) y else x
                  r@U <- if(upper) x else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
              } else lu(.M2gen(x), ...)
          })

setMethod("lu", signature(x = "dgRMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2C(x), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", signature(x = "dsRMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2gen(.tCRT(x)), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", signature(x = "dtRMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@Dim <- d
                  r@Dimnames <- x@Dimnames
                  r@L <- if(upper) y else .M2C(x)
                  r@U <- if(upper) .M2C(x) else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
              } else lu(.M2gen(.M2C(x)), ...)
          })

setMethod("lu", signature(x = "dgTMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2C(x), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", signature(x = "dsTMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2gen(.M2C(x)), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", signature(x = "dtTMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@Dim <- d
                  r@Dimnames <- x@Dimnames
                  r@L <- if(upper) y else .M2C(x)
                  r@U <- if(upper) .M2C(x) else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
                  r
              } else lu(.M2gen(.M2C(x)), ...)
          })

setMethod("lu", "diagonalMatrix",
          function(x, ...) {
              x <- .M2kind(x, ",")
              n <- (d <- x@Dim)[1L]
              L <- new("dtCMatrix")
              r <- new("sparseLU")
              L@Dim <- d
              L@uplo <- "L"
              L@diag <- "U"
              L@p <- integer(n + 1L)
              r@L <- L
              if(x@diag == "N") {
                  L@diag <- "N"
                  L@p <- seq.int(from = 0L, length.out = n + 1L)
                  L@x <- x@x
              }
              r@U <- L
              r@Dim <- d
              r@Dimnames <- x@Dimnames
              r@p <- r@q <- seq.int(from = 0L, length.out = n)
              r
          })


## METHODS FOR CLASS: denseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## TODO: we have C-level functions that do this, but no R interface

.mkL <- function(x, m, n) {
    if(m == n)
        x
    else if(m < n)
        x[seq_len(prod(m, m))]
    else {
        x[seq.int(from = 1L, by = m + 1, length.out = n)] <- 1
        if(n > 1L) {
            nvec <- seq_len(n - 1L)
            from <- seq.int(from = m + 1, by = m, length.out = n - 1L)
            x[sequence.default(nvec, from)] <- 0
        }
        x
    }
}

.mkU <- function(x, m, n) {
    if(m == n)
        x
    else if(m > n) {
        nvec <- rep.int(n, n)
        from <- seq.int(from = 1L, by = m, length.out = n)
        x[sequence.default(nvec, from)]
    } else {
        if(m > 1L) {
            nvec <- seq.int(to = 1L, by = -1L, length.out = m - 1L)
            from <- seq.int(from = 2L, by = m + 1, length.out = m - 1L)
            x[sequence.default(nvec, from)] <- 0
        }
        x
    }
}

setAs("denseLU", "dgeMatrix",
      function(from) {
          to <- new("dgeMatrix")
          to@Dim <- from@Dim
          to@x <- from@x
          to
      })

setMethod("expand1", signature(x = "denseLU"),
          function(x, which, ...) {
              d <- x@Dim
              m <- d[1L]
              n <- d[2L]
              switch(which,
                     "P1" =, "P1." = {
                         r <- new("pMatrix")
                         r@Dim <- c(m, m)
                         r@perm <- asPerm(x@perm, n = m)
                         if(which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "L" = {
                         if(m <= n) {
                             r <- new("dtrMatrix")
                             r@Dim <- c(m, m)
                             r@uplo <- "L"
                             r@diag <- "U"
                         } else {
                             r<- new("dgeMatrix")
                             r@Dim <- d
                         }
                         r@x <- .mkL(x@x, m, n)
                         r
                     },
                     "U" = {
                         if (m >= n) {
                             r <- new("dtrMatrix")
                             r@Dim <- c(n, n)
                         } else {
                             r <- new("dgeMatrix")
                             r@Dim <- d
                         }
                         r@x <- .mkU(x@x, m, n)
                         r
                     },
                     stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%3$s\", or \"%4$s\"",
                                   "which", "P", "L", "U"),
                          domain = NA))
          })

## returning list(P1', L, U), where A = P1' L U
setMethod("expand2", signature(x = "denseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              m <- d[1L]
              n <- d[2L]

              P1. <- new("pMatrix")
              P1.@Dim <- c(m, m)
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@perm <- invertPerm(asPerm(x@perm, n = m))

              if(m <= n) {
                  L <- new("dtrMatrix")
                  L@Dim <- c(m, m)
                  L@uplo <- "L"
                  L@diag <- "U"
              } else {
                  L <- new("dgeMatrix")
                  L@Dim <- d
              }
              L@x <- .mkL(x@x, m, n)

              if (m >= n) {
                  U <- new("dtrMatrix")
                  U@Dim <- c(n, n)
              } else {
                  U <- new("dgeMatrix")
                  U@Dim <- d
              }
              U@Dimnames <- c(list(NULL), dn[2L])
              U@x <- .mkU(x@x, m, n)

              list(P1. = P1., L = L, U = U)
          })

## returning list(L, U, P), where A = P L U
## MJ: for backwards compatibility
setMethod("expand", signature(x = "denseLU"),
          function(x, ...) {
              r <- expand2(x)[c(2L, 3L, 1L)]
              names(r) <- c("L", "U", "P")
              NN <- list(NULL, NULL)
              for(i in 1:3)
                  r[[i]]@Dimnames <- NN
              r
          })


## METHODS FOR CLASS: sparseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("expand1", signature(x = "sparseLU"),
          function(x, which, ...) {
              switch(which,
                     "P1" =, "P1." = {
                         r <- new("pMatrix")
                         r@Dim <- x@Dim
                         r@perm <- x@p + 1L
                         if(which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "P2" =, "P2." = {
                         r <- new("pMatrix")
                         r@Dim <- d <- x@Dim
                         r@perm <- if(length(x@q)) x@q + 1L else seq_len(d[1L])
                         if(which == "P2")
                             r@margin <- 2L
                         r
                     },
                     "L" = x@L,
                     "U" = x@U,
                     stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%2$s2\", \"%2$s2.\", \"%3$s\", or \"%4$s\"",
                                   "which", "P", "L", "U"),
                          domain = NA))
          })

## returning list(P1', L, U, P2'), where A = P1' L U P2'
setMethod("expand2", signature(x = "sparseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              p1 <- x@p
              p2 <- x@q

              P1. <- new("pMatrix")
              P1.@Dim <- d
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@perm <- invertPerm(p1, 0L, 1L)

              P2. <- new("pMatrix")
              P2.@Dim <- d
              P2.@Dimnames <- c(list(NULL), dn[2L])
              P2.@margin <- 2L
              P2.@perm <-
                  if(length(p2))
                      invertPerm(p2, 0L, 1L)
                  else seq_len(d[1L])

              list(P1. = P1., L = x@L, U = x@U, P2. = P2.)
          })

## returning list(P, L, U, Q), where A = P' L U Q
## MJ: for backwards compatibility
setMethod("expand", signature(x = "sparseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              p <- x@p + 1L
              q <- x@q + 1L
              if(length(p) && !is.null(rn <- dn[[1L]]))
                  dn[[1L]] <- rn[p]
              if(length(q) && !is.null(cn <- dn[[2L]]))
                  dn[[2L]] <- cn[q]
              P <- new("pMatrix")
              P@Dim <- d
              P@perm <- p
              Q <- new("pMatrix")
              Q@Dim <- d
              Q@perm <- if(length(q)) q else seq_len(d[1L])
              L <- x@L
              L@Dimnames <- c(dn[1L], list(NULL))
              U <- x@U
              U@Dimnames <- c(list(NULL), dn[2L])
              list(P = P, L = L, U = U, Q = Q)
          })
