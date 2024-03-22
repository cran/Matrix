## METHODS FOR GENERIC: qr
## pivoted QR factorization of dense and sparse matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: Well, I'd very much like to _not_ truncate by default ...
##     not least because qr.qy and qr.qty should remain inverses ...
.qr.rank.def.truncating <- TRUE

.qr.rank.def.warn <- function(qr) {
    if(m0 <- qr@V@Dim[1L] - qr@Dim[1L])
        warning(gettextf("matrix is structurally rank deficient; using augmented matrix with additional %d row(s) of zeros",
                         m0),
                domain = NA)
    m0
}

setMethod("qr", c(x = "sparseMatrix"),
          function(x, ...)
              qr(.M2gen(.M2C(x), ","), ...))

setMethod("qr", c(x = "dgCMatrix"),
          function(x, order = 3L, ...) {
              r <- .Call(dgCMatrix_orf, x, order, TRUE)
              .qr.rank.def.warn(r)
              r
          })


## METHODS FOR CLASS: sparseQR
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("expand1", c(x = "sparseQR"),
          function(x, which, ...) {
              .qr.rank.def.warn(x)
              R <- x@R
              d <- R@Dim
              m <- d[1L]
              n <- d[2L]
              switch(which,
                     "P1" =, "P1." = {
                         r <- new("pMatrix")
                         r@Dim <- c(m, m)
                         r@perm <- x@p + 1L
                         if(which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "P2" =, "P2." = {
                         r <- new("pMatrix")
                         r@Dim <- c(n, n)
                         r@perm <- if(length(x@q)) x@q + 1L else seq_len(n)
                         if(which == "P2")
                             r@margin <- 2L
                         r
                     },
                     "Q"  = .Call(sparseQR_matmult, x, NULL, 6L,  TRUE, NULL),
                     "Q1" = .Call(sparseQR_matmult, x, NULL, 6L, FALSE, NULL),
                     "R"  = R,
                     "R1" = triu(if(m == n) R else R[seq_len(n), , drop = FALSE]),
                     stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%2$s2\", \"%2$s2.\", \"%3$s\", \"%3$s1\", \"%4$s\", or \"%4$s1\"",
                                   "which", "P", "Q", "R"),
                          domain = NA))
          })

## returning list(P1', Q, R, P2'), where A = P1' Q R P2'
setMethod("expand2", c(x = "sparseQR"),
          function(x, complete = FALSE, ...) {
              m0 <- .qr.rank.def.warn(x)
              R <- x@R
              d <- R@Dim
              m <- d[1L]
              n <- d[2L]
              dn <- x@Dimnames
              if(m0 && !is.null(dn[[1L]]))
                  length(dn[[1L]]) <- m
              Q <- .Call(sparseQR_matmult, x, NULL, 6L, complete, NULL)
              if(!complete && n < m)
                  R <- R[seq_len(n), , drop = FALSE]
              p1 <- x@p
              p2 <- x@q
              P1. <- new("pMatrix",
                         Dim = c(m, m),
                         Dimnames = c(dn[1L], list(NULL)),
                         margin = 1L,
                         perm = invertPerm(p1, 0L, 1L))
              P2. <- new("pMatrix",
                         Dim = c(n, n),
                         Dimnames = c(list(NULL), dn[2L]),
                         margin = 2L,
                         perm = if(length(p2)) invertPerm(p2, 0L, 1L) else seq_len(n))
              if(complete)
                  list(P1. = P1., Q = Q, R = R, P2. = P2.)
              else
                  list(P1. = P1., Q1 = Q, R1 = triu(R), P2. = P2.)
          })

setMethod("qr.Q", c(qr = "sparseQR"),
          function(qr, complete = FALSE, Dvec) {
              m0 <- .qr.rank.def.warn(qr)
              if(missing(Dvec))
                  Dvec <- NULL
              else {
                  storage.mode(Dvec) <- "double"
                  if(length(Dvec) != qr@V@Dim[if(complete) 1L else 2L])
                      stop(gettextf("'%s' has the wrong length", "Dvec"),
                           domain = NA)
              }
              Q <- .Call(sparseQR_matmult, qr, NULL, 4L, complete, Dvec)
              dn <- c(qr@Dimnames[1L], list(NULL))
              if(!is.null(rn <- dn[[1L]])) {
                  if(m0)
                      length(rn) <- length(rn) + m0
                  if(is.unsorted(p1 <- qr@p, strictly = TRUE))
                      rn <- rn[p1 + 1L]
                  dn[[1L]] <- rn
              }
              Q@Dimnames <- dn
              if(m0 && .qr.rank.def.truncating) {
                  i <- seq_len(Q@Dim[1L] - m0)
                  Q <- if(complete)
                           Q[i, i, drop = FALSE]
                       else Q[i, , drop = FALSE]
              }
              Q
          })

qrR <- function(qr, complete = FALSE, backPermute = TRUE, row.names = TRUE) {
    m0 <- .qr.rank.def.warn(qr)
    R <- qr@R
    d <- R@Dim
    m <- d[1L]
    n <- d[2L]
    dn <- qr@Dimnames
    p2 <- qr@q + 1L
    p2.uns <- is.unsorted(p2, strictly = TRUE) # FALSE if length is 0
    if(!row.names)
        dn <- c(list(NULL), dn[2L])
    else if(m0 && !is.null(rn <- dn[[1L]]))
        length(dn[[1L]]) <- length(rn) + m0
    if(p2.uns && !is.null(cn <- dn[[2L]]))
        dn[[2L]] <- cn[p2]
    R@Dimnames <- dn
    R <-
        if(!complete && n < m) {
            if(backPermute && p2.uns)
                R[seq_len(n), invertPerm(p2), drop = FALSE]
            else R[seq_len(n), , drop = FALSE]
        } else {
            if(backPermute && p2.uns)
                R[, invertPerm(p2), drop = FALSE]
            else R
        }
    if(m0 && .qr.rank.def.truncating && complete)
        R <- R[seq_len(m - m0), , drop = FALSE]
    if(complete || backPermute) R else triu(R)
}

setMethod("qr.R", c(qr = "sparseQR"),
          function(qr, complete = FALSE, backPermute = FALSE, ...)
              qrR(qr, complete = complete, backPermute = backPermute,
                  row.names = FALSE))

## https://stat.ethz.ch/pipermail/r-devel/2023-June/082649.html
setMethod("qr.X", c(qr = "sparseQR"),
          function(qr, complete = FALSE, ncol) {
              m0 <- .qr.rank.def.warn(qr)
              R <- qr@R
              d <- R@Dim
              m <- d[1L]
              n <- d[2L]
              if(missing(ncol))
                  ncol <- if(complete) m else min(m, n)
              else {
                  ncol <- as.integer(ncol)
                  if(ncol < 0L || ncol > m)
                      stop(gettextf("invalid '%s': not in %d:%d",
                                    "ncol", 0L, m),
                           domain = NA)
              }
              p2 <- qr@q + 1L
              p2.uns <- is.unsorted(p2, strictly = TRUE) # FALSE if length is 0
              if(p2.uns && ncol < n)
                  stop(gettextf("need greater '%s' as pivoting occurred",
                                "ncol"),
                       domain = NA)
              else if(ncol < n)
                  R <- R[, seq_len(ncol), drop = FALSE]
              else if(ncol > n) {
                  Rp <- R@p
                  Ri <- R@i
                  Rx <- R@x
                  Rnnz <- Rp[length(Rp)]
                  R@Dim[2L] <- ncol
                  R@p <- c(Rp, Rnnz + seq_len(ncol - n))
                  R@i <- c(if(length(Ri) == Rnnz) Ri else Ri[seq_len(Rnnz)],
                           n:(ncol - 1L))
                  R@x <- c(if(length(Rx) == Rnnz) Rx else Rx[seq_len(Rnnz)],
                           rep.int(1, ncol - n))
              }
              r <- .Call(sparseQR_matmult, qr, .sparse2dense(R), 4L, NA, NULL)
              if(p2.uns) {
                  j <- invertPerm(p2)
                  if(ncol > n)
                      j <- c(j, (n + 1L):ncol)
                  r <- r[, j, drop = FALSE]
              }
              dn <- qr@Dimnames
              if(m0 && !is.null(rn <- dn[[1L]]))
                  length(dn[[1L]]) <- length(rn) + m0
              if(!is.null(cn <- dn[[2L]]) && length(cn) != ncol)
                  length(dn[[2L]]) <- ncol
              r@Dimnames <- dn
              if(m0 && .qr.rank.def.truncating) {
                  i <- seq_len(r@Dim[1L] - m0)
                  r <- if(ncol > length(j))
                           r[i, i, drop = FALSE]
                       else r[i, , drop = FALSE]
              }
              r
          })

.qr.y0 <- function(y, m0) {
    d <- y@Dim
    d[1L] <- (m <- d[1L]) + m0
    dn <- y@Dimnames
    if(!is.null(dn[[1L]]))
        length(dn[[1L]]) <- d[1L]
    y0 <- new("dgeMatrix")
    y0@Dim <- d
    y0@Dimnames <- dn
    y0@x <- as.double(`[<-`(array(0, d), seq_len(m), , y@x))
    y0
}

setMethod("qr.coef", c(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y) {
              if(m0 <- .qr.rank.def.warn(qr))
                  y <- .qr.y0(y, m0)
              r <- .Call(sparseQR_matmult, qr, y, 0L, NA, NULL)
              r@Dimnames <- c(qr@Dimnames[2L], y@Dimnames[2L])
              r
          })

setMethod("qr.coef", c(qr = "sparseQR", y = "vector"),
          function(qr, y)
              drop(qr.coef(qr, .m2dense(y, ",ge"))))

setMethod("qr.coef", c(qr = "sparseQR", y = "matrix"),
          function(qr, y)
              qr.coef(qr, .m2dense(y, ",ge")))

setMethod("qr.coef", c(qr = "sparseQR", y = "Matrix"),
          function(qr, y)
              qr.coef(qr, .m2dense(.M2m(y), ",ge")))

setMethod("qr.fitted", c(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y, k = qr$rank) {
              if(m0 <- .qr.rank.def.warn(qr))
                  y <- .qr.y0(y, m0)
              r <- .Call(sparseQR_matmult, qr, y, 1L, NA, NULL)
              r@Dimnames <- y@Dimnames
              if(m0 && .qr.rank.def.truncating)
                  r <- r[seq_len(r@Dim[1L] - m0), , drop = FALSE]
              r
          })

setMethod("qr.fitted", c(qr = "sparseQR", y = "vector"),
          function(qr, y, k = qr$rank)
              drop(qr.fitted(qr, .m2dense(y, ",ge"))))

setMethod("qr.fitted", c(qr = "sparseQR", y = "matrix"),
          function(qr, y, k = qr$rank)
              qr.fitted(qr, .m2dense(y, ",ge")))

setMethod("qr.fitted", c(qr = "sparseQR", y = "Matrix"),
          function(qr, y, k = qr$rank)
              qr.fitted(qr, .m2dense(.M2m(y), ",ge")))

setMethod("qr.resid", c(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y) {
              if(m0 <- .qr.rank.def.warn(qr))
                  y <- .qr.y0(y, m0)
              r <- .Call(sparseQR_matmult, qr, y, 2L, NA, NULL)
              r@Dimnames <- y@Dimnames
              if(m0 && .qr.rank.def.truncating)
                  r <- r[seq_len(r@Dim[1L] - m0), , drop = FALSE]
              r
          })

setMethod("qr.resid", c(qr = "sparseQR", y = "vector"),
          function(qr, y)
              drop(qr.resid(qr, .m2dense(y, ",ge"))))

setMethod("qr.resid", c(qr = "sparseQR", y = "matrix"),
          function(qr, y)
              qr.resid(qr, .m2dense(y, ",ge")))

setMethod("qr.resid", c(qr = "sparseQR", y = "Matrix"),
          function(qr, y)
              qr.resid(qr, .m2dense(.M2m(y), ",ge")))

setMethod("qr.qty", c(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y) {
              if(m0 <- .qr.rank.def.warn(qr))
                  y <- .qr.y0(y, m0)
              r <- .Call(sparseQR_matmult, qr, y, 3L, NA, NULL)
              r@Dimnames <- c(list(NULL), y@Dimnames[2L])
              if(m0 && .qr.rank.def.truncating)
                  r <- r[seq_len(r@Dim[1L] - m0), , drop = FALSE]
              r
          })

setMethod("qr.qty", c(qr = "sparseQR", y = "vector"),
          function(qr, y)
              drop(qr.qty(qr, .m2dense(y, ",ge"))))

setMethod("qr.qty", c(qr = "sparseQR", y = "matrix"),
          function(qr, y)
              qr.qty(qr, .m2dense(y, ",ge")))

setMethod("qr.qty", c(qr = "sparseQR", y = "Matrix"),
          function(qr, y)
              qr.qty(qr, .m2dense(.M2m(y), ",ge")))

setMethod("qr.qy", c(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y) {
              if(m0 <- .qr.rank.def.warn(qr))
                  y <- .qr.y0(y, m0)
              r <- .Call(sparseQR_matmult, qr, y, 4L, NA, NULL)
              dn <- c(qr@Dimnames[1L], y@Dimnames[2L])
              if(!is.null(rn <- dn[[1L]])) {
                  if(m0)
                      length(rn) <- length(rn) + m0
                  if(is.unsorted(p1 <- qr@p, strictly = TRUE))
                      rn <- rn[p1 + 1L]
                  dn[[1L]] <- rn
              }
              r@Dimnames <- dn
              if(m0 && .qr.rank.def.truncating)
                  r <- r[seq_len(r@Dim[1L] - m0), , drop = FALSE]
              r
          })

setMethod("qr.qy", c(qr = "sparseQR", y = "vector"),
          function(qr, y)
              drop(qr.qy(qr, .m2dense(y, ",ge"))))

setMethod("qr.qy", c(qr = "sparseQR", y = "matrix"),
          function(qr, y)
              qr.qy(qr, .m2dense(y, ",ge")))

setMethod("qr.qy", c(qr = "sparseQR", y = "Matrix"),
          function(qr, y)
              qr.qy(qr, .m2dense(.M2m(y), ",ge")))
