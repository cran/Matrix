## METHODS FOR GENERIC: rcond
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("rcond", signature(x = "ANY", norm = "missing"),
          function(x, norm, ...)
              rcond(x, norm = "O", ...))

setMethod("rcond", signature(x = "sparseMatrix", norm = "character"),
          function(x, norm, useInv = FALSE, ...) {
              d <- x@Dim
              if((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  stop("rcond(x) is undefined: 'x' has length 0")
              if(m == n) {
                  if(isS4(useInv) || useInv) {
                      if(!isS4(useInv))
                          useInv <- solve(x)
                      1 / (norm(x, type = norm) * norm(useInv, type = norm))
                  } else {
                      warning("'rcond' via sparse -> dense coercion")
                      rcond(as(x, "denseMatrix"), norm = norm, ...)
                  }
              } else {
                  ## MJ: norm(A = P1' Q R P2') = norm(R) holds in general
                  ##     only for norm == "2", but La_rcond_type() disallows
                  ##     norm == "2" ... FIXME ??
                  if(m < n) {
                      x <- t(x)
                      n <- m
                  }
                  R <- triu(qr(x)@R[seq_len(n), , drop = FALSE])
                  rcond(R, norm = norm, ...)
              }
           })

setMethod("rcond", signature(x = "diagonalMatrix", norm = "character"),
          function(x, norm, ...) {
              if((n <- x@Dim[1L]) == 0L)
                  stop("rcond(x) is undefined: 'x' has length 0")
              switch(EXPR = norm[1L],
                     "O" = , "o" = , "1" = ,
                     "I" = , "i" = ,
                     "2" = ,
                     "M" = , "m" =
                         if(x@diag == "N") {
                             rx <- range(abs(x@x))
                             rx[1L] / rx[2L]
                         } else 1,
                     "F" = , "f" = , "E" = , "e" =
                         if(x@diag == "N") {
                             xx <- x@x * x@x
                             1 / sqrt(sum(xx) * sum(1 / xx))
                         } else 1 / n,
                     stop("invalid 'norm'"))
          })

setMethod("rcond", signature(x = "indMatrix", norm = "character"),
          function(x, norm, ...) {
              d <- x@Dim
              if((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  stop("rcond(x) is undefined: 'x' has length 0")
              if (m == n) {
                  if(anyDuplicated.default(x@perm))
                      return(0)
                  switch(EXPR = norm[1L],
                         "O" = , "o" = , "1" = ,
                         "I" = , "i" = ,
                         "2" = ,
                         "M" = , "m" =
                             1,
                         "F" = , "f" = , "E" = , "e" =
                             1 / n,
                         stop("invalid 'norm'"))
              } else {
                  if(m < n) {
                      x <- t(x)
                      n <- m
                  }
                  R <- triu(qr(x)@R[seq_len(n), , drop = FALSE])
                  rcond(R, norm = norm, ...)
              }
          })

setMethod("rcond", signature(x = "pMatrix", norm = "character"),
          function(x, norm, ...) {
              if((n <- x@Dim[1L]) == 0L)
                  stop("rcond(x) is undefined: 'x' has length 0")
              switch(EXPR = norm[1L],
                     "O" = , "o" = , "1" = ,
                     "I" = , "i" = ,
                     "2" = ,
                     "M" = , "m" =
                         1,
                     "F" = , "f" = , "E" = , "e" =
                         1 / n,
                     stop("invalid 'norm'"))
          })

setMethod("rcond", signature(x = "denseMatrix", norm = "character"),
          function(x, norm, ...)
              rcond(..dense2d(x), norm = norm, ...))

setMethod("rcond", signature(x = "dgeMatrix", norm = "character"),
          function(x, norm, ...) {
              d <- x@Dim
              m <- d[1L]
              n <- d[2L]
              if(m == n) {
                  trf <- lu(x, warnSing = FALSE)
                  .Call(dgeMatrix_rcond, x, trf, norm)
              } else {
                  ## MJ: norm(A = P1' Q R P2') = norm(R) holds in general
                  ##     only for norm == "2", but La_rcond_type() disallows
                  ##     norm == "2" ... FIXME ??
                  if(m < n) {
                      x <- t(x)
                      n <- m
                  }
                  R <- triu(qr(x)[["qr"]][seq_len(n), , drop = FALSE])
                  rcond(R, norm = norm, ...)
              }
          })

setMethod("rcond", signature(x = "dtrMatrix", norm = "character"),
          function(x, norm, ...)
              .Call(dtrMatrix_rcond, x, norm))

setMethod("rcond", signature(x = "dtpMatrix", norm = "character"),
          function(x, norm, ...)
              .Call(dtpMatrix_rcond, x, norm))

setMethod("rcond", signature(x = "dsyMatrix", norm = "character"),
          function(x, norm, ...) {
              trf <- BunchKaufman(x, warnSing = FALSE)
              .Call(dsyMatrix_rcond, x, trf, norm)
          })

setMethod("rcond", signature(x = "dspMatrix", norm = "character"),
          function(x, norm, ...) {
              trf <- BunchKaufman(x, warnSing = FALSE)
              .Call(dspMatrix_rcond, x, trf, norm)
          })

setMethod("rcond", signature(x = "dpoMatrix", norm = "character"),
          function(x, norm, ...) {
              ok <- TRUE
              trf <- tryCatch(
                  Cholesky(x, perm = FALSE),
                  error = function(e) {
                      ok <<- FALSE
                      BunchKaufman(x, warnSing = FALSE)
                  })
              if(ok)
                  .Call(dpoMatrix_rcond, x, trf, norm)
              else .Call(dsyMatrix_rcond, x, trf, norm)
          })

setMethod("rcond", signature(x = "dppMatrix", norm = "character"),
          function(x, norm, ...) {
              ok <- TRUE
              trf <- tryCatch(
                  Cholesky(x, perm = FALSE),
                  error = function(e) {
                      ok <<- FALSE
                      BunchKaufman(x, warnSing = FALSE)
                  })
              if(ok)
                  .Call(dppMatrix_rcond, x, trf, norm)
              else .Call(dspMatrix_rcond, x, trf, norm)
          })
