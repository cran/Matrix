## METHODS FOR GENERIC: norm
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("norm", c(x = "ANY", type = "missing"),
          function(x, type, ...) norm(x, type = "O", ...))

setMethod("norm", c(x = "denseMatrix", type = "character"),
          function(x, type, ...) {
              if(identical(type, "2"))
                  return(base::norm(.M2m(x), type = "2"))
              x <- .M2kind(x, ",")
              switch(substr(.M.nonvirtual(x), 2L, 3L),
                     "ge" = .Call(dgeMatrix_norm, x, type),
                     "sy" = .Call(dsyMatrix_norm, x, type),
                     "sp" = .Call(dspMatrix_norm, x, type),
                     "tr" = .Call(dtrMatrix_norm, x, type),
                     "tp" = .Call(dtpMatrix_norm, x, type))
          })

setMethod("norm", c(x = "sparseMatrix", type = "character"),
          function(x, type, ...) {
              if(any(x@Dim == 0L))
                  return(0)
              switch(EXPR = type[1L],
                     "O" = , "o" = , "1" =
                         max(colSums(abs(x))),
                     "I" = , "i" =
                         max(rowSums(abs(x))),
                     "2" =
                         {
                             warning(gettextf("'%s' via sparse -> dense coercion",
                                              "norm"),
                                     domain = NA)
                             base::norm(.M2m(x), type = "2")
                         },
                     "M" = , "m" =
                         max(abs(x)),
                     "F" = , "f" = , "E" = , "e" =
                         {
                             if(.M.kind(x) == "z")
                                 x <- abs(x)
                             sqrt(sum(x * x))
                         },
                     stop(gettextf("invalid %s=\"%s\"", "type", type[1L]),
                          domain = NA))
          })

setMethod("norm", c(x = "diagonalMatrix", type = "character"),
          function(x, type, ...) {
              if((n <- x@Dim[1L]) == 0L)
                  return(0)
              if(nonunit <- x@diag == "N") {
                  y <- x@x
                  if(.M.kind(x) == "n" && anyNA(y))
                      y <- y | is.na(y)
              }
              switch(EXPR = type[1L],
                     "O" = , "o" = , "1" = ,
                     "I" = , "i" = ,
                     "2" = ,
                     "M" = , "m" =
                         if(nonunit) max(abs(y)) else 1,
                     "F" = , "f" = , "E" = , "e" =
                         if(nonunit) {
                             if(is.complex(y))
                                 y <- abs(y)
                             sqrt(sum(y * y))
                         } else sqrt(n),
                     stop(gettextf("invalid %s=\"%s\"", "type", type[1L]),
                          domain = NA))
          })

setMethod("norm", c(x = "indMatrix", type = "character"),
          function(x, type, ...) {
              d <- x@Dim
              if((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  return(0)
              switch(EXPR = type[1L],
                     "O" = , "o" = , "1" =
                         if(x@margin == 1L) max(tabulate(x@perm, n)) else 1,
                     "I" = , "i" =
                         if(x@margin == 1L) 1 else max(tabulate(x@perm, m)),
                     "2" =
                         sqrt(max(tabulate(x@perm, if(x@margin == 1L) n else m))),
                     "M" = , "m" =
                         1,
                     "F" = , "f" = , "E" = , "e" =
                         if(x@margin == 1L) sqrt(m) else sqrt(n),
                     stop(gettextf("invalid %s=\"%s\"", "type", type[1L]),
                          domain = NA))
          })

setMethod("norm", c(x = "pMatrix", type = "character"),
          function(x, type, ...) {
              if((n <- x@Dim[1L]) == 0L)
                  return(0)
              switch(EXPR = type[1L],
                     "O" = , "o" = , "1" = ,
                     "I" = , "i" = ,
                     "2" = ,
                     "M" = , "m" =
                         1,
                     "F" = , "f" = , "E" = , "e" =
                         sqrt(n),
                     stop(gettextf("invalid %s=\"%s\"", "type", type[1L]),
                          domain = NA))
          })


## METHODS FOR GENERIC: rcond
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("rcond", c(x = "ANY", norm = "missing"),
          function(x, norm, ...) rcond(x, norm = "O", ...))

setMethod("rcond", c(x = "denseMatrix", norm = "character"),
          function(x, norm, ...) {
              x <- .M2kind(x, ",")
              switch(substr(.M.nonvirtual(x, strict = TRUE), 2L, 3L),
                     "ge" =
                         {
                             d <- x@Dim
                             m <- d[1L]
                             n <- d[2L]
                             if(m == n) {
                                 trf <- lu(x, warnSing = FALSE)
                                 .Call(dgeMatrix_rcond, x, trf, norm)
                             } else {
                                 ## MJ: norm(A = P1' Q R P2') = norm(R) holds
                                 ##     in general only for norm == "2", but
                                 ##     La_rcond_type() disallows norm == "2"
                                 ##     ... FIXME ??
                                 if(m < n) {
                                     x <- t(x)
                                     n <- m
                                 }
                                 R <- triu(qr(x)[["qr"]][seq_len(n), , drop = FALSE])
                                 rcond(R, norm = norm, ...)
                             }
                         },
                     "sy" =
                         {
                             trf <- BunchKaufman(x, warnSing = FALSE)
                             .Call(dsyMatrix_rcond, x, trf, norm)
                         },
                     "sp" =
                         {
                             trf <- BunchKaufman(x, warnSing = FALSE)
                             .Call(dspMatrix_rcond, x, trf, norm)
                         },
                     "po" = ,
                     "or" = # corMatrix
                         {
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
                         },
                     "pp" = ,
                     "op" = # copMatrix
                         {
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
                         },
                     "tr" = .Call(dtrMatrix_rcond, x, norm),
                     "tp" = .Call(dtpMatrix_rcond, x, norm))
          })

setMethod("rcond", c(x = "sparseMatrix", norm = "character"),
          function(x, norm, useInv = FALSE, ...) {
              d <- x@Dim
              if((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  return(Inf)
              if(m == n) {
                  if(isS4(useInv) || useInv) {
                      if(!isS4(useInv))
                          useInv <- solve(x)
                      1 / (norm(x, type = norm) * norm(useInv, type = norm))
                  } else {
                      warning(gettextf("'%s' via sparse -> dense coercion",
                                       "rcond"),
                              domain = NA)
                      rcond(.M2unpacked(x), norm = norm, ...)
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

setMethod("rcond", c(x = "diagonalMatrix", norm = "character"),
          function(x, norm, ...) {
              if((n <- x@Dim[1L]) == 0L)
                  return(Inf)
              if(nonunit <- x@diag == "N") {
                  y <- x@x
                  if(.M.kind(x) == "n" && anyNA(y))
                      y <- y | is.na(y)
              }
              switch(EXPR = norm[1L],
                     "O" = , "o" = , "1" = ,
                     "I" = , "i" = ,
                     "2" = ,
                     "M" = , "m" =
                         if(nonunit) {
                             ry <- range(abs(y))
                             ry[1L] / ry[2L]
                         } else 1,
                     "F" = , "f" = , "E" = , "e" =
                         if(nonunit) {
                             if(is.complex(y))
                                 y <- abs(y)
                             yy <- y * y
                             1 / sqrt(sum(yy) * sum(1 / yy))
                         } else 1 / n,
                     stop(gettext("invalid %s=\"%s\"", "norm", norm[1L]),
                          domain = NA))
          })

setMethod("rcond", c(x = "indMatrix", norm = "character"),
          function(x, norm, ...) {
              d <- x@Dim
              if((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  return(Inf)
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
                         stop(gettext("invalid %s=\"%s\"", "norm", norm[1L]),
                              domain = NA))
              } else {
                  if(m < n) {
                      x <- t(x)
                      n <- m
                  }
                  R <- triu(qr(x)@R[seq_len(n), , drop = FALSE])
                  rcond(R, norm = norm, ...)
              }
          })

setMethod("rcond", c(x = "pMatrix", norm = "character"),
          function(x, norm, ...) {
              if((n <- x@Dim[1L]) == 0L)
                  return(Inf)
              switch(EXPR = norm[1L],
                     "O" = , "o" = , "1" = ,
                     "I" = , "i" = ,
                     "2" = ,
                     "M" = , "m" =
                         1,
                     "F" = , "f" = , "E" = , "e" =
                         1 / n,
                     stop(gettext("invalid %s=\"%s\"", "norm", norm[1L]),
                          domain = NA))
          })
