## METHODS FOR GENERIC: Schur
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Schur", c(x = "dgeMatrix"),
          function(x, vectors = TRUE, ...) {
              if(length(x.x <- x@x) && !all(is.finite(range(x.x))))
                  stop(gettextf("'%s' has non-finite values", "x"), domain = NA)
              cl <- .Call(dgeMatrix_sch, x, vectors, TRUE)
              if(all(cl$WI == 0)) {
                  vals <- cl$WR
                  T <- triu(cl$T)
              } else {
                  vals <- complex(real = cl$WR, imaginary = cl$WI)
                  T <- .m2dense(cl$T, ",ge")
              }
              if(vectors)
                  new("Schur", Dim = x@Dim, Dimnames = x@Dimnames,
                      Q = .m2dense(cl$Z, ",ge"), T = T, EValues = vals)
              else list(T = T, EValues = vals)
          })

setMethod("Schur", c(x = "dsyMatrix"),
          function(x, vectors = TRUE, ...) {
              e <- eigen(x, symmetric = TRUE, only.values = !vectors)
              vals <- as.double(e$values)
              T <- new("ddiMatrix", Dim = x@Dim, x = vals)
              if(vectors)
                  new("Schur", Dim = x@Dim, Dimnames = symDN(x@Dimnames),
                      Q = .m2dense(e$vectors, ",ge"), T = T, EValues = vals)
              else list(T = T, EValues = vals)
          })

setMethod("Schur", c(x = "matrix"),
          function(x, vectors = TRUE, ...) {
              ## FIXME: wrong for complex, but package 'control' seems to
              ##        rely on the complex->double coercion (!?)
              storage.mode(x) <- "double"
              if(length(x) && !all(is.finite(range(x))))
                  stop(gettextf("'%s' has non-finite values", "x"), domain = NA)
              cl <- .Call(dgeMatrix_sch, x, vectors, FALSE)
              vals <-
                  if(all(cl$WI == 0))
                      cl$WR
                  else complex(real = cl$WR, imaginary = cl$WI)
              if(vectors)
                  list(Q = cl$Z, T = cl$T, EValues = vals)
              else list(T = cl$T, EValues = vals)
          })

## FIXME: don't coerce from sparse to dense
setMethod("Schur", c(x = "generalMatrix"),
          function(x, vectors = TRUE, ...)
              Schur(.M2unpacked(.M2kind(x, ",")), vectors, ...))

## FIXME: don't coerce from sparse to dense
setMethod("Schur", c(x = "symmetricMatrix"),
          function(x, vectors = TRUE, ...)
              Schur(.M2unpacked(.M2kind(x, ",")), vectors, ...))

setMethod("Schur", c(x = "triangularMatrix"),
          function(x, vectors = TRUE, ...) {
              x <- .M2kind(x, ",")
              n <- (d <- x@Dim)[1L]
              if(n == 0L)
                  x@uplo <- "U"
              else if(.M.kind(x) != "n" && !all(is.finite(range(x))))
                  stop(gettextf("'%s' has non-finite values", "x"), domain = NA)
              vals <- diag(x, names = FALSE)
              if(x@uplo == "U") {
                  if(vectors) {
                      Q <- new("ddiMatrix", Dim = d, diag = "U")
                      new("Schur", Dim = d, Dimnames = x@Dimnames,
                          Q = Q, T = x, EValues = vals)
                  } else list(T = x, EValues = vals)
              } else {
                  perm <- n:1L
                  vals <- vals[perm]
                  T <- triu(x[perm, perm, drop = FALSE])
                  if(vectors) {
                      Q <- new("pMatrix", Dim = d, perm = perm)
                      new("Schur", Dim = d, Dimnames = x@Dimnames,
                          Q = Q, T = T, EValues = vals)
                  } else list(T = x, EValues = vals)
              }
          })

setMethod("Schur", c(x = "diagonalMatrix"),
          function(x, vectors = TRUE, ...) {
              x <- .M2kind(x, ",")
              d <- x@Dim
              if(x@diag != "N") {
                  vals <- rep.int(1, d[1L])
                  T <- new("ddiMatrix", Dim = d, diag = "U")
              } else {
                  vals <- x@x
                  if(length(vals) && !all(is.finite(range(vals))))
                      stop(gettextf("'%s' has non-finite values", "x"), domain = NA)
                  T <- new("ddiMatrix", Dim = d, x = vals)
              }
              if(vectors) {
                  Q <- new("ddiMatrix", Dim = d, diag = "U")
                  new("Schur", Dim = d, Dimnames = x@Dimnames,
                      Q = Q, T = T, EValues = vals)
              } else list(T = T, EValues = vals)
          })


## METHODS FOR CLASS: Schur
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("expand1", c(x = "Schur"),
          function(x, which, ...)
              switch(which, "Q" = x@Q, "T" = x@T, "Q." = t(x@Q),
                     stop(gettextf("'%1$s' is not \"%2$s\", \"%3$s\", or \"%2$s.\"",
                                   "which", "Q", "T"),
                          domain = NA)))

setMethod("expand2", c(x = "Schur"),
          function(x, ...) {
              Q  <- x@Q
              Q. <- t(Q)
              dn <- x@Dimnames
              Q @Dimnames <- c(dn[1L], list(NULL))
              Q.@Dimnames <- c(list(NULL), dn[2L])
              list(Q = Q, T = x@T, Q. = Q.)
          })
