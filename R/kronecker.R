## METHODS FOR GENERIC: kronecker
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.kroneckerDimnames <- function(dnX, dX = lengths(dnX, FALSE),
                               dnY, dY = lengths(dnX, FALSE),
                               sep = ":") {
    dnr <- list(NULL, NULL)
    if(identical(dnX, dnr) && identical(dnY, dnr))
        return(NULL)
    for(i in 1:2) {
        if(have.sX <- (nY <- dY[i]) && !is.null(sX <- dnX[[i]]))
            sX <- rep    (sX,  each = nY)
        if(have.sY <- (nX <- dX[i]) && !is.null(sY <- dnY[[i]]))
            sY <- rep.int(sY, times = nX)
        dnr[[i]] <-
            if(have.sX && have.sY)
                paste0(sX, sep, sY)
            else if(have.sX)
                paste0(sX, sep    )
            else if(have.sY)
                paste0(    sep, sY)
            else if(nX && nY)
                rep.int(sep, nX * nY)
    }
    dnr
}

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "diagonalMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              r <- new("ddiMatrix")
              r@Dim <- dX * dY
              if((uX <- X@diag != "N") & (uY <- Y@diag != "N"))
                  r@diag <- "U"
              else if(uX)
                  r@x <- rep.int(as.double(Y@x), dX[1L])
              else if(uY)
                  r@x <- rep(as.double(X@x), each = dY[1L])
              else r@x <- rep(X@x, each = dY[1L]) * Y@x
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY)))
                  r@Dimnames <- dnr
              r
          })

if(FALSE) { # --NOT YET--
setMethod("kronecker", signature(X = "diagonalMatrix", Y = "denseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              uX <- X@diag != "N"
              uY <- FALSE
              shape <- .M.shape(Y)
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  uplo <- r@uplo <- Y@uplo
                  if(shape == "t")
                      uY <- Y@diag != "N"
              }
              if(!all(dr)) {
                  r@p <- integer(dr[2L] + 1)
                  if(uX && uY)
                      r@diag <- "U"
              } else {
                  m <- dY[1L]
                  nX <- dX[1L]
                  nY <- length(y <- Y@x)
                  if(shape != "g" && nY > 1L && nY == prod(dY)) {
                      Y <- pack(Y)
                      nY <- length(y <- Y@x)
                  }
                  if(as.double(nX) * nY > .Machine$integer.max)
                      stop("number of nonzero entries cannot exceed 2^31-1")
                  if(!uX && uY) {
                      diag(Y) <- TRUE
                      nY <- length(y <- Y@x)
                  }
                  if(is.logical(y) && .M.kind(Y) == "n")
                      y <- y | is.na(y)
                  if(shape == "g") {
                      r@p <- seq.int(0L, by = m, length.out = dr[2L] + 1)
                      r@i <-
                          rep(seq.int(0L, by = m, length.out = nX),
                              each = nY) +
                          0:(m-1L)
                  } else if(uplo == "U") {
                      r@p <- c(0L, cumsum(rep.int(1:m, nX)))
                      r@i <-
                          rep(seq.int(0L, by = m, length.out = nX),
                              each = nY) +
                          sequence.default(nvec = 1:m, from = 0L)
                  } else {
                      r@p <- c(0L, cumsum(rep.int(m:1, nX)))
                      r@i <-
                          rep(seq.int(0L, by = m, length.out = nX),
                              each = nY) +
                          sequence.default(nvec = m:1, from = 0:(m-1L))
                  }
                  r@x <-
                      if(uX)
                          rep.int(as.double(y), nX)
                      else as.double(y) * rep(X@x, each = nY)
                  if(uX && uY)
                      r <- ..diagN2U(r, sparse = TRUE)
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "denseMatrix", Y = "diagonalMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              shape <- .M.shape(X)
              uX <- FALSE
              uY <- Y@diag != "N"
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  uplo <- r@uplo <- X@uplo
                  if(shape == "t")
                      uX <- X@diag != "N"
              }
              if(!all(dr)) {
                  r@p <- integer(dr[2L] + 1)
                  if(uX && uY)
                      r@diag <- "U"
              } else {
                  m <- dX[1L]
                  n <- dX[2L]
                  nX <- length(x <- X@x)
                  nY <- dY[1L]
                  if(shape != "g" && nX > 1L && nX == prod(dX)) {
                      X <- pack(X)
                      nX <- length(x <- X@x)
                  }
                  if(as.double(nX) * nY > .Machine$integer.max)
                      stop("number of nonzero entries cannot exceed 2^31-1")
                  if(uX && !uY) {
                      diag(X) <- TRUE
                      nX <- length(x <- X@x)
                  }
                  if(is.logical(x) && .M.kind(X) == "n")
                      x <- x | is.na(x)
                  if(shape == "g") {
                      x. <- function() {
                          j <- rep(1:n, each = nY)
                          as.vector(`dim<-`(as.double(x), dX)[, j])
                      }
                      r@p <- seq.int(0L, by = m, length.out = dr[2L] + 1)
                      r@i <- rep.int(sequence.default(nvec = rep.int(m, nY),
                                                      from = 0:(nY-1L),
                                                      by = nY),
                                     n)
                      r@x <-
                          if(uY)
                              x.()
                          else x.() * rep(Y@x, each = m)
                  } else if(uplo == "U") {
                      rep.1.n <- rep(1:n, each = nY)
                      s <- sequence.default(
                          nvec = rep.1.n,
                          from = rep(1L + cumsum(0:(n-1L)), each = nY),
                          by = 1L)

                      r@p <- c(0L, cumsum(rep.1.n))
                      r@i <- sequence.default(nvec = rep.1.n,
                                              from = rep.int(0:(nY-1L), n),
                                              by = nY)
                      r@x <-
                          if(uY)
                              as.double(x)[s]
                          else
                              as.double(x)[s] *
                                  rep.int(rep.int(Y@x, n), rep.1.n)
                  } else {
                      rep.n.1 <- rep(n:1, each = nY)
                      s <- sequence.default(
                          nvec = rep.n.1,
                          from = rep(1L + cumsum(c(0L, if(n > 1L) n:2)),
                                     each = nY),
                          by = 1L)
                      r@p <- c(0L, cumsum(rep.n.1))
                      r@i <- sequence.default(nvec = rep.n.1,
                                              from = 0:(n*nY-1L),
                                              by = nY)
                      r@x <-
                          if(uY)
                              as.double(x)[s]
                          else
                              as.double(x)[s] *
                                  rep.int(rep.int(Y@x, n), rep.n.1)
                  }
                  if(uX && uY)
                      r <- ..diagN2U(r, sparse = TRUE)
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })
} # --NOT YET--

setMethod("kronecker", signature(X = "denseMatrix", Y = "denseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              shape <- switch(.M.shape(X),
                              g = "g",
                              t = if(.M.shape(Y) == "t" && X@uplo == Y@uplo)
                                      "t"
                                  else "g",
                              s = if(.M.shape(Y) == "s")
                                      "s"
                                  else "g")
              r <- new(switch(shape,
                              g = "dgeMatrix",
                              t = "dtrMatrix",
                              s = "dsyMatrix"))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- X@uplo
                  if(shape == "t" && X@diag != "N" && Y@diag != "N")
                      r@diag <- "U"
              }
              if(all(dr)) {
                  X <- .dense2m(X)
                  Y <- .dense2m(Y)
                  storage.mode(X) <- "double"
                  storage.mode(Y) <- "double"
                  r@x <-
                      as.vector(X[rep(seq_len(dX[1L]), each = dY[1L]),
                                  rep(seq_len(dX[2L]), each = dY[2L])]) *
                      as.vector(Y[rep.int(seq_len(dY[1L]), dX[1L]), ])
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "CsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              uX <- X@diag != "N"
              uY <- FALSE
              shape <- .M.shape(Y)
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- Y@uplo
                  if(shape == "t" && (uX & (uY <- Y@diag != "N")))
                      r@diag <- "U"
              }
              if(!all(dr))
                  r@p <- integer(dr[2L] + 1)
              else {
                  if(!uX && uY)
                      Y <- ..diagU2N(Y)
                  if((nY <- (p <- Y@p)[length(p)]) == 0L)
                      r@p <- integer(dr[2L] + 1)
                  else if(as.double(nX <- dX[1L]) * nY > .Machine$integer.max)
                      stop("number of nonzero entries cannot exceed 2^31-1")
                  else {
                      head. <-
                          if(length(Y@i) > nY)
                              function(x) x[seq_len(nY)]
                          else identity
                      r@p <- c(0L, cumsum(rep.int(p[-1L] - p[-length(p)], nX)))
                      r@i <-
                          rep(seq.int(0L, by = dY[1L], length.out = nX),
                              each = nY) +
                          head.(Y@i)
                      r@x <-
                          if(uX) {
                              if(.M.kind(Y) != "n")
                                  rep.int(as.double(head.(Y@x)), nX)
                              else rep.int(1, nX * nY)
                          } else {
                              if(.M.kind(Y) != "n")
                                  rep(as.double(X@x), each = nY) *
                                      as.double(head.(Y@x))
                              else rep(as.double(X@x), each = nY)
                          }
                  }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "CsparseMatrix", Y = "diagonalMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              uX <- FALSE
              uY <- Y@diag != "N"
              shape <- .M.shape(X)
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- X@uplo
                  if(shape == "t" && (uX <- X@diag != "N") && uY)
                      r@diag <- "U"
              }
              if(!all(dr))
                  r@p <- integer(dr[2L] + 1)
              else {
                  if(uX && !uY)
                      X <- ..diagU2N(X)
                  if((nX <- (p <- X@p)[length(p)]) == 0L)
                      r@p <- integer(dr[2L] + 1)
                  else if(as.double(nY <- dY[1L]) * nX > .Machine$integer.max)
                      stop("number of nonzero entries cannot exceed 2^31-1")
                  else {
                      dp <- p[-1L] - p[-length(p)]
                      j. <- which(dp > 0L)
                      nj. <- length(j.)

                      rep.dp <- rep(dp[j.], each = nY)
                      s. <- sequence.default(
                          nvec = rep.dp,
                          from = rep(1L + p[j.], each = nY),
                          by = 1L)

                      r@p <- c(0L, cumsum(rep(dp, each = nY)))
                      r@i <- (nY * X@i)[s.] +
                          rep.int(rep.int(0:(nY-1L), nj.), rep.dp)
                      r@x <-
                          if(uY) {
                              if(.M.kind(X) != "n")
                                 as.double(X@x)[s.]
                              else
                                 rep.int(1, nX * nY)
                          } else {
                              if(.M.kind(X) != "n")
                                 rep.int(rep.int(as.double(Y@x), nj.), rep.dp) *
                                     as.double(X@x)[s.]
                              else
                                 rep.int(rep.int(as.double(Y@x), nj.), rep.dp)
                          }
                  }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "CsparseMatrix", Y = "CsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              uX <- uY <- FALSE
              if((sX <- .M.shape(X)) == "t")
                  uX <- X@diag != "N"
              if((sY <- .M.shape(Y)) == "t")
                  uY <- Y@diag != "N"
              shape <-
                  switch(sX,
                         g = "g",
                         t = if(sY == "t" && X@uplo == Y@uplo) "t" else "g",
                         s = if(sY == "s") "s" else "g")
              cl <- "d.CMatrix"
              if(!all(dr <- dX * dY)) {
                  substr(cl, 2L, 2L) <- shape
                  r <- new(cl)
                  r@Dim <- dr
                  r@p <- integer(dr[2L] + 1)
                  if(shape != "g") {
                      r@uplo <- X@uplo
                      if(shape == "t" && uX && uY)
                          r@diag <- "U"
                  }
              } else {
                  substr(cl, 2L, 2L) <- if(shape == "s") "g" else shape
                  r <- new(cl)
                  r@Dim <- dr
                  if(uX)
                      X <- ..diagU2N(X)
                  if(uY)
                      Y <- ..diagU2N(Y)
                  if(sY == "s")
                      Y <- .sparse2g(Y)
                  else if(sX == "s")
                      X <- .sparse2g(X)
                  if((nX <- (pX <- X@p)[length(pX)]) == 0L ||
                     (nY <- (pY <- Y@p)[length(pY)]) == 0L)
                      r@p <- integer(dr[2L] + 1)
                  else if(as.double(nX) * nY > .Machine$integer.max)
                      stop("number of nonzero entries cannot exceed 2^31-1")
                  else {
                      dpX <- pX[-1L] - (pX. <- pX[-length(pX)])
                      dpY <- pY[-1L] - (pY. <- pY[-length(pY)])

                      rep.dpX <- rep(dpX, each = dY[2L])
                      rep.dpY <- rep.int(dpY, dX[2L])

                      t1 <- rep.int(rep.dpY, rep.dpX)

                      s1 <- sequence.default(
                          nvec = rep.dpX,
                          from = rep(1L + pX., each = dY[2L]),
                          by = 1L)

                      s2 <- sequence.default(
                          nvec = t1,
                          from = rep.int(rep.int(1L + pY., dX[2L]), rep.dpX),
                          by = 1L)

                      r@p <- c(0L, cumsum(rep.dpX * dpY))
                      r@i <- rep.int((dY[1L] * X@i)[s1], t1) + Y@i[s2]
                      r@x <-
                          if(.M.kind(X) != "n") {
                              if(.M.kind(Y) != "n")
                                  rep.int(as.double(X@x)[s1], t1) *
                                      as.double(Y@x)[s2]
                              else
                                  rep.int(as.double(X@x)[s1], t1)
                          } else {
                              if(.M.kind(Y) != "n")
                                  as.double(Y@x)[s2]
                              else
                                  rep.int(1, nX * nY)
                          }
                  }
                  if(shape == "t") {
                      r@uplo <- X@uplo
                      if(uX && uY)
                          r <- ..diagN2U(r, sparse = TRUE)
                  } else if(shape == "s")
                      r <- .Call(R_sparse_force_symmetric, r, X@uplo)
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "RsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              .tCR2RC(kronecker(t(X), .tCR2RC(Y), FUN, make.dimnames, ...)))

setMethod("kronecker", signature(X = "RsparseMatrix", Y = "diagonalMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              .tCR2RC(kronecker(.tCR2RC(X), t(Y), FUN, make.dimnames, ...)))

setMethod("kronecker", signature(X = "RsparseMatrix", Y = "RsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              .tCR2RC(kronecker(.tCR2RC(X), .tCR2RC(Y), FUN, make.dimnames, ...)))

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "TsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              uX <- X@diag != "N"
              uY <- FALSE
              shape <- .M.shape(Y)
              r <- new(`substr<-`("d.TMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- Y@uplo
                  if(shape == "t" && (uX & (uY <- Y@diag != "N")))
                      r@diag <- "U"
              }
              if(all(dr)) {
                  if(any((kind <- .M.kind(Y)) == c("n", "l")))
                      Y <- .Call(Tsparse_aggregate, Y)
                  if(!uX && uY)
                      Y <- ..diagU2N(Y)
                  nX <- dX[1L]
                  nY <- length(Y@i)
                  r@i <-
                      rep(seq.int(0L, by = dY[1L], length.out = nX),
                          each = nY) +
                      Y@i
                  r@j <-
                      rep(seq.int(0L, by = dY[2L], length.out = nX),
                          each = nY) +
                      Y@j
                  r@x <-
                      if(uX) {
                          if(kind != "n")
                              rep.int(as.double(Y@x), nX)
                          else rep.int(1, as.double(nX) * nY)
                      } else {
                          if(kind != "n")
                              rep(as.double(X@x), each = nY) * as.double(Y@x)
                          else rep(as.double(X@x), each = nY)
                      }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "TsparseMatrix", Y = "diagonalMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              uX <- FALSE
              uY <- Y@diag != "N"
              shape <- .M.shape(X)
              r <- new(`substr<-`("d.TMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- X@uplo
                  if(shape == "t" && (uX <- X@diag != "N") && uY)
                      r@diag <- "U"
              }
              if(all(dr)) {
                  if(any((kind <- .M.kind(X)) == c("n", "l")))
                      X <- .Call(Tsparse_aggregate, X)
                  if(uX && !uY)
                      X <- ..diagU2N(X)
                  nX <- length(X@i)
                  nY <- dY[1L]
                  r@i <- rep(nY * X@i, each = nY) + 0:(nY-1L)
                  r@j <- rep(nY * X@j, each = nY) + 0:(nY-1L)
                  r@x <-
                      if(uY) {
                          if(kind != "n")
                              rep.int(as.double(Y@x), nY)
                          else rep.int(1, as.double(nX) * nY)
                      } else {
                          if(kind != "n")
                              rep(as.double(X@x), each = nY) * as.double(Y@x)
                          else rep(as.double(X@x), each = nY)
                      }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "TsparseMatrix", Y = "TsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              uX <- uY <- FALSE
              if((sX <- .M.shape(X)) == "t")
                  uX <- X@diag != "N"
              if((sY <- .M.shape(Y)) == "t")
                  uY <- Y@diag != "N"
              shape <-
                  switch(sX,
                         g = "g",
                         t = if(sY == "t" && X@uplo == Y@uplo) "t" else "g",
                         s = if(sY == "s") "s" else "g")
              cl <- "d.TMatrix"
              if(!all(dr <- dX * dY)) {
                  substr(cl, 2L, 2L) <- shape
                  r <- new(cl)
                  r@Dim <- dr
                  if(shape != "g") {
                      r@uplo <- X@uplo
                      if(shape == "t" && uX && uY)
                          r@diag <- "U"
                  }
              } else {
                  substr(cl, 2L, 2L) <- if(shape == "s") "g" else shape
                  r <- new(cl)
                  r@Dim <- dr
                  if(uX)
                      X <- ..diagU2N(X)
                  if(uY)
                      Y <- ..diagU2N(Y)
                  if(sY == "s")
                      Y <- .sparse2g(Y)
                  else if(sX == "s")
                      X <- .sparse2g(X)
                  nY <- length(Y@i)
                  r@i <- i. <- rep(dY[1L] * X@i, each = nY) + Y@i
                  r@j <-       rep(dY[2L] * X@j, each = nY) + Y@j
                  r@x <-
                      if(.M.kind(X) != "n") {
                          if(.M.kind(Y) != "n")
                              rep(as.double(X@x), each = nY) * as.double(Y@x)
                          else
                              rep(as.double(X@x), each = nY)
                      } else {
                          if(.M.kind(Y) != "n")
                              rep.int(as.double(X@x), nY)
                          else
                              rep.int(1, length(i.))
                      }

                  if(shape == "t") {
                      r@uplo <- X@uplo
                      if(uX && uY)
                          r <- ..diagN2U(r, sparse = TRUE)
                  } else if(shape == "s")
                      r <- .Call(R_sparse_force_symmetric, r, X@uplo)
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "indMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, as(Y, "nsparseMatrix"), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "indMatrix", Y = "diagonalMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(as(X, "nsparseMatrix"), Y, FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "indMatrix", Y = "indMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if((margin <- X@margin) != Y@margin)
                  kronecker(as(X, "CsparseMatrix"), as(Y, "CsparseMatrix"),
                            FUN, make.dimnames, ...)
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
                  stop("dimensions cannot exceed 2^31-1")
              r <- new("indMatrix")
              r@Dim <- dX * dY
              r@perm <-
                  if(margin == 1L)
                      rep(dY[2L] * (X@perm - 1L), each = dY[1L]) +
                          rep.int(Y@perm, dX[1L])
                  else {
                      r@margin <- 1L
                      rep(dY[1L] * (X@perm - 1L), each = dY[2L]) +
                          rep.int(Y@perm, dX[2L])
                  }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY)))
                  r@Dimnames <- dnr
              r
          })

## Catch everything else with these:

setMethod("kronecker", signature(X = "Matrix", Y = "matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, .m2ge(Y, "d"), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "Matrix", Y = "vector"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, .m2ge(Y, "d"), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "matrix", Y = "Matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(.m2ge(X, "d"), Y, FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "vector", Y = "Matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(.m2ge(X, "d"), Y, FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "denseMatrix", Y = "Matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(as(X, "CsparseMatrix"), Y, FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "CsparseMatrix", Y = "Matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, as(Y, "CsparseMatrix"), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "RsparseMatrix", Y = "Matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, as(Y, "RsparseMatrix"), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "TsparseMatrix", Y = "Matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, as(Y, "TsparseMatrix"), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "Matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, as(Y, "CsparseMatrix"), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "indMatrix", Y = "Matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(as(X, "CsparseMatrix"), Y, FUN, make.dimnames, ...))
