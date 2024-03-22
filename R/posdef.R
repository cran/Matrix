## METHODS FOR CLASS: dpoMatrix, dppMatrix
## dense symmetric positive semidefinite matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Operations such as rounding can lose positive semidefiniteness
## but not symmetry, hence:
.indefinite <- function(x) {
    cl <- .M.nonvirtual(x, TRUE)
    if(any(cl == c("dpoMatrix", "corMatrix")))
        as(x, "dsyMatrix")
    else if(any(cl == c("dppMatrix", "copMatrix")))
        as(x, "dspMatrix")
    else x
}

.dsy2dpo <- .dsp2dpp <- function(from) {
    if(is.null(tryCatch(Cholesky(from, perm = FALSE),
                        error = function(e) NULL)))
        stop("not a positive definite matrix (and positive semidefiniteness is not checked)")
    to <- new(.CLASS)
    to@Dim <- from@Dim
    to@Dimnames <- from@Dimnames
    to@uplo <- from@uplo
    to@x <- from@x
    to@factors <- from@factors
    to
}
body(.dsy2dpo)[[3L]][[3L]][[2L]] <- "dpoMatrix"
body(.dsp2dpp)[[3L]][[3L]][[2L]] <- "dppMatrix"

setAs("dppMatrix", "dpoMatrix", function(from) unpack(from))
setAs("dpoMatrix", "dppMatrix", function(from)   pack(from))

setAs("dsyMatrix", "dpoMatrix", .dsy2dpo)
setAs("dspMatrix", "dppMatrix", .dsp2dpp)

setAs("Matrix", "dpoMatrix",
      function(from) .dsy2dpo(.M2unpacked(.M2sym(.M2kind(from, "d")))))
setAs("Matrix", "dppMatrix",
      function(from) .dsp2dpp(.M2packed  (.M2sym(.M2kind(from, "d")))))

setAs("matrix", "dpoMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .dsy2dpo(.M2sym(from))
      })
setAs("matrix", "dppMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .dsp2dpp(pack(from, symmetric = TRUE))
      })


## METHODS FOR CLASS: corMatrix, copMatrix
## dense correlation matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dpo2cor <- function(from) {
    if(!is.null(to <- from@factors$correlation))
        return(to)
    sd <- sqrt(diag(from, names = FALSE))

    to <- new("corMatrix")
    to@Dim <- d <- from@Dim
    to@Dimnames <- from@Dimnames
    to@uplo <- from@uplo
    to@sd <- sd

    n <- d[1L]
    x <- from@x / sd / rep(sd, each = n)
    x[indDiag(n)] <- 1
    to@x <- x

    .set.factor(from, "correlation", to)
}

.dpp2cop <- function(from) {
    if(!is.null(to <- from@factors$correlation))
        return(to)
    sd <- sqrt(diag(from, names = FALSE))

    to <- new("copMatrix")
    to@Dim <- d <- from@Dim
    to@Dimnames <- from@Dimnames
    to@uplo <- uplo <- from@uplo
    to@sd <- sd

    n <- d[1L]
    u <- uplo == "U"
    if(u) {
        r <- seq_len(n)
        s <- 1L
    } else {
        r <- seq.int(to = 1L, by = -1L, length.out = n)
        s <- seq_len(n)
    }
    x <-  from@x / rep.int(sd, r) / sd[sequence.default(r, s)]
    x[indDiag(n, upper = u, packed = TRUE)] <- 1
    to@x <- x

    .set.factor(from, "correlation", to)
}

setAs("copMatrix", "corMatrix", function(from) unpack(from))
setAs("corMatrix", "copMatrix", function(from)   pack(from))

setAs("dpoMatrix", "corMatrix", .dpo2cor)
setAs("dppMatrix", "copMatrix", .dpp2cop)

setAs("dsyMatrix", "corMatrix",
      function(from) .dpo2cor(.dsy2dpo(from)))
setAs("dspMatrix", "copMatrix",
      function(from) .dpp2cop(.dsp2dpp(from)))

setAs("Matrix", "corMatrix",
      function(from) .dpo2cor(.dsy2dpo(.M2unpacked(.M2sym(.M2kind(from, "d"))))))
setAs("Matrix", "copMatrix",
      function(from) .dpp2cop(.dsp2dpp(.M2packed  (.M2sym(.M2kind(from, "d"))))))

setAs("matrix", "corMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .dpo2cor (.dsy2dpo(.M2sym(from)))
      })
setAs("matrix", "copMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .dpp2cop(.dsp2dpp(pack(from, symmetric = TRUE)))
      })


## METHODS FOR GENERIC: cov2cor
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("cov2cor", c(V = "unpackedMatrix"),
          function(V) {
              d <- V@Dim
              if(d[1L] != d[2L] || .M.kind(V) == "z")
                  stop(gettextf("'%s' is not a square numeric matrix", "V"),
                       domain = NA)
              as(forceSymmetric(V), "corMatrix")
          })

setMethod("cov2cor", c(V = "packedMatrix"),
          function(V) {
              d <- V@Dim
              if(d[1L] != d[2L] || .M.kind(V) == "z")
                  stop(gettextf("'%s' is not a square numeric matrix", "V"),
                       domain = NA)
              as(forceSymmetric(V), "copMatrix")
          })

setMethod("cov2cor", c(V = "sparseMatrix"),
          function(V) {
              d <- V@Dim
              if(d[1L] != d[2L] || .M.kind(V) == "z")
                  stop(gettextf("'%s' is not a square numeric matrix", "V"),
                       domain = NA)
              dn <- symDN(V@Dimnames)
              V <- .M2kind(V, "d")
              V.ii <- diag(V, names = FALSE)
              if(length(V.ii) > 0L && is.na(m <- min(V.ii)) || m <= 0)
                  warning(gettextf("diag(%s) has non-positive or non-finite entries; finite result is doubtful",
                                   "V"),
                          domain = NA)
              D <- Diagonal(x = sqrt(1/V.ii))
              r <- forceSymmetric(D %*% V %*% D)
              diag(r) <- 1
              r@Dimnames <- dn
              r
          })
