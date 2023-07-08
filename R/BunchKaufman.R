## METHODS FOR GENERIC: BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("BunchKaufman", signature(x = "dsyMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(dsyMatrix_trf, x, as.logical(warnSing)))

setMethod("BunchKaufman", signature(x = "dspMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(dspMatrix_trf, x, as.logical(warnSing)))

setMethod("BunchKaufman", signature(x = "matrix"),
          function(x, uplo = "U", ...)
              BunchKaufman(.m2dense(x, "dsy", uplo), ...))


## METHODS FOR CLASS: p?BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("BunchKaufman", "dtrMatrix",
      function(from) {
          to <- new("dtrMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          to
      })

setAs("pBunchKaufman", "dtpMatrix",
      function(from) {
          to <- new("dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          to
      })

.def.unpacked <- .def.packed <- function(x, which, ...) {
    r <- .Call(BunchKaufman_expand, x, .PACKED)
    b <- length(r) - 1L
    switch(which,
           "DU" =, "DL" = {
               if(!endsWith(which, x@uplo))
                   stop(gettextf("which=\"%s\" invalid for x@uplo=\"%s\"",
                                 which, x@uplo),
                        domain = NA)
               r[[b + 1L]]
           },
           "U" =, "U." =, "L" =, "L." = {
               if(!startsWith(which, x@uplo))
                   stop(gettextf("which=\"%s\" invalid for x@uplo=\"%s\"",
                                 which, x@uplo),
                        domain = NA)
               if(b > 0L) {
                   m <- r[[b]]
                   if(b > 1L)
                       for(i in (b - 1L):1L)
                           m <- r[[i]] %*% m
                   if(endsWith(which, ".")) t(m) else m
               } else {
                   m <- new("ddiMatrix")
                   m@Dim <- x@Dim
                   m@diag <- "U"
                   m
               }
           },
           stop(gettextf("'which' is not \"%1$s\", \"D%1$s\", or \"%1$s.\"",
                         x@uplo),
                domain = NA))
}
body(.def.unpacked) <-
    do.call(substitute, list(body(.def.unpacked), list(.PACKED = FALSE)))
body(.def.packed) <-
    do.call(substitute, list(body(.def.packed  ), list(.PACKED =  TRUE)))

setMethod("expand1", signature(x =  "BunchKaufman"), .def.unpacked)
setMethod("expand1", signature(x = "pBunchKaufman"), .def.packed)
rm(.def.unpacked, .def.packed)

.def.unpacked <- .def.packed <- function(x, complete = FALSE, ...) {
    r <- .Call(BunchKaufman_expand, x, .PACKED)
    b <- length(r) - 1L
    if(complete) {
        if(b > 0L)
            r <- c(r, lapply(r[b:1L], t))
    } else {
        if(b > 0L) {
            m <- r[[b]]
            if(b > 1L)
                for(i in (b - 1L):1L)
                    m <- r[[i]] %*% m
            r <- list(m, r[[b + 1L]], t(m))
        } else {
            m <- new("ddiMatrix")
            m@Dim <- x@Dim
            m@diag <- "U"
            r <- list(m, r[[1L]], m)
        }
        names(r) <- if(x@uplo == "U") c("U", "DU", "U.") else c("L", "DL", "L.")
    }
    dn <- x@Dimnames
    if(length(r) == 1L)
        r[[1L]]@Dimnames <- dn
    else {
        r[[1L]]@Dimnames <- c(dn[1L], list(NULL))
        r[[length(r)]]@Dimnames <- c(list(NULL), dn[2L])
    }
    r
}
body(.def.unpacked) <-
    do.call(substitute, list(body(.def.unpacked), list(.PACKED = FALSE)))
body(.def.packed) <-
    do.call(substitute, list(body(.def.packed  ), list(.PACKED =  TRUE)))

## returning
## list(U, DU, U') where A = U DU U' and U = P[b] U[b] ... P[1] U[1]
## OR
## list(L, DL, L') where A = L DL L' and L = P[1] L[1] ... P[b] L[b]
setMethod("expand2", signature(x =  "BunchKaufman"), .def.unpacked)
setMethod("expand2", signature(x = "pBunchKaufman"), .def.packed)
rm(.def.unpacked, .def.packed)
