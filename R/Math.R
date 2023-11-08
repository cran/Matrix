## METHODS FOR GENERIC: Math (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Math")
##  [1] "abs"      "sign"     "sqrt"     "ceiling"  "floor"    "trunc"
##  [7] "cummax"   "cummin"   "cumprod"  "cumsum"   "exp"      "expm1"
## [13] "log"      "log10"    "log2"     "log1p"    "cos"      "cosh"
## [19] "sin"      "sinh"     "tan"      "tanh"     "acos"     "acosh"
## [25] "asin"     "asinh"    "atan"     "atanh"    "cospi"    "sinpi"
## [31] "tanpi"    "gamma"    "lgamma"   "digamma"  "trigamma"

setMethod("Math", signature(x = "denseMatrix"),
          function(x) {
              g <- get(.Generic, mode = "function")
              if(startsWith(.Generic, "cum"))
                  return(g(.M2v(x)))
              cl <- .M.nonvirtual(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              if (kind == "z") {
                  zero <- 0+0i; one <- 1+0i; a <- as.complex
              } else {
                  zero <- 0   ; one <- 1   ; a <- as.double
                  substr(cl, 1L, 1L) <- "d"
              }
              if(shape == "t") {
                  stay0 <- is0(a(g(zero)))
                  if(!stay0) {
                      x <- .M2gen(x)
                      substr(cl, 2L, 3L) <- "ge"
                  }
              }
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if(shape == "s" || (shape == "t" && stay0))
                  r@uplo <- x@uplo
              r@x <- a(g({ y <- x@x; if(kind == "n") y | is.na(y) else y }))
              if(shape == "t" && stay0 && x@diag != "N") {
                  if(is1(g1 <- a(g(one))))
                      r@diag <- "U"
                  else diag(r) <- g1
              }
              r
          })

setMethod("log", signature(x = "denseMatrix"),
          function(x, ...) {
              cl <- .M.nonvirtual(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              if(kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              if(shape == "t") {
                  x <- .M2gen(x)
                  substr(cl, 2L, 3L) <- "ge"
              }
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if(shape == "s")
                  r@uplo <- x@uplo
              r@x <- log({ y <- x@x; if(kind == "n") y | is.na(y) else y }, ...)
              r
          })

setMethod("Math", signature(x = "sparseMatrix"),
          function(x) {
              g <- get(.Generic, mode = "function")
              if(startsWith(.Generic, "cum"))
                  return(g(.M2v(x)))
              cl <- .M.nonvirtual(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              if (kind == "z") {
                  zero <- 0+0i; one <- 1+0i; a <- as.complex
              } else {
                  zero <- 0   ; one <- 1   ; a <- as.double
                  substr(cl, 1L, 1L) <- "d"
              }
              stay0 <- is0(g0 <- a(g(zero)))
              if(!stay0)
                  substr(cl, 2L, 3L) <- if(shape == "s") "sy" else "ge"
              r <- new(cl)
              r@Dim      <- x@Dim
              r@Dimnames <- x@Dimnames
              if(shape == "s" || (shape == "t" && stay0))
                  r@uplo <- x@uplo
              if(!stay0) {
                  y <- .Call(CR2spV, if(repr == "T") .M2C(x) else x)
                  tmp <- rep.int(g0, y@length)
                  tmp[y@i] <- a(g(if(kind == "n") one else y@x))
                  r@x <- tmp
              } else {
                  if(shape == "t" && x@diag != "N") {
                      if(is1(a(g(one))))
                          r@diag <- "U"
                      else diag(x) <- TRUE
                  }
                  nnz <- length(
                      switch(repr,
                             "C" = { r@p <- x@p; r@i <- x@i },
                             "R" = { r@p <- x@p; r@j <- x@j },
                             "T" = { r@i <- x@i; r@j <- x@j }))
                  r@x <- if(kind == "n") rep.int(a(g(one)), nnz) else a(g(x@x))
              }
              r
          })

setMethod("log", signature(x = "sparseMatrix"),
          function(x, ...) {
              cl <- .M.nonvirtual(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              if(kind == "z") {
                  zero <- 0+0i; one <- 1+0i
              } else {
                  zero <- 0   ; one <- 1
                  substr(cl, 1L, 1L) <- "d"
              }
              substr(cl, 2L, 3L) <- if(shape == "s") "sy" else "ge"
              r <- new(cl)
              r@Dim      <- x@Dim
              r@Dimnames <- x@Dimnames
              if(shape == "s")
                  r@uplo <- x@uplo
              y <- .Call(CR2spV, if(repr == "T") .M2C(x) else x)
              tmp <- rep.int(log(zero, ...), y@length)
              tmp[y@i] <- log(if(kind == "n") one else y@x, ...)
              r@x <- tmp
              r
          })

setMethod("Math", signature(x = "diagonalMatrix"),
          function(x) {
              g <- get(.Generic, mode = "function")
              if(startsWith(.Generic, "cum"))
                  return(g(.M2v(x)))
              cl <- .M.nonvirtual(x)
              kind <- substr(cl, 1L, 1L)
              if (kind == "z") {
                  zero <- 0+0i; one <- 1+0i; a <- as.complex
              } else {
                  zero <- 0   ; one <- 1   ; a <- as.double
                  substr(cl, 1L, 1L) <- "d"
              }
              stay0 <- is0(g0 <- a(g(zero)))
              if(!stay0)
                  substr(cl, 2L, 3L) <- "ge"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if(!stay0) {
                  if((n <- d[2L]) > 0L) {
                      tmp <- matrix(g0, n, n)
                      diag(tmp) <- a(g(if(x@diag != "N") one else { y <- x@x; if(kind == "n" && anyNA(y)) y | is.na(y) else y }))
                      dim(tmp) <- NULL
                      r@x <- tmp
                  }
              } else {
                  if(x@diag != "N") {
                      if(is1(g1 <- a(g(one))))
                          r@diag <- "U"
                      else r@x <- rep.int(g1, d[1L])
                  } else r@x <- a(g({ y <- x@x; if(kind == "n" && anyNA(y)) y | is.na(y) else y }))
              }
              r
          })

setMethod("log", signature(x = "diagonalMatrix"),
          function(x, ...) {
              cl <- .M.nonvirtual(x)
              kind <- substr(cl, 1L, 1L)
              if(kind == "z") {
                  zero <- 0+0i; one <- 1+0i
              } else {
                  zero <- 0   ; one <- 1
                  substr(cl, 1L, 1L) <- "d"
              }
              substr(cl, 2L, 3L) <- "ge"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if((n <- d[2L]) > 0L) {
                  tmp <- matrix(log(zero, ...), n, n)
                  diag(tmp) <- log(if(x@diag != "N") one else { y <- x@x; if(kind == "n" && anyNA(y)) y | is.na(y) else y }, ...)
                  dim(tmp) <- NULL
                  r@x <- tmp
              }
              r
          })

setMethod("Math", signature(x = "indMatrix"),
          function(x)
              get(.Generic, mode = "function")(.M2kind(x, "n")))

setMethod("log", signature(x = "indMatrix"),
          function(x, ...)
              log(.M2kind(x, "n"), ...))

setMethod("Math", signature(x = "sparseVector"),
          function(x) {
              g <- get(.Generic, mode = "function")
              if(startsWith(.Generic, "cum"))
                  return(g(.V2v(x)))
              kind <- .M.kind(x)
              if(kind == "z") {
                  zero <- 0+0i; one <- 1+0i; l <- "z"
              } else if(kind == "d" || .Generic != "abs") {
                  zero <- 0   ; one <- 1   ; l <- "d"
              } else {
                  zero <- 0L  ; one <- 1L  ; l <- "i"
              }
              if(isN0(g0 <- g(zero))) {
                  r <- rep.int(g0, x@length)
                  if((nnz <- length(x@i)) > 0L)
                      r[x@i] <- if(kind == "n") rep.int(g(one), nnz) else g(x@x)
              } else {
                  r <- new(paste0(l, "sparseVector"))
                  r@length <- x@length
                  r@i <- x@i
                  if((nnz <- length(x@i)) > 0L)
                      r@x <- if(kind == "n") rep.int(g(one), nnz) else g(x@x)
              }
              r
          })

setMethod("log", signature(x = "sparseVector"),
          function(x, ...) {
              kind <- .M.kind(x)
              if(kind == "z") {
                  zero <- 0+0i; one <- 1+0i
              } else {
                  zero <- 0   ; one <- 1
              }
              r <- rep.int(log(zero, ...), x@length)
              if(length(x@i) > 0L)
                  r[x@i] <- log(if(kind == "n") one else x@x, ...)
              r
          })


## METHODS FOR GENERIC: Math2 (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Math2")
## [1] "round"  "signif"

setMethod("Math2", signature(x = "Matrix"),
          function(x, digits) {
              x <- .indefinite(.M2kind(x, ","))
              x@x <- get(.Generic, mode = "function")(x@x, digits = digits)
              if(.hasSlot(x, "factors") && length(x@factors) > 0L)
                  x@factors <- list()
              x
          })

setMethod("Math2", signature(x = "sparseVector"),
          function(x, digits) {
              x <- .V2kind(x, ",")
              x@x <- get(.Generic, mode = "function")(x@x, digits = digits)
              x
          })


## METHODS FOR GENERIC: zapsmall
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("zapsmall", signature(x = "Matrix"),
          function(x, digits = getOption("digits")) {
              x <- .indefinite(.M2kind(x, ","))
              x@x <- zapsmall(x@x, digits = digits)
              if(.hasSlot(x, "factors") && length(x@factors) > 0L)
                  x@factors <- list()
              x
          })

setMethod("zapsmall", signature(x = "sparseVector"),
          function(x, digits = getOption("digits")) {
              x <- .V2kind(x, ",")
              x@x <- zapsmall(x@x, digits = digits)
              x
          })
