## METHODS FOR GENERIC: lu
## partially pivoted LU factorization, returning denseLU or sparseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lu", signature(x = "matrix"),
          function(x, ...) lu(.m2ge(x, "d"), ...))

setMethod("lu", signature(x = "denseMatrix"),
	  function(x, ...) lu(..dense2d(x), ...))

setMethod("lu", signature(x = "dgeMatrix"),
	  function(x, warnSing = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              .Call(dgeMatrix_LU, x, warnSing)
          })

setMethod("lu", signature(x = "dsyMatrix"),
	  function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.dense2g(x), ...)
              if(cache) .set.factors(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dspMatrix"),
	  function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.dense2g(x), ...)
              if(cache) .set.factors(x, "LU", r) else r
          })

for(.cl in c("dtrMatrix", "dtpMatrix"))
setMethod("lu", signature(x = .cl),
	  function(x, ...) {
              if(x@uplo == "U" || x@diag == "U") {
                  r <- new("denseLU")
                  r@Dim <- d <- x@Dim
                  r@perm <- seq_len(d[1L])
                  r@x <- .dense2g(x)@x
                  r
              } else lu(.dense2g(x), ...)
          })
rm(.cl)

setMethod("lu", signature(x = "sparseMatrix"),
	  function(x, ...) lu(..sparse2d(as(x, "CsparseMatrix")), ...))

setMethod("lu", signature(x = "dgCMatrix"),
          function(x, errSing = TRUE, order = TRUE,
                   tol = 1.0, keep.dimnames = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              .Call(dgCMatrix_LU, x, order, tol, errSing, keep.dimnames)
          })

setMethod("lu", signature(x = "dsCMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.sparse2g(x), ...)
              if(cache) .set.factors(x, "LU", r) else r
          })

setMethod("lu", "dtCMatrix",
	  function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- r@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@L <- if(upper) y else x
                  r@U <- if(upper) x else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
              } else lu(.sparse2g(x), ...)
	  })

setMethod("lu", signature(x = "dgRMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.CR2RC(x), ...)
              if(cache) .set.factors(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dsRMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.sparse2g(.tCR2RC(x)), ...)
              if(cache) .set.factors(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dtRMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- r@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@L <- if(upper) y else .CR2RC(x)
                  r@U <- if(upper) .CR2RC(x) else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
              } else lu(.sparse2g(.CR2RC(x)), ...)
          })

setMethod("lu", signature(x = "dgTMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.T2C(x), ...)
              if(cache) .set.factors(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dsTMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.sparse2g(.T2C(x)), ...)
              if(cache) .set.factors(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dtTMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- r@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@L <- if(upper) y else .T2C(x)
                  r@U <- if(upper) .T2C(x) else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
                  r
              } else lu(.sparse2g(.T2C(x)), ...)
          })

setMethod("lu", "diagonalMatrix",
	  function(x, ...) {
              n <- (d <- x@Dim)[1L]
              L <- new("dtCMatrix")
              r <- new("sparseLU")
              L@Dim <- r@Dim <- d
              L@uplo <- "L"
              L@diag <- "U"
              L@p <- integer(n + 1L)
              r@L <- L
              if(x@diag == "N") {
                  L@diag <- "N"
                  L@p <- seq.int(from = 0L, length.out = n + 1L)
                  L@x <- as.double(x@x)
              }
              r@U <- L
              r@p <- r@q <- seq.int(from = 0L, length.out = n)
              r
	  })


## METHODS FOR CLASS: denseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## returning list(L, U, P), where A = P L U
setMethod("expand", signature(x = "denseLU"),
          function(x, ...) .Call(LU_expand, x))


## METHODS FOR CLASS: sparseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## returning list(P, L, U, Q), where A = P L U Q
setMethod("expand", signature(x = "sparseLU"),
          function(x, ...) list(P = as(x@p + 1L, "pMatrix"),
                                L = x@L,
                                U = x@U,
                                Q = as(x@q + 1L, "pMatrix")))
