## METHODS FOR GENERIC: norm
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("norm", signature(x = "ANY", type = "missing"),
          function(x, type, ...) norm(x, type = "O", ...))

setMethod("norm", signature(x = "sparseMatrix", type = "character"),
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
                             warning("'norm' via sparse -> dense coercion")
                             base::norm(.sparse2m(x), type = "2")
                         },
                     "M" = , "m" =
                         max(abs(x)),
                     "F" = , "f" = , "E" = , "e" =
                         sqrt(sum(x * x)),
                     stop("invalid 'type'"))
          })

setMethod("norm", signature(x = "diagonalMatrix", type = "character"),
          function(x, type, ...) {
              if((n <- x@Dim[1L]) == 0L)
                  return(0)
              switch(EXPR = type[1L],
                     "O" = , "o" = , "1" = ,
                     "I" = , "i" = ,
                     "2" = ,
                     "M" = , "m" =
                         if(x@diag == "N") max(abs(x@x)) else 1,
                     "F" = , "f" = , "E" = , "e" =
                         if(x@diag == "N") sqrt(sum(x@x * x@x)) else sqrt(n),
                     stop("invalid 'type'"))
          })

setMethod("norm", signature(x = "indMatrix", type = "character"),
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
                     stop("invalid 'type'"))
          })

setMethod("norm", signature(x = "pMatrix", type = "character"),
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
                     stop("invalid 'type'"))
          })

setMethod("norm", signature(x = "denseMatrix", type = "character"),
          function(x, type, ...) norm(..dense2d(x), type = type, ...))

setMethod("norm", signature(x = "dgeMatrix", type = "character"),
          function(x, type, ...)
              if(identical(type, "2"))
                  base::norm(.dense2m(x), type = "2")
              else .Call(dgeMatrix_norm, x, type))

setMethod("norm", signature(x = "dtrMatrix", type = "character"),
          function(x, type, ...) {
              if(identical(type, "2"))
                  base::norm(.dense2m(x), type = "2")
              else .Call(dtrMatrix_norm, x, type)
          })

setMethod("norm", signature(x = "dtpMatrix", type = "character"),
          function(x, type, ...)
              if(identical(type, "2"))
                  base::norm(.dense2m(x), type = "2")
              else .Call(dtpMatrix_norm, x, type))

setMethod("norm", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...)
              if(identical(type, "2"))
                  base::norm(.dense2m(x), type = "2")
              else .Call(dsyMatrix_norm, x, type))

setMethod("norm", signature(x = "dspMatrix", type = "character"),
          function(x, type, ...)
              if(identical(type, "2"))
                  base::norm(.dense2m(x), type = "2")
              else .Call(dspMatrix_norm, x, type))
