## METHODS FOR GENERIC: nnzero
## * used to retrieve number of nonzero elements,
##   i.e., number of elements excl. both structural and non-structural zeros
## * like MATLAB's nnz() but more sophisticated due to handling of NA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## na.counted:
## FALSE ... NA is treated as    zero and so excluded from count
##  TRUE ... NA is treated as nonzero and so included   in count
##    NA ... NA is indeterminate (could be zero or nonzero) hence count is NA

 sparseDefault <- function(x) length(x) > 2 *  nnzero(x, na.counted = TRUE)
.sparseDefault <- function(x) length(x) > 2 * .nnzero(x, na.counted = TRUE)

## For logical, integer, double, and complex vectors
.nnzero <- function(x, na.counted = NA, nnzmax = length(x))
    .Call(R_nnz, x, na.counted, nnzmax)

## For any class with methods for 'is.na' and '!='
.nnzero.fallback <- function(x, na.counted = NA)
    sum(if(is.na(na.counted))
            x != 0
        else if(na.counted)
            is.na(x) | x != 0
        else !is.na(x) & x != 0)

.nnzero.dispatching <- function(x, na.counted = NA)
    switch(typeof(x), logical =, integer =, double =, complex = .nnzero,
           .nnzero.fallback)(x, na.counted)

setMethod("nnzero",    "ANY", .nnzero.fallback)
setMethod("nnzero", "vector", .nnzero.dispatching)

setMethod("nnzero", "denseMatrix",
          function(x, na.counted = NA) {
              d <- x@Dim
              if(any(d == 0L))
                  return(0L)
              if(.M.kind(x) == "n")
                  na.counted <- TRUE
              if((shape <- .M.shape(x)) != "g")
                  x <- .M2packed(x)
              N <- .nnzero(x@x, na.counted)
              switch(shape,
                     "g" = N,
                     "s" = N + N - .nnzero(diag(x, names = FALSE), na.counted),
                     "t" = if(x@diag == "N") N else N + d[1L] - .nnzero(x@x[indDiag(d[1L], upper = x@uplo == "U", packed = TRUE)], na.counted))
          })

setMethod("nnzero", "sparseMatrix",
          function(x, na.counted = NA) {
              d <- x@Dim
              if(any(d == 0L))
                  return(0L)
              N <- switch(.M.repr(x),
                          "C" = x@p[d[2L]+1L],
                          "R" = x@p[d[1L]+1L],
                          "T" = length((x <- aggregateT(x))@i))
              if(.M.kind(x) != "n")
                  N <- .nnzero(x@x, na.counted, N)
              switch(.M.shape(x),
                     "g" = N,
                     "s" = N + N - .nnzero(diag(x, names = FALSE), na.counted),
                     "t" = if(x@diag == "N") N else N + d[1L])
          })

setMethod("nnzero", "diagonalMatrix",
          function(x, na.counted = NA) {
              if(x@diag != "N")
                  x@Dim[1L]
              else {
                  y <- x@x
                  if(.M.kind(x) == "n" && anyNA(y))
                      y <- y | is.na(y)
                  .nnzero(y, na.counted)
              }
          })

setMethod("nnzero", "indMatrix",
          function(x, na.counted = NA)
              length(x@perm))

setMethod("nnzero", "CHMfactor",
          function(x, na.counted = NA)
              nnzero(as(x, "CsparseMatrix"), na.counted))

rm(.nnzero.dispatching)
