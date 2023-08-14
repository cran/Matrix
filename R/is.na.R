## METHODS FOR GENERIC: anyNA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("anyNA", signature(x = "diagonalMatrix"),
          function(x) anyNA(x@x))

setMethod("anyNA", signature(x = "indMatrix"),
          function(x) FALSE)

setMethod("anyNA", signature(x = "nMatrix"),
          function(x) FALSE)

for(.kind in c("d", "l")) {
setMethod("anyNA", signature(x = paste0(.kind, "sparseMatrix")),
          function(x) anyNA(x@x))

setMethod("anyNA", signature(x = paste0(.kind, "denseMatrix")),
          function(x) {
              if(!.hasSlot(x, "uplo"))
                  return(anyNA(x@x))
              packed <- .isPacked(x)
              nonunit <- !.hasSlot(x, "diag") || x@diag == "N"
              if(packed && nonunit)
                  return(anyNA(x@x))
              k <- indTri(n = x@Dim[1L], upper = x@uplo == "U",
                          diag = nonunit, packed = packed)
              anyNA(x@x[k])
          })
}
rm(.kind)

setMethod("anyNA", signature(x = "nsparseVector"),
          function(x) FALSE)

setMethod("anyNA", signature(x = "sparseVector"),
          function(x) anyNA(x@x))


## METHODS FOR GENERIC: is.na
## [[ one more in ./abIndex.R ]]
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Returning an all-FALSE matrix preserving the structure of 'x',
## which can be any R object inheriting from Matrix
allFalseMatrix <- function(x) {
    r <- new(if(.hasSlot(x, "diag")) # triangularMatrix, diagonalMatrix
                 "ntCMatrix"
             else if(.hasSlot(x, "uplo")) # symmetricMatrix
                 "nsCMatrix"
             else "ngCMatrix") # generalMatrix
    r@Dim <- d <- x@Dim
    r@Dimnames <- x@Dimnames
    r@p <- integer(d[2L] + 1)
    r
}

setMethod("is.na", signature(x = "diagonalMatrix"),
          function(x) {
              r <- new("ldiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if(x@diag == "N") is.na(x@x) else logical(d[1L])
              r
          })

setMethod("is.na", signature(x = "indMatrix"),
          allFalseMatrix)

setMethod("is.na", signature(x = "nMatrix"),
          allFalseMatrix)

setMethod("is.na", signature(x = "dsparseMatrix"),
          function(x) {
              if(anyNA(x@x)) { # don't allocate in FALSE case
                  r <- .M2kind(diagU2N(x), "l")
                  r@x <- is.na(r@x)
                  .M2kind(.drop0(r), "n")
              } else allFalseMatrix(x)
          })

setMethod("is.na", signature(x = "lsparseMatrix"),
          function(x) {
              if(anyNA(x@x)) { # don't allocate in FALSE case
                  r <- diagU2N(x)
                  r@x <- is.na(r@x)
                  .M2kind(.drop0(r), "n")
              } else allFalseMatrix(x)
          })

.is.na.ge <- function(x) {
    if(anyNA(x@x)) # don't allocate in FALSE case
        new("ngeMatrix", Dim = x@Dim, Dimnames = x@Dimnames, x = is.na(x@x))
    else allFalseMatrix(x)
}

.is.na.tr <- function(x) {
    if(anyNA(x@x)) { # don't allocate in FALSE case
        d <- x@Dim
        i <- is.na(x@x)
        k <- indTri(n = d[1L], upper = x@uplo != "U",
                    diag = x@diag != "N", packed = FALSE)
        i[k] <- FALSE
        if(any(i))
            new("ntrMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.na.tp <- function(x) {
    if(anyNA(x@x)) { # don't allocate in FALSE case
        d <- x@Dim
        i <- is.na(x@x)
        if(x@diag != "N") {
            k <- indDiag(n = d[1L], upper = x@uplo == "U",
                         packed = TRUE)
            i[k] <- FALSE
        }
        if(any(i))
            new("ntpMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.na.sy <- function(x) {
    if(anyNA(x@x)) { # don't allocate in FALSE case
        d <- x@Dim
        i <- is.na(x@x)
        k <- indTri(n = d[1L], upper = x@uplo != "U",
                    diag = FALSE, packed = FALSE)
        i[k] <- FALSE
        if(any(i))
            new("nsyMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.na.sp <- function(x) {
    if(anyNA(x@x)) # don't allocate in FALSE case
        new("nspMatrix", Dim = x@Dim, Dimnames = x@Dimnames,
            x = is.na(x@x), uplo = x@uplo)
    else allFalseMatrix(x)
}

for(.kind in c("d", "l"))
    for(.xx in c("ge", "tr", "tp", "sy", "sp"))
        setMethod("is.na", signature(x = paste0(.kind, .xx, "Matrix")),
                  get(paste0(".is.na.", .xx),
                      mode = "function", inherits = FALSE))
rm(.is.na.ge, .is.na.tr, .is.na.tp,
   .is.na.sy, .is.na.sp,
   .kind, .xx)

setMethod("is.na", signature(x = "sparseVector"),
          function(x) new("nsparseVector", length = x@length,
                          i = x@i[is.na(x@x)]))

setMethod("is.na", signature(x = "nsparseVector"),
          function(x) new("nsparseVector", length = x@length))


## METHODS FOR GENERIC: is.finite
## [[ one more in ./abIndex.R ]]
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

allTrueMatrix <- function(x, symmetric = NA, packed = TRUE) {
    if(is.na(symmetric))
        ##  TRUE for  symmetricMatrix, diagonalMatrix
        ## FALSE for triangularMatrix,  generalMatrix
        symmetric <- if(.hasSlot(x, "uplo"))
                         !.hasSlot(x, "diag")
                     else .hasSlot(x, "diag")
    r <- new(if(!symmetric)
                 "ngeMatrix"
             else if(packed)
                 "nspMatrix"
             else "nsyMatrix")
    r@Dim <- d <- x@Dim
    r@Dimnames <- x@Dimnames
    r@x <- rep.int(TRUE, if(symmetric && packed) {
                             n <- d[1L]
                             0.5 * n * (n + 1)
                         } else prod(d))
    if(symmetric && .hasSlot(x, "uplo"))
        r@uplo <- x@uplo
    r
}

..allTrueMatrix <- function(x) allTrueMatrix(x)

setMethod("is.finite", signature(x = "diagonalMatrix"),
          function(x) {
              r <- allTrueMatrix(x, symmetric = TRUE, packed = TRUE)
              if(x@diag == "N") {
                  k <- indDiag(n = x@Dim[1L], upper = r@uplo == "U",
                               packed = TRUE)
                  r@x[k] <- is.finite(x@x)
              }
              r
          })

setMethod("is.finite", signature(x = "indMatrix"),
          ..allTrueMatrix)

setMethod("is.finite", signature(x = "nMatrix"),
          ..allTrueMatrix)

setMethod("is.finite", signature(x = "dsparseMatrix"),
          function(x) {
              if(!all(is.finite(x@x))) {
                  ## FIXME: use packed=TRUE once [<- is fast for packedMatrix
                  r <- allTrueMatrix(x, symmetric = NA, packed = FALSE)
                  if(.hasSlot(x, "p"))
                      x <- .M2T(x)
                  n <- x@Dim[1L]
                  w <- which(!is.finite(x@x))
                  r@x[as.double(n) * x@j[w] + x@i[w] + 1] <- FALSE
                  r
              } else allTrueMatrix(x, symmetric = NA, packed = TRUE)
          })

setMethod("is.finite", signature(x = "lsparseMatrix"),
          function(x) {
              if(anyNA(x@x)) { # don't allocate in FALSE case
                  ## FIXME: use packed=TRUE once [<- is fast for packedMatrix
                  r <- allTrueMatrix(x, symmetric = NA, packed = FALSE)
                  if(.hasSlot(x, "p"))
                      x <- .M2T(x)
                  n <- x@Dim[1L]
                  w <- which(is.na(x@x))
                  r@x[as.double(n) * x@j[w] + x@i[w] + 1] <- FALSE
                  r
              } else allTrueMatrix(x, symmetric = NA, packed = TRUE)
          })

.is.finite.ge <- function(x)
    new("ngeMatrix", Dim = x@Dim, Dimnames = x@Dimnames, x = is.finite(x@x))

.is.finite.tr <- function(x) {
    d <- x@Dim
    i <- is.finite(x@x)
    k <- indTri(n = d[1L], upper = x@uplo != "U",
                diag = x@diag != "N", packed = FALSE)
    i[k] <- TRUE
    new("ngeMatrix", Dim = d, Dimnames = x@Dimnames, x = i)
}

.is.finite.tp <- function(x) {
    d <- x@Dim
    i <- rep.int(TRUE, prod(d))
    k <- indTri(n = d[1L], upper = x@uplo == "U",
                diag = TRUE, packed = FALSE)
    i[k] <- is.finite(x@x)
    if(x@diag != "N") {
        k <- indDiag(n = d[1L], packed = FALSE)
        i[k] <- TRUE
    }
    new("ngeMatrix", Dim = d, Dimnames = x@Dimnames, x = i)
}

.is.finite.sy <- function(x)
    new("nsyMatrix", Dim = x@Dim, Dimnames = x@Dimnames,
        x = is.finite(x@x), uplo = x@uplo)

.is.finite.sp <- function(x)
    new("nspMatrix", Dim = x@Dim, Dimnames = x@Dimnames,
        x = is.finite(x@x), uplo = x@uplo)

for(.kind in c("d", "l"))
    for(.xx in c("ge", "tr", "tp", "sy", "sp"))
        setMethod("is.finite", signature(x = paste0(.kind, .xx, "Matrix")),
                  get(paste0(".is.finite.", .xx),
                      mode = "function", inherits = FALSE))
rm(.is.finite.ge, .is.finite.tr, .is.finite.tp,
   .is.finite.sy, .is.finite.sp, ..allTrueMatrix,
   .kind, .xx)

setMethod("is.finite", signature(x = "sparseVector"),
          function(x)  {
              r <- rep.int(TRUE, x@length)
              r[x@i[!is.finite(x@x)]] <- FALSE
              r
          })

setMethod("is.finite", signature(x = "nsparseVector"),
          function(x) rep.int(TRUE, x@length))


## METHODS FOR GENERIC: is.infinite
## NB: completely (!) parallel to 'is.infinite'
## [[ one more in ./abIndex.R ]]
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.infinite", signature(x = "ddiMatrix"),
          function(x) {
              r <- new("ldiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if(x@diag == "N") is.infinite(x@x) else logical(d[1L])
              r
          })

setMethod("is.infinite", signature(x = "indMatrix"),
          allFalseMatrix)

setMethod("is.infinite", signature(x = "nMatrix"),
          allFalseMatrix)

setMethod("is.infinite", signature(x = "lMatrix"),
          allFalseMatrix)

setMethod("is.infinite", signature(x = "dsparseMatrix"),
          function(x) {
              if(any(is.infinite(x@x))) {
                  r <- .M2kind(x <- diagU2N(x), "l")
                  r@x <- is.infinite(x@x)
                  .M2kind(.drop0(r), "n")
              } else allFalseMatrix(x)
          })

.is.infinite.ge <- function(x) {
    if(any(i <- is.infinite(x@x)))
        new("ngeMatrix", Dim = x@Dim, Dimnames = x@Dimnames, x = i)
    else allFalseMatrix(x)
}

.is.infinite.tr <- function(x) {
    if(any(i <- is.infinite(x@x))) {
        d <- x@Dim
        k <- indTri(n = d[1L], upper = x@uplo != "U",
                    diag = x@diag != "N", packed = FALSE)
        i[k] <- FALSE
        if(any(i))
            new("ntrMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.infinite.tp <- function(x) {
    if(any(i <- is.infinite(x@x))) {
        d <- x@Dim
        if(x@diag != "N") {
            k <- indDiag(n = d[1L], upper = x@uplo == "U",
                         packed = TRUE)
            i[k] <- FALSE
        }
        if(any(i))
            new("ntpMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.infinite.sy <- function(x) {
    if(any(i <- is.infinite(x@x))) {
        d <- x@Dim
        k <- indTri(n = d[1L], upper = x@uplo != "U",
                    diag = FALSE, packed = FALSE)
        i[k] <- FALSE
        if(any(i))
            new("nsyMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.infinite.sp <- function(x) {
    if(any(i <- is.infinite(x@x)))
        new("nspMatrix", Dim = x@Dim, Dimnames = x@Dimnames,
            x = i, uplo = x@uplo)
    else allFalseMatrix(x)
}

for(.xx in c("ge", "tr", "tp", "sy", "sp"))
    setMethod("is.infinite", signature(x = paste0("d", .xx, "Matrix")),
              get(paste0(".is.infinite.", .xx),
                  mode = "function", inherits = FALSE))
rm(.is.infinite.ge, .is.infinite.tr, .is.infinite.tp,
   .is.infinite.sy, .is.infinite.sp,
   .xx)

setMethod("is.infinite", signature(x = "sparseVector"),
          function(x) new("nsparseVector", length = x@length,
                          i = x@i[is.infinite(x@x)]))

setMethod("is.infinite", signature(x = "nsparseVector"),
          function(x) new("nsparseVector", length = x@length))


## METHODS FOR GENERIC: is.nan
## NB: completely (!) parallel to 'is.infinite'
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.nan", signature(x = "ddiMatrix"),
          function(x) {
              r <- new("ldiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if(x@diag == "N") is.nan(x@x) else logical(d[1L])
              r
          })

setMethod("is.nan", signature(x = "indMatrix"),
          allFalseMatrix)

setMethod("is.nan", signature(x = "nMatrix"),
          allFalseMatrix)

setMethod("is.nan", signature(x = "lMatrix"),
          allFalseMatrix)

setMethod("is.nan", signature(x = "dsparseMatrix"),
          function(x) {
              if(any(is.nan(x@x))) {
                  r <- .M2kind(x <- diagU2N(x), "l")
                  r@x <- is.nan(x@x)
                  .M2kind(.drop0(r), "n")
              } else allFalseMatrix(x)
          })

.is.nan.ge <- function(x) {
    if(any(i <- is.nan(x@x)))
        new("ngeMatrix", Dim = x@Dim, Dimnames = x@Dimnames, x = i)
    else allFalseMatrix(x)
}

.is.nan.tr <- function(x) {
    if(any(i <- is.nan(x@x))) {
        d <- x@Dim
        k <- indTri(n = d[1L], upper = x@uplo != "U",
                    diag = x@diag != "N", packed = FALSE)
        i[k] <- FALSE
        if(any(i))
            new("ntrMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.nan.tp <- function(x) {
    if(any(i <- is.nan(x@x))) {
        d <- x@Dim
        if(x@diag != "N") {
            k <- indDiag(n = d[1L], upper = x@uplo == "U",
                         packed = TRUE)
            i[k] <- FALSE
        }
        if(any(i))
            new("ntpMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.nan.sy <- function(x) {
    if(any(i <- is.nan(x@x))) {
        d <- x@Dim
        k <- indTri(n = d[1L], upper = x@uplo != "U",
                    diag = FALSE, packed = FALSE)
        i[k] <- FALSE
        if(any(i))
            new("nsyMatrix", Dim = d, Dimnames = x@Dimnames,
                x = i, uplo = x@uplo)
        else allFalseMatrix(x)
    } else allFalseMatrix(x)
}

.is.nan.sp <- function(x) {
    if(any(i <- is.nan(x@x)))
        new("nspMatrix", Dim = x@Dim, Dimnames = x@Dimnames,
            x = i, uplo = x@uplo)
    else allFalseMatrix(x)
}

for(.xx in c("ge", "tr", "tp", "sy", "sp"))
    setMethod("is.nan", signature(x = paste0("d", .xx, "Matrix")),
              get(paste0(".is.nan.", .xx),
                  mode = "function", inherits = FALSE))
rm(.is.nan.ge, .is.nan.tr, .is.nan.tp,
   .is.nan.sy, .is.nan.sp,
   .xx)

setMethod("is.nan", signature(x = "sparseVector"),
          function(x) new("nsparseVector", length = x@length,
                          i = x@i[is.nan(x@x)]))

setMethod("is.nan", signature(x = "nsparseVector"),
          function(x) new("nsparseVector", length = x@length))
