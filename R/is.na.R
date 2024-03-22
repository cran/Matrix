## METHODS FOR GENERIC: anyNA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("anyNA", c(x = "denseMatrix"),
          function(x, recursive = FALSE) {
              cl <- .M.nonvirtual(x)
              if(substr(cl, 1L, 1L)  == "n")
                  return(FALSE)
              if((shape <- substr(cl, 2L, 2L)) == "g")
                  anyNA(x@x)
              else {
                  if(shape == "t" && x@diag != "N") {
                      x@diag <- "N"
                      if(anyNA(diag(x, names = FALSE)))
                          diag(x) <- TRUE
                  }
                  anyNA(pack(x)@x)
              }
          })

setMethod("anyNA", c(x = "sparseMatrix"),
          function(x, recursive = FALSE)
              .M.kind(x) != "n" && anyNA(x@x))

setMethod("anyNA", c(x = "diagonalMatrix"),
          function(x, recursive = FALSE)
              .M.kind(x) != "n" && length(y <- x@x) > 0L && anyNA(y))

setMethod("anyNA", c(x = "indMatrix"),
          function(x, recursive = FALSE)
              FALSE)

setMethod("anyNA", c(x = "sparseVector"),
          function(x, recursive = FALSE)
              .M.kind(x) != "n" && anyNA(x@x))


## METHODS FOR GENERIC: is.na
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.na", c(x = "denseMatrix"),
          function(x) {
              cl <- .M.nonvirtual(x)
              never <- substr(cl, 1L, 1L) == "n"
              substr(cl, 1L, 1L) <- "n"
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if((shape <- substr(cl, 2L, 2L)) != "g") {
                  r@uplo <- x@uplo
                  if(!never && shape == "t" && x@diag != "N") {
                      x@diag <- "N"
                      if(anyNA(diag(x, names = FALSE)))
                          diag(x) <- TRUE
                  }
              }
              r@x <- if(never)
                         logical(length(x@x))
                     else is.na(x@x)
              r
          })

setMethod("is.na", c(x = "sparseMatrix"),
          function(x) {
              cl <- .M.nonvirtual(x)
              never <- substr(cl, 1L, 1L) == "n"
              substr(cl, 1L, 1L) <- if(never) "n" else "l"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if(substr(cl, 2L, 2L) != "g")
                  r@uplo <- x@uplo
              if(never) {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- integer(d[2L] + 1) },
                         "R" = { r@p <- integer(d[1L] + 1) })
                  r
              } else {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- x@p; r@i <- x@i },
                         "R" = { r@p <- x@p; r@j <- x@j },
                         "T" = { r@i <- x@i; r@j <- x@j })
                  r@x <- is.na(x@x)
                  .M2kind(.drop0(r), "n")
              }
          })

setMethod("is.na", c(x = "diagonalMatrix"),
          function(x) {
              r <- new("ndiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if(x@diag != "N" || .M.kind(x) == "n")
                         logical(d[1L])
                     else is.na(x@x)
              r
          })

setMethod("is.na", c(x = "indMatrix"),
          function(x) {
              m <- x@margin
              r <- new(if(m == 1L) "ngRMatrix" else "ngCMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@p <- integer(d[m] + 1)
              r
          })

setMethod("is.na", c(x = "sparseVector"),
          function(x) {
              r <- new("nsparseVector")
              r@length <- x@length
              if(.M.kind(x) != "n")
                  r@i <- x@i[is.na(x@x)]
              r
          })


## METHODS FOR GENERIC: is.nan
## NB: mostly parallel to is.na, completely parallel to is.infinite
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.nan", c(x = "denseMatrix"),
          function(x) {
              cl <- .M.nonvirtual(x)
              never <- switch(substr(cl, 1L, 1L), "d" = , "z" = FALSE, TRUE)
              substr(cl, 1L, 1L) <- "n"
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if((shape <- substr(cl, 2L, 2L)) != "g") {
                  r@uplo <- x@uplo
                  if(!never && shape == "t" && x@diag != "N") {
                      x@diag <- "N"
                      if(any(is.nan(diag(x, names = FALSE))))
                          diag(x) <- TRUE
                  }
              }
              r@x <- if(never)
                         logical(length(x@x))
                     else is.nan(x@x)
              r
          })

setMethod("is.nan", c(x = "sparseMatrix"),
          function(x) {
              cl <- .M.nonvirtual(x)
              never <- switch(substr(cl, 1L, 1L), "d" = , "z" = FALSE, TRUE)
              substr(cl, 1L, 1L) <- if(never) "n" else "l"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if(substr(cl, 2L, 2L) != "g")
                  r@uplo <- x@uplo
              if(never) {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- integer(d[2L] + 1) },
                         "R" = { r@p <- integer(d[1L] + 1) })
                  r
              } else {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- x@p; r@i <- x@i },
                         "R" = { r@p <- x@p; r@j <- x@j },
                         "T" = { r@i <- x@i; r@j <- x@j })
                  r@x <- is.nan(x@x)
                  .M2kind(.drop0(r), "n")
              }
          })

setMethod("is.nan", c(x = "diagonalMatrix"),
          function(x) {
              r <- new("ndiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if(x@diag != "N")
                         logical(d[1L])
                     else switch(.M.kind(x), "d" = , "z" = is.nan(x@x),
                                 logical(d[1L]))
              r
          })

setMethod("is.nan", c(x = "indMatrix"),
          function(x) {
              m <- x@margin
              r <- new(if(m == 1L) "ngRMatrix" else "ngCMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@p <- integer(d[m] + 1)
              r
          })

setMethod("is.nan", c(x = "sparseVector"),
          function(x) {
              r <- new("nsparseVector")
              r@length <- x@length
              switch(.M.kind(x), "d" = , "z" = { r@i <- x@i[is.nan(x@x)] })
              r
          })


## METHODS FOR GENERIC: is.infinite
## NB: mostly parallel to is.na, completely parallel to is.nan
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.infinite", c(x = "denseMatrix"),
          function(x) {
              cl <- .M.nonvirtual(x)
              never <- switch(substr(cl, 1L, 1L), "d" = , "z" = FALSE, TRUE)
              substr(cl, 1L, 1L) <- "n"
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if((shape <- substr(cl, 2L, 2L)) != "g") {
                  r@uplo <- x@uplo
                  if(!never && shape == "t" && x@diag != "N") {
                      x@diag <- "N"
                      if(any(is.infinite(diag(x, names = FALSE))))
                          diag(x) <- TRUE
                  }
              }
              r@x <- if(never)
                         logical(length(x@x))
                     else is.infinite(x@x)
              r
          })

setMethod("is.infinite", c(x = "sparseMatrix"),
          function(x) {
              cl <- .M.nonvirtual(x)
              never <- switch(substr(cl, 1L, 1L), "d" = , "z" = FALSE, TRUE)
              substr(cl, 1L, 1L) <- if(never) "n" else "l"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if(substr(cl, 2L, 2L) != "g")
                  r@uplo <- x@uplo
              if(never) {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- integer(d[2L] + 1) },
                         "R" = { r@p <- integer(d[1L] + 1) })
                  r
              } else {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- x@p; r@i <- x@i },
                         "R" = { r@p <- x@p; r@j <- x@j },
                         "T" = { r@i <- x@i; r@j <- x@j })
                  r@x <- is.infinite(x@x)
                  .M2kind(.drop0(r), "n")
              }
          })

setMethod("is.infinite", c(x = "diagonalMatrix"),
          function(x) {
              r <- new("ndiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if(x@diag != "N")
                         logical(d[1L])
                     else switch(.M.kind(x), "d" = , "z" = is.infinite(x@x),
                                 logical(d[1L]))
              r
          })

setMethod("is.infinite", c(x = "indMatrix"),
          function(x) {
              m <- x@margin
              r <- new(if(m == 1L) "ngRMatrix" else "ngCMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@p <- integer(d[m] + 1)
              r
          })

setMethod("is.infinite", c(x = "sparseVector"),
          function(x) {
              r <- new("nsparseVector")
              r@length <- x@length
              switch(.M.kind(x), "d" = , "z" = { r@i <- x@i[is.infinite(x@x)] })
              r
          })


## METHODS FOR GENERIC: is.finite
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.finite", c(x = "denseMatrix"),
          function(x) {
              cl <- .M.nonvirtual(x)
              always <- substr(cl, 1L, 1L) == "n"
              packed <- substr(cl, 3L, 3L) == "p"
              if((shape <- substr(cl, 2L, 2L)) != "s")
                  r <- new("ngeMatrix")
              else {
                  r <- new(if(!packed) "nsyMatrix" else "nspMatrix")
                  r@uplo <- x@uplo
              }
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <-
              if(shape != "t") {
                  if(always)
                      rep.int(TRUE, length(x@x))
                  else is.finite(x@x)
              } else {
                  if(always)
                      rep.int(TRUE, prod(d))
                  else if(!packed) {
                      tmp <- is.finite(x@x)
                      tmp[indTri(d[1L], x@uplo != "U", x@diag != "N", FALSE)] <-
                          TRUE
                      tmp
                  } else {
                      tmp <- rep.int(TRUE, prod(d))
                      tmp[indTri(d[1L], x@uplo == "U",          TRUE, FALSE)] <-
                          is.finite(x@x)
                      if(x@diag != "N") {
                          dim(tmp) <- d
                          diag(tmp) <- TRUE
                          dim(tmp) <- NULL
                      }
                      tmp
                  }
              }
              r
          })

setMethod("is.finite", c(x = "sparseMatrix"),
          function(x) {
              cl <- .M.nonvirtual(x)
              always <- substr(cl, 1L, 1L) == "n"
              if(substr(cl, 2L, 2L) != "s")
                  r <- new("ngeMatrix")
              else {
                  r <- new("nsyMatrix")
                  r@uplo <- x@uplo
              }
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              tmp <- rep.int(TRUE, prod(d))
              if(!always && !all(k <- is.finite(x@x))) {
                  if(substr(cl, 3L, 3L) != "T") {
                      x <- .M2T(x)
                      if(length(k) > length(x@x)) # was overallocated
                          k <- is.finite(x@x)
                  }
                  i <- c(x@i, x@j) + 1L
                  dim(i) <- c(length(k), 2L)
                  dim(tmp) <- d
                  tmp[i] <- k
                  dim(tmp) <- NULL
              }
              r@x <- tmp
              r
          })

setMethod("is.finite", c(x = "diagonalMatrix"),
          function(x) {
              r <- new("nsyMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              tmp <- rep.int(TRUE, prod(d))
              if(x@diag == "N" && .M.kind(x) != "n" && !all(k <- is.finite(x@x))) {
                  dim(tmp) <- d
                  diag(tmp) <- k
                  dim(tmp) <- NULL
              }
              r@x <- tmp
              r
          })

setMethod("is.finite", c(x = "indMatrix"),
          function(x)  {
              r <- new("ngeMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- rep.int(TRUE, prod(d))
              r
          })

setMethod("is.finite", c(x = "sparseVector"),
          function(x)  {
              r <- rep.int(TRUE, x@length)
              if(.M.kind(x) != "n")
                  r[x@i[!is.finite(x@x)]] <- FALSE
              r
          })
