matmultDim <- function(d.a, d.b, type = 1L) {
    ## Return the 'dim' of the product indicated by 'type':
    ##     type 1:    a  %*%   b
    ##          2:  t(a) %*%   b   { crossprod}
    ##          3:    a  %*% t(b)  {tcrossprod}
    ## after asserting that ncol(<left operand>) == nrow(<right operand>)
    i.a <- 1L + (type != 2L)
    i.b <- 1L + (type == 3L)
    if(d.a[i.a] != d.b[i.b])
        stop("non-conformable arguments")
    c(d.a[-i.a], d.b[-i.b])
}

matmultDN <- function(dn.a, dn.b, type = 1L)
    ## Return the 'dimnames' of the product indicated by 'type':
    ##     type 1:    a  %*%   b
    ##          2:  t(a) %*%   b   { crossprod}
    ##          3:    a  %*% t(b)  {tcrossprod}
    c(if(is.null(dn.a)) list(NULL) else dn.a[2L - (type != 2L)],
      if(is.null(dn.b)) list(NULL) else dn.b[2L - (type == 3L)])


## METHODS FOR GENERIC: %*%
## NB: x %*% y == t(t(y) %*% t(x))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(.cl in c("Matrix", "sparseVector")) {
setMethod("%*%", signature(x = .cl, y = "ANY"),
          function(x, y)
              x %*% (if(length(dim(y)) == 2L) as.matrix else as.vector)(y))

setMethod("%*%", signature(x = "ANY", y = .cl),
          function(x, y)
              (if(length(dim(x)) == 2L) as.matrix else as.vector)(x) %*% y)
}


## .... denseMatrix ....................................................

setMethod("%*%", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y)
              .Call(R_dense_matmult, x, y, FALSE, FALSE))

for(.cl in c("matrix", "vector")) {
setMethod("%*%", signature(x = "denseMatrix", y = .cl),
          function(x, y)
              .Call(R_dense_matmult, x, y, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "denseMatrix"),
          function(x, y)
              .Call(R_dense_matmult, x, y, FALSE, FALSE))
}


## .... CsparseMatrix ..................................................

setMethod("%*%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = "CsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = "CsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE, FALSE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%*%", signature(x = "CsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE, FALSE))
}


## .... RsparseMatrix ..................................................

setMethod("%*%", signature(x = "RsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE, FALSE))

setMethod("%*%", signature(x = "RsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE, FALSE))

setMethod("%*%", signature(x = "RsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE, FALSE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%*%", signature(x = "RsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE, FALSE))
}


## .... TsparseMatrix ..................................................

setMethod("%*%", signature(x = "TsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = "TsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE, FALSE))

setMethod("%*%", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE, FALSE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%*%", signature(x = "TsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE, FALSE))
}


## .... diagonalMatrix .................................................

setMethod("%*%", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
          function(x, y) {
              r <- new("ddiMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              xu <- x@diag != "N"
              yu <- y@diag != "N"
              if(xu && yu)
                  r@diag <- "U"
              else {
                  if(!xu) {
                      xii <- x@x
                      if(.M.kind(x) == "n" && anyNA(xii))
                          xii <- xii | is.na(xii)
                  }
                  if(!yu) {
                      yii <- y@x
                      if(.M.kind(y) == "n" && anyNA(yii))
                          yii <- yii | is.na(yii)
                  }
                  r@x <-
                      if(xu)
                          as.double(yii)
                      else if(yu)
                          as.double(xii)
                      else as.double(xii * yii)
              }
              r
          })

for(.cl in c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix",
             "denseMatrix", "matrix", "vector")) {
setMethod("%*%", signature(x = "diagonalMatrix", y = .cl),
          function(x, y)
              .Call(R_diagonal_matmult, x, y, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "diagonalMatrix"),
          function(x, y)
              .Call(R_diagonal_matmult, x, y, FALSE, FALSE, FALSE))
}


## .... indMatrix ......................................................

setMethod("%*%", signature(x = "indMatrix", y = "indMatrix"),
          function(x, y) {
              mx <- x@margin
              my <- y@margin
              px <- x@perm
              py <- y@perm
              r <- new(if(mx == my)
                           "indMatrix"
                       else if(mx == 1L)
                           "dgeMatrix"
                       else "dgTMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              if(mx == my)
                  r@perm <- if(mx == 1L) py[px] else { r@margin <- 2L; px[py] }
              else if(mx == 1L)
                  r@x <- as.double(px == rep(py, each = length(px)))
              else {
                  r@i <- px - 1L
                  r@j <- py - 1L
                  r@x <- rep.int(1, length(px))
              }
              r
          })

setMethod("%*%", signature(x = "indMatrix", y = "Matrix"),
          function(x, y) {
              if(x@margin != 1L)
                  return(.M2kind(x, "d") %*% y)
              matmultDim(x@Dim, y@Dim, type = 1L)
              r <- .M2kind(y[x@perm, , drop = FALSE], ",")
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%*%", signature(x = "Matrix", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %*% .M2kind(y, "d"))
              matmultDim(x@Dim, y@Dim, type = 1L)
              r <- .M2kind(x[, y@perm, drop = FALSE], ",")
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%*%", signature(x = "indMatrix", y = "matrix"),
          function(x, y) {
              if(x@margin != 1L)
                  return(.M2kind(x, "d") %*% y)
              matmultDim(x@Dim, dim(y), type = 1L)
              r <- .m2dense(y[x@perm, , drop = FALSE], ",ge")
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%*%", signature(x = "matrix", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %*% .M2kind(y, "d"))
              matmultDim(dim(x), y@Dim, type = 1L)
              r <- .m2dense(x[, y@perm, drop = FALSE], ",ge")
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%*%", signature(x = "indMatrix", y = "vector"),
          function(x, y) {
              if(x@margin != 1L)
                  return(.M2kind(x, "d") %*% y)
              k <- (d <- x@Dim)[2L]
              r <-
              if(k == length(y))
                  .m2dense(y[x@perm], ",ge")
              else if(k == 1L)
                  .m2dense(matrix(y, d[1L], length(y), byrow = TRUE), ",ge")
              else stop("non-conformable arguments")
              r@Dimnames <- c(x@Dimnames[1L], list(NULL))
              r
          })

setMethod("%*%", signature(x = "vector", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %*% .M2kind(y, "d"))
              k <- (d <- y@Dim)[1L]
              r <-
              if(k == length(x))
                  .m2dense(x[y@perm], ",ge", trans = TRUE)
              else if(k == 1L)
                  .m2dense(matrix(x, length(x), d[2L]), ",ge")
              else stop("non-conformable arguments")
              r@Dimnames <- c(list(NULL), y@Dimnames[2L])
              r
          })


## .... pMatrix ........................................................

setMethod("%*%", signature(x = "pMatrix", y = "pMatrix"),
          function(x, y) {
              r <- new("pMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
              if(y@margin == 1L)
                  y@perm[if(x@margin == 1L) x@perm else invertPerm(x@perm)]
              else {
                  r@margin <- 2L
                  (if(x@margin == 1L) invertPerm(x@perm) else x@perm)[y@perm]
              }
              r
          })

setMethod("%*%", signature(x = "pMatrix", y = "indMatrix"),
          function(x, y) {
              r <- new("indMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
              if(y@margin == 1L)
                  y@perm[if(x@margin == 1L) x@perm else invertPerm(x@perm)]
              else {
                  r@margin <- 2L
                  (if(x@margin == 1L) invertPerm(x@perm) else x@perm)[y@perm]
              }
              r
          })

setMethod("%*%", signature(x = "indMatrix", y = "pMatrix"),
          function(x, y) {
              r <- new("indMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
              if(x@margin == 1L)
                  (if(y@margin == 1L) y@perm else invertPerm(y@perm))[x@perm]
              else {
                  r@margin <- 2L
                  x@perm[if(y@margin == 1L) invertPerm(x@perm) else y@perm]
              }
              r
          })

setMethod("%*%", signature(x = "pMatrix", y = "Matrix"),
          function(x, y) {
              matmultDim(x@Dim, y@Dim, type = 1L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .M2kind(y[perm, , drop = FALSE], ",")
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%*%", signature(x = "Matrix", y = "pMatrix"),
          function(x, y) {
              matmultDim(x@Dim, y@Dim, type = 1L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .M2kind(x[, perm, drop = FALSE], ",")
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%*%", signature(x = "pMatrix", y = "matrix"),
          function(x, y) {
              matmultDim(x@Dim, dim(y), type = 1L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .m2dense(y[perm, , drop = FALSE], ",ge")
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%*%", signature(x = "matrix", y = "pMatrix"),
          function(x, y) {
              matmultDim(dim(x), y@Dim, type = 1L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .m2dense(x[, perm, drop = FALSE], ",ge")
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%*%", signature(x = "pMatrix", y = "vector"),
          function(x, y) {
              k <- x@Dim[2L]
              r <-
              if(k == length(y)) {
                  perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
                  .m2dense(y[perm], ",ge")
              }
              else if(k == 1L)
                  .m2dense(y, ",ge", trans = TRUE)
              else stop("non-conformable arguments")
              r@Dimnames <- c(x@Dimnames[1L], list(NULL))
              r
          })

setMethod("%*%", signature(x = "vector", y = "pMatrix"),
          function(x, y) {
              k <- y@Dim[1L]
              r <-
              if(k == length(x)) {
                  perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
                  .m2dense(x[perm], ",ge", trans = TRUE)
              }
              else if(k == 1L)
                  .m2dense(x, ",ge")
              else stop("non-conformable arguments")
              r@Dimnames <- c(list(NULL), y@Dimnames[2L])
              r
          })


## .... sparseVector ...................................................

setMethod("%*%", signature(x = "sparseVector", y = "sparseVector"),
          function(x, y)
              if((nx <- length(x)) == 1L)
                  .tCRT(.V2C(x * y))
              else if(nx == length(y))
           ## else if((ny <- length(y)) == 1L)
           ##           .V2C(x * y)
           ## else if(nx == ny)
                  .m2sparse(sum(x * y), ",gR")
              else stop("non-conformable arguments"))

for(.cl in c("Matrix", "matrix")) {
setMethod("%*%", signature(x = "sparseVector", y = .cl),
          function(x, y)
              if((k <- dim(y)[1L]) == length(x))
                  .tCRT(.V2C(x)) %*% y
              else if(k == 1L)
                        .V2C(x)  %*% y
              else stop("non-conformable arguments"))

setMethod("%*%", signature(x = .cl, y = "sparseVector"),
          function(x, y)
              if((k <- dim(x)[2L]) == length(y))
                  x %*%       .V2C(y)
              else if(k == 1L)
                  x %*% .tCRT(.V2C(y))
              else stop("non-conformable arguments"))
}

setMethod("%*%", signature(x = "sparseVector", y = "vector"),
          function(x, y)
              if((nx <- length(x)) == 1L)
                  .m2dense(.V2v(x * y), ",ge", trans = TRUE)
              else if(nx == length(y))
           ## else if((ny <- length(y)) == 1L)
           ##     .m2dense(.V2v(x * y), ",ge")
           ## else if(nx == ny)
                  .m2dense(sum(x * y), ",ge")
              else stop("non-conformable arguments"))

setMethod("%*%", signature(x = "vector", y = "sparseVector"),
          function(x, y)
              if((nx <- length(x)) == 1L)
                  .m2dense(.V2v(x * y), ",ge", trans = TRUE)
              else if(nx == length(y))
           ## else if((ny <- length(y)) == 1L)
           ##     .m2dense(.V2v(x * y), ",ge")
           ## else if(nx == ny)
                  .m2dense(sum(x * y), ",ge")
              else stop("non-conformable arguments"))


## METHODS FOR GENERIC: %&%
## NB: x %*% y == t(t(y) %*% t(x))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("%&%", signature(x = "ANY", y = "ANY"),
          function(x, y)
              (if(length(dim(x)) == 2L) as.matrix else as.vector)(x) %&%
              (if(length(dim(y)) == 2L) as.matrix else as.vector)(y))

for(.cl in c("Matrix", "sparseVector", "matrix", "vector")) {
setMethod("%&%", signature(x = .cl, y = "ANY"),
          function(x, y)
              x %&% (if(length(dim(y)) == 2L) as.matrix else as.vector)(y))

setMethod("%&%", signature(x = "ANY", y = .cl),
          function(x, y)
              (if(length(dim(x)) == 2L) as.matrix else as.vector)(x) %&% y)
}


setMethod("%&%", signature(x = "matrix", y = "matrix"),
          function(x, y)
              .m2sparse(x, "ngC") %&% .m2sparse(y, "ngC"))

setMethod("%&%", signature(x = "matrix", y = "vector"),
          function(x, y)
              .m2sparse(x, "ngC") %&%           y        )

setMethod("%&%", signature(x = "vector", y = "matrix"),
          function(x, y)
                        x         %&% .m2sparse(y, "ngC"))

setMethod("%&%", signature(x = "vector", y = "vector"),
          function(x, y) {
              r <-
              if((nx <- length(x)) == 1L)
                  .m2sparse(x, "ngC") %&% .m2sparse(y, "ngR", trans = TRUE)
              else if(nx == length(y))
           ## else if((ny <- length(y)) == 1L || nx == ny)
                  .m2sparse(x, "ngR", trans = TRUE) %&% .m2sparse(y, "ngC")
              else stop("non-conformable arguments")
              r@Dimnames <- list(NULL, NULL)
              r
          })


## .... denseMatrix ....................................................

setMethod("%&%", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

for(.cl in c("matrix", "vector")) {
setMethod("%&%", signature(x = "denseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "denseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))
}


## .... CsparseMatrix ..................................................

setMethod("%&%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = "CsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = "CsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%&%", signature(x = "CsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))
}


## .... RsparseMatrix ..................................................

setMethod("%&%", signature(x = "RsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

setMethod("%&%", signature(x = "RsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

setMethod("%&%", signature(x = "RsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%&%", signature(x = "RsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE,  TRUE))
}


## .... TsparseMatrix ..................................................

setMethod("%&%", signature(x = "TsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = "TsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

setMethod("%&%", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%&%", signature(x = "TsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_matmult, x, y, FALSE, FALSE, FALSE,  TRUE))
}


## .... diagonalMatrix .................................................

setMethod("%&%", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
          function(x, y) {
              r <- new("ndiMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              xu <- x@diag != "N"
              yu <- y@diag != "N"
              if(xu && yu)
                  r@diag <- "U"
              else
                  r@x <-
                      if(xu)
                          as.logical(y@x)
                      else if(yu)
                          as.logical(x@x)
                      else x@x & y@x
              r
          })

for(.cl in c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix",
             "denseMatrix", "matrix", "vector")) {
setMethod("%&%", signature(x = "diagonalMatrix", y = .cl),
          function(x, y)
              .Call(R_diagonal_matmult, x, y, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "diagonalMatrix"),
          function(x, y)
              .Call(R_diagonal_matmult, x, y, FALSE, FALSE,  TRUE))
}


## .... indMatrix ......................................................

setMethod("%&%", signature(x = "indMatrix", y = "indMatrix"),
          function(x, y) {
              mx <- x@margin
              my <- y@margin
              px <- x@perm
              py <- y@perm
              r <- new(if(mx == my)
                           "indMatrix"
                       else if(mx == 1L)
                           "ngeMatrix"
                       else "ngTMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              if(mx == my)
                  r@perm <- if(mx == 1L) py[px] else { r@margin <- 2L; px[py] }
              else if(mx == 1L)
                  r@x <- px == rep(py, each = length(px))
              else {
                  r@i <- px - 1L
                  r@j <- py - 1L
              }
              r
          })

setMethod("%&%", signature(x = "indMatrix", y = "Matrix"),
          function(x, y) {
              if(x@margin != 1L)
                  return(.M2kind(x, "n") %&% y)
              matmultDim(x@Dim, y@Dim, type = 1L)
              r <- .M2kind(y[x@perm, , drop = FALSE], "n")
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "Matrix", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %&% .M2kind(y, "n"))
              matmultDim(x@Dim, y@Dim, type = 1L)
              r <- .M2kind(x[, y@perm, drop = FALSE], "n")
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "indMatrix", y = "matrix"),
          function(x, y) {
              if(x@margin != 1L)
                  return(.M2kind(x, "n") %&% y)
              matmultDim(x@Dim, dim(y), type = 1L)
              r <- .m2dense(y[x@perm, , drop = FALSE], "nge")
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "matrix", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %&% .M2kind(y, "n"))
              matmultDim(dim(x), y@Dim, type = 1L)
              r <- .m2dense(x[, y@perm, drop = FALSE], "nge")
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "indMatrix", y = "vector"),
          function(x, y) {
              if(x@margin != 1L)
                  return(.M2kind(x, "n") %&% y)
              k <- (d <- x@Dim)[2L]
              r <-
              if(k == length(y))
                  .m2dense(y[x@perm], "nge")
              else if(k == 1L)
                  .m2dense(matrix(y, d[1L], length(y), byrow = TRUE), "nge")
              else stop("non-conformable arguments")
              r@Dimnames <- c(x@Dimnames[1L], list(NULL))
              r
          })

setMethod("%&%", signature(x = "vector", y = "indMatrix"),
          function(x, y) {
              if(y@margin == 1L)
                  return(x %&% .M2kind(y, "n"))
              k <- (d <- y@Dim)[1L]
              r <-
              if(k == length(x))
                  .m2dense(x[y@perm], "nge", trans = TRUE)
              else if(k == 1L)
                  .m2dense(matrix(x, length(x), d[2L]), "nge")
              else stop("non-conformable arguments")
              r@Dimnames <- c(list(NULL), y@Dimnames[2L])
              r
          })


## .... pMatrix ........................................................

setMethod("%&%", signature(x = "pMatrix", y = "pMatrix"),
          function(x, y) {
              r <- new("pMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
              if(y@margin == 1L)
                  y@perm[if(x@margin == 1L) x@perm else invertPerm(x@perm)]
              else {
                  r@margin <- 2L
                  (if(x@margin == 1L) invertPerm(x@perm) else x@perm)[y@perm]
              }
              r
          })

setMethod("%&%", signature(x = "pMatrix", y = "indMatrix"),
          function(x, y) {
              r <- new("indMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
              if(y@margin == 1L)
                  y@perm[if(x@margin == 1L) x@perm else invertPerm(x@perm)]
              else {
                  r@margin <- 2L
                  (if(x@margin == 1L) invertPerm(x@perm) else x@perm)[y@perm]
              }
              r
          })

setMethod("%&%", signature(x = "indMatrix", y = "pMatrix"),
          function(x, y) {
              r <- new("indMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 1L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 1L)
              r@perm <-
              if(x@margin == 1L)
                  (if(y@margin == 1L) y@perm else invertPerm(y@perm))[x@perm]
              else {
                  r@margin <- 2L
                  x@perm[if(y@margin == 1L) invertPerm(x@perm) else y@perm]
              }
              r
          })

setMethod("%&%", signature(x = "pMatrix", y = "Matrix"),
          function(x, y) {
              matmultDim(x@Dim, y@Dim, type = 1L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .M2kind(y[perm, , drop = FALSE], "n")
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "Matrix", y = "pMatrix"),
          function(x, y) {
              matmultDim(x@Dim, y@Dim, type = 1L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .M2kind(x[, perm, drop = FALSE], "n")
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "pMatrix", y = "matrix"),
          function(x, y) {
              matmultDim(x@Dim, dim(y), type = 1L)
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .m2dense(y[perm, , drop = FALSE], "nge")
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 1L)
              r
          })

setMethod("%&%", signature(x = "matrix", y = "pMatrix"),
          function(x, y) {
              matmultDim(dim(x), y@Dim, type = 1L)
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .m2dense(x[, perm, drop = FALSE], "nge")
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 1L)
              r
          })

setMethod("%&%", signature(x = "pMatrix", y = "vector"),
          function(x, y) {
              k <- x@Dim[2L]
              r <-
              if(k == length(y)) {
                  perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
                  .m2dense(y[perm], "nge")
              }
              else if(k == 1L)
                  .m2dense(y, "nge", trans = TRUE)
              else stop("non-conformable arguments")
              r@Dimnames <- c(x@Dimnames[1L], list(NULL))
              r
          })

setMethod("%&%", signature(x = "vector", y = "pMatrix"),
          function(x, y) {
              k <- y@Dim[1L]
              r <-
              if(k == length(x)) {
                  perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
                  .m2dense(x[perm], "nge", trans = TRUE)
              }
              else if(k == 1L)
                  .m2dense(x, "nge")
              else stop("non-conformable arguments")
              r@Dimnames <- c(list(NULL), y@Dimnames[2L])
              r
          })


## .... sparseVector ...................................................

setMethod("%&%", signature(x = "sparseVector", y = "sparseVector"),
          function(x, y) {
              x <- .V2kind(.drop0(x, isM = FALSE), "n")
              y <- .V2kind(.drop0(y, isM = FALSE), "n")
              if((nx <- length(x)) == 1L) {
                  if(length(x@i) == 0L)
                      y@i <- integer(0L)
                  .tCRT(.V2C(y))
              } else if(nx == length(y))
           ## } else if((ny <- length(y)) == 1L) {
           ##     if(length(y@i))
           ##         x@i <- integer(0L)
           ##           .V2C(x)
           ## } else if(nx == ny)
                  .m2sparse(any(match(x@i, y@i, 0L)), "ngR")
              else stop("non-conformable arguments")
          })

for(.cl in c("Matrix", "matrix")) {
setMethod("%&%", signature(x = "sparseVector", y = .cl),
          function(x, y) {
              x <- .V2kind(.drop0(x, isM = FALSE), "n")
              if((k <- dim(y)[1L]) == length(x))
                  .tCRT(.V2C(x)) %&% y
              else if(k == 1L)
                        .V2C(x)  %&% y
              else stop("non-conformable arguments")
          })

setMethod("%&%", signature(x = .cl, y = "sparseVector"),
          function(x, y) {
              y <- .V2kind(.drop0(y, isM = FALSE), "n")
              if((k <- dim(x)[2L]) == length(y))
                  x %&%       .V2C(y)
              else if(k == 1L)
                  x %&% .tCRT(.V2C(y))
              else stop("non-conformable arguments")
          })
}

setMethod("%&%", signature(x = "sparseVector", y = "vector"),
          function(x, y)
              .V2kind(.drop0(x, isM = FALSE), "n") %&% .m2V(y, "n"))

setMethod("%&%", signature(x = "vector", y = "sparseVector"),
          function(x, y)
              .m2V(x, "n") %&% .V2kind(.drop0(y, isM = FALSE), "n"))


## METHODS FOR GENERIC: crossprod
## NB: t(x) %*% y == t(t(y) %*% x)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(.cl in c("Matrix", "sparseVector")) {
setMethod("crossprod", signature(x = .cl, y = "ANY"),
          function(x, y = NULL, ...)
              crossprod(x, (if(length(dim(y)) == 2L) as.matrix else as.vector)(y), ...))

setMethod("crossprod", signature(x = "ANY", y = .cl),
          function(x, y = NULL, ...)
              crossprod((if(length(dim(x)) == 2L) as.matrix else as.vector)(x), y, ...))
}


## .... denseMatrix ....................................................

setMethod("crossprod", signature(x = "denseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" else boolArith)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE,  TRUE)
              else
              .Call(R_dense_matmult, x, y,  TRUE, FALSE))

setMethod("crossprod", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" && .M.kind(y) == "n" else boolArith)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE,  TRUE)
              else
              .Call(R_dense_matmult, x, y,  TRUE, FALSE))

for(.cl in c("matrix", "vector")) {
setMethod("crossprod", signature(x = "denseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE,  TRUE)
              else
              .Call(R_dense_matmult, x, y,  TRUE, FALSE))

setMethod("crossprod", signature(x = .cl, y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE,  TRUE)
              else
              .Call(R_dense_matmult, x, y,  TRUE, FALSE))
}


## .... CsparseMatrix ..................................................

setMethod("crossprod", signature(x = "CsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("crossprod", signature(x = "CsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = .cl, y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  TRUE, FALSE,  TRUE, boolArith))
}


## .... RsparseMatrix ..................................................

setMethod("crossprod", signature(x = "RsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE,  TRUE, boolArith))

setMethod("crossprod", signature(x = "RsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  TRUE, FALSE,  TRUE, boolArith))

setMethod("crossprod", signature(x = "RsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  TRUE, FALSE,  TRUE, boolArith))

setMethod("crossprod", signature(x = "RsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  TRUE, FALSE,  TRUE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("crossprod", signature(x = "RsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = .cl, y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  TRUE, FALSE,  TRUE, boolArith))
}


## .... TsparseMatrix ..................................................

setMethod("crossprod", signature(x = "TsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "TsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "TsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("crossprod", signature(x = "TsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = .cl, y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  TRUE, FALSE,  TRUE, boolArith))
}


## .... diagonalMatrix .................................................

setMethod("crossprod", signature(x = "diagonalMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              boolArith <- !is.na(boolArith) && boolArith
              r <- new(if(boolArith) "ndiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames[c(2L, 2L)]
              if(x@diag != "N")
                  r@diag <- "U"
              else {
                  xii <- x@x
                  r@x <-
                      if(boolArith)
                          as.logical(xii)
                      else {
                          if(.M.kind(x) == "n" && anyNA(xii))
                              xii <- xii | is.na(xii)
                          as.double(xii * xii)
                      }
              }
              r
          })

setMethod("crossprod", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(t(x), y))

for(.cl in c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix",
             "denseMatrix", "matrix", "vector")) {
setMethod("crossprod", signature(x = "diagonalMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_diagonal_matmult, x, y,  TRUE, FALSE, boolArith))

setMethod("crossprod", signature(x = .cl, y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_diagonal_matmult, x, y,  TRUE, FALSE, boolArith))
}


## .... indMatrix ......................................................

setMethod("crossprod", signature(x = "indMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              if(x@margin != 1L)
                  return(tcrossprod(t(x), boolArith = boolArith, ...))
              boolArith <- !is.na(boolArith) && boolArith
              tt <- tabulate(x@perm, x@Dim[2L])
              r <- new(if(boolArith) "ndiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim[c(2L, 2L)]
              r@Dimnames <- x@Dimnames[c(2L, 2L)]
              r@x <- if(boolArith) as.logical(tt) else as.double(tt)
              r
          })

for(.cl in c("Matrix", "matrix", "vector"))
setMethod("crossprod", signature(x = "indMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(t(x), y))

setMethod("crossprod", signature(x = "Matrix", y = "indMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              matmultDim(x@Dim, y@Dim, type = 2L)
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              if(y@margin == 1L)
                  r <- crossprod(x, .M2kind(y, l), boolArith = boolArith, ...)
              else {
                  r <- .M2kind(t(x)[, y@perm, drop = FALSE], l)
                  r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 2L)
              }
              r
          })

setMethod("crossprod", signature(x = "matrix", y = "indMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              matmultDim(dim(x), y@Dim, type = 2L)
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              if(y@margin == 1L)
                  r <- crossprod(x, .M2kind(y, l), boolArith = boolArith, ...)
              else {
                  r <- .m2dense(t(x)[, y@perm, drop = FALSE], paste0(l, "ge"))
                  r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 2L)
              }
              r
          })

setMethod("crossprod", signature(x = "vector", y = "indMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              if(y@Dim[1L] != length(x))
                  stop("non-conformable arguments")
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              if(y@margin == 1L)
                  r <- crossprod(x, .M2kind(y, l), boolArith = boolArith, ...)
              else {
                  r <- .m2dense(x[y@perm], paste0(l, "ge"), trans = TRUE)
                  r@Dimnames <- c(list(NULL), y@Dimnames[2L])
              }
              r
          })


## .... pMatrix ........................................................

setMethod("crossprod", signature(x = "pMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              boolArith <- !is.na(boolArith) && boolArith
              r <- new(if(boolArith) "ndiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames[c(2L, 2L)]
              r@diag <- "U"
              r
          })

setMethod("crossprod", signature(x = "pMatrix", y = "pMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <- new("pMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 2L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 2L)
              r@perm <-
              if(y@margin == 1L)
                  y@perm[if(x@margin == 1L) invertPerm(x@perm) else x@perm]
              else {
                  r@margin <- 2L
                  (if(x@margin == 1L) x@perm else invertPerm(x@perm))[y@perm]
              }
              r
          })

setMethod("crossprod", signature(x = "Matrix", y = "pMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              matmultDim(x@Dim, y@Dim, type = 2L)
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .M2kind(t(x)[, perm, drop = FALSE], l)
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("crossprod", signature(x = "matrix", y = "pMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              matmultDim(dim(x), y@Dim, type = 2L)
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .m2dense(t(x)[, perm, drop = FALSE], paste0(l, "ge"))
              r@Dimnames <- matmultDN(dimnames(x), y@Dimnames, type = 2L)
              r
          })

setMethod("crossprod", signature(x = "vector", y = "pMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              if(y@Dim[1L] != length(x))
                  stop("non-conformable arguments")
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              perm <- if(y@margin == 1L) invertPerm(y@perm) else y@perm
              r <- .m2dense(x[perm], paste0(l, "ge"), trans = TRUE)
              r@Dimnames <- c(list(NULL), y@Dimnames[2L])
              r
          })


## .... sparseVector ...................................................

setMethod("crossprod", signature(x = "sparseVector", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" else boolArith) {
                  if(!is.na(boolArith))
                      x <- .V2kind(.drop0(x, isM = FALSE), "n")
                  .m2sparse(length(x@i) > 0L, "nsR")
              } else .m2sparse(sum(x * x), ",sR"))

setMethod("crossprod", signature(x = "sparseVector", y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" && .M.kind(y) == "n" else boolArith) {
                  if(!is.na(boolArith)) {
                      x <- .V2kind(.drop0(x, isM = FALSE), "n")
                      y <- .V2kind(.drop0(y, isM = FALSE), "n")
                  }
                  if((nx <- length(x)) == 1L) {
                      if(length(x@i) == 0L)
                          y@i <- integer(0L)
                      .tCRT(.V2C(y))
                  } else if(nx == length(y))
                      .m2sparse(any(match(x@i, y@i, 0L)), "ngR")
                  else stop("non-conformable arguments")
              } else {
                  if((nx <- length(x)) == 1L)
                      .tCRT(.V2C(x * y))
                  else if(nx == length(y))
                      .m2sparse(sum(x * y), ",gR")
                  else stop("non-conformable arguments")
              })

for(.cl in c("Matrix", "matrix")) {
setMethod("crossprod", signature(x = "sparseVector", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              crossprod(.tCRT(.V2C(x)), y, boolArith = boolArith, ...))

setMethod("crossprod", signature(x = .cl, y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(dim(x)[1L] == length(y))
              crossprod(x,       .V2C(y) , boolArith = boolArith, ...)
              else
              crossprod(x, .tCRT(.V2C(y)), boolArith = boolArith, ...))
}

setMethod("crossprod", signature(x = "sparseVector", y = "vector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith) {
                  x <- .V2kind(.drop0(x, isM = FALSE), "n")
                  y <- .m2V(y, "n")
                  if((nx <- length(x)) == 1L) {
                      if(length(x@i) == 0L)
                          y@i <- integer(0L)
                      .tCRT(.V2C(y))
                  } else if(nx == length(y))
                      .m2sparse(any(match(x@i, y@i, 0L)), "ngR")
                  else stop("non-conformable arguments")
              } else {
                  if((nx <- length(x)) == 1L)
                      .m2dense(.V2v(x * y), ",ge", trans = TRUE)
                  else if(nx == length(y))
                      .m2dense(sum(x * y), ",ge")
                  else stop("non-conformable arguments")
              })

setMethod("crossprod", signature(x = "vector", y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith) {
                  x <- .m2V(x, "n")
                  y <- .V2kind(.drop0(y, isM = FALSE), "n")
                  if((nx <- length(x)) == 1L) {
                      if(length(x@i) == 0L)
                          y@i <- integer(0L)
                      .tCRT(.V2C(y))
                  } else if(nx == length(y))
                      .m2sparse(any(match(x@i, y@i, 0L)), "ngR")
                  else stop("non-conformable arguments")
              } else {
                  if((nx <- length(x)) == 1L)
                      .m2dense(.V2v(x * y), ",ge", trans = TRUE)
                  else if(nx == length(y))
                      .m2dense(sum(x * y), ",ge")
                  else stop("non-conformable arguments")
              })


## METHODS FOR GENERIC: tcrossprod
## NB: x %*% t(y) == t(y %*% t(x))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(.cl in c("Matrix", "sparseVector")) {
setMethod("tcrossprod", signature(x = .cl, y = "ANY"),
          function(x, y = NULL, ...)
              tcrossprod(x, (if(length(dim(y)) == 2L) as.matrix else as.vector)(y), ...))

setMethod("tcrossprod", signature(x = "ANY", y = .cl),
          function(x, y = NULL, ...)
              tcrossprod((if(length(dim(x)) == 2L) as.matrix else as.vector)(x), y, ...))
}


## .... denseMatrix ....................................................

setMethod("tcrossprod", signature(x = "denseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" else boolArith)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE,  TRUE)
              else
              .Call(R_dense_matmult, x, y, FALSE,  TRUE))

setMethod("tcrossprod", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" && .M.kind(y) == "n" else boolArith)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE,  TRUE)
              else
              .Call(R_dense_matmult, x, y, FALSE,  TRUE))

for(.cl in c("matrix", "vector")) {
setMethod("tcrossprod", signature(x = "denseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE,  TRUE)
              else
              .Call(R_dense_matmult, x, y, FALSE,  TRUE))

setMethod("tcrossprod", signature(x = .cl, y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE,  TRUE)
              else
              .Call(R_dense_matmult, x, y, FALSE,  TRUE))
}


## .... CsparseMatrix ..................................................

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("tcrossprod", signature(x = "CsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = .cl, y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x, FALSE,  TRUE,  TRUE, boolArith))
}


## .... RsparseMatrix ..................................................

setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y,  FALSE,  TRUE,  TRUE, boolArith))

setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  FALSE,  TRUE,  TRUE, boolArith))

setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  FALSE,  TRUE,  TRUE, boolArith))

setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x,  FALSE,  TRUE,  TRUE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("tcrossprod", signature(x = "RsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = .cl, y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x, FALSE,  TRUE,  TRUE, boolArith))
}


## .... TsparseMatrix ..................................................

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("tcrossprod", signature(x = "TsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = .cl, y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_matmult, y, x, FALSE,  TRUE,  TRUE, boolArith))
}


## .... diagonalMatrix .................................................

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              boolArith <- !is.na(boolArith) && boolArith
              r <- new(if(boolArith) "ndiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames[c(1L, 1L)]
              if(x@diag != "N")
                  r@diag <- "U"
              else {
                  xii <- x@x
                  r@x <-
                      if(boolArith)
                          as.logical(xii)
                      else {
                          if(.M.kind(x) == "n" && anyNA(xii))
                              xii <- xii | is.na(xii)
                          as.double(xii * xii)
                      }
              }
              r
          })

setMethod("tcrossprod", signature(x = "diagonalMatrix", y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(x, t(y)))

for(.cl in c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix",
             "denseMatrix", "matrix", "vector")) {
setMethod("tcrossprod", signature(x = "diagonalMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_diagonal_matmult, x, y, FALSE,  TRUE, boolArith))

setMethod("tcrossprod", signature(x = .cl, y = "diagonalMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_diagonal_matmult, x, y, FALSE,  TRUE, boolArith))
}


## .... indMatrix ......................................................

setMethod("tcrossprod", signature(x = "indMatrix", y = "missing"),
          function(x, y = NULL, boolArith = TRUE, ...) {
              if(x@margin != 1L)
                  return(crossprod(t(x), boolArith = boolArith, ...))
              boolArith <- !is.na(boolArith) && boolArith
              r <- new(if(boolArith) "ngeMatrix" else "dgeMatrix")
              r@Dim <- x@Dim[c(1L, 1L)]
              r@Dimnames <- x@Dimnames[c(1L, 1L)]
              r@x <- as.vector(`storage.mode<-`(
                  .M2m(x), if(boolArith) "logical" else "double")[, x@perm])
              r
          })

for(.cl in c("Matrix", "matrix", "vector"))
setMethod("tcrossprod", signature(x = .cl, y = "indMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(x, t(y)))

setMethod("tcrossprod", signature(x = "indMatrix", y = "Matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              matmultDim(x@Dim, y@Dim, type = 3L)
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              if(y@margin != 1L)
                  r <- tcrossprod(.M2kind(x, l), y, boolArith = boolArith, ...)
              else {
                  r <- .M2kind(t(y)[x@perm, , drop = FALSE], l)
                  r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 3L)
              }
              r
          })

setMethod("tcrossprod", signature(x = "indMatrix", y = "matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              matmultDim(x@Dim, dim(y), type = 3L)
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              if(y@margin != 1L)
                  r <- tcrossprod(.M2kind(x, l), y, boolArith = boolArith, ...)
              else {
                  r <- .m2dense(t(y)[x@perm, , drop = FALSE], paste0(l, "ge"))
                  r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 3L)
              }
              r
          })

setMethod("tcrossprod", signature(x = "indMatrix", y = "vector"),
          function(x, y = NULL, boolArith = NA, ...) {
              d <- x@Dim
              m <- d[1L]
              k <- d[2L]
              if(k != (if(m == 1L) length(y) else 1L))
                  stop("non-conformable arguments")
              boolArith <- !is.na(boolArith) && boolArith
              h <- if(boolArith) isN0 else as.double
              r <- new(if(boolArith) "ngeMatrix" else "dgeMatrix")
              r@Dim <- d <- c(m, if(m == 1L) 1L else length(y))
              r@Dimnames <- c(x@Dimnames[1L], list(NULL))
              r@x <-
              if(m == 1L) {
                  if(x@margin == 1L)
                      h(y[x@perm])
                  else (if(boolArith) any else sum)(h(y))
              } else {
                  if(x@margin == 1L)
                      as.vector(matrix(h(y), m, length(y), byrow = TRUE))
                  else {
                      tmp <- array(if(boolArith) FALSE else 0, d)
                      tmp[y@perm, ] <- h(y)
                      dim(tmp) <- NULL
                      tmp
                  }
              }
              r
          })


## .... pMatrix ........................................................

setMethod("tcrossprod", signature(x = "pMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...) {
              boolArith <- !is.na(boolArith) && boolArith
              r <- new(if(boolArith) "ndiMatrix" else "ddiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames[c(1L, 1L)]
              r@diag <- "U"
              r
          })

setMethod("tcrossprod", signature(x = "pMatrix", y = "pMatrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <- new("pMatrix")
              r@Dim <- matmultDim(x@Dim, y@Dim, type = 2L)
              r@Dimnames <- matmultDN(x@Dimnames, y@Dimnames, type = 2L)
              r@perm <-
              if(y@margin != 1L)
                  y@perm[if(x@margin == 1L) x@perm else invertPerm(x@perm)]
              else {
                  r@margin <- 2L
                  (if(x@margin == 1L) invertPerm(x@perm) else x@perm)[y@perm]
              }
              r
          })

setMethod("tcrossprod", signature(x = "pMatrix", y = "Matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              matmultDim(x@Dim, y@Dim, type = 3L)
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .M2kind(t(y)[perm, , drop = FALSE], l)
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 3L)
              r
          })

setMethod("tcrossprod", signature(x = "pMatrix", y = "matrix"),
          function(x, y = NULL, boolArith = NA, ...) {
              matmultDim(x@Dim, dim(y), type = 3L)
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              perm <- if(x@margin == 1L) x@perm else invertPerm(x@perm)
              r <- .m2dense(t(y)[perm, , drop = FALSE], paste0(l, "ge"))
              r@Dimnames <- matmultDN(x@Dimnames, dimnames(y), type = 3L)
              r
          })

setMethod("tcrossprod", signature(x = "pMatrix", y = "vector"),
          function(x, y = NULL, boolArith = NA, ...) {
              if(x@Dim[2L] != 1L || length(y) != 1L)
                  stop("non-conformable arguments")
              l <- if(!is.na(boolArith) && boolArith) "n" else ","
              r <- .m2dense(y, paste0(l, "ge"), trans = TRUE)
              r@Dimnames <- c(x@Dimnames[1L], list(NULL))
              r
          })


## .... sparseVector ...................................................

setMethod("tcrossprod", signature(x = "sparseVector", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              tcrossprod(.V2C(x), boolArith = boolArith, ...))

setMethod("tcrossprod", signature(x = "sparseVector", y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(if(is.na(boolArith)) .M.kind(x) == "n" && .M.kind(y) == "n" else boolArith) `%&%` else `%*%`)(x, .tCRT(.V2C(y))))

for(.cl in c("Matrix", "matrix")) {
setMethod("tcrossprod", signature(x = "sparseVector", y = .cl),
          function(x, y = NULL, boolArith = NA, ...) {
              x <- if(dim(y)[2L] == length(x)) .tCRT(.V2C(x)) else .V2C(x)
              tcrossprod(x, y, boolArith = boolArith, ...)
          })

setMethod("tcrossprod", signature(x = .cl, y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...) {
              y <- if(dim(x)[1L] == 1L) .tCRT(.V2C(y)) else .V2C(y)
              tcrossprod(x, y, boolArith = boolArith, ...)
          })
}

setMethod("tcrossprod", signature(x = "sparseVector", y = "vector"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <-
              if(!is.na(boolArith) && boolArith)
                  x %&% .m2sparse(y, "ngR", trans = TRUE)
              else
                  x %*% .m2dense (y, ",ge", trans = TRUE)
              r@Dimnames <- list(NULL, NULL)
              r
          })

setMethod("tcrossprod", signature(x = "vector", y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(x, .tCRT(.V2C(y))))

rm(.cl)
