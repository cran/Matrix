## METHODS FOR GENERIC: coerce, as.*
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.M2kind <- function(from, kind = ".", sparse = NA)
    .Call(R_Matrix_as_kind, from, kind, sparse)

.M2gen <- function(from, kind = ".")
    .Call(R_Matrix_as_general, from, kind)
..M2gen <- function(from) # for setAs()
    .Call(R_Matrix_as_general, from, ".")

.M2sym <- function(from, ...) {
    if(isSymmetric(from, ...))
        forceSymmetric(from)
    else
        stop("matrix is not symmetric; consider forceSymmetric(.) or symmpart(.)")
}
..M2sym <- .M2sym # for setAs()
formals(..M2sym) <- formals(..M2sym)[-2L]
body(..M2sym)[[2L]][[2L]] <-
    body(..M2sym)[[2L]][[2L]][-3L]

.M2tri <- function(from, ...) {
    if(!(it <- isTriangular(from, ...)))
        stop("matrix is not triangular; consider triu(.) or tril(.)")
    else if(attr(it, "kind") == "U")
        triu(from)
    else
        tril(from)
}
..M2tri <- .M2tri # for setAs()
formals(..M2tri) <- formals(..M2tri)[-2L]
body(..M2tri)[[2L]][[2L]][[2L]][[2L]][[3L]] <-
    body(..M2tri)[[2L]][[2L]][[2L]][[2L]][[3L]][-3L]

.M2diag <- function(from) {
    if (!isDiagonal(from))
        stop("matrix is not diagonal; consider Diagonal(x=diag(.))")
    forceDiagonal(from)
}

.M2v <- function(from)
    .Call(R_Matrix_as_vector, from)

.M2m <- function(from)
    .Call(R_Matrix_as_matrix, from)

.M2unpacked <- function(from)
    .Call(R_Matrix_as_unpacked, from)

.M2packed <- function(from)
    .Call(R_Matrix_as_packed, from)

.M2C <- function(from)
    .Call(R_Matrix_as_Csparse, from)

.M2R <- function(from)
    .Call(R_Matrix_as_Rsparse, from)

.M2T <- function(from)
    .Call(R_Matrix_as_Tsparse, from)

.sparse2dense <- function(from, packed = FALSE)
    .Call(R_sparse_as_dense, from, packed)

.diag2dense <- function(from, shape = "t", packed = FALSE, uplo = "U")
    .Call(R_diagonal_as_dense, from, shape, packed, uplo)

.ind2dense <- function(from, kind = "n")
    .Call(R_index_as_dense, from, kind)

.dense2sparse <- function(from, repr = "C")
    .Call(R_dense_as_sparse, from, repr)

.diag2sparse <- function(from, shape = "t", repr = "C", uplo = "U")
    .Call(R_diagonal_as_sparse, from, shape, repr, uplo)

.ind2sparse <- function(from, kind = "n", repr = ".")
    .Call(R_index_as_sparse, from, kind, repr)

.m2dense <- function(from, class, uplo = "U", diag = "N")
    .Call(R_matrix_as_dense, from, class, uplo, diag)

.m2dense.checking <- function(from, kind = ".", ...) {
    switch(typeof(from), logical =, integer =, double = NULL,
           stop(gettextf("matrix of invalid type \"%s\" to .m2dense.checking()",
                         typeof(from)),
                domain = NA))
    if(kind != ".") {
        ## These must happen before isSymmetric() call
        storage.mode(from) <-
            switch(kind, n =, l = "logical", d = "double",
                   stop(gettextf("invalid kind \"%s\" to .m2dense.checking()",
                                 kind),
                        domain = NA))
        if(kind == "n" && anyNA(from))
            from[is.na(from)] <- TRUE
    }
    if(isSymmetric(from, ...))
        .m2dense(from, paste0(kind, "sy"), "U", NULL)
    else if(it <- isTriangular(from))
        .m2dense(from, paste0(kind, "tr"), attr(it, "kind"), "N")
    else
        .m2dense(from, paste0(kind, "ge"), NULL, NULL)
}

.m2sparse <- function(from, class, uplo = "U", diag = "N")
    .Call(R_matrix_as_sparse, from, class, uplo, diag)

.m2sparse.checking <- function(from, kind = ".", repr = "C", ...) {
    switch(typeof(from), logical =, integer =, double = NULL,
           stop(gettextf("matrix of invalid type \"%s\" to .m2sparse.checking()",
                         typeof(from)),
                domain = NA))
    if(kind != ".") {
        ## These must happen before isSymmetric() call
        storage.mode(from) <-
            switch(kind, n =, l = "logical", d = "double",
                   stop(gettextf("invalid kind \"%s\" to .m2sparse.checking()",
                                 kind),
                        domain = NA))
        if(kind == "n" && anyNA(from))
            from[is.na(from)] <- TRUE
    }
    if(isSymmetric(from, ...))
        .m2sparse(from, paste0(kind, "s", repr), "U", NULL)
    else if(it <- isTriangular(from))
        .m2sparse(from, paste0(kind, "t", repr), attr(it, "kind"), "N")
    else
        .m2sparse(from, paste0(kind, "g", repr), NULL, NULL)
}


## ==== From Matrix to vector ==========================================

## Need 'base' functions calling as.*() to dispatch to our S4 methods:
if (FALSE) {
## 2023-08-10: breaks iGraphMatch, mcmcsae, mcompanion
## which define proper subclasses of Matrix not extending
## any of _our_ proper subclasses of Matrix
as.vector.Matrix <- function(x, mode = "any") as.vector(.M2v(x), mode)
as.matrix.Matrix <- function(x, ...)                    .M2m(x)
 as.array.Matrix <- function(x, ...)                    .M2m(x)
} else {
## Hence, for now ...
as.vector.Matrix <- function(x, mode = "any") as.vector(as(x, "matrix"), mode)
as.matrix.Matrix <- function(x, ...)                    as(x, "matrix")
 as.array.Matrix <- function(x, ...)                    as(x, "matrix")
}

setMethod("as.vector" , signature(x = "Matrix"),
          as.vector.Matrix)
setMethod("as.matrix" , signature(x = "Matrix"),
          as.matrix.Matrix)
setMethod("as.array"  , signature(x = "Matrix"),
           as.array.Matrix)
setMethod("as.logical", signature(x = "Matrix"),
          function(x, ...) as.logical(.M2v(x)))
setMethod("as.integer", signature(x = "Matrix"),
          function(x, ...) as.integer(.M2v(x)))
setMethod("as.double" , signature(x = "Matrix"),
          function(x, ...) as.double (.M2v(x)))
setMethod("as.numeric", signature(x = "Matrix"),
          function(x, ...) as.numeric(.M2v(x)))
setMethod("as.complex", signature(x = "Matrix"),
          function(x, ...) as.complex(.M2v(x)))

setAs("Matrix",  "vector", .M2v)
setAs("Matrix",  "matrix", .M2m)
setAs("Matrix",   "array", .M2m)
setAs("Matrix", "logical", function(from) as.logical(.M2v(from)))
setAs("Matrix", "integer", function(from) as.integer(.M2v(from)))
setAs("Matrix",  "double", function(from) as.double (.M2v(from)))
setAs("Matrix", "numeric", function(from) as.numeric(.M2v(from)))
setAs("Matrix", "complex", function(from) as.complex(.M2v(from)))


## ==== From vector to Matrix ==========================================

setAs("matrix", "Matrix",
      function(from) {
          if(isDiagonal(from))
              forceDiagonal(from)
          else if(.sparseDefault(from))
              .m2sparse.checking(from, ".", "C")
          else .m2dense.checking(from, ".")
      })
setAs("vector", "Matrix",
      function(from) {
          if(is.object(from)) # e.g., data.frame
              as(as.matrix(from), "Matrix")
          else if(.sparseDefault(from))
              .m2sparse(from, ".gC")
          else .m2dense(from, ".ge")
      })
setAs("ANY", "Matrix",
      function(from) as(as(from, "matrix"), "Matrix"))


## ==== To "kind" ======================================================

setAs("Matrix", "nMatrix",
      function(from) .M2kind(from, "n",    NA))
setAs("Matrix", "lMatrix",
      function(from) .M2kind(from, "l",    NA))
setAs("Matrix", "dMatrix",
      function(from) .M2kind(from, "d",    NA))

setAs("Matrix", "ndenseMatrix",
      function(from) .M2kind(from, "n", FALSE))
setAs("Matrix", "ldenseMatrix",
      function(from) .M2kind(from, "l", FALSE))
setAs("Matrix", "ddenseMatrix",
      function(from) .M2kind(from, "d", FALSE))

setAs("Matrix", "nsparseMatrix",
      function(from) .M2kind(from, "n",  TRUE))
setAs("Matrix", "lsparseMatrix",
      function(from) .M2kind(from, "l",  TRUE))
setAs("Matrix", "dsparseMatrix",
      function(from) .M2kind(from, "d",  TRUE))

setAs("matrix", "nMatrix",
      function(from) {
          if(.sparseDefault(from))
              .m2sparse.checking(from, "n", "C")
          else .m2dense.checking(from, "n")
      })
setAs("matrix", "lMatrix",
      function(from) {
          if(isDiagonal(from))
              forceDiagonal(`storage.mode<-`(from, "logical"))
          else if(.sparseDefault(from))
              .m2sparse.checking(from, "l", "C")
          else .m2dense.checking(from, "l")
      })
setAs("matrix", "dMatrix",
      function(from) {
          if(isDiagonal(from))
              forceDiagonal(`storage.mode<-`(from, "double"))
          else if(.sparseDefault(from))
              .m2sparse.checking(from, "d", "C")
          else .m2dense.checking(from, "d")
      })

setAs("matrix", "ndenseMatrix",
      function(from) .m2dense.checking(from, "n"))
setAs("matrix", "ldenseMatrix",
      function(from) .m2dense.checking(from, "l"))
setAs("matrix", "ddenseMatrix",
      function(from) .m2dense.checking(from, "d"))

setAs("matrix", "nsparseMatrix",
      function(from) .m2sparse.checking(from, "n", "C"))
setAs("matrix", "lsparseMatrix",
      function(from) .m2sparse.checking(from, "l", "C"))
setAs("matrix", "dsparseMatrix",
      function(from) .m2sparse.checking(from, "d", "C"))

setAs("vector", "nMatrix",
      function(from) {
          if(.sparseDefault(from))
              .m2sparse(from, "ngC")
          else .m2dense(from, "nge")
      })
setAs("vector", "lMatrix",
      function(from) {
          if(.sparseDefault(from))
              .m2sparse(from, "lgC")
          else .m2dense(from, "lge")
      })
setAs("vector", "dMatrix",
      function(from) {
          if(.sparseDefault(from))
              .m2sparse(from, "dgC")
          else .m2dense(from, "dge")
      })

setAs("vector", "ndenseMatrix",
      function(from) .m2dense(from, "nge"))
setAs("vector", "ldenseMatrix",
      function(from) .m2dense(from, "lge"))
setAs("vector", "ddenseMatrix",
      function(from) .m2dense(from, "dge"))

setAs("vector", "nsparseMatrix",
      function(from) .m2sparse(from, "ngC"))
setAs("vector", "lsparseMatrix",
      function(from) .m2sparse(from, "lgC"))
setAs("vector", "dsparseMatrix",
      function(from) .m2sparse(from, "dgC"))


## ==== To "shape" =====================================================

..m2gen <- function(from) .Call(R_matrix_as_dense, from, ".ge", NULL, NULL)

setAs("Matrix",    "generalMatrix", ..M2gen)
setAs("matrix",    "generalMatrix", ..m2gen)
setAs("vector",    "generalMatrix", ..m2gen)

setAs("Matrix",  "symmetricMatrix", ..M2sym)
setAs("matrix",  "symmetricMatrix", ..M2sym)

setAs("Matrix", "triangularMatrix", ..M2tri)
setAs("matrix", "triangularMatrix", ..M2tri)

setAs("Matrix",   "diagonalMatrix",  .M2diag)
setAs("matrix",   "diagonalMatrix",  .M2diag)

rm(..m2gen)

setAs("diagonalMatrix",  "symmetricMatrix",
      function(from) {
          if(!isSymmetricDN(from@Dimnames))
              stop("matrix is not symmetric; consider forceSymmetric(.) or symmpart(.)")
          .diag2sparse(from, "s", "C", "U")
      })

setAs("diagonalMatrix", "triangularMatrix",
      function(from)
          .diag2sparse(from, "t", "C", "U"))


## ==== To "representation" ============================================

setAs("Matrix",    "denseMatrix", .M2unpacked)
setAs("Matrix", "unpackedMatrix", .M2unpacked)
setAs("Matrix",   "packedMatrix", .M2packed)
setAs("Matrix",   "sparseMatrix", .M2C)
setAs("Matrix",  "CsparseMatrix", .M2C)
setAs("Matrix",  "RsparseMatrix", .M2R)
setAs("Matrix",  "TsparseMatrix", .M2T)

## Do test for structure:
setAs("generalMatrix", "packedMatrix", function(from) pack(from))

setAs("matrix",    "denseMatrix",
      function(from)  .m2dense.checking(from, "."))
setAs("matrix", "unpackedMatrix",
      function(from)  .m2dense.checking(from, "."))
setAs("matrix",   "packedMatrix",
      function(from) pack(from))
setAs("matrix",   "sparseMatrix",
      function(from) .m2sparse.checking(from, ".", "C"))
setAs("matrix",  "CsparseMatrix",
      function(from) .m2sparse.checking(from, ".", "C"))
setAs("matrix",  "RsparseMatrix",
      function(from) .m2sparse.checking(from, ".", "R"))
setAs("matrix",  "TsparseMatrix",
      function(from) .m2sparse.checking(from, ".", "T"))

## Many people want this coercion to be available and fast:
setAs("matrix", "dgCMatrix",
      function(from) .m2sparse(from, "dgC"))

setAs("vector",    "denseMatrix",
      function(from)
          if(is.object(from)) # e.g., data.frame
              as(as.matrix(from),  "denseMatrix")
          else .m2dense(from, ".ge"))
setAs("vector", "unpackedMatrix",
      function(from)  .m2dense(from, ".ge"))
setAs("vector",   "sparseMatrix",
      function(from)
          if(is.object(from)) # e.g., data.frame
              as(as.matrix(from), "sparseMatrix")
          else .m2sparse(from, ".gC"))
setAs("vector",  "CsparseMatrix",
      function(from) .m2sparse(from, ".gC"))
setAs("vector",  "RsparseMatrix",
      function(from) .m2sparse(from, ".gR"))
setAs("vector",  "TsparseMatrix",
      function(from) .m2sparse(from, ".gT"))

setAs("ANY",  "denseMatrix",
      function(from) as(as(from, "matrix"),  "denseMatrix"))
setAs("ANY", "sparseMatrix",
      function(from) as(as(from, "matrix"), "sparseMatrix"))

setAs("Matrix", "indMatrix",
      function(from) as(as(from, "nsparseMatrix"), "indMatrix"))
setAs("matrix", "indMatrix",
      function(from) as(as(from, "nsparseMatrix"), "indMatrix"))

setAs("Matrix",   "pMatrix",
      function(from) as(as(from, "nsparseMatrix"),   "pMatrix"))
setAs("matrix",   "pMatrix",
      function(from) as(as(from, "nsparseMatrix"),   "pMatrix"))
