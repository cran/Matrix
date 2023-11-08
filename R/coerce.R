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

.diag2dense <- function(from, kind = ".", shape = "t", packed = FALSE, uplo = "U")
    .Call(R_diagonal_as_dense, from, kind, shape, packed, uplo)

.ind2dense <- function(from, kind = "n")
    .Call(R_index_as_dense, from, kind)

.dense2sparse <- function(from, repr = "C")
    .Call(R_dense_as_sparse, from, repr)

.diag2sparse <- function(from, kind = ".", shape = "t", repr = "C", uplo = "U")
    .Call(R_diagonal_as_sparse, from, kind, shape, repr, uplo)

.ind2sparse <- function(from, kind = "n", repr = ".")
    .Call(R_index_as_sparse, from, kind, repr)

.m2dense <- function(from, class = ".ge", uplo = "U", diag = "N",
                     trans = FALSE)
    .Call(R_matrix_as_dense, from, class, uplo, diag, trans)

.m2dense.checking <- function(from, kind = ".", ...) {
    switch(typeof(from), logical =, integer =, double = NULL,
           stop(gettextf("invalid type \"%s\" in '%s'",
                         typeof(from), ".m2dense.checking"),
                domain = NA))
    if(kind != ".") {
        ## These must happen before isSymmetric() call
        storage.mode(from) <-
            switch(kind, n =, l = "logical", d = "double",
                   stop(gettextf("invalid %s=\"%s\" to '%s'",
                                 "kind", kind, ".m2dense.checking"),
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

.m2sparse <- function(from, class = ".gC", uplo = "U", diag = "N",
                      trans = FALSE)
    .Call(R_matrix_as_sparse, from, class, uplo, diag, trans)

.m2sparse.checking <- function(from, kind = ".", repr = "C", ...) {
    switch(typeof(from), logical =, integer =, double = NULL,
           stop(gettextf("invalid type \"%s\" in '%s'",
                         typeof(from), ".m2sparse.checking"),
                domain = NA))
    if(kind != ".") {
        ## These must happen before isSymmetric() call
        storage.mode(from) <-
            switch(kind, n =, l = "logical", d = "double",
                   stop(gettextf("invalid %s=\"%s\" to '%s'",
                                 "kind", kind, ".m2sparse.checking"),
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

.V2kind <- function(from, kind = ".") {
    if(kind == ".")
        return(from)
    kind. <- .M.kind(from)
    if(kind == ",")
        kind <- if(kind. == "z") "z" else "d"
    if(kind == kind.)
        return(from)
    to <- new(paste0(kind, "sparseVector"))
    to@length <- from@length
    to@i <- from@i
    if(kind != "n")
        to@x <-
            if(kind. == "n")
                rep.int(switch(kind, "l" = TRUE, "i" = 1L, "d" = 1, "z" = 1+0i), length(from@i))
            else as.vector(from@x, typeof(to@x))
    to
}

.V2v <- function(from) {
    if(.M.kind(from) != "n") {
        to <- vector(typeof(from@x), from@length)
        to[from@i] <- from@x
    } else {
        to <- logical(from@length)
        to[from@i] <- TRUE
    }
    to
}

.V2m <- function(from) {
    if(is.double(m <- from@length)) {
        if(m > .Machine$integer.max)
            stop(gettextf("dimensions cannot exceed %s", "2^31-1"), domain = NA)
        m <- as.integer(m)
    }
    to <- .V2v(from)
    dim(to) <- c(m, 1L)
    to
}

.V2a <- function(from) {
    if(is.double(m <- from@length)) {
        if(m > .Machine$integer.max)
            stop(gettextf("dimensions cannot exceed %s", "2^31-1"), domain = NA)
        m <- as.integer(m)
    }
    to <- .V2v(from)
    dim(to) <- m
    to
}

.V2unpacked <- function(from) {
    if(is.double(m <- from@length)) {
        if(m > .Machine$integer.max)
            stop(gettextf("dimensions cannot exceed %s", "2^31-1"), domain = NA)
        m <- as.integer(m)
    }
    kind <- .M.kind(from)
    to <- new(paste0(if(kind == "i") "d" else kind, "geMatrix"))
    to@Dim <- c(m, 1L)
    to@x <- replace(vector(typeof(to@x), m), from@i,
                    if(kind == "n") TRUE else from@x)
    to
}

.V2C <- function(from) {
    if(is.double(m <- from@length)) {
        if(m > .Machine$integer.max)
            stop(gettextf("dimensions cannot exceed %s", "2^31-1"), domain = NA)
        m <- as.integer(m)
    }
    kind <- .M.kind(from)
    to <- new(paste0(if(kind == "i") "d" else kind, "gCMatrix"))
    to@Dim <- c(m, 1L)
    to@p <- c(0L, length(from@i))
    to@i <- as.integer(from@i) - 1L
    if(kind != "n")
        to@x <- if(kind == "i") as.double(from@x) else from@x
    to
}

.V2R <- function(from) {
    if(is.double(m <- from@length)) {
        if(m > .Machine$integer.max)
            stop(gettextf("dimensions cannot exceed %s", "2^31-1"), domain = NA)
        m <- as.integer(m)
    }
    kind <- .M.kind(from)
    to <- new(paste0(if(kind == "i") "d" else kind, "gRMatrix"))
    to@Dim <- c(m, 1L)
    to@p <- c(0L, cumsum(replace(logical(m), from@i, TRUE)))
    to@j <- integer(length(from@i))
    if(kind != "n")
        to@x <- if(kind == "i") as.double(from@x) else from@x
    to
}

.V2T <- function(from) {
    if(is.double(m <- from@length)) {
        if(m > .Machine$integer.max)
            stop(gettextf("dimensions cannot exceed %s", "2^31-1"), domain = NA)
        m <- as.integer(m)
    }
    kind <- .M.kind(from)
    to <- new(paste0(if(kind == "i") "d" else kind, "gTMatrix"))
    to@Dim <- c(m, 1L)
    to@i <- as.integer(from@i) - 1L
    to@j <- integer(length(from@i))
    if(kind != "n")
        to@x <- if(kind == "i") as.double(from@x) else from@x
    to
}

## FIXME: define R_Matrix_as_sparseVector in ../src/coerce.c and use here
.M2V <- function(from) {
    repr <- .M.repr(from)
    if(repr == "u" || repr == "p")
        return(.Call( v2spV, .M2v(from)))
    if(repr == "C" || repr == "R")
        return(.Call(CR2spV,      from ))
    if(repr == "T")
        return(.Call(CR2spV, .M2C(from)))
    if(repr != "d" && repr != "i") {
        if(is.object(from))
            stop(gettextf("invalid class \"%s\" in '%s'",
                          class(from)[1L], ".M2V"),
                 domain = NA)
        else
            stop(gettextf("invalid type \"%s\" in '%s'",
                          typeof(from), ".M2V"),
                 domain = NA)
    }
    d <- from@Dim
    m <- d[1L]
    n <- d[2L]
    mn <- prod(d)
    if(mn <= .Machine$integer.max)
        mn <- as.integer(mn)
    else if(mn > 0x1p+53)
        stop(gettextf("%s length cannot exceed %s", "sparseVector", "2^53"),
             domain = NA)
    kind <- .M.kind(from)
    to <- new(paste0(kind, "sparseVector"))
    to@length <- mn
    to@i <-
        if(repr == "d") {
            if(kind == "n" && from@diag == "N")
                indDiag(n)[from@x | is.na(from@x)]
            else
                indDiag(n)
        } else if(is.integer(mn)) {
            if(from@margin == 1L)
                seq.int(to =   0L, by = 1L, length.out = m) +
                    from@perm * m
            else
                seq.int(from = 0L, by =  m, length.out = n) +
                    from@perm
        } else {
            if(from@margin == 1L)
                seq.int(  to = 0, by =            1, length.out = m) +
                    from@perm * as.double(m)
            else
                seq.int(from = 0, by = as.double(m), length.out = n) +
                    as.double(from@perm)
        }
    if(kind != "n")
        to@x <-
            if(from@diag == "N")
                from@x
            else rep.int(switch(kind, "l" = TRUE, "i" = 1L, "d" = 1, "z" = 1+0i), n)
    to
}

.m2V <- function(from, kind = ".") {
    to <- .Call(v2spV, from)
    if(kind == ".")
        to
    else {
        to. <- new(paste0(kind, "sparseVector"))
        to.@length <- to@length
        to.@i <- to@i
        if(kind != "n")
            to.@x <- as.vector(to@x, typeof(to.@x))
        to.
    }
}


## ==== To vector ======================================================

## Need 'base' functions calling as.*() to dispatch to our S4 methods:
if (FALSE) {
## 2023-08-10: breaks iGraphMatch, mcmcsae, mcompanion
## which define proper subclasses of Matrix not extending
## any of _our_ proper subclasses of Matrix
as.matrix.Matrix <- function(x, ...) .M2m(x)
 as.array.Matrix <- function(x, ...) .M2m(x)
} else {
as.matrix.Matrix <- function(x, ...) as(x, "matrix")
 as.array.Matrix <- function(x, ...) as(x, "matrix")
setAs("Matrix", "matrix", .M2m)
}
as.matrix.sparseVector <- function(x, ...) .V2m(x)
 as.array.sparseVector <- function(x, ...) .V2a(x)

setMethod("as.vector" , signature(x = "Matrix"),
          function(x, mode = "any") as.vector(.M2v(x), mode))
setMethod("as.matrix" , signature(x = "Matrix"),
          as.matrix.Matrix)
setMethod("as.array"  , signature(x = "Matrix"),
          as.array.Matrix)
setMethod("as.logical", signature(x = "Matrix"),
          function(x, ...) as.logical(.M2v(x)))
setMethod("as.integer", signature(x = "Matrix"),
          function(x, ...) as.integer(.M2v(x)))
setMethod("as.numeric", signature(x = "Matrix"),
          function(x, ...) as.numeric(.M2v(x)))
setMethod("as.complex", signature(x = "Matrix"),
          function(x, ...) as.complex(.M2v(x)))

setMethod("as.vector" , signature(x = "sparseVector"),
          function(x, mode = "any") as.vector(.V2v(x), mode))
setMethod("as.matrix" , signature(x = "sparseVector"),
          as.matrix.sparseVector)
setMethod("as.array"  , signature(x = "sparseVector"),
          as.array.sparseVector)
setMethod("as.logical", signature(x = "sparseVector"),
          function(x, ...) as.logical(.V2v(x)))
setMethod("as.integer", signature(x = "sparseVector"),
          function(x, ...) as.integer(.V2v(x)))
setMethod("as.numeric", signature(x = "sparseVector"),
          function(x, ...) as.numeric(.V2v(x)))
setMethod("as.complex", signature(x = "sparseVector"),
          function(x, ...) as.complex(.V2v(x)))


## ==== To Matrix ======================================================

setAs("sparseVector", "Matrix",
      .V2C)
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
          if(is.object(from) && length(dim(from)) == 2L) # e.g., data.frame
              as(as.matrix(from), "Matrix")
          else if(.sparseDefault(from))
              .m2sparse(from, ".gC")
          else .m2dense(from, ".ge")
      })
setAs(   "ANY", "Matrix",
      function(from) as(as(from, "matrix"), "Matrix"))

if(FALSE) {
## MJ: not yet ... existing as(<CHMfactor>, "Matrix") must become defunct first
setAs("MatrixFactorization", "Matrix",
      function(from) {
          n <- length(x <- expand2(from))
          to <- x[[1L]]
          if(n >= 2L) for(i in 2L:n) to <- to %*% x[[i]]
          to
      })
}


## ==== To sparseVector ================================================

setAs("Matrix", "sparseVector",
      function(from) .M2V(from))
setAs("matrix", "sparseVector",
      function(from) .m2V(from))
setAs("vector", "sparseVector",
      function(from) .m2V(from))
setAs(   "ANY", "sparseVector",
      function(from) as(as.vector(from), "sparseVector"))


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

setAs("sparseVector", "nsparseVector",
      function(from) .V2kind(from, "n"))
setAs("sparseVector", "lsparseVector",
      function(from) .V2kind(from, "l"))
setAs("sparseVector", "isparseVector",
      function(from) .V2kind(from, "i"))
setAs("sparseVector", "dsparseVector",
      function(from) .V2kind(from, "d"))
setAs("sparseVector", "zsparseVector",
      function(from) .V2kind(from, "z"))

setAs("vector", "nsparseVector",
      function(from) .m2V(from, "n"))
setAs("vector", "lsparseVector",
      function(from) .m2V(from, "l"))
setAs("vector", "isparseVector",
      function(from) .m2V(from, "i"))
setAs("vector", "dsparseVector",
      function(from) .m2V(from, "d"))
setAs("vector", "zsparseVector",
      function(from) .m2V(from, "z"))


## ==== To "shape" =====================================================

..m2gen <- function(from) .m2dense(from, ".ge")

setAs(      "Matrix", "generalMatrix", ..M2gen)
setAs(      "matrix", "generalMatrix", ..m2gen)
setAs(      "vector", "generalMatrix", ..m2gen)
setAs("sparseVector", "generalMatrix",  .V2C)

setAs("Matrix",  "symmetricMatrix", ..M2sym)
setAs("matrix",  "symmetricMatrix", ..M2sym)

setAs("Matrix", "triangularMatrix", ..M2tri)
setAs("matrix", "triangularMatrix", ..M2tri)

rm(..m2gen)

setAs("diagonalMatrix",  "symmetricMatrix",
      function(from) {
          if(!isSymmetricDN(from@Dimnames))
              stop("matrix is not symmetric; consider forceSymmetric(.) or symmpart(.)")
          .diag2sparse(from, ".", "s", "C", "U")
      })

setAs("diagonalMatrix", "triangularMatrix",
      function(from)
          .diag2sparse(from, ".", "t", "C", "U"))


## ==== To "representation" ============================================

setAs("Matrix",    "denseMatrix", .M2unpacked)
setAs("Matrix", "unpackedMatrix", .M2unpacked)
setAs("Matrix",   "packedMatrix", .M2packed)
setAs("Matrix",   "sparseMatrix", .M2C)
setAs("Matrix",  "CsparseMatrix", .M2C)
setAs("Matrix",  "RsparseMatrix", .M2R)
setAs("Matrix",  "TsparseMatrix", .M2T)

## Do test for structure:
## FIXME: wrongly assumes that methods are defined for pack(<sparseMatrix>) ...
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
          if(is.object(from) && length(dim(from)) == 2L) # e.g., data.frame
              as(as.matrix(from),  "denseMatrix")
          else .m2dense(from, ".ge"))
setAs("vector", "unpackedMatrix",
      function(from)  .m2dense(from, ".ge"))
setAs("vector",   "sparseMatrix",
      function(from)
          if(is.object(from) && length(dim(from)) == 2L) # e.g., data.frame
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

setAs("sparseVector",    "denseMatrix", .V2unpacked)
setAs("sparseVector", "unpackedMatrix", .V2unpacked)
setAs("sparseVector",   "sparseMatrix", .V2C)
setAs("sparseVector",  "CsparseMatrix", .V2C)
setAs("sparseVector",  "RsparseMatrix", .V2R)
setAs("sparseVector",  "TsparseMatrix", .V2T)

setAs("Matrix", "diagonalMatrix", .M2diag)
setAs("matrix", "diagonalMatrix", .M2diag)

setAs("Matrix", "indMatrix",
      function(from) as(as(from, "nsparseMatrix"), "indMatrix"))
setAs("matrix", "indMatrix",
      function(from) as(as(from, "nsparseMatrix"), "indMatrix"))

setAs("Matrix",   "pMatrix",
      function(from) as(as(from, "nsparseMatrix"),   "pMatrix"))
setAs("matrix",   "pMatrix",
      function(from) as(as(from, "nsparseMatrix"),   "pMatrix"))
